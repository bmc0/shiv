/*
 * Copyright (c) 2013-2017 Michael Barbour <barbour.michael.0@gmail.com>
 * 
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

/*
 * This file implements macros for dealing with doubly-linked lists.
 *
 * Each node structure must include "next" and "prev" pointers:
 *     struct node *next, *prev;
 *
 * Each list structure must include "head" and "tail" pointers:
 *     struct node *head, *tail;
 * The "head" and "tail" pointers must be initialized to NULL before inserting
 * elements into the list.
 *
 * All macros expect pointers.
 */

#ifndef _LIST_H
#define _LIST_H

#define LIST_FOREACH(_l, _n) \
	for \
	( \
		(_n) = (_l)->head; \
		(_n) != NULL; \
		(_n) = (_n)->next \
	)

#define LIST_ADD_HEAD(_l, _n) do { \
	if ((_l)->head == NULL) { \
		(_l)->head = (_l)->tail = (_n); \
		(_n)->next = (_n)->prev = NULL; \
	} \
	else { \
		(_n)->prev = NULL; \
		(_n)->next = (_l)->head; \
		(_l)->head->prev = (_n); \
		(_l)->head = (_n); \
	} \
} while(0)

#define LIST_ADD_TAIL(_l, _n) do { \
	if ((_l)->tail == NULL) { \
		(_l)->head = (_l)->tail = (_n); \
		(_n)->next = (_n)->prev = NULL; \
	} \
	else { \
		(_n)->next = NULL; \
		(_n)->prev = (_l)->tail; \
		(_l)->tail->next = (_n); \
		(_l)->tail = (_n); \
	} \
} while(0)

/* insert _n after _p */
#define LIST_INSERT(_l, _n, _p) do { \
	if ((_p) == (_l)->tail) LIST_ADD_TAIL(_l, _n); \
	else { \
		(_n)->prev = (_p); \
		(_n)->next = (_p)->next; \
		(_p)->next->prev = (_n); \
		(_p)->next = (_n); \
	} \
} while(0)

#define LIST_REMOVE_HEAD(_l) do { \
	if ((_l)->head != NULL) { \
		if ((_l)->head->next == NULL) { \
			(_l)->head = NULL; \
			(_l)->tail = NULL; \
		} \
		else { \
			(_l)->head->next->prev = NULL; \
			(_l)->head = (_l)->head->next; \
		} \
	} \
} while(0)

#define LIST_REMOVE_TAIL(_l) do { \
	if ((_l)->tail != NULL) { \
		if ((_l)->tail->prev == NULL) { \
			(_l)->head = NULL; \
			(_l)->tail = NULL; \
		} \
		else { \
			(_l)->tail->prev->next = NULL; \
			(_l)->tail = (_l)->tail->prev; \
		} \
	} \
} while(0)

#define LIST_REMOVE(_l, _n) do { \
	if ((_n) == (_l)->head) LIST_REMOVE_HEAD(_l); \
	else if ((_n) == (_l)->tail) LIST_REMOVE_TAIL(_l); \
	else { \
		(_n)->prev->next = (_n)->next; \
		(_n)->next->prev = (_n)->prev; \
	} \
} while(0)

#endif
