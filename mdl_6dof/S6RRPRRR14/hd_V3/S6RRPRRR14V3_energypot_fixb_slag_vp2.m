% Calculate potential energy for
% S6RRPRRR14V3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRR14V3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(1,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR14V3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_energypot_fixb_slag_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14V3_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR14V3_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:03:18
% EndTime: 2019-04-12 15:03:18
% DurationCPUTime: 0.21s
% Computational Cost: add. (103->47), mult. (232->52), div. (0->0), fcn. (240->10), ass. (0->26)
t74 = -mrSges(3,2) + mrSges(4,3) + mrSges(5,3) + (m(4) + m(5) + m(6) + m(7)) * qJ(3);
t53 = sin(qJ(2));
t54 = sin(qJ(1));
t73 = t53 * t54;
t57 = cos(qJ(4));
t72 = t53 * t57;
t59 = cos(qJ(1));
t71 = t53 * t59;
t58 = cos(qJ(2));
t70 = t54 * t58;
t52 = sin(qJ(4));
t69 = t59 * t52;
t68 = t59 * t57;
t67 = -mrSges(3,1) - mrSges(4,1);
t66 = mrSges(6,2) - mrSges(7,3);
t65 = mrSges(2,2) - mrSges(3,3) - mrSges(4,2);
t50 = sin(qJ(6));
t55 = cos(qJ(6));
t62 = -t55 * mrSges(7,1) + t50 * mrSges(7,2) - mrSges(6,1);
t61 = -t50 * mrSges(7,1) - t55 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t60 = -t74 * t53 + t67 * t58 - mrSges(2,1);
t56 = cos(qJ(5));
t51 = sin(qJ(5));
t49 = t54 * t52 + t58 * t68;
t47 = t57 * t70 - t69;
t1 = (-mrSges(1,3) - mrSges(2,3) + t66 * (t51 * t72 + t58 * t56) + t62 * (-t58 * t51 + t56 * t72) + (-mrSges(5,1) * t57 + t61 * t52 + t67) * t53 + t74 * t58) * g(3) + (-t47 * mrSges(5,1) - mrSges(1,2) + t62 * (t47 * t56 + t51 * t73) + t61 * (t52 * t70 + t68) + t66 * (t47 * t51 - t56 * t73) - t65 * t59 + t60 * t54) * g(2) + (-t49 * mrSges(5,1) - mrSges(1,1) + t62 * (t49 * t56 + t51 * t71) + t61 * (-t54 * t57 + t58 * t69) + t66 * (t49 * t51 - t56 * t71) + t65 * t54 + t60 * t59) * g(1);
U  = t1;
