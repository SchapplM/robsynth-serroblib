% Calculate Gravitation load on the joints for
% S6RRRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRP10_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP10_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP10_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP10_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:29:30
% EndTime: 2019-03-09 17:29:34
% DurationCPUTime: 1.60s
% Computational Cost: add. (837->133), mult. (1723->183), div. (0->0), fcn. (2050->12), ass. (0->68)
t150 = -mrSges(5,3) - mrSges(6,3) - mrSges(7,2);
t149 = -m(5) * qJ(4) + mrSges(4,2);
t148 = -t149 - t150;
t70 = sin(pkin(11));
t72 = cos(pkin(11));
t147 = m(5) * pkin(3) + t72 * mrSges(5,1) - t70 * mrSges(5,2) + mrSges(4,1);
t100 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t69 = pkin(11) + qJ(5);
t67 = cos(t69);
t145 = -t100 * t67 - t147;
t143 = -t70 * mrSges(5,1) - t72 * mrSges(5,2) + mrSges(3,2) - mrSges(4,3);
t142 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t134 = m(4) + m(5);
t133 = m(6) + m(7);
t65 = pkin(4) * t72 + pkin(3);
t73 = -pkin(10) - qJ(4);
t74 = sin(qJ(3));
t76 = cos(qJ(3));
t138 = -t65 * t76 + t73 * t74;
t137 = -m(5) * pkin(9) + t143;
t128 = cos(qJ(1));
t71 = sin(pkin(6));
t108 = t71 * t128;
t126 = sin(qJ(1));
t127 = cos(qJ(2));
t75 = sin(qJ(2));
t111 = cos(pkin(6));
t95 = t111 * t128;
t49 = t126 * t127 + t75 * t95;
t22 = -t74 * t108 + t49 * t76;
t48 = t126 * t75 - t127 * t95;
t66 = sin(t69);
t1 = t22 * t66 - t48 * t67;
t136 = t22 * t67 + t48 * t66;
t135 = t147 * t76 + t148 * t74 + mrSges(3,1);
t132 = m(5) + t133;
t131 = (m(4) + t132) * pkin(9) - t143;
t130 = t142 * t66 - t145;
t129 = pkin(4) * t70;
t122 = t66 * t76;
t121 = t67 * t76;
t120 = t71 * t75;
t107 = t71 * t127;
t113 = pkin(2) * t107 + pkin(9) * t120;
t106 = t71 * t126;
t112 = t128 * pkin(1) + pkin(8) * t106;
t94 = t111 * t126;
t51 = t128 * t127 - t75 * t94;
t110 = t51 * pkin(2) + t112;
t105 = t74 * t127;
t42 = t48 * pkin(2);
t103 = t49 * pkin(9) - t42;
t50 = t127 * t94 + t128 * t75;
t44 = t50 * pkin(2);
t102 = t51 * pkin(9) - t44;
t99 = t71 * t105;
t98 = t66 * t107;
t96 = -t126 * pkin(1) + pkin(8) * t108;
t90 = t49 * pkin(2) - t96;
t21 = t76 * t108 + t49 * t74;
t47 = t111 * t74 + t76 * t120;
t46 = -t111 * t76 + t74 * t120;
t26 = t74 * t106 + t51 * t76;
t25 = -t76 * t106 + t51 * t74;
t15 = t67 * t107 + t47 * t66;
t6 = t26 * t67 + t50 * t66;
t5 = t26 * t66 - t50 * t67;
t2 = [(-t128 * mrSges(2,1) + t126 * mrSges(2,2) - m(3) * t112 - t51 * mrSges(3,1) - mrSges(3,3) * t106 - t100 * t6 - t147 * t26 - t131 * t50 - t142 * t5 - t148 * t25 - t133 * (t50 * t129 - t25 * t73 + t26 * t65 + t110) - t134 * t110) * g(2) + (t126 * mrSges(2,1) + t128 * mrSges(2,2) - m(3) * t96 + t49 * mrSges(3,1) - mrSges(3,3) * t108 + t147 * t22 + t131 * t48 + t100 * t136 + t142 * t1 + t148 * t21 + t134 * t90 + t133 * (t48 * t129 - t21 * t73 + t22 * t65 + t90)) * g(1) (-t133 * (t120 * t129 - t73 * t99 + t113) - t142 * (-t67 * t120 + t76 * t98) - t134 * t113 + t150 * t99 + (t149 * t105 + (-t100 * t66 + t143) * t75 + (-mrSges(3,1) + (-t133 * t65 + t145) * t76) * t127) * t71) * g(3) + (-m(4) * t103 + m(5) * t42 - t133 * (t49 * t129 + t138 * t48 + t103) - t100 * (-t48 * t121 + t49 * t66) - t142 * (-t48 * t122 - t49 * t67) + t137 * t49 + t135 * t48) * g(2) + (-m(4) * t102 + m(5) * t44 - t133 * (t51 * t129 + t138 * t50 + t102) - t142 * (-t50 * t122 - t51 * t67) - t100 * (-t50 * t121 + t51 * t66) + t137 * t51 + t135 * t50) * g(1) (-t133 * (-t46 * t65 - t47 * t73) - t148 * t47 + t130 * t46) * g(3) + (-t133 * (-t21 * t65 - t22 * t73) - t148 * t22 + t130 * t21) * g(2) + (-t133 * (-t25 * t65 - t26 * t73) - t148 * t26 + t130 * t25) * g(1), t132 * (-g(1) * t25 - g(2) * t21 - g(3) * t46) (-t142 * (t47 * t67 - t98) + t100 * t15) * g(3) + (t100 * t1 - t136 * t142) * g(2) + (t100 * t5 - t142 * t6) * g(1) (-g(1) * t5 - g(2) * t1 - g(3) * t15) * m(7)];
taug  = t2(:);
