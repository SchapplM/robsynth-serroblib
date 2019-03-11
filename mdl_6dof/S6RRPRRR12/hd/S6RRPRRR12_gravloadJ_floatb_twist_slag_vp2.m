% Calculate Gravitation load on the joints for
% S6RRPRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR12_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR12_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR12_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR12_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:35:16
% EndTime: 2019-03-09 14:35:19
% DurationCPUTime: 1.31s
% Computational Cost: add. (665->140), mult. (1231->182), div. (0->0), fcn. (1403->12), ass. (0->75)
t146 = mrSges(6,2) - mrSges(7,3);
t94 = -m(7) * pkin(11) + t146;
t138 = m(6) + m(7);
t139 = m(4) + m(5);
t147 = -t139 - t138;
t150 = t147 * qJ(3);
t67 = sin(qJ(6));
t71 = cos(qJ(6));
t145 = mrSges(7,1) * t71 - mrSges(7,2) * t67 + mrSges(6,1);
t116 = mrSges(3,2) - mrSges(4,3);
t148 = -m(7) * pkin(5) - t145;
t65 = qJ(4) + qJ(5);
t62 = sin(t65);
t63 = cos(t65);
t68 = sin(qJ(4));
t72 = cos(qJ(4));
t92 = t68 * mrSges(5,1) + t72 * mrSges(5,2);
t149 = t148 * t62 - t94 * t63 + t116 - t92;
t84 = -m(5) * pkin(9) - mrSges(3,1) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t143 = mrSges(7,1) * t67 + mrSges(7,2) * t71 - t84;
t66 = sin(pkin(6));
t74 = cos(qJ(1));
t117 = t66 * t74;
t69 = sin(qJ(2));
t70 = sin(qJ(1));
t73 = cos(qJ(2));
t111 = cos(pkin(6));
t96 = t74 * t111;
t45 = t69 * t70 - t73 * t96;
t142 = t68 * t117 + t45 * t72;
t119 = t66 * t70;
t97 = t70 * t111;
t47 = t74 * t69 + t73 * t97;
t19 = -t68 * t119 + t47 * t72;
t118 = t66 * t73;
t29 = -t111 * t62 - t118 * t63;
t30 = t111 * t63 - t118 * t62;
t136 = -t145 * t29 + t146 * t30;
t15 = t117 * t62 + t45 * t63;
t17 = t117 * t63 - t45 * t62;
t135 = -t145 * t15 - t146 * t17;
t13 = t119 * t62 - t47 * t63;
t14 = t119 * t63 + t47 * t62;
t134 = t145 * t13 + t146 * t14;
t132 = t116 + t150;
t130 = t149 + t150;
t127 = pkin(4) * t68;
t124 = t45 * t68;
t122 = t47 * t68;
t120 = t66 * t69;
t114 = t142 * pkin(4);
t113 = pkin(2) * t118 + qJ(3) * t120;
t112 = t74 * pkin(1) + pkin(8) * t119;
t48 = -t69 * t97 + t73 * t74;
t106 = t48 * pkin(2) + t112;
t101 = -pkin(1) * t70 + pkin(8) * t117;
t100 = -t13 * pkin(5) + pkin(11) * t14;
t99 = t15 * pkin(5) - pkin(11) * t17;
t98 = t29 * pkin(5) + t30 * pkin(11);
t95 = -m(5) * pkin(3) - mrSges(4,1) - mrSges(3,3);
t46 = t69 * t96 + t70 * t73;
t93 = t46 * pkin(2) - t101;
t88 = t19 * pkin(4);
t61 = pkin(4) * t72 + pkin(3);
t75 = -pkin(10) - pkin(9);
t85 = pkin(4) * t122 + t61 * t119 - t48 * t75 + t106;
t79 = -t111 * t68 - t118 * t72;
t77 = t79 * pkin(4);
t55 = t72 * t117;
t43 = t47 * pkin(2);
t41 = t45 * pkin(2);
t20 = t119 * t72 + t122;
t2 = t14 * t71 + t48 * t67;
t1 = -t14 * t67 + t48 * t71;
t3 = [(-t74 * mrSges(2,1) - m(3) * t112 - t20 * mrSges(5,1) - t19 * mrSges(5,2) - m(6) * t85 - t14 * mrSges(6,1) - m(7) * (pkin(5) * t14 + t85) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t132 * t47 + t94 * t13 + (t66 * t95 + mrSges(2,2)) * t70 + t84 * t48 - t139 * t106) * g(2) + (t70 * mrSges(2,1) - m(3) * t101 - t55 * mrSges(5,1) + t94 * t15 + (-t132 + t92) * t45 + (mrSges(2,2) + (mrSges(5,2) * t68 + t95) * t66) * t74 + t148 * t17 + t143 * t46 + t139 * t93 + t138 * (pkin(4) * t124 - t61 * t117 - t46 * t75 + t93)) * g(1) (-t139 * t113 - t138 * (t120 * t127 + t113) + ((t138 * t75 - t143) * t73 + t149 * t69) * t66) * g(3) + (t139 * t41 - t138 * (t46 * t127 + t45 * t75 - t41) + t130 * t46 + t143 * t45) * g(2) + (t139 * t43 - t138 * (t48 * t127 + t47 * t75 - t43) + t130 * t48 + t143 * t47) * g(1) -(-g(1) * t47 - g(2) * t45 + g(3) * t118) * t147 (-t79 * mrSges(5,1) - (-t111 * t72 + t118 * t68) * mrSges(5,2) - m(6) * t77 - m(7) * (t77 + t98) + t136) * g(3) + (-t142 * mrSges(5,1) - (t55 - t124) * mrSges(5,2) - m(6) * t114 - m(7) * (t99 + t114) + t135) * g(2) + (-t19 * mrSges(5,1) + t20 * mrSges(5,2) - m(6) * t88 - m(7) * (t100 + t88) + t134) * g(1) (-m(7) * t98 + t136) * g(3) + (-m(7) * t99 + t135) * g(2) + (-m(7) * t100 + t134) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t17 * t67 + t46 * t71) * mrSges(7,1) + (t17 * t71 - t46 * t67) * mrSges(7,2)) - g(3) * ((t120 * t71 - t30 * t67) * mrSges(7,1) + (-t120 * t67 - t30 * t71) * mrSges(7,2))];
taug  = t3(:);
