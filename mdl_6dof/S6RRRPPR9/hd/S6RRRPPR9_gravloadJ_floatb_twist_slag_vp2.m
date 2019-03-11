% Calculate Gravitation load on the joints for
% S6RRRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 16:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPPR9_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR9_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR9_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:11:50
% EndTime: 2019-03-09 16:11:53
% DurationCPUTime: 1.64s
% Computational Cost: add. (813->161), mult. (2049->219), div. (0->0), fcn. (2494->12), ass. (0->88)
t143 = mrSges(5,2) - mrSges(6,3);
t73 = sin(qJ(6));
t76 = cos(qJ(6));
t156 = t73 * mrSges(7,1) + t76 * mrSges(7,2) - t143;
t136 = -m(7) * pkin(5) - mrSges(5,1) - mrSges(6,1);
t155 = t76 * mrSges(7,1) - t73 * mrSges(7,2) - t136;
t70 = sin(pkin(11));
t72 = cos(pkin(11));
t154 = -pkin(4) * t72 - qJ(5) * t70;
t149 = mrSges(7,3) - mrSges(6,2) - mrSges(5,3);
t153 = m(7) * pkin(10) + t149;
t129 = cos(qJ(3));
t74 = sin(qJ(3));
t148 = -t129 * pkin(3) - qJ(4) * t74;
t130 = cos(qJ(1));
t71 = sin(pkin(6));
t110 = t71 * t130;
t128 = sin(qJ(1));
t75 = sin(qJ(2));
t77 = cos(qJ(2));
t114 = cos(pkin(6));
t96 = t114 * t130;
t54 = t128 * t77 + t75 * t96;
t28 = -t74 * t110 + t129 * t54;
t53 = t128 * t75 - t77 * t96;
t5 = t28 * t70 - t53 * t72;
t6 = t28 * t72 + t53 * t70;
t146 = m(6) + m(7);
t145 = -mrSges(4,1) * t129 + mrSges(4,2) * t74 - mrSges(3,1);
t144 = mrSges(3,2) - mrSges(4,3);
t109 = t71 * t129;
t27 = t109 * t130 + t54 * t74;
t141 = t154 * t27;
t108 = t71 * t128;
t95 = t114 * t128;
t56 = t130 * t77 - t75 * t95;
t31 = -t108 * t129 + t56 * t74;
t140 = t154 * t31;
t125 = t71 * t75;
t51 = -t114 * t129 + t125 * t74;
t139 = t154 * t51;
t138 = t155 * t72 + t156 * t70 + mrSges(4,1);
t137 = mrSges(4,2) - m(7) * (-pkin(10) + qJ(4)) + t149;
t134 = -t153 * t74 - t145;
t124 = t71 * t77;
t120 = pkin(2) * t124 + pkin(9) * t125;
t119 = t130 * pkin(1) + pkin(8) * t108;
t118 = qJ(4) * t31;
t115 = t27 * qJ(4);
t113 = t74 * t124;
t111 = t70 * t129;
t107 = t72 * t129;
t106 = t77 * t129;
t105 = -t53 * pkin(2) + pkin(9) * t54;
t55 = t130 * t75 + t77 * t95;
t104 = -t55 * pkin(2) + pkin(9) * t56;
t21 = t27 * pkin(3);
t103 = qJ(4) * t28 - t21;
t23 = t31 * pkin(3);
t32 = t108 * t74 + t129 * t56;
t102 = qJ(4) * t32 - t23;
t46 = t51 * pkin(3);
t52 = t109 * t75 + t114 * t74;
t101 = qJ(4) * t52 - t46;
t99 = t71 * t106;
t100 = pkin(3) * t99 + qJ(4) * t113 + t120;
t98 = -pkin(1) * t128 + pkin(8) * t110;
t92 = t148 * t53 + t105;
t91 = t148 * t55 + t104;
t90 = t56 * pkin(2) + pkin(9) * t55 + t119;
t89 = t32 * pkin(3) + t90;
t84 = -t54 * pkin(2) - t53 * pkin(9) + t98;
t81 = -pkin(3) * t28 + t84;
t10 = t32 * t72 + t55 * t70;
t9 = t32 * t70 - t55 * t72;
t80 = t10 * pkin(4) + qJ(5) * t9 + t89;
t79 = -pkin(4) * t6 - qJ(5) * t5 + t81;
t35 = (t106 * t72 + t70 * t75) * t71;
t34 = -t125 * t72 + t70 * t99;
t26 = -t124 * t70 + t52 * t72;
t25 = t124 * t72 + t52 * t70;
t16 = -t107 * t55 + t56 * t70;
t15 = -t111 * t55 - t56 * t72;
t14 = -t107 * t53 + t54 * t70;
t13 = -t111 * t53 - t54 * t72;
t2 = t10 * t76 + t73 * t9;
t1 = -t10 * t73 + t76 * t9;
t3 = [(-t130 * mrSges(2,1) + t128 * mrSges(2,2) - m(3) * t119 - t56 * mrSges(3,1) - mrSges(3,3) * t108 - m(4) * t90 - t32 * mrSges(4,1) - m(5) * (t89 + t118) - m(6) * (t80 + t118) - m(7) * t80 - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t143 * t9 + t144 * t55 + t136 * t10 + t137 * t31) * g(2) + (t128 * mrSges(2,1) + t130 * mrSges(2,2) - m(3) * t98 + t54 * mrSges(3,1) - mrSges(3,3) * t110 - m(4) * t84 + t28 * mrSges(4,1) - m(5) * (t81 - t115) - m(6) * (t79 - t115) - m(7) * t79 + t155 * t6 + t156 * t5 - t144 * t53 - t137 * t27) * g(1) (-m(4) * t120 - m(5) * t100 - t146 * (t35 * pkin(4) + t34 * qJ(5) + t100) + (t144 * t75 + t145 * t77) * t71 - t155 * t35 - t156 * t34 + t153 * t113) * g(3) + (-m(4) * t105 - m(5) * t92 - t146 * (t14 * pkin(4) + qJ(5) * t13 + t92) + t144 * t54 - t155 * t14 - t156 * t13 + t134 * t53) * g(2) + (-m(4) * t104 - m(5) * t91 - t146 * (t16 * pkin(4) + qJ(5) * t15 + t91) + t144 * t56 - t155 * t16 - t156 * t15 + t134 * t55) * g(1) (-m(5) * t101 - m(6) * (t101 + t139) - m(7) * (-t46 + t139) + t137 * t52 + t138 * t51) * g(3) + (-m(5) * t103 - m(6) * (t103 + t141) - m(7) * (-t21 + t141) + t137 * t28 + t138 * t27) * g(2) + (-m(5) * t102 - m(6) * (t102 + t140) - m(7) * (-t23 + t140) + t137 * t32 + t138 * t31) * g(1) (m(5) + t146) * (-g(1) * t31 - g(2) * t27 - g(3) * t51) t146 * (-g(1) * t9 - g(2) * t5 - g(3) * t25) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t5 * t76 - t6 * t73) * mrSges(7,1) + (-t5 * t73 - t6 * t76) * mrSges(7,2)) - g(3) * ((t25 * t76 - t26 * t73) * mrSges(7,1) + (-t25 * t73 - t26 * t76) * mrSges(7,2))];
taug  = t3(:);
