% Calculate Gravitation load on the joints for
% S6RRRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-10 00:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR13_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR13_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR13_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR13_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:51:08
% EndTime: 2019-03-09 23:51:13
% DurationCPUTime: 1.70s
% Computational Cost: add. (958->169), mult. (2410->230), div. (0->0), fcn. (2970->12), ass. (0->91)
t163 = -m(7) * pkin(5) - mrSges(5,1) - mrSges(6,1);
t154 = mrSges(5,2) - mrSges(6,3);
t75 = sin(qJ(6));
t79 = cos(qJ(6));
t144 = -t75 * mrSges(7,1) - t79 * mrSges(7,2) + t154;
t143 = -t79 * mrSges(7,1) + t75 * mrSges(7,2) + t163;
t76 = sin(qJ(4));
t80 = cos(qJ(4));
t162 = pkin(4) * t80 + qJ(5) * t76;
t159 = mrSges(7,3) - mrSges(6,2) - mrSges(5,3);
t125 = sin(pkin(6));
t137 = cos(qJ(1));
t105 = t137 * t125;
t136 = cos(qJ(3));
t126 = cos(pkin(6));
t107 = t126 * t137;
t135 = sin(qJ(1));
t78 = sin(qJ(2));
t81 = cos(qJ(2));
t58 = t107 * t78 + t135 * t81;
t77 = sin(qJ(3));
t31 = -t77 * t105 + t136 * t58;
t57 = -t107 * t81 + t135 * t78;
t7 = t31 * t76 - t57 * t80;
t8 = t31 * t80 + t57 * t76;
t158 = m(6) + m(7);
t130 = t77 * mrSges(4,2);
t157 = mrSges(4,1) * t136 + mrSges(3,1) - t130;
t155 = mrSges(3,2) - mrSges(4,3);
t104 = t136 * t125;
t114 = -t137 * t104 - t58 * t77;
t152 = t162 * t114;
t103 = t125 * t135;
t106 = t126 * t135;
t60 = -t106 * t78 + t137 * t81;
t34 = -t103 * t136 + t60 * t77;
t151 = t162 * t34;
t113 = t78 * t125;
t55 = -t113 * t77 + t126 * t136;
t150 = t162 * t55;
t149 = t158 * pkin(4) - t143;
t148 = mrSges(4,2) - m(7) * (pkin(10) - pkin(11)) + t159;
t147 = t143 * t80 + t144 * t76 - mrSges(4,1);
t146 = m(7) * pkin(11) + t159;
t139 = pkin(10) * t34;
t138 = t114 * pkin(10);
t133 = t57 * t77;
t59 = t106 * t81 + t137 * t78;
t131 = t59 * t77;
t112 = t81 * t125;
t129 = pkin(2) * t112 + pkin(9) * t113;
t128 = t137 * pkin(1) + pkin(8) * t103;
t122 = t136 * pkin(3);
t121 = t76 * t136;
t120 = t80 * t136;
t119 = -t57 * pkin(2) + pkin(9) * t58;
t118 = -t59 * pkin(2) + pkin(9) * t60;
t24 = t114 * pkin(3);
t117 = pkin(10) * t31 + t24;
t26 = t34 * pkin(3);
t35 = t103 * t77 + t136 * t60;
t116 = pkin(10) * t35 - t26;
t50 = t55 * pkin(3);
t56 = t104 * t78 + t126 * t77;
t115 = pkin(10) * t56 + t50;
t110 = t77 * t112;
t95 = t81 * t104;
t111 = pkin(3) * t95 + pkin(10) * t110 + t129;
t108 = -pkin(1) * t135 + pkin(8) * t105;
t98 = -pkin(10) * t133 - t57 * t122 + t119;
t97 = -pkin(10) * t131 - t59 * t122 + t118;
t96 = t60 * pkin(2) + pkin(9) * t59 + t128;
t94 = t35 * pkin(3) + t96;
t90 = -t58 * pkin(2) - t57 * pkin(9) + t108;
t87 = -pkin(3) * t31 + t90;
t11 = t35 * t76 - t59 * t80;
t12 = t35 * t80 + t59 * t76;
t85 = t12 * pkin(4) + qJ(5) * t11 + t94;
t84 = -t158 * qJ(5) + t144;
t83 = -pkin(4) * t8 - qJ(5) * t7 + t87;
t38 = t113 * t76 + t80 * t95;
t37 = -t113 * t80 + t76 * t95;
t29 = -t112 * t76 + t56 * t80;
t28 = t112 * t80 + t56 * t76;
t18 = -t120 * t59 + t60 * t76;
t17 = -t121 * t59 - t60 * t80;
t16 = -t120 * t57 + t58 * t76;
t15 = -t121 * t57 - t58 * t80;
t2 = t11 * t75 + t12 * t79;
t1 = t11 * t79 - t12 * t75;
t3 = [(-t137 * mrSges(2,1) + t135 * mrSges(2,2) - m(3) * t128 - t60 * mrSges(3,1) - mrSges(3,3) * t103 - m(4) * t96 - t35 * mrSges(4,1) - m(5) * (t94 + t139) - m(6) * (t85 + t139) - m(7) * t85 - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t155 * t59 + t163 * t12 + t154 * t11 + t148 * t34) * g(2) + (t135 * mrSges(2,1) + t137 * mrSges(2,2) - m(3) * t108 + t58 * mrSges(3,1) - mrSges(3,3) * t105 - m(4) * t90 + t31 * mrSges(4,1) - m(5) * (t87 + t138) - m(6) * (t83 + t138) - m(7) * t83 - t155 * t57 - t144 * t7 - t143 * t8 + t148 * t114) * g(1) (-mrSges(3,1) * t112 - m(4) * t129 - (mrSges(4,1) * t104 - t125 * t130) * t81 - m(5) * t111 - t158 * (t38 * pkin(4) + t37 * qJ(5) + t111) + t143 * t38 + t144 * t37 + t155 * t113 + t146 * t110) * g(3) + (-m(4) * t119 - m(5) * t98 - t158 * (t16 * pkin(4) + qJ(5) * t15 + t98) + t155 * t58 + t157 * t57 + t143 * t16 + t144 * t15 - t146 * t133) * g(2) + (-m(4) * t118 - m(5) * t97 - t158 * (t18 * pkin(4) + qJ(5) * t17 + t97) + t155 * t60 + t157 * t59 + t143 * t18 + t144 * t17 - t146 * t131) * g(1) (-m(5) * t115 - m(6) * (t115 + t150) - m(7) * (t50 + t150) + t148 * t56 + t147 * t55) * g(3) + (-m(5) * t117 - m(6) * (t117 + t152) - m(7) * (t24 + t152) + t148 * t31 + t147 * t114) * g(2) + (-m(5) * t116 - m(6) * (t116 - t151) - m(7) * (-t26 - t151) + t148 * t35 - t147 * t34) * g(1) (t149 * t28 + t29 * t84) * g(3) + (t149 * t7 + t8 * t84) * g(2) + (t149 * t11 + t12 * t84) * g(1), t158 * (-g(1) * t11 - g(2) * t7 - g(3) * t28) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t7 * t79 - t75 * t8) * mrSges(7,1) + (-t7 * t75 - t79 * t8) * mrSges(7,2)) - g(3) * ((t28 * t79 - t29 * t75) * mrSges(7,1) + (-t28 * t75 - t29 * t79) * mrSges(7,2))];
taug  = t3(:);
