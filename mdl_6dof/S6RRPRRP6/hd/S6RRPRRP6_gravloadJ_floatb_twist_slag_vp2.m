% Calculate Gravitation load on the joints for
% S6RRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 17:14
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP6_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:14:19
% EndTime: 2018-11-23 17:14:20
% DurationCPUTime: 1.12s
% Computational Cost: add. (1995->159), mult. (1702->212), div. (0->0), fcn. (1657->20), ass. (0->93)
t119 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t117 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t85 = sin(pkin(6));
t94 = cos(qJ(1));
t139 = t85 * t94;
t82 = qJ(2) + pkin(11);
t121 = pkin(6) - t82;
t110 = sin(t121);
t120 = pkin(6) + t82;
t109 = sin(t120);
t67 = t109 / 0.2e1;
t129 = t67 - t110 / 0.2e1;
t79 = cos(t82);
t90 = sin(qJ(1));
t136 = t90 * t79;
t39 = t129 * t94 + t136;
t88 = sin(qJ(4));
t92 = cos(qJ(4));
t18 = -t139 * t88 + t39 * t92;
t76 = sin(t82);
t106 = cos(t121) / 0.2e1;
t111 = cos(t120);
t97 = t111 / 0.2e1 + t106;
t38 = t76 * t90 - t94 * t97;
t87 = sin(qJ(5));
t91 = cos(qJ(5));
t1 = t18 * t87 - t38 * t91;
t156 = t18 * t91 + t38 * t87;
t148 = -m(6) - m(7);
t132 = mrSges(5,3) - mrSges(4,2);
t155 = mrSges(6,3) + mrSges(7,2);
t154 = t92 * mrSges(5,1) - mrSges(5,2) * t88 + mrSges(4,1);
t83 = pkin(6) + qJ(2);
t74 = cos(t83) / 0.2e1;
t84 = pkin(6) - qJ(2);
t81 = cos(t84);
t62 = t81 / 0.2e1 + t74;
t153 = mrSges(5,2) - t155;
t152 = t117 * t87 - t119 * t91 - mrSges(5,1);
t127 = m(5) - t148;
t151 = pkin(9) * t127 + t132;
t150 = sin(t83) / 0.2e1;
t78 = sin(t84);
t147 = pkin(2) * t78;
t146 = pkin(4) * t92;
t144 = t38 * t88;
t41 = t76 * t94 + t90 * t97;
t142 = t41 * t88;
t105 = t110 / 0.2e1;
t58 = t105 + t67;
t141 = t58 * t88;
t140 = t85 * t90;
t138 = t87 * t92;
t89 = sin(qJ(2));
t137 = t89 * t94;
t135 = t90 * t89;
t134 = t91 * t92;
t133 = t94 * t79;
t42 = -t129 * t90 + t133;
t93 = cos(qJ(2));
t75 = pkin(2) * t93 + pkin(1);
t63 = t94 * t75;
t130 = pkin(3) * t42 + t63;
t70 = pkin(2) * t150;
t128 = t147 / 0.2e1 + t70;
t125 = m(4) + t127;
t17 = -t139 * t92 - t39 * t88;
t60 = t62 * pkin(2);
t118 = -pkin(2) * t135 + t60 * t94;
t59 = t106 - t111 / 0.2e1;
t112 = pkin(3) * t58 + pkin(9) * t59 + t128;
t108 = m(3) * pkin(1) + mrSges(3,1) * t93 + mrSges(2,1);
t107 = -pkin(2) * t137 - t60 * t90;
t104 = pkin(10) * t148 + t153;
t96 = -t109 / 0.2e1 + t105;
t40 = -t94 * t96 + t136;
t102 = -pkin(3) * t38 + pkin(9) * t40 + t118;
t43 = t90 * t96 + t133;
t101 = -pkin(3) * t41 + pkin(9) * t43 + t107;
t61 = t150 - t78 / 0.2e1;
t95 = mrSges(3,1) * t61 + mrSges(2,2) + (-m(3) * pkin(8) - mrSges(3,3) - mrSges(4,3)) * t85 + t125 * (-t147 / 0.2e1 + t70 - t85 * (pkin(8) + qJ(3)));
t86 = cos(pkin(6));
t48 = -t62 * t90 - t137;
t47 = -t62 * t94 + t135;
t46 = t59 * t92 + t86 * t88;
t45 = -t59 * t88 + t86 * t92;
t35 = t39 * pkin(3);
t22 = t140 * t88 + t42 * t92;
t21 = -t140 * t92 + t42 * t88;
t11 = t46 * t87 + t58 * t91;
t6 = t22 * t91 + t41 * t87;
t5 = t22 * t87 - t41 * t91;
t2 = [(-t48 * mrSges(3,2) - m(4) * t63 - t42 * mrSges(4,1) - m(5) * t130 - t22 * mrSges(5,1) - t119 * t6 + t117 * t5 - t151 * t41 - t108 * t94 + t95 * t90 + t104 * t21 + t148 * (pkin(4) * t22 + t130)) * g(2) + (-t47 * mrSges(3,2) + t39 * mrSges(4,1) + m(5) * t35 + t18 * mrSges(5,1) + t119 * t156 + t151 * t38 - t117 * t1 + (t125 * t75 + t108) * t90 + t95 * t94 + t104 * t17 + t148 * (-pkin(4) * t18 - t35)) * g(1) (-(t150 + t78 / 0.2e1) * mrSges(3,1) - (t74 - t81 / 0.2e1) * mrSges(3,2) - m(4) * t128 - m(5) * t112 - t132 * t59 - t154 * t58 + t148 * (pkin(10) * t141 + t146 * t58 + t112) - t119 * (t134 * t58 + t59 * t87) + t117 * (t138 * t58 - t59 * t91) - t155 * t141) * g(3) + (t47 * mrSges(3,1) - (-t61 * t94 - t90 * t93) * mrSges(3,2) - m(4) * t118 - m(5) * t102 + t148 * (-pkin(10) * t144 - t146 * t38 + t102) - t119 * (-t134 * t38 + t40 * t87) + t117 * (-t138 * t38 - t40 * t91) - t132 * t40 + t154 * t38 + t155 * t144) * g(2) + (-t48 * mrSges(3,1) - (t61 * t90 - t93 * t94) * mrSges(3,2) - m(4) * t107 - m(5) * t101 + t148 * (-pkin(10) * t142 - t146 * t41 + t101) + t117 * (-t138 * t41 - t43 * t91) - t132 * t43 + t154 * t41 + t155 * t142 - t119 * (-t134 * t41 + t43 * t87)) * g(1) (-g(3) * t86 + (-g(1) * t90 + g(2) * t94) * t85) * t125 (t148 * (pkin(4) * t45 + pkin(10) * t46) + t153 * t46 + t152 * t45) * g(3) + (t148 * (pkin(4) * t17 + pkin(10) * t18) + t153 * t18 + t152 * t17) * g(2) + (t148 * (-pkin(4) * t21 + pkin(10) * t22) + t153 * t22 - t152 * t21) * g(1) (t117 * (t46 * t91 - t58 * t87) + t119 * t11) * g(3) + (t119 * t1 + t117 * t156) * g(2) + (t117 * t6 + t119 * t5) * g(1) (-g(1) * t5 - g(2) * t1 - g(3) * t11) * m(7)];
taug  = t2(:);
