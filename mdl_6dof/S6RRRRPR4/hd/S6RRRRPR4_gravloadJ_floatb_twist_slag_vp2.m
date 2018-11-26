% Calculate Gravitation load on the joints for
% S6RRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2018-11-23 18:14
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:13:58
% EndTime: 2018-11-23 18:13:59
% DurationCPUTime: 1.00s
% Computational Cost: add. (620->140), mult. (619->161), div. (0->0), fcn. (553->12), ass. (0->80)
t59 = qJ(4) + pkin(11);
t53 = qJ(6) + t59;
t45 = sin(t53);
t51 = sin(t59);
t62 = sin(qJ(4));
t158 = -t62 * mrSges(5,2) - t51 * mrSges(6,2) - t45 * mrSges(7,2);
t157 = -mrSges(5,3) - mrSges(6,3) - mrSges(7,3);
t46 = cos(t53);
t52 = cos(t59);
t65 = cos(qJ(4));
t156 = t65 * mrSges(5,1) + t52 * mrSges(6,1) + t46 * mrSges(7,1);
t60 = qJ(2) + qJ(3);
t54 = sin(t60);
t155 = t158 * t54;
t154 = -m(5) * pkin(9) + t157;
t64 = sin(qJ(1));
t67 = cos(qJ(1));
t142 = g(1) * t67 + g(2) * t64;
t152 = t156 * t54;
t55 = cos(t60);
t151 = -t55 * mrSges(4,1) + (mrSges(4,2) + t157) * t54;
t56 = t65 * pkin(4);
t49 = t56 + pkin(3);
t61 = -qJ(5) - pkin(9);
t144 = t55 * t49 - t54 * t61;
t27 = pkin(5) * t52 + t56;
t23 = pkin(3) + t27;
t58 = -pkin(10) + t61;
t145 = t55 * t23 - t54 * t58;
t148 = t55 * pkin(3) + t54 * pkin(9);
t150 = -m(5) * t148 - m(6) * t144 - m(7) * t145;
t132 = m(6) * pkin(4);
t149 = m(6) + m(7);
t146 = mrSges(5,1) + t132;
t141 = -m(3) * pkin(7) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t128 = pkin(3) * t54;
t63 = sin(qJ(2));
t129 = pkin(2) * t63;
t79 = -t49 * t54 - t55 * t61;
t81 = -t23 * t54 - t55 * t58;
t140 = -m(7) * (t81 - t129) - m(6) * (t79 - t129) - m(5) * (-t128 - t129) + t152;
t139 = m(4) + m(5) + t149;
t138 = t151 + (-t156 - t158) * t55;
t127 = pkin(4) * t62;
t26 = pkin(5) * t51 + t127;
t137 = -m(6) * t127 - m(7) * t26;
t136 = m(5) * t128 - m(6) * t79 - m(7) * t81 + t152;
t66 = cos(qJ(2));
t87 = t66 * mrSges(3,1) - t63 * mrSges(3,2);
t135 = m(3) * pkin(1) + mrSges(2,1) - t151 + t87;
t115 = t55 * t67;
t134 = t154 * t115 + t155 * t67;
t116 = t55 * t64;
t133 = t154 * t116 + t155 * t64;
t5 = t45 * t116 + t46 * t67;
t6 = -t46 * t116 + t45 * t67;
t131 = -t5 * mrSges(7,1) + t6 * mrSges(7,2);
t7 = -t45 * t115 + t46 * t64;
t8 = t46 * t115 + t45 * t64;
t130 = t7 * mrSges(7,1) - t8 * mrSges(7,2);
t124 = g(3) * t54;
t57 = t66 * pkin(2);
t113 = t62 * t64;
t112 = t62 * t67;
t111 = t64 * t26;
t110 = t64 * t65;
t108 = t65 * t67;
t84 = mrSges(4,1) * t54 + mrSges(4,2) * t55;
t83 = -mrSges(7,1) * t45 - mrSges(7,2) * t46;
t16 = -t55 * t112 + t110;
t14 = t55 * t113 + t108;
t68 = -pkin(8) - pkin(7);
t50 = t57 + pkin(1);
t17 = t55 * t108 + t113;
t15 = -t55 * t110 + t112;
t12 = t52 * t115 + t51 * t64;
t11 = -t51 * t115 + t52 * t64;
t10 = -t52 * t116 + t51 * t67;
t9 = t51 * t116 + t52 * t67;
t1 = [(-t113 * t132 - m(7) * t111 - t17 * mrSges(5,1) - t12 * mrSges(6,1) - t8 * mrSges(7,1) - t16 * mrSges(5,2) - t11 * mrSges(6,2) - t7 * mrSges(7,2) - t139 * (t67 * t50 - t64 * t68) + t141 * t64 + (-t135 + t150) * t67) * g(2) + (-t15 * mrSges(5,1) - t10 * mrSges(6,1) - t6 * mrSges(7,1) - t14 * mrSges(5,2) - t9 * mrSges(6,2) - t5 * mrSges(7,2) + (t139 * t68 + t137 + t141) * t67 + (m(4) * t50 - m(5) * (-t50 - t148) - m(6) * (-t50 - t144) - m(7) * (-t50 - t145) + t135) * t64) * g(1) (t140 * t64 + t133) * g(2) + (t140 * t67 + t134) * g(1) + (-t87 - m(4) * t57 - m(5) * (t57 + t148) - m(6) * (t57 + t144) - m(7) * (t57 + t145) + t138) * g(3) + t142 * (m(4) * t129 + mrSges(3,1) * t63 + mrSges(3,2) * t66 + t84) t142 * t84 + (t136 * t64 + t133) * g(2) + (t136 * t67 + t134) * g(1) + (t138 + t150) * g(3) (mrSges(5,1) * t62 + mrSges(6,1) * t51 + mrSges(5,2) * t65 + mrSges(6,2) * t52 - t137 - t83) * t124 + (-t15 * mrSges(5,2) + t9 * mrSges(6,1) - t10 * mrSges(6,2) - m(7) * (-t55 * t111 - t27 * t67) - t131 + t146 * t14) * g(2) + (t17 * mrSges(5,2) - t11 * mrSges(6,1) + t12 * mrSges(6,2) - m(7) * (-t26 * t115 + t27 * t64) - t130 - t146 * t16) * g(1) (g(3) * t55 - t142 * t54) * t149, -g(1) * t130 - g(2) * t131 - t83 * t124];
taug  = t1(:);
