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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 17:48
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

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
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP10_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP10_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:47:44
% EndTime: 2018-11-23 17:47:45
% DurationCPUTime: 1.21s
% Computational Cost: add. (1761->125), mult. (2031->159), div. (0->0), fcn. (2050->16), ass. (0->71)
t107 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t155 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t76 = sin(pkin(11));
t77 = cos(pkin(11));
t157 = m(5) * pkin(3) + t77 * mrSges(5,1) - t76 * mrSges(5,2) + mrSges(4,1);
t75 = pkin(11) + qJ(5);
t72 = sin(t75);
t73 = cos(t75);
t141 = t107 * t73 + t155 * t72 + t157;
t97 = -m(5) * qJ(4) + mrSges(4,2) - mrSges(7,2) - mrSges(5,3) - mrSges(6,3);
t146 = m(6) + m(7);
t71 = pkin(4) * t77 + pkin(3);
t78 = -pkin(10) - qJ(4);
t79 = sin(qJ(3));
t82 = cos(qJ(3));
t156 = t146 * (-t71 * t82 + t78 * t79) + t97 * t79 - mrSges(3,1) - t141 * t82;
t154 = -m(4) - t146;
t153 = t76 * mrSges(5,1) + t77 * mrSges(5,2) - mrSges(3,2) + mrSges(4,3);
t122 = sin(pkin(6));
t138 = cos(qJ(1));
t104 = t138 * t122;
t121 = pkin(6) - qJ(2);
t109 = sin(t121);
t120 = pkin(6) + qJ(2);
t108 = sin(t120);
t67 = t108 / 0.2e1;
t125 = t67 - t109 / 0.2e1;
t81 = sin(qJ(1));
t83 = cos(qJ(2));
t129 = t81 * t83;
t49 = t138 * t125 + t129;
t24 = -t79 * t104 + t49 * t82;
t80 = sin(qJ(2));
t103 = cos(t120) / 0.2e1;
t110 = cos(t121);
t87 = t110 / 0.2e1 + t103;
t48 = -t138 * t87 + t80 * t81;
t1 = t24 * t72 - t48 * t73;
t151 = t24 * t73 + t48 * t72;
t139 = pkin(4) * t76;
t149 = m(5) * pkin(9) + t107 * t72 + t139 * t146 - t155 * t73 + t153;
t147 = m(4) + m(5);
t144 = m(5) + t146;
t142 = (m(4) + t144) * pkin(9) + t153;
t111 = t81 * t122;
t124 = t138 * pkin(1) + pkin(8) * t111;
t123 = cos(pkin(6));
t117 = t138 * t83;
t52 = -t81 * t125 + t117;
t119 = t52 * pkin(2) + t124;
t116 = -pkin(1) * t81 + pkin(8) * t104;
t106 = t49 * pkin(2) - t116;
t102 = t109 / 0.2e1;
t23 = t82 * t104 + t49 * t79;
t86 = t102 - t108 / 0.2e1;
t62 = t103 - t110 / 0.2e1;
t61 = t67 + t102;
t60 = t61 * pkin(2);
t53 = t81 * t86 + t117;
t51 = t138 * t80 + t81 * t87;
t50 = -t138 * t86 + t129;
t47 = t123 * t79 - t62 * t82;
t46 = -t123 * t82 - t62 * t79;
t44 = t51 * pkin(2);
t42 = t48 * pkin(2);
t28 = t79 * t111 + t52 * t82;
t27 = -t82 * t111 + t52 * t79;
t11 = t47 * t72 + t61 * t73;
t6 = t28 * t73 + t51 * t72;
t5 = t28 * t72 - t51 * t73;
t2 = [(-t138 * mrSges(2,1) - m(3) * t124 - t52 * mrSges(3,1) + (-t122 * mrSges(3,3) + mrSges(2,2)) * t81 - t107 * t6 - t157 * t28 - t142 * t51 - t155 * t5 + t97 * t27 - t146 * (t51 * t139 - t27 * t78 + t28 * t71 + t119) - t147 * t119) * g(2) + (t81 * mrSges(2,1) + t138 * mrSges(2,2) - m(3) * t116 + t49 * mrSges(3,1) - mrSges(3,3) * t104 + t157 * t24 + t142 * t48 + t107 * t151 + t155 * t1 - t97 * t23 + t146 * (t48 * t139 - t23 * t78 + t24 * t71 + t106) + t147 * t106) * g(1) (-m(5) * t60 + t154 * (-pkin(9) * t62 + t60) + t149 * t62 + t156 * t61) * g(3) + (m(5) * t42 + t154 * (pkin(9) * t50 - t42) - t149 * t50 - t156 * t48) * g(2) + (m(5) * t44 + t154 * (pkin(9) * t53 - t44) - t149 * t53 - t156 * t51) * g(1) (-t146 * (-t46 * t71 - t47 * t78) + t97 * t47 + t141 * t46) * g(3) + (-t146 * (-t23 * t71 - t24 * t78) + t97 * t24 + t141 * t23) * g(2) + (-t146 * (-t27 * t71 - t28 * t78) + t97 * t28 + t141 * t27) * g(1), t144 * (-g(1) * t27 - g(2) * t23 - g(3) * t46) (-t155 * (t47 * t73 - t61 * t72) + t107 * t11) * g(3) + (t107 * t1 - t151 * t155) * g(2) + (t107 * t5 - t155 * t6) * g(1) (-g(1) * t5 - g(2) * t1 - g(3) * t11) * m(7)];
taug  = t2(:);
