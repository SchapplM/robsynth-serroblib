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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 17:39
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

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
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR9_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:39:19
% EndTime: 2018-11-23 17:39:20
% DurationCPUTime: 1.62s
% Computational Cost: add. (1908->166), mult. (2414->215), div. (0->0), fcn. (2494->16), ass. (0->95)
t152 = mrSges(5,2) - mrSges(6,3);
t77 = sin(qJ(6));
t80 = cos(qJ(6));
t162 = t77 * mrSges(7,1) + t80 * mrSges(7,2) - t152;
t146 = -m(7) * pkin(5) - mrSges(5,1) - mrSges(6,1);
t161 = t80 * mrSges(7,1) - t77 * mrSges(7,2) - t146;
t75 = sin(pkin(11));
t76 = cos(pkin(11));
t160 = -pkin(4) * t76 - qJ(5) * t75;
t159 = mrSges(4,2) + mrSges(7,3) - mrSges(6,2) - mrSges(5,3);
t139 = cos(qJ(3));
t78 = sin(qJ(3));
t156 = t139 * pkin(3) + qJ(4) * t78;
t126 = sin(pkin(6));
t140 = cos(qJ(1));
t107 = t140 * t126;
t138 = sin(qJ(1));
t81 = cos(qJ(2));
t119 = t138 * t81;
t124 = pkin(6) + qJ(2);
t110 = sin(t124);
t103 = t110 / 0.2e1;
t125 = pkin(6) - qJ(2);
t111 = sin(t125);
t96 = t103 - t111 / 0.2e1;
t54 = t140 * t96 + t119;
t29 = -t78 * t107 + t139 * t54;
t79 = sin(qJ(2));
t105 = cos(t124) / 0.2e1;
t112 = cos(t125);
t85 = t112 / 0.2e1 + t105;
t53 = t138 * t79 - t140 * t85;
t5 = t29 * t75 - t53 * t76;
t6 = t29 * t76 + t53 * t75;
t154 = m(6) + m(7);
t153 = mrSges(3,2) - mrSges(4,3);
t28 = t107 * t139 + t54 * t78;
t151 = t160 * t28;
t106 = t126 * t138;
t120 = t140 * t81;
t57 = -t138 * t96 + t120;
t32 = -t106 * t139 + t57 * t78;
t150 = t160 * t32;
t127 = cos(pkin(6));
t66 = t105 - t112 / 0.2e1;
t51 = -t127 * t139 - t66 * t78;
t149 = t160 * t51;
t148 = -m(7) * (-pkin(10) + qJ(4)) + t159;
t147 = t161 * t76 + t162 * t75 + mrSges(4,1);
t144 = -mrSges(4,1) * t139 - mrSges(3,1) + (m(7) * pkin(10) + t159) * t78;
t132 = t140 * pkin(1) + pkin(8) * t106;
t131 = qJ(4) * t32;
t128 = t28 * qJ(4);
t122 = t75 * t139;
t121 = t76 * t139;
t104 = t111 / 0.2e1;
t84 = t104 - t110 / 0.2e1;
t55 = -t140 * t84 + t119;
t118 = -t53 * pkin(2) + pkin(9) * t55;
t56 = t138 * t85 + t140 * t79;
t58 = t138 * t84 + t120;
t117 = -t56 * pkin(2) + pkin(9) * t58;
t65 = t103 + t104;
t116 = t65 * pkin(2) - pkin(9) * t66;
t23 = t28 * pkin(3);
t115 = qJ(4) * t29 - t23;
t25 = t32 * pkin(3);
t33 = t106 * t78 + t139 * t57;
t114 = qJ(4) * t33 - t25;
t46 = t51 * pkin(3);
t52 = t127 * t78 - t139 * t66;
t113 = qJ(4) * t52 - t46;
t109 = -pkin(1) * t138 + pkin(8) * t107;
t100 = -t156 * t53 + t118;
t99 = -t156 * t56 + t117;
t98 = t57 * pkin(2) + pkin(9) * t56 + t132;
t97 = t156 * t65 + t116;
t95 = t33 * pkin(3) + t98;
t91 = -t54 * pkin(2) - t53 * pkin(9) + t109;
t87 = -pkin(3) * t29 + t91;
t10 = t33 * t76 + t56 * t75;
t9 = t33 * t75 - t56 * t76;
t86 = t10 * pkin(4) + qJ(5) * t9 + t95;
t83 = -pkin(4) * t6 - qJ(5) * t5 + t87;
t35 = t121 * t65 - t66 * t75;
t34 = t122 * t65 + t66 * t76;
t18 = t52 * t76 - t65 * t75;
t17 = t52 * t75 + t65 * t76;
t16 = -t121 * t56 + t58 * t75;
t15 = -t122 * t56 - t58 * t76;
t14 = -t121 * t53 + t55 * t75;
t13 = -t122 * t53 - t55 * t76;
t2 = t10 * t80 + t77 * t9;
t1 = -t10 * t77 + t80 * t9;
t3 = [(-t140 * mrSges(2,1) + t138 * mrSges(2,2) - m(3) * t132 - t57 * mrSges(3,1) - mrSges(3,3) * t106 - m(4) * t98 - t33 * mrSges(4,1) - m(5) * (t95 + t131) - m(6) * (t86 + t131) - m(7) * t86 - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t152 * t9 + t153 * t56 + t146 * t10 + t148 * t32) * g(2) + (t138 * mrSges(2,1) + t140 * mrSges(2,2) - m(3) * t109 + t54 * mrSges(3,1) - mrSges(3,3) * t107 - m(4) * t91 + t29 * mrSges(4,1) - m(5) * (t87 - t128) - m(6) * (t83 - t128) - m(7) * t83 + t161 * t6 + t162 * t5 - t153 * t53 - t148 * t28) * g(1) (-m(4) * t116 - m(5) * t97 - t154 * (t35 * pkin(4) + qJ(5) * t34 + t97) - t153 * t66 - t161 * t35 - t162 * t34 + t144 * t65) * g(3) + (-m(4) * t118 - m(5) * t100 - t154 * (t14 * pkin(4) + qJ(5) * t13 + t100) + t153 * t55 - t161 * t14 - t162 * t13 - t144 * t53) * g(2) + (-m(4) * t117 - m(5) * t99 - t154 * (t16 * pkin(4) + qJ(5) * t15 + t99) + t153 * t58 - t161 * t16 - t162 * t15 - t144 * t56) * g(1) (-m(5) * t113 - m(6) * (t113 + t149) - m(7) * (-t46 + t149) + t148 * t52 + t147 * t51) * g(3) + (-m(5) * t115 - m(6) * (t115 + t151) - m(7) * (-t23 + t151) + t148 * t29 + t147 * t28) * g(2) + (-m(5) * t114 - m(6) * (t114 + t150) - m(7) * (-t25 + t150) + t148 * t33 + t147 * t32) * g(1) (m(5) + t154) * (-g(1) * t32 - g(2) * t28 - g(3) * t51) t154 * (-g(1) * t9 - g(2) * t5 - g(3) * t17) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t5 * t80 - t6 * t77) * mrSges(7,1) + (-t5 * t77 - t6 * t80) * mrSges(7,2)) - g(3) * ((t17 * t80 - t18 * t77) * mrSges(7,1) + (-t17 * t77 - t18 * t80) * mrSges(7,2))];
taug  = t3(:);
