% Calculate joint inertia matrix for
% S6RRPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR4_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR4_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR4_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:01:15
% EndTime: 2019-03-09 09:01:18
% DurationCPUTime: 1.25s
% Computational Cost: add. (2103->297), mult. (4966->411), div. (0->0), fcn. (5228->10), ass. (0->122)
t161 = -Ifges(5,4) + Ifges(4,5);
t160 = Ifges(5,1) + Ifges(3,3) + Ifges(4,3);
t104 = sin(pkin(11));
t106 = cos(pkin(11));
t105 = sin(pkin(6));
t113 = cos(qJ(2));
t132 = t105 * t113;
t110 = sin(qJ(2));
t133 = t105 * t110;
t57 = t104 * t133 - t106 * t132;
t159 = (Ifges(5,5) - Ifges(4,6)) * t57;
t90 = -pkin(2) * t106 - pkin(3);
t86 = -pkin(9) + t90;
t158 = -0.2e1 * t86;
t108 = sin(qJ(6));
t111 = cos(qJ(6));
t107 = cos(pkin(6));
t109 = sin(qJ(5));
t112 = cos(qJ(5));
t40 = t107 * t112 + t109 * t57;
t58 = (t104 * t113 + t106 * t110) * t105;
t28 = -t108 * t40 + t111 * t58;
t39 = t107 * t109 - t112 * t57;
t13 = -mrSges(7,2) * t39 + mrSges(7,3) * t28;
t29 = t108 * t58 + t111 * t40;
t14 = mrSges(7,1) * t39 - mrSges(7,3) * t29;
t157 = -t108 * t14 + t111 * t13;
t72 = -mrSges(7,1) * t111 + mrSges(7,2) * t108;
t156 = m(7) * pkin(5) + mrSges(6,1) - t72;
t88 = pkin(2) * t104 + qJ(4);
t154 = t88 ^ 2;
t153 = t28 / 0.2e1;
t152 = t29 / 0.2e1;
t74 = Ifges(7,5) * t108 + Ifges(7,6) * t111;
t151 = t74 / 0.2e1;
t150 = pkin(3) + pkin(9);
t149 = -t108 / 0.2e1;
t148 = t108 / 0.2e1;
t147 = t111 / 0.2e1;
t101 = t109 ^ 2;
t103 = t112 ^ 2;
t128 = t103 + t101;
t146 = m(6) * t128 + m(5);
t145 = pkin(1) * t107;
t144 = pkin(5) * t109;
t85 = t113 * t145;
t41 = pkin(2) * t107 + t85 + (-pkin(8) - qJ(3)) * t133;
t63 = pkin(8) * t132 + t110 * t145;
t48 = qJ(3) * t132 + t63;
t26 = -t104 * t48 + t106 * t41;
t15 = pkin(4) * t58 - t107 * t150 - t26;
t66 = (-pkin(2) * t113 - pkin(1)) * t105;
t117 = -qJ(4) * t58 + t66;
t20 = t150 * t57 + t117;
t5 = -t109 * t20 + t112 * t15;
t3 = -pkin(5) * t58 - t5;
t143 = t112 * t3;
t142 = t58 * Ifges(6,6);
t62 = -pkin(8) * t133 + t85;
t141 = t62 * mrSges(3,1);
t140 = t63 * mrSges(3,2);
t11 = -mrSges(7,1) * t28 + mrSges(7,2) * t29;
t32 = mrSges(6,1) * t58 - mrSges(6,3) * t40;
t139 = t11 - t32;
t6 = t109 * t15 + t112 * t20;
t27 = t104 * t41 + t106 * t48;
t138 = Ifges(7,4) * t108;
t137 = Ifges(7,4) * t111;
t135 = t109 * t86;
t131 = t108 * t112;
t130 = t111 * t112;
t129 = t108 ^ 2 + t111 ^ 2;
t7 = Ifges(7,5) * t29 + Ifges(7,6) * t28 + Ifges(7,3) * t39;
t127 = Ifges(6,5) * t40 - Ifges(6,6) * t39 + Ifges(6,3) * t58;
t22 = -qJ(4) * t107 - t27;
t43 = t58 * mrSges(5,1) + mrSges(5,2) * t107;
t126 = t128 * mrSges(6,3);
t124 = t129 * t112;
t18 = -pkin(4) * t57 - t22;
t10 = pkin(5) * t39 - pkin(10) * t40 + t18;
t4 = pkin(10) * t58 + t6;
t1 = t10 * t111 - t108 * t4;
t2 = t10 * t108 + t111 * t4;
t123 = -t1 * t108 + t111 * t2;
t122 = t109 * t6 + t112 * t5;
t73 = t109 * mrSges(6,1) + t112 * mrSges(6,2);
t121 = mrSges(7,1) * t108 + mrSges(7,2) * t111;
t64 = -pkin(10) * t112 + t144 + t88;
t33 = -t108 * t135 + t111 * t64;
t34 = t108 * t64 + t111 * t135;
t120 = -t108 * t33 + t111 * t34;
t68 = -mrSges(7,2) * t109 - mrSges(7,3) * t131;
t69 = mrSges(7,1) * t109 - mrSges(7,3) * t130;
t119 = -t108 * t69 + t111 * t68;
t59 = Ifges(7,5) * t130 - Ifges(7,6) * t131 + Ifges(7,3) * t109;
t31 = -mrSges(6,2) * t58 - mrSges(6,3) * t39;
t118 = t31 + t157;
t116 = Ifges(3,5) * t133 + Ifges(3,6) * t132 + t107 * t160 + t161 * t58 + t159;
t99 = Ifges(6,5) * t112;
t80 = t86 ^ 2;
t78 = Ifges(6,1) * t112 - Ifges(6,4) * t109;
t77 = Ifges(7,1) * t108 + t137;
t76 = Ifges(6,4) * t112 - Ifges(6,2) * t109;
t75 = Ifges(7,2) * t111 + t138;
t71 = t103 * t86;
t70 = t103 * t80;
t65 = t121 * t112;
t61 = Ifges(7,5) * t109 + (Ifges(7,1) * t111 - t138) * t112;
t60 = Ifges(7,6) * t109 + (-Ifges(7,2) * t108 + t137) * t112;
t51 = t58 * mrSges(5,3);
t50 = t58 * mrSges(4,2);
t45 = mrSges(4,1) * t107 - mrSges(4,3) * t58;
t44 = -mrSges(4,2) * t107 - mrSges(4,3) * t57;
t42 = mrSges(5,1) * t57 - mrSges(5,3) * t107;
t30 = pkin(3) * t57 + t117;
t23 = -pkin(3) * t107 - t26;
t21 = mrSges(6,1) * t39 + mrSges(6,2) * t40;
t17 = Ifges(6,1) * t40 - Ifges(6,4) * t39 + Ifges(6,5) * t58;
t16 = Ifges(6,4) * t40 - Ifges(6,2) * t39 + t142;
t9 = Ifges(7,1) * t29 + Ifges(7,4) * t28 + Ifges(7,5) * t39;
t8 = Ifges(7,4) * t29 + Ifges(7,2) * t28 + Ifges(7,6) * t39;
t12 = [m(3) * (t62 ^ 2 + t63 ^ 2) + 0.2e1 * t66 * (t57 * mrSges(4,1) + t50) + t40 * t17 + 0.2e1 * t22 * t42 + 0.2e1 * t23 * t43 + 0.2e1 * t27 * t44 + 0.2e1 * t26 * t45 + t28 * t8 + t29 * t9 + 0.2e1 * t6 * t31 + 0.2e1 * t5 * t32 + 0.2e1 * t2 * t13 + 0.2e1 * t1 * t14 + 0.2e1 * t18 * t21 + 0.2e1 * t3 * t11 + (Ifges(5,3) + Ifges(4,2)) * t57 ^ 2 + ((Ifges(3,5) * t110 + Ifges(3,6) * t113) * t107 + 0.2e1 * (-t110 * t62 + t113 * t63) * mrSges(3,3) + (t113 * (Ifges(3,4) * t110 + Ifges(3,2) * t113) + t110 * (Ifges(3,1) * t110 + Ifges(3,4) * t113) - 0.2e1 * pkin(1) * (-mrSges(3,1) * t113 + mrSges(3,2) * t110) + m(3) * pkin(1) ^ 2) * t105) * t105 + 0.2e1 * t30 * (-t57 * mrSges(5,2) - t51) + ((Ifges(5,2) + Ifges(4,1)) * t58 + 0.2e1 * (-Ifges(4,4) - Ifges(5,6)) * t57 + t161 * t107 + t127) * t58 + (0.2e1 * t141 - 0.2e1 * t140 + t159 + t116) * t107 + Ifges(2,3) + (t7 - t16) * t39 + m(4) * (t26 ^ 2 + t27 ^ 2 + t66 ^ 2) + m(5) * (t22 ^ 2 + t23 ^ 2 + t30 ^ 2) + m(6) * (t18 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2); m(5) * (-t22 * t88 + t23 * t90) + (-t76 / 0.2e1 + t59 / 0.2e1) * t39 + t58 * t99 / 0.2e1 + t116 + t18 * t73 + t40 * t78 / 0.2e1 + t90 * t43 + t3 * t65 + t2 * t68 + t1 * t69 - t140 + t141 + m(6) * (t122 * t86 + t18 * t88) + t33 * t14 + t34 * t13 + t23 * mrSges(5,2) + t26 * mrSges(4,1) - t27 * mrSges(4,2) - t22 * mrSges(5,3) + t61 * t152 + t60 * t153 + (t21 - t42) * t88 + (t17 / 0.2e1 + t9 * t147 + t8 * t149 - t5 * mrSges(6,3) - t139 * t86) * t112 + m(7) * (t1 * t33 - t143 * t86 + t2 * t34) + (-t142 / 0.2e1 - t16 / 0.2e1 + t7 / 0.2e1 + t86 * t31 - t6 * mrSges(6,3)) * t109 + (m(4) * (t104 * t27 + t106 * t26) + t104 * t44 + t106 * t45) * pkin(2); 0.2e1 * t90 * mrSges(5,2) + 0.2e1 * t33 * t69 + 0.2e1 * t34 * t68 + (t59 - t76) * t109 + (-t108 * t60 + t111 * t61 + t158 * t65 + t78) * t112 + m(7) * (t33 ^ 2 + t34 ^ 2 + t70) + m(6) * (t101 * t80 + t154 + t70) + m(5) * (t90 ^ 2 + t154) + m(4) * (t104 ^ 2 + t106 ^ 2) * pkin(2) ^ 2 + 0.2e1 * (t73 + mrSges(5,3)) * t88 + 0.2e1 * (mrSges(4,1) * t106 - mrSges(4,2) * t104) * pkin(2) + t126 * t158 + t160; t50 - t51 + (-mrSges(5,2) + mrSges(4,1)) * t57 + t139 * t109 + t118 * t112 + m(7) * (t109 * t3 + t112 * t123) + m(6) * (-t109 * t5 + t112 * t6) + m(5) * t30 + m(4) * t66; t109 * t65 + (m(7) * (t120 - t135) + t119) * t112; m(4) + m(7) * (t103 * t129 + t101) + t146; -t139 * t112 + t118 * t109 + m(7) * (t109 * t123 - t143) + m(6) * t122 + m(5) * t23 + t43; -t112 * t65 + mrSges(5,2) + t119 * t109 - t126 + m(7) * (t109 * t120 + t71) + m(6) * (t101 * t86 + t71) + m(5) * t90; m(7) * (-0.1e1 + t129) * t112 * t109; m(7) * (t101 * t129 + t103) + t146; t5 * mrSges(6,1) - t6 * mrSges(6,2) + t123 * mrSges(7,3) + t8 * t147 + t9 * t148 + t39 * t151 + t77 * t152 + t75 * t153 + t3 * t72 + t127 + (-m(7) * t3 - t11) * pkin(5) + (m(7) * t123 + t157) * pkin(10); -pkin(5) * t65 + t61 * t148 + t60 * t147 + t99 + (m(7) * t120 + t119) * pkin(10) + (t77 * t147 + t75 * t149 + t156 * t86) * t112 + t120 * mrSges(7,3) + (-t86 * mrSges(6,2) - Ifges(6,6) + t151) * t109; t109 * t72 + m(7) * (pkin(10) * t124 - t144) + mrSges(7,3) * t124 - t73; t156 * t112 + (-mrSges(6,2) + (m(7) * pkin(10) + mrSges(7,3)) * t129) * t109; Ifges(6,3) + m(7) * (pkin(10) ^ 2 * t129 + pkin(5) ^ 2) + t108 * t77 + t111 * t75 - 0.2e1 * pkin(5) * t72 + 0.2e1 * t129 * pkin(10) * mrSges(7,3); mrSges(7,1) * t1 - mrSges(7,2) * t2 + t7; mrSges(7,1) * t33 - mrSges(7,2) * t34 + t59; -t65; -t121 * t109; -pkin(10) * t121 + t74; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t12(1) t12(2) t12(4) t12(7) t12(11) t12(16); t12(2) t12(3) t12(5) t12(8) t12(12) t12(17); t12(4) t12(5) t12(6) t12(9) t12(13) t12(18); t12(7) t12(8) t12(9) t12(10) t12(14) t12(19); t12(11) t12(12) t12(13) t12(14) t12(15) t12(20); t12(16) t12(17) t12(18) t12(19) t12(20) t12(21);];
Mq  = res;
