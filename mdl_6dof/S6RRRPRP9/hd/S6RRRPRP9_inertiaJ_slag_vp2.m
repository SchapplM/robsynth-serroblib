% Calculate joint inertia matrix for
% S6RRRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 17:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP9_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP9_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP9_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP9_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP9_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP9_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:21:54
% EndTime: 2019-03-09 17:21:56
% DurationCPUTime: 1.23s
% Computational Cost: add. (1333->330), mult. (2474->433), div. (0->0), fcn. (2201->6), ass. (0->118)
t136 = Ifges(7,2) + Ifges(6,3);
t106 = sin(qJ(3));
t109 = cos(qJ(3));
t152 = t106 ^ 2 + t109 ^ 2;
t151 = 2 * pkin(7);
t139 = mrSges(6,2) - mrSges(7,3);
t138 = mrSges(7,2) + mrSges(6,3);
t107 = sin(qJ(2));
t125 = t107 * t109;
t127 = t106 * t107;
t150 = -Ifges(5,6) * t127 + (-Ifges(5,4) - Ifges(4,5)) * t125;
t149 = -m(7) * pkin(5) - mrSges(7,1);
t148 = mrSges(6,1) - t149;
t147 = m(7) * qJ(6) - t139;
t146 = 2 * mrSges(7,1);
t111 = -pkin(3) - pkin(4);
t145 = pkin(8) - pkin(9);
t144 = pkin(3) * t106;
t110 = cos(qJ(2));
t143 = pkin(7) * t110;
t105 = sin(qJ(5));
t108 = cos(qJ(5));
t66 = -qJ(4) * t105 + t108 * t111;
t142 = t66 * mrSges(6,1);
t67 = t108 * qJ(4) + t105 * t111;
t141 = t67 * mrSges(6,2);
t140 = mrSges(6,1) + mrSges(7,1);
t137 = Ifges(5,2) + Ifges(4,3);
t100 = t110 * pkin(3);
t69 = -pkin(2) * t110 - pkin(8) * t107 - pkin(1);
t84 = t106 * t143;
t14 = pkin(4) * t110 + t100 + t84 + (-pkin(9) * t107 - t69) * t109;
t36 = t106 * t69 + t109 * t143;
t29 = -qJ(4) * t110 + t36;
t16 = pkin(9) * t127 + t29;
t4 = t105 * t14 + t108 * t16;
t126 = t106 * t108;
t47 = t105 * t125 - t107 * t126;
t31 = -mrSges(7,2) * t47 + mrSges(7,3) * t110;
t32 = -mrSges(6,2) * t110 - mrSges(6,3) * t47;
t135 = t31 + t32;
t57 = t105 * t106 + t108 * t109;
t48 = t57 * t107;
t33 = mrSges(6,1) * t110 - mrSges(6,3) * t48;
t34 = -t110 * mrSges(7,1) + t48 * mrSges(7,2);
t134 = t33 - t34;
t63 = t110 * mrSges(5,1) + mrSges(5,2) * t125;
t133 = t152 * pkin(8) ^ 2;
t132 = Ifges(4,4) * t106;
t131 = Ifges(4,4) * t109;
t130 = Ifges(5,5) * t106;
t129 = Ifges(5,5) * t109;
t128 = t110 * Ifges(5,6);
t123 = t145 * t106;
t76 = t145 * t109;
t26 = t105 * t76 - t108 * t123;
t28 = t105 * t123 + t108 * t76;
t124 = t26 ^ 2 + t28 ^ 2;
t68 = -t109 * pkin(3) - t106 * qJ(4) - pkin(2);
t51 = t109 * pkin(4) - t68;
t35 = t109 * t69 - t84;
t52 = Ifges(7,6) * t57;
t53 = Ifges(6,6) * t57;
t58 = -t105 * t109 + t126;
t54 = Ifges(6,5) * t58;
t55 = Ifges(7,4) * t58;
t121 = t52 + t54 + t55 - t53;
t119 = t106 * mrSges(4,1) + t109 * mrSges(4,2);
t118 = t106 * mrSges(5,1) - t109 * mrSges(5,3);
t117 = qJ(4) * t109 - t144;
t3 = -t105 * t16 + t108 * t14;
t115 = (Ifges(7,4) + Ifges(6,5)) * t48 + (-Ifges(6,6) + Ifges(7,6)) * t47 + t136 * t110;
t79 = qJ(4) * t125;
t25 = t79 + (t111 * t106 - pkin(7)) * t107;
t1 = qJ(6) * t110 + t4;
t2 = -pkin(5) * t110 - t3;
t114 = t3 * mrSges(6,1) - t2 * mrSges(7,1) - t4 * mrSges(6,2) + t1 * mrSges(7,3) + t115;
t113 = pkin(7) ^ 2;
t104 = t110 ^ 2;
t102 = t107 ^ 2;
t96 = t102 * t113;
t94 = Ifges(5,4) * t106;
t93 = Ifges(4,5) * t106;
t91 = Ifges(4,6) * t109;
t75 = Ifges(4,1) * t106 + t131;
t74 = Ifges(5,1) * t106 - t129;
t73 = Ifges(4,2) * t109 + t132;
t72 = -Ifges(5,3) * t109 + t130;
t71 = -mrSges(4,1) * t109 + mrSges(4,2) * t106;
t70 = -mrSges(5,1) * t109 - mrSges(5,3) * t106;
t65 = pkin(5) - t66;
t64 = -mrSges(5,2) * t127 - mrSges(5,3) * t110;
t62 = -mrSges(4,1) * t110 - mrSges(4,3) * t125;
t61 = mrSges(4,2) * t110 - mrSges(4,3) * t127;
t60 = -qJ(6) + t67;
t50 = t119 * t107;
t49 = t118 * t107;
t46 = -t79 + (pkin(7) + t144) * t107;
t45 = -Ifges(4,5) * t110 + (Ifges(4,1) * t109 - t132) * t107;
t44 = -Ifges(5,4) * t110 + (Ifges(5,1) * t109 + t130) * t107;
t43 = -Ifges(4,6) * t110 + (-Ifges(4,2) * t106 + t131) * t107;
t42 = -t128 + (Ifges(5,3) * t106 + t129) * t107;
t30 = t100 - t35;
t22 = Ifges(6,1) * t58 - Ifges(6,4) * t57;
t21 = Ifges(7,1) * t58 + Ifges(7,5) * t57;
t20 = Ifges(6,4) * t58 - Ifges(6,2) * t57;
t19 = Ifges(7,5) * t58 + Ifges(7,3) * t57;
t18 = mrSges(6,1) * t57 + mrSges(6,2) * t58;
t17 = mrSges(7,1) * t57 - mrSges(7,3) * t58;
t12 = mrSges(6,1) * t47 + mrSges(6,2) * t48;
t11 = mrSges(7,1) * t47 - mrSges(7,3) * t48;
t10 = pkin(5) * t57 - qJ(6) * t58 + t51;
t9 = Ifges(6,1) * t48 - Ifges(6,4) * t47 + Ifges(6,5) * t110;
t8 = Ifges(7,1) * t48 + Ifges(7,4) * t110 + Ifges(7,5) * t47;
t7 = Ifges(6,4) * t48 - Ifges(6,2) * t47 + Ifges(6,6) * t110;
t6 = Ifges(7,5) * t48 + Ifges(7,6) * t110 + Ifges(7,3) * t47;
t5 = pkin(5) * t47 - qJ(6) * t48 + t25;
t13 = [0.2e1 * t1 * t31 + 0.2e1 * t5 * t11 + 0.2e1 * t25 * t12 + 0.2e1 * t2 * t34 + 0.2e1 * t29 * t64 + 0.2e1 * t3 * t33 + 0.2e1 * t30 * t63 + 0.2e1 * t4 * t32 + 0.2e1 * t35 * t62 + 0.2e1 * t36 * t61 + 0.2e1 * t46 * t49 + Ifges(2,3) + (t8 + t9) * t48 + (t6 - t7) * t47 + (t102 + t104) * mrSges(3,3) * t151 + (-0.2e1 * pkin(1) * mrSges(3,2) + Ifges(3,1) * t107 + t50 * t151 + (t44 + t45) * t109 + (t42 - t43) * t106) * t107 + (0.2e1 * pkin(1) * mrSges(3,1) + (Ifges(3,2) + t137) * t110 + (Ifges(4,6) * t106 + (2 * Ifges(3,4))) * t107 + t115 + t150) * t110 + m(3) * (pkin(1) ^ 2 + t104 * t113 + t96) + m(4) * (t35 ^ 2 + t36 ^ 2 + t96) + m(5) * (t29 ^ 2 + t30 ^ 2 + t46 ^ 2) + m(6) * (t25 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2); Ifges(3,5) * t107 - pkin(2) * t50 + t10 * t11 + t51 * t12 + t5 * t17 + t25 * t18 + t46 * t70 + t68 * t49 + (t21 / 0.2e1 + t22 / 0.2e1) * t48 + (t19 / 0.2e1 - t20 / 0.2e1) * t47 + t135 * t28 - t134 * t26 + (-t93 / 0.2e1 - t91 / 0.2e1 + t54 / 0.2e1 - t53 / 0.2e1 + t55 / 0.2e1 + t52 / 0.2e1 + Ifges(3,6) - t94 / 0.2e1) * t110 + (-t110 * mrSges(3,2) + (-mrSges(3,1) + t71) * t107) * pkin(7) + (t2 * mrSges(7,2) - t3 * mrSges(6,3) + t8 / 0.2e1 + t9 / 0.2e1) * t58 + (-t1 * mrSges(7,2) - t4 * mrSges(6,3) + t6 / 0.2e1 - t7 / 0.2e1) * t57 + (-t42 / 0.2e1 + t43 / 0.2e1 + t128 / 0.2e1 + t36 * mrSges(4,3) + t29 * mrSges(5,2) + (t74 / 0.2e1 + t75 / 0.2e1) * t107 + (t61 + t64) * pkin(8)) * t109 + (t44 / 0.2e1 + t45 / 0.2e1 - t35 * mrSges(4,3) + t30 * mrSges(5,2) + (t72 / 0.2e1 - t73 / 0.2e1) * t107 + (-t62 + t63) * pkin(8)) * t106 + m(4) * (-pkin(2) * pkin(7) * t107 + (-t106 * t35 + t109 * t36) * pkin(8)) + m(5) * (t46 * t68 + (t106 * t30 + t109 * t29) * pkin(8)) + m(6) * (t25 * t51 - t26 * t3 + t28 * t4) + m(7) * (t1 * t28 + t10 * t5 + t2 * t26); -0.2e1 * pkin(2) * t71 + 0.2e1 * t10 * t17 + 0.2e1 * t51 * t18 + 0.2e1 * t68 * t70 + Ifges(3,3) + (-t72 + t73) * t109 + (t74 + t75) * t106 + (0.2e1 * t138 * t26 + t21 + t22) * t58 + (-0.2e1 * t138 * t28 + t19 - t20) * t57 + m(5) * (t68 ^ 2 + t133) + m(4) * (pkin(2) ^ 2 + t133) + m(6) * (t51 ^ 2 + t124) + m(7) * (t10 ^ 2 + t124) + 0.2e1 * (mrSges(5,2) + mrSges(4,3)) * pkin(8) * t152; m(7) * (t1 * t60 + t2 * t65) + m(6) * (t3 * t66 + t4 * t67) + m(5) * (-pkin(3) * t30 + qJ(4) * t29) - t137 * t110 - Ifges(4,6) * t127 - t114 + t66 * t33 + t67 * t32 + t60 * t31 - pkin(3) * t63 + qJ(4) * t64 + t65 * t34 + t29 * mrSges(5,3) - t30 * mrSges(5,1) + t35 * mrSges(4,1) - t36 * mrSges(4,2) - t150; -Ifges(5,6) * t109 + t91 + t93 + t94 + t139 * t28 + t140 * t26 + (-t57 * t67 - t58 * t66) * mrSges(6,3) + (-t57 * t60 + t58 * t65) * mrSges(7,2) + t117 * mrSges(5,2) + m(7) * (t26 * t65 + t28 * t60) + m(6) * (-t26 * t66 + t28 * t67) + (m(5) * t117 - t118 - t119) * pkin(8) - t121; 0.2e1 * pkin(3) * mrSges(5,1) - 0.2e1 * t142 + t65 * t146 + 0.2e1 * t141 + 0.2e1 * qJ(4) * mrSges(5,3) - 0.2e1 * t60 * mrSges(7,3) + m(6) * (t66 ^ 2 + t67 ^ 2) + m(7) * (t60 ^ 2 + t65 ^ 2) + m(5) * (pkin(3) ^ 2 + qJ(4) ^ 2) + t136 + t137; t134 * t108 + t135 * t105 + m(7) * (t1 * t105 - t108 * t2) + m(6) * (t105 * t4 + t108 * t3) + m(5) * t30 + t63; (m(5) * pkin(8) + mrSges(5,2)) * t106 + (m(6) + m(7)) * (t105 * t28 - t108 * t26) + t138 * (-t105 * t57 - t108 * t58); -m(5) * pkin(3) - mrSges(5,1) - t140 * t108 + t139 * t105 + m(6) * (t105 * t67 + t108 * t66) + m(7) * (t105 * t60 - t108 * t65); m(5) + 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (t105 ^ 2 + t108 ^ 2); m(7) * (-pkin(5) * t2 + qJ(6) * t1) + qJ(6) * t31 - pkin(5) * t34 + t114; (-pkin(5) * t58 - qJ(6) * t57) * mrSges(7,2) + t121 + t147 * t28 - t148 * t26; t142 - t141 + m(7) * (-pkin(5) * t65 + qJ(6) * t60) + (t60 - qJ(6)) * mrSges(7,3) + (-t65 - pkin(5)) * mrSges(7,1) - t136; t147 * t105 + t148 * t108; pkin(5) * t146 + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) + t136; m(7) * t2 + t34; m(7) * t26 + t58 * mrSges(7,2); m(7) * t65 + mrSges(7,1); -m(7) * t108; t149; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t13(1) t13(2) t13(4) t13(7) t13(11) t13(16); t13(2) t13(3) t13(5) t13(8) t13(12) t13(17); t13(4) t13(5) t13(6) t13(9) t13(13) t13(18); t13(7) t13(8) t13(9) t13(10) t13(14) t13(19); t13(11) t13(12) t13(13) t13(14) t13(15) t13(20); t13(16) t13(17) t13(18) t13(19) t13(20) t13(21);];
Mq  = res;
