% Calculate time derivative of joint inertia matrix for
% S6PRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRP4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP4_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP4_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP4_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP4_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:40:19
% EndTime: 2019-03-08 21:40:23
% DurationCPUTime: 2.40s
% Computational Cost: add. (1446->328), mult. (3645->460), div. (0->0), fcn. (2785->8), ass. (0->142)
t182 = Ifges(6,5) + Ifges(7,5);
t97 = sin(qJ(3));
t184 = t182 * t97;
t183 = Ifges(6,6) + Ifges(7,6);
t168 = pkin(4) + pkin(8);
t181 = Ifges(6,3) + Ifges(7,3);
t100 = cos(qJ(3));
t99 = cos(qJ(5));
t164 = Ifges(7,4) * t99;
t96 = sin(qJ(5));
t119 = Ifges(7,1) * t96 + t164;
t166 = Ifges(6,4) * t99;
t120 = Ifges(6,1) * t96 + t166;
t154 = t184 + (-t119 - t120) * t100;
t180 = t154 * t96;
t165 = Ifges(7,4) * t96;
t117 = Ifges(7,2) * t99 + t165;
t167 = Ifges(6,4) * t96;
t118 = Ifges(6,2) * t99 + t167;
t155 = t183 * t97 + (-t117 - t118) * t100;
t179 = t155 * t99;
t140 = qJD(5) * t100;
t127 = t96 * t140;
t146 = qJD(3) * t97;
t108 = t99 * t146 + t127;
t128 = t99 * t140;
t135 = t96 * t146;
t107 = -t128 + t135;
t178 = -t182 * t96 - t183 * t99;
t173 = -2 * mrSges(7,3);
t177 = m(7) + m(6);
t102 = -pkin(3) - pkin(9);
t126 = -qJ(4) * t97 - pkin(2);
t46 = t100 * t102 + t126;
t75 = t168 * t97;
t19 = t99 * t46 + t96 * t75;
t101 = cos(qJ(2));
t94 = sin(pkin(6));
t149 = t101 * t94;
t98 = sin(qJ(2));
t163 = t94 * t98;
t138 = t97 * t163;
t95 = cos(pkin(6));
t40 = -t95 * t100 + t138;
t111 = t149 * t99 - t40 * t96;
t22 = t149 * t96 + t40 * t99;
t147 = qJD(2) * t98;
t136 = t94 * t147;
t142 = qJD(2) * t101;
t129 = t94 * t142;
t41 = t100 * t163 + t95 * t97;
t20 = qJD(3) * t41 + t129 * t97;
t7 = qJD(5) * t111 - t136 * t96 + t20 * t99;
t8 = qJD(5) * t22 + t136 * t99 + t20 * t96;
t175 = qJD(5) * (t111 * t99 + t22 * t96) - t7 * t99 - t8 * t96;
t174 = 2 * m(7);
t172 = m(5) / 0.2e1;
t21 = -qJD(3) * t138 + (qJD(3) * t95 + t129) * t100;
t9 = t41 * t21;
t162 = -mrSges(4,1) + mrSges(5,2);
t161 = -mrSges(6,2) - mrSges(7,2);
t160 = mrSges(5,3) - mrSges(4,2);
t141 = qJD(3) * t100;
t26 = -mrSges(7,2) * t141 + mrSges(7,3) * t108;
t27 = -mrSges(6,2) * t141 + mrSges(6,3) * t108;
t157 = t26 + t27;
t28 = mrSges(7,1) * t141 - mrSges(7,3) * t107;
t29 = mrSges(6,1) * t141 - mrSges(6,3) * t107;
t156 = t28 + t29;
t151 = t100 * t96;
t58 = mrSges(7,1) * t97 + mrSges(7,3) * t151;
t59 = mrSges(6,1) * t97 + mrSges(6,3) * t151;
t153 = t58 + t59;
t150 = t100 * t99;
t60 = -mrSges(7,2) * t97 - mrSges(7,3) * t150;
t61 = -mrSges(6,2) * t97 - mrSges(6,3) * t150;
t152 = t60 + t61;
t76 = t168 * t100;
t148 = qJ(6) * t100;
t145 = qJD(5) * t96;
t144 = qJD(5) * t99;
t143 = qJ(6) - t102;
t139 = qJD(6) * t100;
t71 = -Ifges(7,2) * t96 + t164;
t72 = -Ifges(6,2) * t96 + t166;
t133 = t71 / 0.2e1 + t72 / 0.2e1;
t73 = Ifges(7,1) * t99 - t165;
t74 = Ifges(6,1) * t99 - t167;
t132 = t73 / 0.2e1 + t74 / 0.2e1;
t131 = m(5) * pkin(8) + mrSges(5,1);
t130 = (m(7) * pkin(5)) + mrSges(7,1);
t65 = t143 * t99;
t125 = -t46 + t148;
t124 = pkin(3) * t146 - qJD(4) * t97;
t123 = mrSges(6,1) + t130;
t50 = mrSges(7,1) * t144 - mrSges(7,2) * t145;
t121 = mrSges(6,1) * t99 - mrSges(6,2) * t96;
t49 = t99 * t75;
t10 = pkin(5) * t97 + t125 * t96 + t49;
t15 = -t148 * t99 + t19;
t116 = t10 * t96 - t15 * t99;
t18 = -t46 * t96 + t49;
t115 = t18 * t96 - t19 * t99;
t113 = t21 * t100 + t20 * t97;
t112 = t21 * qJ(4) + t41 * qJD(4);
t30 = (pkin(9) * t97 - qJ(4) * t100) * qJD(3) + t124;
t63 = t168 * t141;
t5 = t75 * t144 - t145 * t46 + t99 * t30 + t96 * t63;
t109 = t108 * t183 + t182 * t135 + t181 * t141;
t16 = -mrSges(7,1) * t108 + mrSges(7,2) * t107;
t86 = pkin(5) * t96 + qJ(4);
t84 = pkin(5) * t144 + qJD(4);
t70 = mrSges(6,1) * t96 + mrSges(6,2) * t99;
t69 = mrSges(7,1) * t96 + mrSges(7,2) * t99;
t68 = t100 * mrSges(5,2) - t97 * mrSges(5,3);
t66 = -pkin(3) * t100 + t126;
t64 = t143 * t96;
t62 = t168 * t146;
t57 = t120 * qJD(5);
t56 = t119 * qJD(5);
t55 = t118 * qJD(5);
t54 = t117 * qJD(5);
t53 = (mrSges(4,1) * t97 + mrSges(4,2) * t100) * qJD(3);
t52 = (-mrSges(5,2) * t97 - mrSges(5,3) * t100) * qJD(3);
t51 = t121 * qJD(5);
t45 = t121 * t100;
t44 = (mrSges(7,1) * t99 - mrSges(7,2) * t96) * t100;
t43 = t99 * t63;
t39 = pkin(5) * t150 + t76;
t38 = -qJ(4) * t141 + t124;
t36 = -qJD(5) * t65 - t96 * qJD(6);
t35 = -t99 * qJD(6) + t143 * t145;
t24 = -pkin(5) * t127 + (-pkin(5) * t99 - t168) * t146;
t17 = -mrSges(6,1) * t108 + mrSges(6,2) * t107;
t14 = -t74 * t140 + (t100 * Ifges(6,5) + t120 * t97) * qJD(3);
t13 = -t73 * t140 + (t100 * Ifges(7,5) + t119 * t97) * qJD(3);
t12 = -t72 * t140 + (t100 * Ifges(6,6) + t118 * t97) * qJD(3);
t11 = -t71 * t140 + (t100 * Ifges(7,6) + t117 * t97) * qJD(3);
t6 = -t19 * qJD(5) - t30 * t96 + t43;
t4 = qJ(6) * t108 - t139 * t99 + t5;
t3 = pkin(5) * t141 + t43 + t125 * t144 + (-qJ(6) * t146 - qJD(5) * t75 + t139 - t30) * t96;
t1 = [0.2e1 * (-t111 * t8 + t22 * t7 + t9) * t177 + 0.2e1 * (m(5) + m(4)) * (-t142 * t94 ^ 2 * t98 + t20 * t40 + t9); t152 * t8 + t153 * t7 + (t16 + t17) * t41 - t157 * t111 + t156 * t22 + (t44 + t45) * t21 + m(6) * (-t111 * t5 + t18 * t7 + t19 * t8 + t21 * t76 + t22 * t6 - t41 * t62) + m(7) * (t10 * t7 - t111 * t4 + t15 * t8 + t21 * t39 + t22 * t3 + t24 * t41) + 0.2e1 * (m(4) / 0.2e1 + t172) * (t141 * t40 - t146 * t41 + t113) * pkin(8) + (mrSges(4,3) + mrSges(5,1)) * ((t100 * t40 - t41 * t97) * qJD(3) + t113) + ((-t52 - t53) * t101 + (-t101 * mrSges(3,2) + (-t100 * mrSges(4,1) + t97 * mrSges(4,2) - mrSges(3,1) + t68) * t98) * qJD(2) - m(4) * pkin(2) * t147 + 0.2e1 * (-t101 * t38 + t147 * t66) * t172) * t94; -0.2e1 * pkin(2) * t53 + 0.2e1 * t10 * t28 + 0.2e1 * t15 * t26 + 0.2e1 * t39 * t16 + 0.2e1 * t76 * t17 + 0.2e1 * t18 * t29 + 0.2e1 * t19 * t27 + 0.2e1 * t24 * t44 + 0.2e1 * t3 * t58 + 0.2e1 * t4 * t60 - 0.2e1 * t62 * t45 + 0.2e1 * t5 * t61 + 0.2e1 * t66 * t52 + 0.2e1 * t6 * t59 + 0.2e1 * (m(5) * t66 + t68) * t38 + 0.2e1 * m(6) * (t18 * t6 + t19 * t5 - t62 * t76) + (t10 * t3 + t15 * t4 + t24 * t39) * t174 + ((0.2e1 * (-Ifges(4,4) - Ifges(5,6)) * t97 + t179 + t180) * qJD(3) + t109) * t97 + ((-t11 - t12) * t99 + (-t13 - t14) * t96 + (t155 * t96 + (-t154 - t184) * t99) * qJD(5) + ((0.2e1 * Ifges(4,4) + 0.2e1 * Ifges(5,6) + t178) * t100 + ((2 * Ifges(4,1)) - (2 * Ifges(4,2)) + (2 * Ifges(5,2)) - (2 * Ifges(5,3)) + t181) * t97) * qJD(3)) * t100; (t50 + t51) * t41 + t162 * t20 + (t69 + t70 + t160) * t21 + m(7) * (-t111 * t36 + t21 * t86 + t22 * t35 + t41 * t84 - t64 * t8 - t65 * t7) + m(5) * (-pkin(3) * t20 + t112) + (mrSges(7,3) + mrSges(6,3)) * t175 + (-t102 * t175 + t112) * m(6); qJ(4) * t17 + qJD(4) * t45 + t86 * t16 + t24 * t69 - t64 * t26 - t65 * t28 + t35 * t58 + t36 * t60 + t39 * t50 + t84 * t44 + t76 * t51 - t62 * t70 + m(7) * (t10 * t35 + t15 * t36 + t24 * t86 - t3 * t65 + t39 * t84 - t4 * t64) + m(6) * (-qJ(4) * t62 + qJD(4) * t76) + (-t3 * mrSges(7,3) - t6 * mrSges(6,3) + t13 / 0.2e1 + t14 / 0.2e1 + (m(6) * t6 + t29) * t102) * t99 + (-t4 * mrSges(7,3) - t5 * mrSges(6,3) - t11 / 0.2e1 - t12 / 0.2e1 + (m(6) * t5 + t27) * t102) * t96 + ((t54 / 0.2e1 + t55 / 0.2e1) * t99 + (t56 / 0.2e1 + t57 / 0.2e1) * t96 + t131 * qJD(4)) * t100 + ((-pkin(3) * mrSges(5,1) - Ifges(5,4) + Ifges(4,5) + (Ifges(6,5) / 0.2e1 + Ifges(7,5) / 0.2e1) * t99 + (-Ifges(6,6) / 0.2e1 - Ifges(7,6) / 0.2e1) * t96 + (-m(5) * pkin(3) + t162) * pkin(8)) * t100 + (-qJ(4) * mrSges(5,1) + Ifges(5,5) - Ifges(4,6) + t133 * t99 + t132 * t96 + (-m(5) * qJ(4) - t160) * pkin(8)) * t97) * qJD(3) + (t116 * mrSges(7,3) + t115 * mrSges(6,3) + (-m(6) * t115 - t96 * t59 + t99 * t61) * t102 + (-t132 * t99 + t133 * t96) * t100 - t180 / 0.2e1 + t178 * t97 / 0.2e1 - t179 / 0.2e1) * qJD(5); (-t35 * t65 - t36 * t64 + t84 * t86) * t174 + 0.2e1 * qJ(4) * t51 + 0.2e1 * t84 * t69 + 0.2e1 * t86 * t50 + (t35 * t173 - t56 - t57) * t99 + (t36 * t173 + t54 + t55) * t96 + 0.2e1 * (mrSges(5,3) + t70 + (m(5) + m(6)) * qJ(4)) * qJD(4) + ((-t173 * t64 - t71 - t72) * t99 + (t65 * t173 - t73 - t74) * t96) * qJD(5); m(5) * t20 - t175 * t177; t156 * t99 + t157 * t96 + t131 * t141 + (t152 * t99 - t153 * t96) * qJD(5) + m(7) * (-t116 * qJD(5) + t3 * t99 + t4 * t96) + m(6) * (-t115 * qJD(5) + t5 * t96 + t6 * t99); m(7) * (t35 * t99 + t36 * t96 + (-t64 * t99 + t65 * t96) * qJD(5)); 0; t123 * t7 + t161 * t8; mrSges(6,1) * t6 + mrSges(7,1) * t3 - mrSges(6,2) * t5 - mrSges(7,2) * t4 - t182 * t128 + (m(7) * t3 + t28) * pkin(5) + t109; -t36 * mrSges(7,2) + t130 * t35 + ((-mrSges(6,2) * t102 - t183) * t99 + (-mrSges(6,1) * t102 + (mrSges(7,3) * pkin(5)) - t182) * t96) * qJD(5); (-t123 * t96 + t161 * t99) * qJD(5); 0; m(7) * t21; m(7) * t24 + t16; m(7) * t84 + t50; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
