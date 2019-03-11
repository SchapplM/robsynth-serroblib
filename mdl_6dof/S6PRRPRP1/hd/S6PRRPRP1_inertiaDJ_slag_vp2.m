% Calculate time derivative of joint inertia matrix for
% S6PRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-03-08 21:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRP1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:22:47
% EndTime: 2019-03-08 21:22:56
% DurationCPUTime: 4.12s
% Computational Cost: add. (2580->329), mult. (6466->484), div. (0->0), fcn. (5977->10), ass. (0->157)
t230 = Ifges(6,4) + Ifges(7,4);
t225 = Ifges(6,1) + Ifges(7,1);
t223 = Ifges(6,2) + Ifges(7,2);
t112 = cos(qJ(3));
t105 = sin(pkin(11));
t109 = sin(qJ(3));
t164 = t105 * t109;
t165 = cos(pkin(11));
t116 = t112 * t165 - t164;
t188 = Ifges(6,5) + Ifges(7,5);
t229 = t116 * t188;
t216 = Ifges(6,6) + Ifges(7,6);
t228 = t216 * t116;
t111 = cos(qJ(5));
t227 = t230 * t111;
t108 = sin(qJ(5));
t226 = t230 * t108;
t222 = -t223 * t108 + t227;
t221 = t225 * t111 - t226;
t220 = t108 * t188 + t216 * t111;
t219 = Ifges(6,3) + Ifges(7,3);
t156 = qJD(5) * t111;
t133 = t165 * t109;
t79 = t105 * t112 + t133;
t146 = t79 * t156;
t73 = t116 * qJD(3);
t172 = t108 * t73;
t119 = t146 + t172;
t157 = qJD(5) * t108;
t168 = t111 * t73;
t118 = t79 * t157 - t168;
t217 = m(7) + m(6);
t72 = t79 * qJD(3);
t215 = -t118 * t230 - t223 * t119 + t216 * t72;
t214 = -t225 * t118 - t119 * t230 + t188 * t72;
t213 = t222 * t79 - t228;
t182 = t221 * t79 - t229;
t212 = t222 * qJD(5);
t211 = t221 * qJD(5);
t210 = t223 * t111 + t226;
t209 = t225 * t108 + t227;
t208 = t188 * t156;
t141 = t165 * pkin(3);
t98 = -t141 - pkin(4);
t88 = -t111 * pkin(5) + t98;
t89 = -mrSges(7,1) * t111 + mrSges(7,2) * t108;
t206 = m(7) * t88 + t89;
t151 = pkin(5) * t157;
t81 = mrSges(7,1) * t157 + mrSges(7,2) * t156;
t205 = m(7) * t151 + t81;
t198 = m(5) * pkin(3);
t204 = t105 * t198 - mrSges(5,2);
t203 = m(6) * t98 - mrSges(6,1) * t111 + mrSges(6,2) * t108 - t165 * t198 - mrSges(5,1);
t202 = 0.2e1 * m(7);
t201 = -2 * mrSges(7,3);
t186 = -qJ(4) - pkin(8);
t134 = qJD(3) * t186;
t115 = -qJD(4) * t109 + t112 * t134;
t71 = qJD(4) * t112 + t109 * t134;
t36 = t105 * t71 - t165 * t115;
t200 = 0.2e1 * t36;
t91 = t186 * t112;
t57 = -t105 * t91 - t186 * t133;
t199 = 0.2e1 * t57;
t197 = m(7) * pkin(5);
t192 = pkin(3) * t105;
t191 = t36 * t57;
t106 = sin(pkin(6));
t113 = cos(qJ(2));
t160 = qJD(2) * t113;
t139 = t106 * t160;
t107 = cos(pkin(6));
t110 = sin(qJ(2));
t163 = t106 * t110;
t75 = t107 * t109 + t112 * t163;
t55 = -qJD(3) * t75 - t109 * t139;
t74 = t107 * t112 - t109 * t163;
t56 = qJD(3) * t74 + t112 * t139;
t20 = t105 * t56 - t165 * t55;
t40 = t105 * t75 - t165 * t74;
t190 = t40 * t20;
t189 = t72 * mrSges(5,3);
t22 = mrSges(7,1) * t72 + mrSges(7,3) * t118;
t23 = mrSges(6,1) * t72 + mrSges(6,3) * t118;
t185 = t22 + t23;
t24 = -mrSges(7,2) * t72 - mrSges(7,3) * t119;
t25 = -mrSges(6,2) * t72 - mrSges(6,3) * t119;
t184 = t24 + t25;
t171 = t108 * t79;
t48 = mrSges(7,2) * t116 - mrSges(7,3) * t171;
t49 = mrSges(6,2) * t116 - mrSges(6,3) * t171;
t181 = t48 + t49;
t167 = t111 * t79;
t50 = -mrSges(7,1) * t116 - mrSges(7,3) * t167;
t51 = -mrSges(6,1) * t116 - mrSges(6,3) * t167;
t180 = t50 + t51;
t99 = -pkin(3) * t112 - pkin(2);
t47 = -pkin(4) * t116 - pkin(9) * t79 + t99;
t58 = t164 * t186 - t165 * t91;
t54 = t111 * t58;
t19 = t108 * t47 + t54;
t179 = mrSges(6,2) * t111;
t174 = t108 * mrSges(6,3);
t173 = t108 * mrSges(7,3);
t170 = t111 * mrSges(6,3);
t169 = t111 * mrSges(7,3);
t97 = pkin(9) + t192;
t166 = qJ(6) + t97;
t162 = t106 * t113;
t161 = qJD(2) * t110;
t159 = qJD(3) * t109;
t158 = qJD(3) * t112;
t155 = 0.2e1 * t109;
t37 = t105 * t115 + t165 * t71;
t152 = pkin(3) * t159;
t38 = pkin(4) * t72 - pkin(9) * t73 + t152;
t153 = t108 * t38 + t111 * t37 + t47 * t156;
t150 = mrSges(7,1) + t197;
t140 = t106 * t161;
t44 = t72 * mrSges(5,1) + t73 * mrSges(5,2);
t138 = t216 * t108;
t135 = -t108 * t37 + t111 * t38;
t18 = -t108 * t58 + t111 * t47;
t132 = t188 * t168 + t219 * t72;
t131 = qJD(5) * t166;
t130 = t106 ^ 2 * t110 * t160;
t127 = t57 * t20 + t36 * t40;
t126 = -(2 * Ifges(5,4)) - t138;
t125 = mrSges(6,1) * t108 + t179;
t120 = -qJ(6) * t73 - qJD(6) * t79;
t41 = t105 * t74 + t165 * t75;
t26 = -t108 * t41 - t111 * t162;
t117 = t108 * t162 - t111 * t41;
t16 = t119 * mrSges(7,1) - t118 * mrSges(7,2);
t114 = -t55 * t109 + t56 * t112 + (-t109 * t75 - t112 * t74) * qJD(3);
t83 = (mrSges(4,1) * t109 + mrSges(4,2) * t112) * qJD(3);
t82 = t125 * qJD(5);
t77 = t166 * t111;
t76 = t166 * t108;
t60 = -qJD(6) * t108 - t111 * t131;
t59 = qJD(6) * t111 - t108 * t131;
t53 = -mrSges(5,1) * t116 + mrSges(5,2) * t79;
t46 = t125 * t79;
t45 = (mrSges(7,1) * t108 + mrSges(7,2) * t111) * t79;
t35 = pkin(5) * t171 + t57;
t21 = t105 * t55 + t165 * t56;
t17 = mrSges(6,1) * t119 - mrSges(6,2) * t118;
t15 = pkin(5) * t119 + t36;
t14 = -qJ(6) * t171 + t19;
t9 = -pkin(5) * t116 - qJ(6) * t167 + t18;
t8 = qJD(5) * t117 - t108 * t21 + t111 * t140;
t7 = qJD(5) * t26 + t108 * t140 + t111 * t21;
t6 = -t19 * qJD(5) + t135;
t5 = -t157 * t58 + t153;
t4 = -qJ(6) * t146 + (-qJD(5) * t58 + t120) * t108 + t153;
t3 = pkin(5) * t72 + t120 * t111 + (-t54 + (qJ(6) * t79 - t47) * t108) * qJD(5) + t135;
t1 = [0.2e1 * m(5) * (t41 * t21 - t130 + t190) + 0.2e1 * m(4) * (t74 * t55 + t75 * t56 - t130) + 0.2e1 * (-t117 * t7 + t26 * t8 + t190) * t217; t180 * t8 + t181 * t7 + (t16 + t17) * t40 - t184 * t117 + t185 * t26 + (t45 + t46) * t20 + (t116 * t21 + t20 * t79 + t40 * t73 - t41 * t72) * mrSges(5,3) + t114 * mrSges(4,3) + ((-t44 - t83) * t113 + (-t113 * mrSges(3,2) + (-mrSges(4,1) * t112 + mrSges(4,2) * t109 - mrSges(3,1) + t53) * t110) * qJD(2)) * t106 + m(5) * (t58 * t21 + t37 * t41 + (-t113 * t152 + t161 * t99) * t106 + t127) + m(6) * (-t117 * t5 + t18 * t8 + t19 * t7 + t26 * t6 + t127) + m(7) * (-t117 * t4 + t14 * t7 + t15 * t40 + t20 * t35 + t26 * t3 + t8 * t9) + (-pkin(2) * t140 + pkin(8) * t114) * m(4); -0.2e1 * t58 * t189 - 0.2e1 * pkin(2) * t83 + 0.2e1 * t14 * t24 + 0.2e1 * t15 * t45 + 0.2e1 * t35 * t16 + t17 * t199 + 0.2e1 * t18 * t23 + 0.2e1 * t19 * t25 + 0.2e1 * t9 * t22 + 0.2e1 * t3 * t50 + t46 * t200 + 0.2e1 * t4 * t48 + 0.2e1 * t99 * t44 + 0.2e1 * t5 * t49 + 0.2e1 * t6 * t51 + (-Ifges(4,4) * t109 + pkin(3) * t53) * qJD(3) * t155 + 0.2e1 * m(5) * (t152 * t99 + t37 * t58 + t191) + 0.2e1 * m(6) * (t18 * t6 + t19 * t5 + t191) + (t14 * t4 + t15 * t35 + t3 * t9) * t202 + (mrSges(5,3) * t199 - t108 * t213 + t111 * t182) * t73 + (0.2e1 * Ifges(4,4) * t112 + (Ifges(4,1) - Ifges(4,2)) * t155) * t158 - (-0.2e1 * t37 * mrSges(5,3) + t126 * t73 + ((2 * Ifges(5,2)) + t219) * t72 + t132) * t116 + (mrSges(5,3) * t200 + 0.2e1 * Ifges(5,1) * t73 + t214 * t111 - t215 * t108 + (t111 * t188 + t126) * t72 + ((-t213 + t228) * t111 + (-t182 + t229) * t108) * qJD(5)) * t79; t7 * t170 - t8 * t174 + t7 * t169 - t8 * t173 + m(7) * (-t117 * t59 + t26 * t60 + t7 * t77 - t76 * t8) + m(6) * (-t8 * t108 + t7 * t111 + (t108 * t117 - t111 * t26) * qJD(5)) * t97 - t56 * mrSges(4,2) + t55 * mrSges(4,1) + (t205 + t82) * t40 + t204 * t21 + (t203 + t206) * t20 + (mrSges(7,3) + mrSges(6,3)) * (t117 * t157 - t26 * t156); (t220 / 0.2e1 - Ifges(5,6)) * t72 + t204 * t37 + t203 * t36 - (-t216 * t157 + t208) * t116 / 0.2e1 + t211 * t167 / 0.2e1 - t212 * t171 / 0.2e1 + t182 * t156 / 0.2e1 + t214 * t108 / 0.2e1 + t215 * t111 / 0.2e1 + t209 * t168 / 0.2e1 + Ifges(4,5) * t158 + t4 * t169 + t5 * t170 + (m(6) * (-t108 * t6 + t111 * t5 + (-t108 * t19 - t111 * t18) * qJD(5)) - t51 * t156 - t49 * t157 + t111 * t25 - t108 * t23) * t97 - t189 * t192 - t6 * t174 - t3 * t173 - (t209 * t79 + t213) * t157 / 0.2e1 - Ifges(4,6) * t159 + (-t172 / 0.2e1 - t146 / 0.2e1) * t210 + (-mrSges(5,3) * t141 + Ifges(5,5)) * t73 + m(7) * (t14 * t59 + t15 * t88 + t151 * t35 - t3 * t76 + t4 * t77 + t60 * t9) + t45 * t151 + t59 * t48 + t60 * t50 - t76 * t22 + t77 * t24 + t35 * t81 + t57 * t82 + t88 * t16 + t15 * t89 + t98 * t17 + (-t156 * t18 - t157 * t19) * mrSges(6,3) + (-t14 * t157 - t156 * t9) * mrSges(7,3) + (-mrSges(4,1) * t158 + mrSges(4,2) * t159) * pkin(8); 0.2e1 * t88 * t81 + (t59 * t77 - t60 * t76) * t202 + 0.2e1 * t98 * t82 + (t60 * t201 + (0.2e1 * t206 * pkin(5) + t77 * t201 - t210) * qJD(5) + t211) * t108 + (0.2e1 * t59 * mrSges(7,3) + (-t201 * t76 + t209) * qJD(5) + t212) * t111; m(5) * t140 + (t108 * t7 + t111 * t8 + (-t108 * t26 - t111 * t117) * qJD(5)) * t217; m(5) * t152 + t185 * t111 + t184 * t108 + (-t108 * t180 + t111 * t181) * qJD(5) + m(7) * (t108 * t4 + t111 * t3 + (-t108 * t9 + t111 * t14) * qJD(5)) + m(6) * (t108 * t5 + t111 * t6 + (-t108 * t18 + t111 * t19) * qJD(5)) + t44; m(7) * (t108 * t59 + t111 * t60 + (t108 * t76 + t111 * t77) * qJD(5)); 0; (-mrSges(6,2) - mrSges(7,2)) * t7 + (mrSges(6,1) + t150) * t8; mrSges(6,1) * t6 + mrSges(7,1) * t3 - mrSges(6,2) * t5 - mrSges(7,2) * t4 - t73 * t138 + (m(7) * t3 + t22) * pkin(5) - t220 * t79 * qJD(5) + t132; -mrSges(7,2) * t59 + t150 * t60 + ((-mrSges(6,1) * t97 - mrSges(7,3) * pkin(5)) * t111 + (mrSges(6,2) * t97 - t216) * t108) * qJD(5) + t208; (-t179 + (-mrSges(6,1) - t197) * t108) * qJD(5) - t81; 0; m(7) * t20; m(7) * t15 + t16; t205; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
