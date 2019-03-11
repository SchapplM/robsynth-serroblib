% Calculate time derivative of joint inertia matrix for
% S6PRRPRP2
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
% Datum: 2019-03-08 21:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRP2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP2_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP2_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP2_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP2_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:28:27
% EndTime: 2019-03-08 21:28:35
% DurationCPUTime: 3.56s
% Computational Cost: add. (2588->335), mult. (6484->490), div. (0->0), fcn. (6000->10), ass. (0->157)
t219 = Ifges(6,1) + Ifges(7,1);
t112 = cos(qJ(3));
t105 = sin(pkin(11));
t109 = sin(qJ(3));
t164 = t105 * t109;
t166 = cos(pkin(11));
t118 = t112 * t166 - t164;
t190 = Ifges(7,4) + Ifges(6,5);
t220 = t118 * t190;
t218 = -Ifges(6,6) + Ifges(7,6);
t111 = cos(qJ(5));
t108 = sin(qJ(5));
t179 = Ifges(7,5) * t108;
t181 = Ifges(6,4) * t108;
t217 = t219 * t111 + t179 - t181;
t216 = Ifges(7,2) + Ifges(6,3);
t155 = qJD(5) * t111;
t136 = t166 * t109;
t78 = t105 * t112 + t136;
t148 = t78 * t155;
t73 = t118 * qJD(3);
t174 = t108 * t73;
t120 = t148 + t174;
t156 = qJD(5) * t108;
t169 = t111 * t73;
t119 = t78 * t156 - t169;
t215 = -Ifges(6,6) * t111 - t108 * t190;
t214 = m(7) + m(6);
t70 = -pkin(5) * t156 + qJ(6) * t155 + qJD(6) * t108;
t213 = m(7) * t70;
t212 = -mrSges(6,1) - mrSges(7,1);
t72 = t78 * qJD(3);
t211 = t190 * t72 + (-Ifges(6,4) + Ifges(7,5)) * t120 - t219 * t119;
t25 = mrSges(6,1) * t72 + mrSges(6,3) * t119;
t26 = -t72 * mrSges(7,1) - mrSges(7,2) * t119;
t188 = t26 - t25;
t27 = -mrSges(6,2) * t72 - mrSges(6,3) * t120;
t28 = -mrSges(7,2) * t120 + mrSges(7,3) * t72;
t187 = t27 + t28;
t185 = t217 * t78 - t220;
t173 = t108 * t78;
t49 = mrSges(6,2) * t118 - mrSges(6,3) * t173;
t52 = -mrSges(7,2) * t173 - mrSges(7,3) * t118;
t184 = t49 + t52;
t168 = t111 * t78;
t50 = -mrSges(6,1) * t118 - mrSges(6,3) * t168;
t51 = mrSges(7,1) * t118 + mrSges(7,2) * t168;
t183 = -t50 + t51;
t210 = t217 * qJD(5);
t178 = Ifges(7,5) * t111;
t180 = Ifges(6,4) * t111;
t209 = t219 * t108 - t178 + t180;
t207 = t190 * t155 + t218 * t156;
t124 = -t111 * pkin(5) - t108 * qJ(6);
t142 = t166 * pkin(3);
t98 = -t142 - pkin(4);
t76 = t124 + t98;
t89 = -t111 * mrSges(7,1) - t108 * mrSges(7,3);
t206 = -m(7) * t76 - t89;
t205 = m(7) * qJ(6) + mrSges(7,3);
t199 = m(5) * pkin(3);
t204 = t105 * t199 - mrSges(5,2);
t99 = -pkin(3) * t112 - pkin(2);
t48 = -pkin(4) * t118 - pkin(9) * t78 + t99;
t189 = -qJ(4) - pkin(8);
t91 = t189 * t112;
t59 = t164 * t189 - t166 * t91;
t182 = t108 * t48 + t111 * t59;
t137 = qJD(3) * t189;
t116 = -qJD(4) * t109 + t112 * t137;
t71 = qJD(4) * t112 + t109 * t137;
t39 = t105 * t116 + t166 * t71;
t158 = qJD(3) * t109;
t152 = pkin(3) * t158;
t40 = pkin(4) * t72 - pkin(9) * t73 + t152;
t5 = -qJD(5) * t182 - t108 * t39 + t111 * t40;
t90 = -t111 * mrSges(6,1) + t108 * mrSges(6,2);
t203 = m(6) * t98 - t166 * t199 - mrSges(5,1) + t90;
t202 = -2 * Ifges(5,4);
t38 = t105 * t71 - t166 * t116;
t201 = 0.2e1 * t38;
t58 = -t105 * t91 - t189 * t136;
t200 = 0.2e1 * t58;
t93 = Ifges(6,2) * t111 + t181;
t196 = -t93 / 0.2e1;
t194 = pkin(3) * t105;
t193 = t38 * t58;
t106 = sin(pkin(6));
t113 = cos(qJ(2));
t159 = qJD(2) * t113;
t140 = t106 * t159;
t107 = cos(pkin(6));
t110 = sin(qJ(2));
t163 = t106 * t110;
t75 = t107 * t109 + t112 * t163;
t56 = -qJD(3) * t75 - t109 * t140;
t74 = t107 * t112 - t109 * t163;
t57 = qJD(3) * t74 + t112 * t140;
t22 = t105 * t57 - t166 * t56;
t42 = t105 * t75 - t166 * t74;
t11 = t42 * t22;
t192 = t72 * mrSges(5,3);
t191 = t118 * Ifges(6,6);
t125 = Ifges(7,3) * t108 + t178;
t32 = -Ifges(7,6) * t118 + t125 * t78;
t126 = -Ifges(6,2) * t108 + t180;
t33 = t126 * t78 - t191;
t186 = t32 - t33;
t176 = t108 * mrSges(7,2);
t175 = t108 * mrSges(6,3);
t171 = t111 * mrSges(7,2);
t170 = t111 * mrSges(6,3);
t97 = pkin(9) + t194;
t167 = t111 * t97;
t165 = qJD(5) * t78;
t162 = t106 * t113;
t160 = qJD(2) * t110;
t157 = qJD(3) * t112;
t154 = qJD(6) * t111;
t153 = 0.2e1 * t109;
t146 = t108 * t162;
t141 = t106 * t160;
t45 = t72 * mrSges(5,1) + t73 * mrSges(5,2);
t134 = t106 ^ 2 * t110 * t159;
t131 = t58 * t22 + t38 * t42;
t130 = mrSges(6,1) * t108 + mrSges(6,2) * t111;
t129 = mrSges(7,1) * t108 - mrSges(7,3) * t111;
t20 = -t108 * t59 + t111 * t48;
t121 = t120 * Ifges(7,6) + t190 * t169 + t216 * t72;
t43 = t105 * t74 + t166 * t75;
t30 = t108 * t43 + t111 * t162;
t4 = t108 * t40 + t111 * t39 + t48 * t155 - t156 * t59;
t115 = -t56 * t109 + t57 * t112 + (-t109 * t75 - t112 * t74) * qJD(3);
t92 = -Ifges(7,3) * t111 + t179;
t84 = t126 * qJD(5);
t83 = t125 * qJD(5);
t82 = (mrSges(4,1) * t109 + mrSges(4,2) * t112) * qJD(3);
t81 = t130 * qJD(5);
t80 = t129 * qJD(5);
t54 = -mrSges(5,1) * t118 + mrSges(5,2) * t78;
t47 = t130 * t78;
t46 = t129 * t78;
t31 = t111 * t43 - t146;
t24 = (pkin(5) * t108 - qJ(6) * t111) * t78 + t58;
t23 = t105 * t56 + t166 * t57;
t19 = mrSges(6,1) * t120 - mrSges(6,2) * t119;
t18 = mrSges(7,1) * t120 + mrSges(7,3) * t119;
t17 = pkin(5) * t118 - t20;
t16 = -qJ(6) * t118 + t182;
t13 = -Ifges(6,4) * t119 - Ifges(6,2) * t120 + t72 * Ifges(6,6);
t12 = -Ifges(7,5) * t119 + t72 * Ifges(7,6) + Ifges(7,3) * t120;
t10 = -qJD(5) * t146 + t108 * t23 - t111 * t141 + t155 * t43;
t9 = -qJD(5) * t30 + t108 * t141 + t111 * t23;
t7 = (pkin(5) * t73 + qJ(6) * t165) * t108 + (-qJ(6) * t73 + (pkin(5) * qJD(5) - qJD(6)) * t78) * t111 + t38;
t3 = -pkin(5) * t72 - t5;
t1 = qJ(6) * t72 - qJD(6) * t118 + t4;
t2 = [0.2e1 * m(5) * (t43 * t23 + t11 - t134) + 0.2e1 * m(4) * (t74 * t56 + t75 * t57 - t134) + 0.2e1 * t214 * (t10 * t30 + t31 * t9 + t11); t184 * t9 + (t18 + t19) * t42 + t187 * t31 + t188 * t30 + (t46 + t47) * t22 + t183 * t10 + (t118 * t23 + t22 * t78 + t42 * t73 - t43 * t72) * mrSges(5,3) + t115 * mrSges(4,3) + ((-t45 - t82) * t113 + (-t113 * mrSges(3,2) + (-mrSges(4,1) * t112 + mrSges(4,2) * t109 - mrSges(3,1) + t54) * t110) * qJD(2)) * t106 + m(7) * (t1 * t31 + t10 * t17 + t16 * t9 + t22 * t24 + t3 * t30 + t42 * t7) + m(6) * (-t10 * t20 + t182 * t9 - t30 * t5 + t31 * t4 + t131) + m(5) * (t59 * t23 + t39 * t43 + (-t113 * t152 + t160 * t99) * t106 + t131) + (-pkin(2) * t141 + pkin(8) * t115) * m(4); -0.2e1 * t59 * t192 - 0.2e1 * pkin(2) * t82 + 0.2e1 * t1 * t52 + 0.2e1 * t16 * t28 + 0.2e1 * t17 * t26 + 0.2e1 * t24 * t18 + t19 * t200 + 0.2e1 * t20 * t25 + 0.2e1 * t182 * t27 + 0.2e1 * t3 * t51 + t47 * t201 + 0.2e1 * t4 * t49 + 0.2e1 * t99 * t45 + 0.2e1 * t7 * t46 + 0.2e1 * t5 * t50 + (-Ifges(4,4) * t109 + pkin(3) * t54) * qJD(3) * t153 + 0.2e1 * m(7) * (t1 * t16 + t17 * t3 + t24 * t7) + 0.2e1 * m(6) * (t182 * t4 + t20 * t5 + t193) + 0.2e1 * m(5) * (t152 * t99 + t39 * t59 + t193) + (mrSges(5,3) * t200 + t186 * t108 + t185 * t111) * t73 + (0.2e1 * Ifges(4,4) * t112 + (Ifges(4,1) - Ifges(4,2)) * t153) * t157 - (-0.2e1 * t39 * mrSges(5,3) + (-Ifges(6,6) * t108 + t202) * t73 + ((2 * Ifges(5,2)) + t216) * t72 + t121) * t118 + (mrSges(5,3) * t201 + 0.2e1 * Ifges(5,1) * t73 + t211 * t111 + (t12 - t13) * t108 + (t218 * t108 + t190 * t111 + t202) * t72 + ((t186 + t191) * t111 + (-t185 + t220) * t108) * qJD(5)) * t78; t56 * mrSges(4,1) - t57 * mrSges(4,2) + (t171 + t170) * t9 + t214 * (t9 * t167 + (t10 * t108 + (-t108 * t31 + t111 * t30) * qJD(5)) * t97) + (t80 + t81 - t213) * t42 + t204 * t23 + (t175 + t176) * t10 + (t203 - t206) * t22 + (mrSges(6,3) + mrSges(7,2)) * (t30 * t155 - t31 * t156); m(7) * (-t24 * t70 + t7 * t76) - (t209 * t78 + t33) * t156 / 0.2e1 - t207 * t118 / 0.2e1 + (t32 / 0.2e1 - t182 * mrSges(6,3) - t16 * mrSges(7,2)) * t156 + t148 * t196 + t4 * t170 + t1 * t171 + t3 * t176 + (t183 * t155 - t184 * t156 + t188 * t108 + m(6) * (-t108 * t5 + t111 * t4 + (-t108 * t182 - t111 * t20) * qJD(5)) + m(7) * (t1 * t111 + t108 * t3 + (-t108 * t16 + t111 * t17) * qJD(5))) * t97 + Ifges(4,5) * t157 - t192 * t194 - t5 * t175 - Ifges(4,6) * t158 + (-mrSges(4,1) * t157 + mrSges(4,2) * t158) * pkin(8) + (-Ifges(7,6) * t111 - t215) * t72 / 0.2e1 + (-mrSges(5,3) * t142 + Ifges(5,5)) * t73 - t70 * t46 - Ifges(5,6) * t72 + t76 * t18 + t24 * t80 + t58 * t81 + t7 * t89 + t98 * t19 + t203 * t38 - t111 * t12 / 0.2e1 + t111 * t13 / 0.2e1 + t204 * t39 + t209 * t169 / 0.2e1 + t210 * t168 / 0.2e1 + (t78 * t92 + t185) * t155 / 0.2e1 + t187 * t167 + t211 * t108 / 0.2e1 + (t17 * mrSges(7,2) - t20 * mrSges(6,3)) * t155 + (t83 / 0.2e1 - t84 / 0.2e1) * t173 + (t196 + t92 / 0.2e1) * t174; 0.2e1 * t76 * t80 + 0.2e1 * t81 * t98 + 0.2e1 * t206 * t70 + (-t83 + t84) * t111 + t210 * t108 + (t209 * t111 + (t92 - t93) * t108) * qJD(5); m(5) * t141 + t214 * (-t10 * t111 + t108 * t9 + t31 * t155 + t156 * t30); m(5) * t152 - t188 * t111 + t187 * t108 + (t183 * t108 + t184 * t111) * qJD(5) + m(7) * (t1 * t108 - t111 * t3 + (t108 * t17 + t111 * t16) * qJD(5)) + m(6) * (t108 * t4 + t111 * t5 + (-t108 * t20 + t111 * t182) * qJD(5)) + t45; 0; 0; m(7) * qJD(6) * t31 + (-mrSges(6,2) + t205) * t9 + (-m(7) * pkin(5) + t212) * t10; -Ifges(6,6) * t174 + m(7) * (-pkin(5) * t3 + qJ(6) * t1 + qJD(6) * t16) + t1 * mrSges(7,3) + qJD(6) * t52 + qJ(6) * t28 + t5 * mrSges(6,1) - t3 * mrSges(7,1) - t4 * mrSges(6,2) - pkin(5) * t26 + t215 * t165 + t121; (qJD(5) * t124 + t154) * mrSges(7,2) + (m(7) * t154 + (m(7) * t124 + t89 + t90) * qJD(5)) * t97 + t207; t213 + ((-mrSges(6,2) + mrSges(7,3)) * t111 + t212 * t108) * qJD(5); 0.2e1 * t205 * qJD(6); m(7) * t10; m(7) * t3 + t26; (m(7) * t97 + mrSges(7,2)) * t155; m(7) * t156; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t2(1) t2(2) t2(4) t2(7) t2(11) t2(16); t2(2) t2(3) t2(5) t2(8) t2(12) t2(17); t2(4) t2(5) t2(6) t2(9) t2(13) t2(18); t2(7) t2(8) t2(9) t2(10) t2(14) t2(19); t2(11) t2(12) t2(13) t2(14) t2(15) t2(20); t2(16) t2(17) t2(18) t2(19) t2(20) t2(21);];
Mq  = res;
