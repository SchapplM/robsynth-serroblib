% Calculate time derivative of joint inertia matrix for
% S5RRPRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR14_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR14_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR14_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR14_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR14_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR14_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR14_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:35:49
% EndTime: 2019-12-31 20:35:59
% DurationCPUTime: 3.49s
% Computational Cost: add. (4703->421), mult. (12444->649), div. (0->0), fcn. (12224->10), ass. (0->182)
t151 = sin(pkin(5));
t225 = 0.2e1 * t151;
t150 = sin(pkin(10));
t152 = cos(pkin(10));
t155 = sin(qJ(4));
t158 = cos(qJ(4));
t125 = t150 * t155 - t158 * t152;
t120 = t125 * qJD(4);
t126 = t150 * t158 + t152 * t155;
t157 = cos(qJ(5));
t154 = sin(qJ(5));
t182 = qJD(5) * t154;
t160 = t157 * t120 + t126 * t182;
t224 = t154 / 0.2e1;
t204 = t157 / 0.2e1;
t153 = cos(pkin(5));
t156 = sin(qJ(2));
t191 = t151 * t156;
t119 = t150 * t153 + t152 * t191;
t159 = cos(qJ(2));
t190 = t151 * t159;
t123 = t153 * t156 * pkin(1) + pkin(7) * t190;
t108 = qJ(3) * t153 + t123;
t109 = (-pkin(2) * t159 - qJ(3) * t156 - pkin(1)) * t151;
t72 = -t150 * t108 + t152 * t109;
t48 = -pkin(3) * t190 - t119 * pkin(8) + t72;
t118 = -t150 * t191 + t152 * t153;
t73 = t152 * t108 + t150 * t109;
t60 = pkin(8) * t118 + t73;
t200 = t155 * t48 + t158 * t60;
t185 = qJD(2) * t151;
t189 = t152 * t159;
t176 = t156 * t185;
t203 = pkin(1) * t159;
t178 = t153 * t203;
t112 = -pkin(7) * t176 + qJD(2) * t178;
t104 = qJD(3) * t153 + t112;
t99 = (-qJD(3) * t156 + (pkin(2) * t156 - qJ(3) * t159) * qJD(2)) * t151;
t63 = -t150 * t104 + t152 * t99;
t46 = (pkin(3) * t156 - pkin(8) * t189) * t185 + t63;
t175 = t159 * t185;
t171 = t150 * t175;
t64 = t152 * t104 + t150 * t99;
t59 = -pkin(8) * t171 + t64;
t9 = -qJD(4) * t200 - t155 * t59 + t158 * t46;
t223 = 2 * m(4);
t222 = 2 * m(5);
t221 = 2 * m(6);
t220 = -2 * mrSges(3,3);
t219 = -2 * mrSges(5,3);
t201 = pkin(8) + qJ(3);
t134 = t201 * t152;
t172 = qJD(4) * t201;
t179 = t158 * qJD(3);
t180 = t155 * qJD(3);
t183 = qJD(4) * t158;
t75 = t152 * t180 + t134 * t183 + (-t155 * t172 + t179) * t150;
t218 = 0.2e1 * t75;
t173 = t201 * t150;
t94 = t134 * t155 + t158 * t173;
t217 = 0.2e1 * t94;
t77 = t118 * t155 + t119 * t158;
t162 = t154 * t190 - t157 * t77;
t164 = t158 * t118 - t119 * t155;
t56 = qJD(4) * t164 - t125 * t175;
t31 = qJD(5) * t162 - t154 * t56 + t157 * t176;
t216 = t31 / 0.2e1;
t121 = t126 * qJD(4);
t181 = qJD(5) * t157;
t161 = -t154 * t120 + t126 * t181;
t38 = -Ifges(6,1) * t160 - Ifges(6,4) * t161 + Ifges(6,5) * t121;
t215 = t38 / 0.2e1;
t65 = -t154 * t77 - t157 * t190;
t214 = t65 / 0.2e1;
t196 = Ifges(6,4) * t154;
t168 = Ifges(6,1) * t157 - t196;
t69 = Ifges(6,5) * t125 + t126 * t168;
t213 = t69 / 0.2e1;
t212 = -t126 / 0.2e1;
t147 = Ifges(6,5) * t181;
t211 = -Ifges(6,6) * t182 / 0.2e1 + t147 / 0.2e1;
t130 = t168 * qJD(5);
t210 = t130 / 0.2e1;
t209 = Ifges(6,5) * t224 + Ifges(6,6) * t204;
t137 = Ifges(6,2) * t157 + t196;
t208 = -t137 / 0.2e1;
t195 = Ifges(6,4) * t157;
t138 = Ifges(6,1) * t154 + t195;
t207 = t138 / 0.2e1;
t206 = t152 / 0.2e1;
t205 = -t154 / 0.2e1;
t202 = t75 * t94;
t199 = mrSges(6,3) * t126;
t198 = Ifges(4,4) * t150;
t197 = Ifges(4,4) * t152;
t194 = Ifges(6,6) * t154;
t146 = -pkin(3) * t152 - pkin(2);
t84 = pkin(4) * t125 - pkin(9) * t126 + t146;
t95 = t158 * t134 - t155 * t173;
t54 = -t154 * t95 + t157 * t84;
t193 = qJD(5) * t54;
t55 = t154 * t84 + t157 * t95;
t192 = qJD(5) * t55;
t187 = -Ifges(5,5) * t120 - Ifges(5,6) * t121;
t100 = t152 * mrSges(4,2) * t175 + mrSges(4,1) * t171;
t184 = qJD(4) * t155;
t30 = qJD(5) * t65 + t154 * t176 + t157 * t56;
t57 = qJD(4) * t77 + t126 * t175;
t3 = Ifges(6,5) * t30 + Ifges(6,6) * t31 + Ifges(6,3) * t57;
t177 = Ifges(5,5) * t56 - Ifges(5,6) * t57 + Ifges(5,3) * t176;
t27 = t57 * mrSges(5,1) + t56 * mrSges(5,2);
t79 = t121 * mrSges(5,1) - t120 * mrSges(5,2);
t19 = -pkin(9) * t190 + t200;
t142 = pkin(7) * t191;
t111 = t142 + (-pkin(2) - t203) * t153;
t78 = -t118 * pkin(3) + t111;
t32 = -pkin(4) * t164 - t77 * pkin(9) + t78;
t11 = -t154 * t19 + t157 * t32;
t113 = t123 * qJD(2);
t93 = pkin(3) * t171 + t113;
t15 = t57 * pkin(4) - t56 * pkin(9) + t93;
t8 = t155 * t46 + t158 * t59 + t48 * t183 - t184 * t60;
t6 = pkin(9) * t176 + t8;
t1 = qJD(5) * t11 + t15 * t154 + t157 * t6;
t12 = t154 * t32 + t157 * t19;
t2 = -qJD(5) * t12 + t15 * t157 - t154 * t6;
t170 = t1 * t157 - t2 * t154;
t135 = -mrSges(6,1) * t157 + mrSges(6,2) * t154;
t169 = mrSges(6,1) * t154 + mrSges(6,2) * t157;
t167 = -Ifges(6,2) * t154 + t195;
t25 = -t155 * t60 + t158 * t48;
t21 = -Ifges(6,4) * t162 + Ifges(6,2) * t65 - Ifges(6,6) * t164;
t22 = -Ifges(6,1) * t162 + Ifges(6,4) * t65 - Ifges(6,5) * t164;
t163 = t204 * t22 + t205 * t21;
t36 = -Ifges(6,5) * t160 - Ifges(6,6) * t161 + Ifges(6,3) * t121;
t140 = Ifges(3,5) * t175;
t129 = t167 * qJD(5);
t127 = t169 * qJD(5);
t122 = -t142 + t178;
t106 = (mrSges(4,1) * t156 - mrSges(4,3) * t189) * t185;
t105 = (-mrSges(4,3) * t150 * t159 - mrSges(4,2) * t156) * t185;
t98 = -mrSges(4,1) * t190 - t119 * mrSges(4,3);
t97 = mrSges(4,2) * t190 + t118 * mrSges(4,3);
t90 = Ifges(5,1) * t126 - Ifges(5,4) * t125;
t89 = Ifges(5,4) * t126 - Ifges(5,2) * t125;
t88 = (t156 * Ifges(4,5) + (t152 * Ifges(4,1) - t198) * t159) * t185;
t87 = (t156 * Ifges(4,6) + (-t150 * Ifges(4,2) + t197) * t159) * t185;
t86 = mrSges(6,1) * t125 - t157 * t199;
t85 = -mrSges(6,2) * t125 - t154 * t199;
t83 = pkin(4) * t121 + pkin(9) * t120;
t82 = t169 * t126;
t81 = -Ifges(5,1) * t120 - Ifges(5,4) * t121;
t80 = -Ifges(5,4) * t120 - Ifges(5,2) * t121;
t74 = t152 * t179 - t134 * t184 + (-t158 * t172 - t180) * t150;
t71 = -mrSges(5,1) * t190 - t77 * mrSges(5,3);
t70 = mrSges(5,2) * t190 + mrSges(5,3) * t164;
t68 = Ifges(6,6) * t125 + t126 * t167;
t67 = Ifges(6,3) * t125 + (Ifges(6,5) * t157 - t194) * t126;
t62 = -mrSges(6,2) * t121 - mrSges(6,3) * t161;
t61 = mrSges(6,1) * t121 + mrSges(6,3) * t160;
t44 = mrSges(6,1) * t161 - mrSges(6,2) * t160;
t42 = -mrSges(5,2) * t176 - mrSges(5,3) * t57;
t41 = mrSges(5,1) * t176 - mrSges(5,3) * t56;
t40 = Ifges(5,1) * t77 + Ifges(5,4) * t164 - Ifges(5,5) * t190;
t39 = Ifges(5,4) * t77 + Ifges(5,2) * t164 - Ifges(5,6) * t190;
t37 = -Ifges(6,4) * t160 - Ifges(6,2) * t161 + Ifges(6,6) * t121;
t35 = -mrSges(6,1) * t164 + mrSges(6,3) * t162;
t34 = mrSges(6,2) * t164 + mrSges(6,3) * t65;
t33 = -mrSges(6,1) * t65 - mrSges(6,2) * t162;
t24 = Ifges(5,1) * t56 - Ifges(5,4) * t57 + Ifges(5,5) * t176;
t23 = Ifges(5,4) * t56 - Ifges(5,2) * t57 + Ifges(5,6) * t176;
t20 = -Ifges(6,5) * t162 + Ifges(6,6) * t65 - Ifges(6,3) * t164;
t18 = pkin(4) * t190 - t25;
t17 = -t154 * t74 + t157 * t83 - t192;
t16 = t154 * t83 + t157 * t74 + t193;
t14 = -mrSges(6,2) * t57 + mrSges(6,3) * t31;
t13 = mrSges(6,1) * t57 - mrSges(6,3) * t30;
t10 = -mrSges(6,1) * t31 + mrSges(6,2) * t30;
t7 = -pkin(4) * t176 - t9;
t5 = Ifges(6,1) * t30 + Ifges(6,4) * t31 + Ifges(6,5) * t57;
t4 = Ifges(6,4) * t30 + Ifges(6,2) * t31 + Ifges(6,6) * t57;
t26 = [((t123 * t220 + Ifges(4,5) * t119 + Ifges(5,5) * t77 - 0.2e1 * Ifges(3,6) * t153 + Ifges(4,6) * t118 + Ifges(5,6) * t164 + (-mrSges(3,1) * pkin(1) - Ifges(3,4) * t156) * t225) * t156 + (t122 * t220 - t150 * (Ifges(4,4) * t119 + Ifges(4,2) * t118) + t152 * (Ifges(4,1) * t119 + Ifges(4,4) * t118) + Ifges(3,5) * t153 + (-pkin(1) * mrSges(3,2) + (-Ifges(4,5) * t152 + Ifges(4,6) * t150 + Ifges(3,4)) * t159) * t225 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - (2 * Ifges(4,3)) - Ifges(5,3)) * t191) * t159) * t185 + 0.2e1 * m(3) * (t112 * t123 - t113 * t122) + (t1 * t12 + t11 * t2 + t18 * t7) * t221 + (t111 * t113 + t63 * t72 + t64 * t73) * t223 + (t20 - t39) * t57 + t118 * t87 + 0.2e1 * t113 * (-mrSges(4,1) * t118 + mrSges(4,2) * t119) + t119 * t88 + 0.2e1 * t111 * t100 + 0.2e1 * t64 * t97 + 0.2e1 * t63 * t98 + 0.2e1 * t73 * t105 + 0.2e1 * t72 * t106 + t77 * t24 + 0.2e1 * t78 * t27 + 0.2e1 * t8 * t70 + 0.2e1 * t9 * t71 + t65 * t4 + t56 * t40 + 0.2e1 * t7 * t33 + 0.2e1 * t1 * t34 + 0.2e1 * t2 * t35 + 0.2e1 * t25 * t41 + t30 * t22 + t31 * t21 + 0.2e1 * t18 * t10 + 0.2e1 * t11 * t13 + 0.2e1 * t12 * t14 - t177 * t190 + 0.2e1 * t112 * (-t153 * mrSges(3,2) + mrSges(3,3) * t190) - 0.2e1 * t113 * (mrSges(3,1) * t153 - mrSges(3,3) * t191) - t162 * t5 - t164 * t3 + t164 * t23 + 0.2e1 * t93 * (-mrSges(5,1) * t164 + mrSges(5,2) * t77) + (t200 * t8 + t25 * t9 + t78 * t93) * t222 + 0.2e1 * t200 * t42 + t153 * t140; (-t89 / 0.2e1 + t67 / 0.2e1) * t57 + t140 - (t40 / 0.2e1 - t25 * mrSges(5,3) + t163) * t120 + t68 * t216 + t30 * t213 + t37 * t214 + (t10 - t41) * t94 + (t33 - t71) * t75 + (-t113 * mrSges(4,1) + t87 / 0.2e1 + t64 * mrSges(4,3) + qJD(3) * t97 + qJ(3) * t105) * t152 + (t88 / 0.2e1 + t113 * mrSges(4,2) - t63 * mrSges(4,3) - qJD(3) * t98 - qJ(3) * t106) * t150 + t146 * t27 + (t3 / 0.2e1 - t23 / 0.2e1 + t93 * mrSges(5,1) - t8 * mrSges(5,3)) * t125 - t112 * mrSges(3,2) - t113 * mrSges(3,1) + t95 * t42 - pkin(2) * t100 + t77 * t81 / 0.2e1 + t7 * t82 + t1 * t85 + t2 * t86 + t56 * t90 / 0.2e1 + t78 * t79 + t74 * t70 + t11 * t61 + t12 * t62 + t54 * t13 + t55 * t14 + t18 * t44 + t16 * t34 + t17 * t35 + (t93 * mrSges(5,2) + t24 / 0.2e1 - t9 * mrSges(5,3) + t4 * t205 + t5 * t204 + (t22 * t205 - t157 * t21 / 0.2e1) * qJD(5)) * t126 + (-t159 * t187 / 0.2e1 + ((-Ifges(3,6) + Ifges(5,5) * t126 / 0.2e1 - Ifges(5,6) * t125 / 0.2e1 + Ifges(4,5) * t150 / 0.2e1 + Ifges(4,6) * t206) * t156 + (-t150 * (Ifges(4,2) * t152 + t198) / 0.2e1 + (Ifges(4,1) * t150 + t197) * t206) * t159) * qJD(2)) * t151 - t162 * t215 - (t36 / 0.2e1 - t80 / 0.2e1) * t164 + m(6) * (t1 * t55 + t11 * t17 + t12 * t16 + t18 * t75 + t2 * t54 + t7 * t94) + m(4) * (-pkin(2) * t113 + (-t150 * t72 + t152 * t73) * qJD(3) + (-t150 * t63 + t152 * t64) * qJ(3)) + (t20 / 0.2e1 - t39 / 0.2e1 - t200 * mrSges(5,3)) * t121 + m(5) * (t146 * t93 + t200 * t74 - t25 * t75 + t8 * t95 - t9 * t94); 0.2e1 * t146 * t79 + 0.2e1 * t16 * t85 + 0.2e1 * t17 * t86 + t44 * t217 + 0.2e1 * t54 * t61 + 0.2e1 * t55 * t62 + t82 * t218 + (t16 * t55 + t17 * t54 + t202) * t221 + (t74 * t95 + t202) * t222 + (t219 * t74 + t36 - t80) * t125 + (t219 * t95 + t67 - t89) * t121 - (mrSges(5,3) * t217 - t154 * t68 + t157 * t69 + t90) * t120 + (mrSges(5,3) * t218 - t154 * t37 + t157 * t38 + t81 + (-t154 * t69 - t157 * t68) * qJD(5)) * t126 + (qJ(3) * t223 + 0.2e1 * mrSges(4,3)) * qJD(3) * (t150 ^ 2 + t152 ^ 2); t157 * t13 + t154 * t14 + (-t154 * t35 + t157 * t34) * qJD(5) + m(6) * (t1 * t154 + t157 * t2 + (-t11 * t154 + t12 * t157) * qJD(5)) + m(5) * t93 + m(4) * t113 + t27 + t100; m(6) * (t154 * t16 + t157 * t17 + (-t154 * t54 + t157 * t55) * qJD(5)) + t85 * t181 + t154 * t62 - t86 * t182 + t157 * t61 + t79; 0; t4 * t204 + t5 * t224 + t7 * t135 + t57 * t209 + t137 * t216 + t30 * t207 + t18 * t127 - t164 * t211 + t129 * t214 - t162 * t210 - t8 * mrSges(5,2) + t9 * mrSges(5,1) + t163 * qJD(5) + (-m(6) * t7 - t10) * pkin(4) + ((-t11 * t157 - t12 * t154) * qJD(5) + t170) * mrSges(6,3) + (-t34 * t182 - t35 * t181 - t154 * t13 + t157 * t14 + m(6) * (-t11 * t181 - t12 * t182 + t170)) * pkin(9) + t177; t121 * t209 + t94 * t127 + t125 * t211 - t74 * mrSges(5,2) - pkin(4) * t44 + (-m(6) * pkin(4) - mrSges(5,1) + t135) * t75 + (t37 / 0.2e1 + t16 * mrSges(6,3) + t126 * t210 - t120 * t207 + (-t54 * mrSges(6,3) + t126 * t208 + t213) * qJD(5) + (m(6) * (t16 - t193) + t62 - qJD(5) * t86) * pkin(9)) * t157 + (t215 - t17 * mrSges(6,3) + t129 * t212 - t120 * t208 + (-t55 * mrSges(6,3) + t138 * t212 - t68 / 0.2e1) * qJD(5) + (-qJD(5) * t85 - t61 + m(6) * (-t17 - t192)) * pkin(9)) * t154 + t187; 0; -0.2e1 * pkin(4) * t127 + t129 * t157 + t130 * t154 + (-t137 * t154 + t138 * t157) * qJD(5); mrSges(6,1) * t2 - mrSges(6,2) * t1 + t3; mrSges(6,1) * t17 - mrSges(6,2) * t16 + t36; -t127; t147 + (pkin(9) * t135 - t194) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t26(1), t26(2), t26(4), t26(7), t26(11); t26(2), t26(3), t26(5), t26(8), t26(12); t26(4), t26(5), t26(6), t26(9), t26(13); t26(7), t26(8), t26(9), t26(10), t26(14); t26(11), t26(12), t26(13), t26(14), t26(15);];
Mq = res;
