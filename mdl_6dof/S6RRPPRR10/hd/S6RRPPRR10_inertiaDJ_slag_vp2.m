% Calculate time derivative of joint inertia matrix for
% S6RRPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-03-09 09:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR10_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR10_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR10_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR10_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR10_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:34:11
% EndTime: 2019-03-09 09:34:19
% DurationCPUTime: 3.37s
% Computational Cost: add. (5199->391), mult. (11002->582), div. (0->0), fcn. (9914->8), ass. (0->168)
t160 = sin(pkin(10));
t161 = cos(pkin(10));
t167 = cos(qJ(5));
t193 = qJD(5) * t167;
t164 = sin(qJ(5));
t194 = qJD(5) * t164;
t119 = -t160 * t194 + t161 * t193;
t174 = t167 * t160 + t164 * t161;
t120 = t174 * qJD(5);
t163 = sin(qJ(6));
t166 = cos(qJ(6));
t173 = t160 * t164 - t161 * t167;
t227 = -t163 * t174 - t166 * t173;
t42 = qJD(6) * t227 + t119 * t166 - t120 * t163;
t88 = t163 * t173 - t166 * t174;
t234 = t42 * t88;
t168 = cos(qJ(2));
t109 = t174 * t168;
t165 = sin(qJ(2));
t162 = -pkin(2) - qJ(4);
t199 = t165 * qJ(3);
t129 = t162 * t168 - pkin(1) - t199;
t216 = pkin(3) + pkin(7);
t143 = t216 * t165;
t134 = t161 * t143;
t77 = t165 * pkin(4) + t134 + (pkin(8) * t168 - t129) * t160;
t200 = t161 * t168;
t92 = t161 * t129 + t160 * t143;
t82 = -pkin(8) * t200 + t92;
t35 = -t164 * t82 + t167 * t77;
t26 = pkin(5) * t165 + pkin(9) * t109 + t35;
t108 = t173 * t168;
t36 = t164 * t77 + t167 * t82;
t27 = pkin(9) * t108 + t36;
t12 = -t163 * t27 + t166 * t26;
t13 = t163 * t26 + t166 * t27;
t169 = qJD(6) * t88 - t119 * t163 - t166 * t120;
t203 = t160 * t165;
t195 = qJD(2) * t168;
t138 = t216 * t195;
t196 = qJD(2) * t165;
t183 = pkin(2) * t196 - t165 * qJD(3);
t99 = -qJD(4) * t168 + (-qJ(3) * t168 + qJ(4) * t165) * qJD(2) + t183;
t64 = t161 * t138 - t160 * t99;
t54 = (pkin(4) * t168 - pkin(8) * t203) * qJD(2) + t64;
t187 = t161 * t196;
t65 = t160 * t138 + t161 * t99;
t56 = pkin(8) * t187 + t65;
t15 = -qJD(5) * t36 - t164 * t56 + t167 * t54;
t79 = qJD(5) * t108 + t174 * t196;
t7 = pkin(5) * t195 - t79 * pkin(9) + t15;
t14 = t164 * t54 + t167 * t56 + t77 * t193 - t194 * t82;
t80 = t120 * t168 - t173 * t196;
t8 = pkin(9) * t80 + t14;
t2 = qJD(6) * t12 + t163 * t7 + t166 * t8;
t3 = -qJD(6) * t13 - t163 * t8 + t166 * t7;
t237 = t12 * t169 + t13 * t42 - t2 * t88 + t227 * t3;
t211 = -pkin(8) + t162;
t139 = t211 * t160;
t140 = t211 * t161;
t94 = -t139 * t164 + t167 * t140;
t68 = pkin(9) * t173 + t94;
t95 = t167 * t139 + t164 * t140;
t69 = -pkin(9) * t174 + t95;
t28 = -t163 * t69 + t166 * t68;
t61 = -qJD(4) * t174 - t139 * t194 + t140 * t193;
t49 = -pkin(9) * t119 + t61;
t62 = t173 * qJD(4) - t95 * qJD(5);
t50 = pkin(9) * t120 + t62;
t10 = qJD(6) * t28 + t163 * t50 + t166 * t49;
t29 = t163 * t68 + t166 * t69;
t11 = -qJD(6) * t29 - t163 * t49 + t166 * t50;
t236 = -t10 * t88 + t11 * t227 + t169 * t28 + t29 * t42;
t235 = (t163 * t227 + t166 * t88) * qJD(6) - t163 * t42 - t166 * t169;
t228 = t169 * t227;
t159 = t161 ^ 2;
t184 = (t160 ^ 2 + t159) * qJD(4);
t224 = 2 * m(5);
t223 = 2 * m(6);
t222 = 2 * m(7);
t221 = -0.2e1 * pkin(1);
t179 = -pkin(2) * t168 - t199;
t141 = -pkin(1) + t179;
t219 = -0.2e1 * t141;
t218 = t88 / 0.2e1;
t217 = t227 / 0.2e1;
t215 = -t174 / 0.2e1;
t214 = -t173 / 0.2e1;
t213 = t161 / 0.2e1;
t156 = t168 * pkin(7);
t210 = Ifges(7,5) * t169 - Ifges(7,6) * t42;
t209 = mrSges(5,2) * t160;
t208 = Ifges(5,4) * t160;
t207 = Ifges(5,4) * t161;
t206 = t160 * Ifges(5,1);
t202 = t160 * t168;
t149 = t160 * pkin(4) + qJ(3);
t198 = -Ifges(6,5) * t120 - Ifges(6,6) * t119;
t144 = t168 * pkin(3) + t156;
t192 = qJD(6) * t163;
t191 = qJD(6) * t166;
t190 = 2 * mrSges(6,3);
t72 = t108 * t166 + t109 * t163;
t24 = qJD(6) * t72 + t163 * t80 + t166 * t79;
t73 = t108 * t163 - t109 * t166;
t25 = -qJD(6) * t73 - t163 * t79 + t166 * t80;
t189 = Ifges(7,5) * t24 + Ifges(7,6) * t25 + Ifges(7,3) * t195;
t188 = Ifges(6,5) * t79 + Ifges(6,6) * t80 + Ifges(6,3) * t195;
t118 = pkin(4) * t200 + t144;
t37 = -t80 * mrSges(6,1) + t79 * mrSges(6,2);
t6 = -t25 * mrSges(7,1) + t24 * mrSges(7,2);
t16 = mrSges(7,1) * t42 + mrSges(7,2) * t169;
t186 = mrSges(7,1) * t169 - t42 * mrSges(7,2);
t84 = t119 * mrSges(6,1) - t120 * mrSges(6,2);
t182 = t11 * mrSges(7,1) - t10 * mrSges(7,2) + t210;
t181 = t168 * mrSges(4,2) - t165 * mrSges(4,3);
t180 = -Ifges(5,5) * t160 - Ifges(5,6) * t161;
t178 = t65 * t160 + t64 * t161;
t177 = -t119 * t174 - t120 * t173;
t111 = -mrSges(5,1) * t187 + t196 * t209;
t172 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + t189;
t105 = (-pkin(4) * t161 - t216) * t196;
t171 = t36 * t119 - t35 * t120 + t14 * t174 - t15 * t173;
t170 = t119 * t95 - t120 * t94 - t173 * t62 + t174 * t61;
t142 = mrSges(5,1) * t160 + mrSges(5,2) * t161;
t137 = t216 * t196;
t136 = -t165 * mrSges(5,2) - mrSges(5,3) * t200;
t135 = t165 * mrSges(5,1) + mrSges(5,3) * t202;
t130 = (-mrSges(7,1) * t163 - mrSges(7,2) * t166) * qJD(6) * pkin(5);
t125 = (mrSges(5,3) * t161 * t165 - mrSges(5,2) * t168) * qJD(2);
t124 = (mrSges(5,1) * t168 - mrSges(5,3) * t203) * qJD(2);
t123 = (mrSges(5,1) * t161 - t209) * t168;
t106 = pkin(5) * t119 + qJD(3);
t104 = pkin(5) * t174 + t149;
t103 = (t168 * Ifges(5,5) + (t206 + t207) * t165) * qJD(2);
t102 = (t168 * Ifges(5,6) + (t161 * Ifges(5,2) + t208) * t165) * qJD(2);
t101 = mrSges(6,1) * t165 + mrSges(6,3) * t109;
t100 = -mrSges(6,2) * t165 + mrSges(6,3) * t108;
t98 = -Ifges(6,1) * t173 - Ifges(6,4) * t174;
t97 = -Ifges(6,4) * t173 - Ifges(6,2) * t174;
t96 = mrSges(6,1) * t174 - mrSges(6,2) * t173;
t91 = -t160 * t129 + t134;
t86 = -Ifges(6,1) * t120 - Ifges(6,4) * t119;
t85 = -Ifges(6,4) * t120 - Ifges(6,2) * t119;
t83 = -pkin(5) * t108 + t118;
t81 = -mrSges(6,1) * t108 - mrSges(6,2) * t109;
t71 = -Ifges(6,1) * t109 + Ifges(6,4) * t108 + Ifges(6,5) * t165;
t70 = -Ifges(6,4) * t109 + Ifges(6,2) * t108 + Ifges(6,6) * t165;
t60 = -mrSges(6,2) * t195 + t80 * mrSges(6,3);
t59 = mrSges(6,1) * t195 - t79 * mrSges(6,3);
t58 = mrSges(7,1) * t165 - mrSges(7,3) * t73;
t57 = -mrSges(7,2) * t165 + mrSges(7,3) * t72;
t53 = -pkin(5) * t80 + t105;
t48 = Ifges(7,1) * t227 + Ifges(7,4) * t88;
t47 = Ifges(7,4) * t227 + Ifges(7,2) * t88;
t46 = -mrSges(7,1) * t88 + mrSges(7,2) * t227;
t34 = -mrSges(7,1) * t72 + mrSges(7,2) * t73;
t33 = Ifges(6,1) * t79 + Ifges(6,4) * t80 + Ifges(6,5) * t195;
t32 = Ifges(6,4) * t79 + Ifges(6,2) * t80 + Ifges(6,6) * t195;
t31 = Ifges(7,1) * t73 + Ifges(7,4) * t72 + Ifges(7,5) * t165;
t30 = Ifges(7,4) * t73 + Ifges(7,2) * t72 + Ifges(7,6) * t165;
t20 = -mrSges(7,2) * t195 + t25 * mrSges(7,3);
t19 = mrSges(7,1) * t195 - t24 * mrSges(7,3);
t18 = Ifges(7,1) * t169 - Ifges(7,4) * t42;
t17 = Ifges(7,4) * t169 - Ifges(7,2) * t42;
t5 = Ifges(7,1) * t24 + Ifges(7,4) * t25 + Ifges(7,5) * t195;
t4 = Ifges(7,4) * t24 + Ifges(7,2) * t25 + Ifges(7,6) * t195;
t1 = [((mrSges(3,2) * t221 + mrSges(4,3) * t219 - Ifges(6,5) * t109 + Ifges(7,5) * t73 + Ifges(6,6) * t108 + Ifges(7,6) * t72 + ((2 * Ifges(3,4)) + (2 * Ifges(4,6)) + t180) * t168) * t168 + (mrSges(3,1) * t221 + mrSges(4,2) * t219 + 0.2e1 * (-Ifges(3,4) - Ifges(4,6) - t180) * t165 + (-t159 * Ifges(5,2) + (2 * Ifges(3,1)) - (2 * Ifges(3,2)) + (2 * Ifges(4,2)) - (2 * Ifges(4,3)) + (2 * Ifges(5,3)) + Ifges(6,3) + Ifges(7,3) + (-t206 - 0.2e1 * t207) * t160) * t168) * t165) * qJD(2) + 0.2e1 * (m(4) * t141 + t181) * (-qJ(3) * t195 + t183) + (t12 * t3 + t13 * t2 + t53 * t83) * t222 + (t105 * t118 + t14 * t36 + t15 * t35) * t223 + 0.2e1 * t64 * t135 + 0.2e1 * t65 * t136 - 0.2e1 * t137 * t123 + 0.2e1 * t144 * t111 + 0.2e1 * t91 * t124 + 0.2e1 * t92 * t125 + 0.2e1 * t118 * t37 + 0.2e1 * t105 * t81 + t108 * t32 - t109 * t33 + 0.2e1 * t14 * t100 + 0.2e1 * t15 * t101 + 0.2e1 * t83 * t6 + (-t137 * t144 + t64 * t91 + t65 * t92) * t224 + t79 * t71 + t80 * t70 + t72 * t4 + t73 * t5 + 0.2e1 * t35 * t59 + 0.2e1 * t36 * t60 + 0.2e1 * t2 * t57 + 0.2e1 * t3 * t58 + 0.2e1 * t53 * t34 + t25 * t30 + t24 * t31 + 0.2e1 * t12 * t19 + 0.2e1 * t13 * t20 - t103 * t202 - t102 * t200 + t165 * t188 + t165 * t189; (m(4) * t156 + m(5) * t144 + m(6) * t118 + mrSges(4,1) * t168 + t123 + t81) * qJD(3) + m(7) * (t10 * t13 + t104 * t53 + t106 * t83 + t11 * t12 + t2 * t29 + t28 * t3) + t33 * t214 + t32 * t215 + t5 * t217 + t4 * t218 - t42 * t30 / 0.2e1 + m(6) * (t105 * t149 + t14 * t95 + t15 * t94 + t35 * t62 + t36 * t61) - t237 * mrSges(7,3) + t149 * t37 - t137 * t142 + t118 * t84 - t119 * t70 / 0.2e1 - t120 * t71 / 0.2e1 + t105 * t96 + t106 * t34 + t108 * t85 / 0.2e1 - t109 * t86 / 0.2e1 + qJ(3) * t111 + t94 * t59 + t95 * t60 + t80 * t97 / 0.2e1 + t79 * t98 / 0.2e1 + t61 * t100 + t62 * t101 + t104 * t6 + t83 * t16 + (t210 + t198) * t165 / 0.2e1 + t169 * t31 / 0.2e1 + t72 * t17 / 0.2e1 + t73 * t18 / 0.2e1 + t10 * t57 + t11 * t58 + t25 * t47 / 0.2e1 + t24 * t48 / 0.2e1 + t53 * t46 + t28 * t19 + t29 * t20 + (t103 / 0.2e1 - t64 * mrSges(5,3) - qJD(4) * t135 + t162 * t124) * t161 + (-t102 / 0.2e1 - t65 * mrSges(5,3) - qJD(4) * t136 + t162 * t125) * t160 + ((-pkin(2) * mrSges(4,1) + Ifges(5,5) * t213 - Ifges(5,6) * t160 / 0.2e1 - Ifges(4,4) + Ifges(3,5) + Ifges(7,5) * t217 + Ifges(7,6) * t218 + Ifges(6,5) * t214 + Ifges(6,6) * t215) * t168 + (-qJ(3) * mrSges(4,1) + t160 * (Ifges(5,1) * t161 - t208) / 0.2e1 + (-Ifges(5,2) * t160 + t207) * t213 + Ifges(4,5) - Ifges(3,6)) * t165 + (m(4) * t179 - t168 * mrSges(3,1) + t165 * mrSges(3,2) + t181) * pkin(7)) * qJD(2) - t171 * mrSges(6,3) + m(5) * (-qJ(3) * t137 + t178 * t162 + (-t160 * t92 - t161 * t91) * qJD(4)); 0.2e1 * t104 * t16 + 0.2e1 * t106 * t46 - t119 * t97 - t120 * t98 - t174 * t85 - t173 * t86 + 0.2e1 * t149 * t84 + t88 * t17 + t227 * t18 + t169 * t48 - t42 * t47 + (qJ(3) * qJD(3) - t162 * t184) * t224 + (t10 * t29 + t104 * t106 + t11 * t28) * t222 + (qJD(3) * t149 + t61 * t95 + t62 * t94) * t223 - t170 * t190 + 0.2e1 * mrSges(5,3) * t184 + 0.2e1 * (m(4) * qJ(3) + mrSges(4,3) + t142 + t96) * qJD(3) - 0.2e1 * t236 * mrSges(7,3); t119 * t100 - t120 * t101 + t161 * t124 + t160 * t125 + t174 * t60 - t173 * t59 + t227 * t19 - t88 * t20 + t42 * t57 + t169 * t58 + (m(4) * pkin(7) + mrSges(4,1)) * t195 + m(7) * t237 + m(6) * t171 + m(5) * t178; t177 * t190 + m(7) * t236 + m(6) * t170 - m(5) * t184 + (-0.2e1 * t228 + 0.2e1 * t234) * mrSges(7,3); -0.2e1 * m(6) * t177 + 0.2e1 * m(7) * (t228 - t234); -m(5) * t137 + m(6) * t105 + m(7) * t53 + t111 + t37 + t6; m(7) * t106 + (m(6) + m(5)) * qJD(3) + t84 + t16; 0; 0; t15 * mrSges(6,1) - t14 * mrSges(6,2) + (m(7) * (-t12 * t192 + t13 * t191 + t163 * t2 + t166 * t3) + t57 * t191 + t163 * t20 - t58 * t192 + t166 * t19) * pkin(5) + t172 + t188; t62 * mrSges(6,1) - t61 * mrSges(6,2) + (m(7) * (t10 * t163 + t11 * t166 + (-t163 * t28 + t166 * t29) * qJD(6)) + t235 * mrSges(7,3)) * pkin(5) + t182 + t198; -m(7) * t235 * pkin(5) - t120 * mrSges(6,1) - t119 * mrSges(6,2) + t186; 0; 0.2e1 * t130; t172; t182; t186; 0; t130; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
