% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 23:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRPRPR6_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR6_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR6_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR6_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:10:53
% EndTime: 2019-05-04 23:11:00
% DurationCPUTime: 2.70s
% Computational Cost: add. (10533->316), mult. (21326->460), div. (0->0), fcn. (14798->12), ass. (0->194)
t164 = sin(pkin(11));
t166 = cos(pkin(11));
t172 = cos(qJ(4));
t200 = qJD(2) * t172;
t132 = -t166 * qJD(4) + t164 * t200;
t134 = t164 * qJD(4) + t166 * t200;
t111 = t134 * t132;
t197 = qJD(2) * qJD(4);
t189 = t172 * t197;
t169 = sin(qJ(4));
t196 = t169 * qJDD(2);
t139 = t189 + t196;
t223 = -t111 + t139;
t230 = t164 * t223;
t229 = t166 * t223;
t168 = sin(qJ(6));
t135 = qJDD(6) + t139;
t171 = cos(qJ(6));
t107 = t171 * t132 + t168 * t134;
t109 = -t168 * t132 + t171 * t134;
t83 = t109 * t107;
t224 = t135 - t83;
t228 = t168 * t224;
t227 = t171 * t224;
t165 = sin(pkin(6));
t167 = cos(pkin(6));
t210 = sin(pkin(10));
t211 = cos(pkin(10));
t178 = t210 * g(1) - t211 * g(2);
t201 = -g(3) + qJDD(1);
t226 = t165 * t201 + t167 * t178;
t154 = t172 * qJDD(2);
t190 = t169 * t197;
t140 = t154 - t190;
t117 = t164 * qJDD(4) + t166 * t140;
t186 = -t166 * qJDD(4) + t164 * t140;
t73 = -t107 * qJD(6) + t171 * t117 - t168 * t186;
t199 = t169 * qJD(2);
t151 = qJD(6) + t199;
t95 = t151 * t107;
t225 = -t95 + t73;
t187 = t168 * t117 + t171 * t186;
t54 = (qJD(6) - t151) * t109 + t187;
t105 = t107 ^ 2;
t106 = t109 ^ 2;
t222 = t132 ^ 2;
t131 = t134 ^ 2;
t150 = t151 ^ 2;
t221 = qJD(4) ^ 2;
t220 = -pkin(8) - pkin(2);
t184 = t169 * pkin(4) - t172 * qJ(5);
t137 = t184 * qJD(2);
t163 = qJDD(2) * pkin(2);
t174 = qJD(2) ^ 2;
t170 = sin(qJ(2));
t173 = cos(qJ(2));
t179 = -t211 * g(1) - t210 * g(2);
t103 = -t170 * t179 + t226 * t173;
t176 = qJDD(3) - t103;
t90 = -t174 * qJ(3) - t163 + t176;
t175 = -qJDD(2) * pkin(8) + t90;
t119 = -t165 * t178 + t167 * t201;
t203 = t172 * t119;
t64 = -t221 * pkin(4) + qJDD(4) * qJ(5) + t203 + (-qJD(2) * t137 + t175) * t169;
t182 = -t140 + t190;
t183 = t139 + t189;
t157 = qJDD(2) * qJ(3);
t104 = t226 * t170 + t173 * t179;
t185 = 0.2e1 * qJD(3) * qJD(2) + t104;
t181 = t157 + t185;
t87 = t220 * t174 + t181;
t67 = t183 * pkin(4) + t182 * qJ(5) + t87;
t36 = 0.2e1 * qJD(5) * t134 + t164 * t64 - t166 * t67;
t192 = t132 * t199;
t99 = -t117 - t192;
t27 = pkin(5) * t223 + t99 * pkin(9) - t36;
t118 = pkin(5) * t199 - t134 * pkin(9);
t37 = -0.2e1 * qJD(5) * t132 + t164 * t67 + t166 * t64;
t30 = -t222 * pkin(5) - t186 * pkin(9) - t118 * t199 + t37;
t11 = t168 * t30 - t171 * t27;
t12 = t168 * t27 + t171 * t30;
t7 = -t171 * t11 + t168 * t12;
t219 = t164 * t7;
t218 = t166 * t7;
t79 = t169 * t119 - t172 * t175;
t63 = -qJDD(4) * pkin(4) - t221 * qJ(5) + t137 * t200 + qJDD(5) + t79;
t217 = t164 * t63;
t216 = t166 * t63;
t38 = t186 * pkin(5) - t222 * pkin(9) + t134 * t118 + t63;
t215 = t168 * t38;
t75 = t135 + t83;
t214 = t168 * t75;
t213 = t171 * t38;
t212 = t171 * t75;
t209 = t151 * t168;
t208 = t151 * t171;
t161 = t172 ^ 2;
t207 = t161 * t174;
t101 = t111 + t139;
t206 = t164 * t101;
t205 = t166 * t101;
t193 = t172 * t174 * t169;
t144 = qJDD(4) + t193;
t204 = t169 * t144;
t145 = qJDD(4) - t193;
t202 = t172 * t145;
t195 = t169 * t83;
t194 = t169 * t111;
t191 = t134 * t199;
t8 = t168 * t11 + t171 * t12;
t20 = t164 * t36 + t166 * t37;
t19 = t164 * t37 - t166 * t36;
t80 = t169 * t175 + t203;
t41 = t169 * t80 - t172 * t79;
t180 = qJ(3) + t184;
t97 = -t186 + t191;
t160 = t169 ^ 2;
t155 = t160 * t174;
t149 = -t207 - t221;
t148 = -t155 - t221;
t143 = -t155 - t207;
t142 = (t160 + t161) * qJDD(2);
t141 = t154 - 0.2e1 * t190;
t138 = 0.2e1 * t189 + t196;
t125 = (-qJDD(2) * t173 + t170 * t174) * t165;
t124 = (qJDD(2) * t170 + t173 * t174) * t165;
t122 = -t131 - t155;
t121 = -t131 + t155;
t120 = -t155 + t222;
t115 = t172 * t149 - t204;
t114 = t169 * t148 + t202;
t113 = t167 * t119;
t110 = -t155 - t222;
t98 = t117 - t192;
t96 = t186 + t191;
t93 = -t131 - t222;
t92 = -t106 + t150;
t91 = t105 - t150;
t89 = -t174 * pkin(2) + t181;
t88 = -t106 - t150;
t85 = -t164 * t122 - t205;
t84 = t166 * t122 - t206;
t82 = t106 - t105;
t81 = -t150 - t105;
t78 = t166 * t110 - t230;
t77 = t164 * t110 + t229;
t72 = -t109 * qJD(6) - t187;
t71 = -t164 * t99 + t166 * t97;
t70 = t164 * t97 + t166 * t99;
t69 = (-t107 * t171 + t109 * t168) * t151;
t68 = (-t107 * t168 - t109 * t171) * t151;
t62 = -t105 - t106;
t60 = t169 * t85 - t172 * t98;
t59 = t169 * t78 - t172 * t96;
t58 = t95 + t73;
t53 = (qJD(6) + t151) * t109 + t187;
t52 = t171 * t91 - t214;
t51 = -t168 * t92 + t227;
t50 = t168 * t91 + t212;
t49 = t171 * t92 + t228;
t48 = -t109 * t209 + t171 * t73;
t47 = t109 * t208 + t168 * t73;
t46 = t107 * t208 - t168 * t72;
t45 = t107 * t209 + t171 * t72;
t44 = t169 * t71 - t172 * t93;
t43 = -t168 * t88 - t212;
t42 = t171 * t88 - t214;
t40 = t171 * t81 - t228;
t39 = t168 * t81 + t227;
t34 = -t168 * t225 - t171 * t53;
t33 = t168 * t58 - t171 * t54;
t32 = -t168 * t53 + t171 * t225;
t31 = -t168 * t54 - t171 * t58;
t29 = -t164 * t42 + t166 * t43;
t28 = t164 * t43 + t166 * t42;
t25 = -pkin(9) * t42 + t213;
t24 = -t164 * t39 + t166 * t40;
t23 = t164 * t40 + t166 * t39;
t22 = -pkin(9) * t39 + t215;
t21 = t169 * t29 - t172 * t225;
t18 = -pkin(5) * t225 + pkin(9) * t43 + t215;
t17 = t169 * t24 - t172 * t53;
t16 = -pkin(5) * t53 + pkin(9) * t40 - t213;
t15 = -t164 * t31 + t166 * t33;
t14 = t164 * t33 + t166 * t31;
t13 = t169 * t20 - t172 * t63;
t9 = t169 * t15 - t172 * t62;
t6 = -pkin(5) * t38 + pkin(9) * t8;
t5 = -pkin(9) * t31 - t7;
t4 = -pkin(5) * t62 + pkin(9) * t33 + t8;
t3 = t166 * t8 - t219;
t2 = t164 * t8 + t218;
t1 = t169 * t3 - t172 * t38;
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t201, 0, 0, 0, 0, 0, 0, -t125, -t124, 0, t113 + (t103 * t173 + t104 * t170) * t165, 0, 0, 0, 0, 0, 0, 0, t125, t124, t113 + (t170 * t89 - t173 * t90) * t165, 0, 0, 0, 0, 0, 0, t167 * (-t169 * t145 + t172 * t148) + (-t173 * t114 + t170 * t138) * t165, t167 * (-t172 * t144 - t169 * t149) + (-t173 * t115 + t170 * t141) * t165, (t142 * t173 + t143 * t170) * t165, t167 * (t169 * t79 + t172 * t80) + (t170 * t87 - t173 * t41) * t165, 0, 0, 0, 0, 0, 0, t167 * (t169 * t96 + t172 * t78) + (t170 * t77 - t173 * t59) * t165, t167 * (t169 * t98 + t172 * t85) + (t170 * t84 - t173 * t60) * t165, t167 * (t169 * t93 + t172 * t71) + (t170 * t70 - t173 * t44) * t165, t167 * (t169 * t63 + t172 * t20) + (-t173 * t13 + t170 * t19) * t165, 0, 0, 0, 0, 0, 0, t167 * (t169 * t53 + t172 * t24) + (-t173 * t17 + t170 * t23) * t165, t167 * (t169 * t225 + t172 * t29) + (t170 * t28 - t173 * t21) * t165, t167 * (t172 * t15 + t169 * t62) + (t170 * t14 - t173 * t9) * t165, t167 * (t169 * t38 + t172 * t3) + (-t173 * t1 + t170 * t2) * t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t103, -t104, 0, 0, qJDD(2), 0, 0, 0, 0, 0, 0, -0.2e1 * t163 + t176, 0.2e1 * t157 + t185, -pkin(2) * t90 + qJ(3) * t89, -t182 * t172, -t172 * t138 - t169 * t141, t202 - t169 * (-t207 + t221), t183 * t169, t172 * (t155 - t221) - t204, 0, qJ(3) * t138 + t220 * t114 + t169 * t87, qJ(3) * t141 + t220 * t115 + t172 * t87, qJ(3) * t143 - t220 * t142 - t41, qJ(3) * t87 + t220 * t41, t172 * (t166 * t117 - t164 * t191) + t194, t172 * (-t164 * t98 - t166 * t96) - t169 * (-t131 + t222), t172 * (-t164 * t121 + t229) - t169 * t99, t172 * (t164 * t186 + t166 * t192) - t194, t172 * (t166 * t120 - t206) + t169 * t97, (t139 + (-t132 * t166 + t134 * t164) * t200) * t169, t172 * (-qJ(5) * t77 + t217) - t169 * (-pkin(4) * t77 + t36) + qJ(3) * t77 + t220 * t59, t172 * (-qJ(5) * t84 + t216) - t169 * (-pkin(4) * t84 + t37) + qJ(3) * t84 + t220 * t60, -t172 * t19 + t180 * t70 + t220 * t44, t220 * t13 + t180 * t19, t172 * (-t164 * t47 + t166 * t48) + t195, t172 * (-t164 * t32 + t166 * t34) + t169 * t82, t172 * (-t164 * t49 + t166 * t51) + t169 * t58, t172 * (-t164 * t45 + t166 * t46) - t195, t172 * (-t164 * t50 + t166 * t52) - t169 * t54, t172 * (-t164 * t68 + t166 * t69) + t169 * t135, t172 * (-qJ(5) * t23 - t164 * t16 + t166 * t22) - t169 * (-pkin(4) * t23 - pkin(5) * t39 + t11) + qJ(3) * t23 + t220 * t17, t172 * (-qJ(5) * t28 - t164 * t18 + t166 * t25) - t169 * (-pkin(4) * t28 - pkin(5) * t42 + t12) + qJ(3) * t28 + t220 * t21, t172 * (-qJ(5) * t14 - t164 * t4 + t166 * t5) - t169 * (-pkin(4) * t14 - pkin(5) * t31) + qJ(3) * t14 + t220 * t9, t172 * (-pkin(9) * t218 - qJ(5) * t2 - t164 * t6) - t169 * (-pkin(4) * t2 - pkin(5) * t7) + qJ(3) * t2 + t220 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t174, t90, 0, 0, 0, 0, 0, 0, t114, t115, -t142, t41, 0, 0, 0, 0, 0, 0, t59, t60, t44, t13, 0, 0, 0, 0, 0, 0, t17, t21, t9, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193, -t155 + t207, t154, -t193, -t196, qJDD(4), -t79, -t80, 0, 0, t164 * t117 + t166 * t191, -t164 * t96 + t166 * t98, t166 * t121 + t230, t164 * t192 - t166 * t186, t164 * t120 + t205, (-t132 * t164 - t134 * t166) * t199, -pkin(4) * t96 + qJ(5) * t78 - t216, -pkin(4) * t98 + qJ(5) * t85 + t217, -pkin(4) * t93 + qJ(5) * t71 + t20, -pkin(4) * t63 + qJ(5) * t20, t164 * t48 + t166 * t47, t164 * t34 + t166 * t32, t164 * t51 + t166 * t49, t164 * t46 + t166 * t45, t164 * t52 + t166 * t50, t164 * t69 + t166 * t68, -pkin(4) * t53 + qJ(5) * t24 + t166 * t16 + t164 * t22, -pkin(4) * t225 + qJ(5) * t29 + t164 * t25 + t166 * t18, -pkin(4) * t62 + qJ(5) * t15 + t164 * t5 + t166 * t4, -pkin(4) * t38 - pkin(9) * t219 + qJ(5) * t3 + t166 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, t98, t93, t63, 0, 0, 0, 0, 0, 0, t53, t225, t62, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, t82, t58, -t83, -t54, t135, -t11, -t12, 0, 0;];
tauJ_reg  = t10;
