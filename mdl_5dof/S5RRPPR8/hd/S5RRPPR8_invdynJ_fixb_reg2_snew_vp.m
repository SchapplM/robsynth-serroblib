% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRPPR8_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR8_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR8_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR8_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:39:16
% EndTime: 2019-12-31 19:39:23
% DurationCPUTime: 2.62s
% Computational Cost: add. (8643->320), mult. (19971->417), div. (0->0), fcn. (12422->8), ass. (0->196)
t168 = qJD(2) ^ 2;
t164 = sin(qJ(2));
t159 = t164 ^ 2;
t169 = qJD(1) ^ 2;
t211 = t159 * t169;
t141 = t168 + t211;
t167 = cos(qJ(2));
t204 = t167 * t169;
t144 = t164 * t204;
t139 = qJDD(2) - t144;
t205 = t167 * t139;
t246 = pkin(6) * (-t164 * t141 + t205);
t161 = sin(pkin(8));
t162 = cos(pkin(8));
t119 = (-t164 * t161 - t167 * t162) * qJD(1);
t198 = qJD(2) * t119;
t194 = qJD(1) * qJD(2);
t190 = t167 * t194;
t193 = t164 * qJDD(1);
t132 = t190 + t193;
t150 = t167 * qJDD(1);
t191 = t164 * t194;
t179 = -t150 + t191;
t98 = t162 * t132 + t161 * t179;
t80 = -t98 - t198;
t163 = sin(qJ(5));
t154 = qJDD(2) - qJDD(5);
t199 = qJD(1) * t167;
t200 = qJD(1) * t164;
t121 = -t161 * t199 + t162 * t200;
t166 = cos(qJ(5));
t89 = -t166 * t119 + t163 * t121;
t91 = t163 * t119 + t166 * t121;
t69 = t91 * t89;
t242 = -t69 - t154;
t245 = t163 * t242;
t244 = t166 * t242;
t185 = t161 * t132 - t162 * t179;
t53 = -t89 * qJD(5) - t163 * t185 + t166 * t98;
t155 = qJD(2) - qJD(5);
t86 = t89 * t155;
t243 = t86 + t53;
t214 = t121 * t119;
t183 = -qJDD(2) + t214;
t241 = t161 * t183;
t240 = t162 * t183;
t165 = sin(qJ(1));
t231 = cos(qJ(1));
t175 = t231 * g(1) + t165 * g(2);
t215 = qJDD(1) * pkin(6);
t124 = -t169 * pkin(1) - t175 + t215;
t210 = t164 * qJ(3);
t227 = t167 * pkin(2);
t181 = -t210 - t227;
t129 = t181 * qJD(1);
t184 = qJD(1) * t129 + t124;
t239 = t164 * t184;
t115 = qJD(2) * t121;
t78 = -t185 - t115;
t235 = t167 ^ 2;
t203 = t235 * t169;
t238 = t205 + t164 * (-t168 + t203);
t187 = t163 * t98 + t166 * t185;
t40 = (qJD(5) + t155) * t91 + t187;
t87 = t89 ^ 2;
t88 = t91 ^ 2;
t236 = t119 ^ 2;
t118 = t121 ^ 2;
t153 = t155 ^ 2;
t234 = 2 * qJD(3);
t233 = 2 * qJD(4);
t232 = pkin(2) + pkin(3);
t137 = -qJD(2) * pkin(3) - qJ(4) * t200;
t113 = t167 * t124;
t182 = qJDD(2) * qJ(3) + qJD(2) * t234 + t129 * t199 + t113;
t228 = t164 * g(3);
t174 = t182 - t228;
t225 = t168 * pkin(2);
t75 = t174 - t225;
t66 = -pkin(3) * t203 + qJ(4) * t179 + qJD(2) * t137 + t75;
t189 = -qJDD(2) * pkin(2) - t168 * qJ(3) + qJDD(3);
t226 = t167 * g(3);
t176 = t189 + t226;
t68 = -qJDD(2) * pkin(3) + (-t132 + t190) * qJ(4) + (-pkin(3) * t204 + t184) * t164 + t176;
t30 = t121 * t233 + t161 * t66 - t162 * t68;
t28 = t183 * pkin(4) + t80 * pkin(7) - t30;
t106 = -qJD(2) * pkin(4) - t121 * pkin(7);
t31 = t119 * t233 + t161 * t68 + t162 * t66;
t29 = -t236 * pkin(4) - pkin(7) * t185 + qJD(2) * t106 + t31;
t11 = t163 * t29 - t166 * t28;
t12 = t163 * t28 + t166 * t29;
t6 = -t166 * t11 + t163 * t12;
t230 = t161 * t6;
t229 = t162 * t6;
t188 = t165 * g(1) - t231 * g(2);
t123 = qJDD(1) * pkin(1) + t169 * pkin(6) + t188;
t171 = -pkin(2) * t191 + t123;
t216 = qJ(3) * t167;
t56 = t150 * pkin(2) + t132 * qJ(3) + qJDD(4) - t179 * pkin(3) - qJ(4) * t203 + (qJD(2) * t216 + (-pkin(2) * qJD(2) + t137 + t234) * t164) * qJD(1) + t171;
t224 = t161 * t56;
t93 = qJDD(2) + t214;
t223 = t161 * t93;
t222 = t162 * t56;
t221 = t162 * t93;
t32 = pkin(4) * t185 - t236 * pkin(7) + t121 * t106 + t56;
t220 = t163 * t32;
t63 = -t69 + t154;
t219 = t163 * t63;
t218 = t166 * t32;
t217 = t166 * t63;
t213 = t155 * t163;
t212 = t155 * t166;
t133 = t150 - 0.2e1 * t191;
t209 = t164 * t133;
t138 = qJDD(2) + t144;
t208 = t164 * t138;
t143 = -t168 - t203;
t202 = pkin(6) * (t167 * t143 - t208) + pkin(1) * t133;
t192 = t159 + t235;
t135 = t192 * t169;
t201 = pkin(1) * t135 + t192 * t215;
t197 = t161 * qJD(2);
t196 = t162 * qJD(2);
t7 = t163 * t11 + t166 * t12;
t104 = t164 * t124 + t226;
t105 = t113 - t228;
t186 = t164 * t104 + t167 * t105;
t180 = t132 + t190;
t15 = t161 * t31 - t162 * t30;
t16 = t161 * t30 + t162 * t31;
t131 = 0.2e1 * t190 + t193;
t177 = t167 * t131 + t209;
t173 = t167 * t232 + pkin(1) + t210;
t172 = t189 + t239;
t170 = -pkin(2) * t179 + t200 * t234 + t171;
t136 = (t159 - t235) * t169;
t109 = -t118 - t168;
t108 = -t118 + t168;
t107 = -t168 + t236;
t103 = t208 + t167 * (t168 - t211);
t102 = t180 * t164;
t101 = t133 * t167;
t92 = -t168 - t236;
t84 = -t88 + t153;
t83 = t87 - t153;
t82 = t172 + t226;
t81 = -t88 - t153;
t79 = -t198 + t98;
t77 = -t115 + t185;
t76 = -t118 - t236;
t73 = -t161 * t109 + t221;
t72 = t162 * t109 + t223;
t71 = t162 * t92 - t241;
t70 = t161 * t92 + t240;
t67 = t88 - t87;
t60 = -t153 - t87;
t58 = (-t163 * t91 + t166 * t89) * t155;
t57 = (t163 * t89 + t166 * t91) * t155;
t55 = -t161 * t80 + t162 * t78;
t54 = t161 * t78 + t162 * t80;
t52 = -t91 * qJD(5) - t187;
t51 = -t87 - t88;
t50 = t166 * t83 + t219;
t49 = -t163 * t84 + t244;
t48 = t163 * t83 - t217;
t47 = t166 * t84 + t245;
t46 = -t163 * t81 + t217;
t45 = t166 * t81 + t219;
t44 = -t86 + t53;
t39 = (qJD(5) - t155) * t91 + t187;
t38 = t166 * t53 + t91 * t213;
t37 = t163 * t53 - t91 * t212;
t36 = -t163 * t52 - t89 * t212;
t35 = t166 * t52 - t89 * t213;
t34 = t166 * t60 - t245;
t33 = t163 * t60 + t244;
t26 = -t161 * t45 + t162 * t46;
t25 = t161 * t46 + t162 * t45;
t24 = t163 * t44 - t166 * t40;
t23 = -t163 * t243 - t166 * t39;
t22 = -t163 * t40 - t166 * t44;
t21 = -t163 * t39 + t166 * t243;
t20 = -pkin(7) * t45 + t218;
t19 = -t161 * t33 + t162 * t34;
t18 = t161 * t34 + t162 * t33;
t17 = -pkin(7) * t33 + t220;
t14 = -pkin(4) * t243 + pkin(7) * t46 + t220;
t13 = -pkin(4) * t39 + pkin(7) * t34 - t218;
t9 = -t161 * t22 + t162 * t24;
t8 = t161 * t24 + t162 * t22;
t5 = -pkin(4) * t32 + pkin(7) * t7;
t4 = -pkin(7) * t22 - t6;
t3 = -pkin(4) * t51 + pkin(7) * t24 + t7;
t2 = t162 * t7 - t230;
t1 = t161 * t7 + t229;
t10 = [0, 0, 0, 0, 0, qJDD(1), t188, t175, 0, 0, t102, t177, t103, t101, t238, 0, t167 * t123 + t202, -pkin(1) * t131 - t164 * t123 - t246, t186 + t201, pkin(1) * t123 + pkin(6) * t186, t102, t103, -t177, 0, -t238, t101, t167 * (pkin(2) * t133 + t170) + (t167 * t180 + t209) * qJ(3) + t202, t167 * (pkin(2) * t135 + t182 - t225) + (qJ(3) * t135 + t172) * t164 + t201, t164 * t170 + t246 + (pkin(1) + t227) * t131 + (t131 + t180) * t210, pkin(6) * (t164 * t82 + t167 * t75) + (pkin(1) - t181) * (qJ(3) * t180 + t170), t164 * (t121 * t197 + t162 * t98) + t167 * (t121 * t196 - t161 * t98), t164 * (-t161 * t79 - t162 * t77) + t167 * (t161 * t77 - t162 * t79), t164 * (-t161 * t108 + t240) + t167 * (-t162 * t108 - t241), t164 * (t119 * t196 + t161 * t185) + t167 * (-t119 * t197 + t162 * t185), t164 * (t162 * t107 + t223) + t167 * (-t161 * t107 + t221), (t164 * (-t162 * t119 - t121 * t161) + t167 * (t161 * t119 - t121 * t162)) * qJD(2), t164 * (-qJ(4) * t70 + t224) + t167 * (-qJ(4) * t71 + t222) + pkin(6) * (t164 * t70 + t167 * t71) + t173 * t77, t164 * (-qJ(4) * t72 + t222) + t167 * (-qJ(4) * t73 - t224) + pkin(6) * (t164 * t72 + t167 * t73) + t173 * t79, t164 * (-qJ(4) * t54 - t15) + t167 * (-qJ(4) * t55 - t16) + pkin(6) * (t164 * t54 + t167 * t55) + t173 * t76, t173 * t56 + (pkin(6) - qJ(4)) * (t164 * t15 + t167 * t16), t164 * (-t161 * t37 + t162 * t38) + t167 * (-t161 * t38 - t162 * t37), t164 * (-t161 * t21 + t162 * t23) + t167 * (-t161 * t23 - t162 * t21), t164 * (-t161 * t47 + t162 * t49) + t167 * (-t161 * t49 - t162 * t47), t164 * (-t161 * t35 + t162 * t36) + t167 * (-t161 * t36 - t162 * t35), t164 * (-t161 * t48 + t162 * t50) + t167 * (-t161 * t50 - t162 * t48), t164 * (-t161 * t57 + t162 * t58) + t167 * (-t161 * t58 - t162 * t57), t164 * (-qJ(4) * t18 - t161 * t13 + t162 * t17) + t167 * (-qJ(4) * t19 - t162 * t13 - t161 * t17) + pkin(6) * (t164 * t18 + t167 * t19) + t173 * t39, t164 * (-qJ(4) * t25 - t161 * t14 + t162 * t20) + t167 * (-qJ(4) * t26 - t162 * t14 - t161 * t20) + pkin(6) * (t164 * t25 + t167 * t26) + t173 * t243, t164 * (-qJ(4) * t8 - t161 * t3 + t162 * t4) + t167 * (-qJ(4) * t9 - t161 * t4 - t162 * t3) + pkin(6) * (t164 * t8 + t167 * t9) + t173 * t51, t164 * (-pkin(7) * t229 - qJ(4) * t1 - t161 * t5) + t167 * (pkin(7) * t230 - qJ(4) * t2 - t162 * t5) + pkin(6) * (t164 * t1 + t167 * t2) + t173 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t144, t136, t193, t144, t150, qJDD(2), -t104, -t105, 0, 0, -t144, t193, -t136, qJDD(2), -t150, t144, pkin(2) * t138 + qJ(3) * t143 - t176 - t239, (-pkin(2) * t164 + t216) * qJDD(1), qJ(3) * t139 + (t141 - t168) * pkin(2) + t174, -pkin(2) * t82 + qJ(3) * t75, t214, -t118 + t236, t80, -t214, -t78, qJDD(2), qJ(3) * t71 - t232 * t70 + t30, qJ(3) * t73 - t232 * t72 + t31, qJ(3) * t55 - t232 * t54, qJ(3) * t16 - t232 * t15, -t69, -t67, -t44, t69, t40, t154, -pkin(4) * t33 + qJ(3) * t19 - t232 * t18 + t11, -pkin(4) * t45 + qJ(3) * t26 - t232 * t25 + t12, -pkin(4) * t22 + qJ(3) * t9 - t232 * t8, -pkin(4) * t6 + qJ(3) * t2 - t232 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t138, t193, -t141, t82, 0, 0, 0, 0, 0, 0, t70, t72, t54, t15, 0, 0, 0, 0, 0, 0, t18, t25, t8, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, t79, t76, t56, 0, 0, 0, 0, 0, 0, t39, t243, t51, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, t67, t44, -t69, -t40, -t154, -t11, -t12, 0, 0;];
tauJ_reg = t10;
