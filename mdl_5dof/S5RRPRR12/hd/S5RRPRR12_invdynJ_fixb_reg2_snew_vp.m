% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRPRR12
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRPRR12_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR12_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR12_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR12_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:30:38
% EndTime: 2019-12-31 20:30:46
% DurationCPUTime: 3.01s
% Computational Cost: add. (8629->328), mult. (18042->402), div. (0->0), fcn. (11236->8), ass. (0->208)
t163 = sin(qJ(4));
t156 = qJDD(2) - qJDD(4);
t164 = sin(qJ(2));
t167 = cos(qJ(4));
t168 = cos(qJ(2));
t119 = (t163 * t164 + t167 * t168) * qJD(1);
t218 = t164 * qJD(1);
t121 = -qJD(1) * t163 * t168 + t167 * t218;
t231 = t121 * t119;
t259 = -t156 - t231;
t261 = t163 * t259;
t260 = t167 * t259;
t170 = qJD(2) ^ 2;
t159 = t164 ^ 2;
t171 = qJD(1) ^ 2;
t230 = t159 * t171;
t140 = t170 + t230;
t144 = t168 * t171 * t164;
t138 = qJDD(2) - t144;
t223 = t168 * t138;
t258 = pkin(6) * (-t140 * t164 + t223);
t212 = qJD(2) - qJD(4);
t199 = t119 * t212;
t215 = qJD(1) * qJD(2);
t206 = t168 * t215;
t214 = t164 * qJDD(1);
t131 = t206 + t214;
t207 = t164 * t215;
t213 = t168 * qJDD(1);
t132 = -t207 + t213;
t83 = -qJD(4) * t119 + t131 * t167 - t132 * t163;
t71 = t83 + t199;
t162 = sin(qJ(5));
t166 = cos(qJ(5));
t101 = t121 * t166 - t162 * t212;
t99 = t121 * t162 + t166 * t212;
t78 = t101 * t99;
t201 = t163 * t131 + t132 * t167;
t82 = -t121 * qJD(4) - t201;
t80 = qJDD(5) - t82;
t253 = -t78 + t80;
t257 = t162 * t253;
t256 = t166 * t253;
t255 = t132 - t207;
t211 = t212 ^ 2;
t254 = pkin(6) - pkin(7);
t160 = t168 ^ 2;
t229 = t160 * t171;
t251 = t223 + t164 * (-t170 + t229);
t116 = qJD(5) + t119;
t202 = t156 * t166 + t162 * t83;
t45 = (qJD(5) - t116) * t101 + t202;
t96 = t99 ^ 2;
t97 = t101 ^ 2;
t114 = t116 ^ 2;
t117 = t119 ^ 2;
t118 = t121 ^ 2;
t249 = 2 * qJD(3);
t248 = pkin(2) + pkin(3);
t247 = pkin(2) * t168;
t246 = t164 * g(3);
t245 = t168 * g(3);
t244 = t170 * pkin(2);
t137 = qJDD(2) + t144;
t165 = sin(qJ(1));
t169 = cos(qJ(1));
t195 = g(1) * t169 + g(2) * t165;
t234 = qJDD(1) * pkin(6);
t123 = -pkin(1) * t171 - t195 + t234;
t228 = t164 * t123;
t182 = -qJDD(2) * pkin(2) - qJ(3) * t170 + qJDD(3) + t228;
t235 = qJ(3) * t164;
t188 = -t235 - t247;
t219 = qJD(2) * t168;
t172 = -t131 * pkin(7) + t245 - t137 * pkin(3) + (pkin(7) * t219 + t188 * t218) * qJD(1) + t182;
t139 = -qJD(2) * pkin(3) - pkin(7) * t218;
t113 = t168 * t123;
t180 = t171 * t188;
t197 = qJDD(2) * qJ(3) + qJD(2) * t249 + t168 * t180 + t113;
t181 = t197 - t246;
t79 = t181 - t244;
t63 = -pkin(3) * t229 - pkin(7) * t132 + qJD(2) * t139 + t79;
t34 = t163 * t172 + t167 * t63;
t33 = t163 * t63 - t167 * t172;
t93 = pkin(4) * t119 - pkin(8) * t121;
t27 = pkin(4) * t156 - pkin(8) * t211 + t121 * t93 + t33;
t243 = t162 * t27;
t52 = t78 + t80;
t242 = t162 * t52;
t204 = t165 * g(1) - g(2) * t169;
t122 = qJDD(1) * pkin(1) + t171 * pkin(6) + t204;
t175 = pkin(2) * t255 + t122;
t54 = t131 * qJ(3) + t132 * pkin(3) - pkin(7) * t229 + (qJ(3) * t219 + (t249 + t139) * t164) * qJD(1) + t175;
t241 = t163 * t54;
t90 = -t231 + t156;
t240 = t163 * t90;
t239 = t166 * t27;
t238 = t166 * t52;
t237 = t167 * t54;
t236 = t167 * t90;
t233 = t116 * t162;
t232 = t116 * t166;
t133 = -0.2e1 * t207 + t213;
t227 = t164 * t133;
t226 = t164 * t137;
t142 = -t170 - t229;
t222 = pkin(6) * (t142 * t168 - t226) + pkin(1) * t133;
t220 = t159 + t160;
t135 = t220 * t171;
t221 = pkin(1) * t135 + t220 * t234;
t216 = qJD(5) + t116;
t210 = t99 * t233;
t209 = t163 * t78;
t208 = t167 * t78;
t198 = t212 * t121;
t26 = -t71 * pkin(8) + (-t82 - t198) * pkin(4) + t54;
t28 = -pkin(4) * t211 - pkin(8) * t156 - t119 * t93 + t34;
t10 = t162 * t26 + t166 * t28;
t9 = t162 * t28 - t166 * t26;
t5 = t10 * t166 + t162 * t9;
t205 = pkin(4) * t163 + qJ(3);
t107 = t228 + t245;
t108 = t113 - t246;
t203 = t107 * t164 + t108 * t168;
t200 = pkin(4) * t167 + t248;
t196 = -pkin(4) * t27 + pkin(8) * t5;
t193 = t10 * t162 - t166 * t9;
t192 = t163 * t199;
t191 = t163 * t198;
t190 = t167 * t199;
t189 = t167 * t198;
t187 = t131 + t206;
t15 = t163 * t34 - t167 * t33;
t16 = t163 * t33 + t167 * t34;
t185 = t162 * t156 - t166 * t83;
t130 = 0.2e1 * t206 + t214;
t184 = t130 * t168 + t227;
t68 = qJD(2) * t121 + t201;
t76 = -t97 - t114;
t32 = -t162 * t76 - t238;
t50 = t216 * t99 + t185;
t179 = pkin(4) * t50 + pkin(8) * t32 + t243;
t66 = -t114 - t96;
t30 = t166 * t66 - t257;
t46 = -t101 * t216 - t202;
t178 = pkin(4) * t46 + pkin(8) * t30 - t239;
t177 = t168 * t248 + pkin(1) + t235;
t58 = -qJD(5) * t99 - t185;
t88 = t116 * t99;
t49 = t58 + t88;
t24 = t162 * t49 - t166 * t45;
t62 = t96 + t97;
t176 = pkin(4) * t62 + pkin(8) * t24 + t5;
t174 = t164 * t180 + t182;
t173 = t218 * t249 + t175;
t81 = t174 + t245;
t136 = (t159 - t160) * t171;
t110 = -t118 + t211;
t109 = t117 - t211;
t106 = -t118 - t211;
t105 = t226 + t168 * (t170 - t230);
t104 = t187 * t164;
t103 = t255 * t168;
t94 = t118 - t117;
t89 = -t211 - t117;
t87 = -t97 + t114;
t86 = t96 - t114;
t85 = t101 * t232;
t77 = t97 - t96;
t74 = -t106 * t163 + t236;
t73 = t106 * t167 + t240;
t72 = -t199 + t83;
t67 = (0.2e1 * qJD(4) - qJD(2)) * t121 + t201;
t65 = t167 * t89 - t261;
t64 = t163 * t89 + t260;
t57 = -qJD(5) * t101 - t202;
t56 = (t101 * t162 - t166 * t99) * t116;
t55 = -t85 - t210;
t48 = t58 - t88;
t44 = -t101 * t233 + t166 * t58;
t43 = -t162 * t58 - t85;
t42 = -t162 * t57 + t232 * t99;
t41 = t166 * t57 + t210;
t40 = t163 * t72 - t167 * t68;
t39 = -t163 * t68 - t167 * t72;
t38 = t166 * t86 - t242;
t37 = -t162 * t87 + t256;
t36 = t162 * t86 + t238;
t35 = t166 * t87 + t257;
t31 = t166 * t76 - t242;
t29 = t162 * t66 + t256;
t23 = -t162 * t48 + t166 * t46;
t22 = -t162 * t45 - t166 * t49;
t21 = t162 * t46 + t166 * t48;
t20 = -t163 * t50 + t167 * t32;
t19 = t163 * t32 + t167 * t50;
t18 = -t163 * t46 + t167 * t30;
t17 = t163 * t30 + t167 * t46;
t14 = -t163 * t62 + t167 * t24;
t13 = t163 * t24 + t167 * t62;
t12 = -pkin(8) * t31 + t239;
t11 = -pkin(8) * t29 + t243;
t7 = -pkin(4) * t31 + t10;
t6 = -pkin(4) * t29 + t9;
t3 = t163 * t27 + t167 * t5;
t2 = t163 * t5 - t167 * t27;
t1 = -pkin(8) * t22 - t193;
t4 = [0, 0, 0, 0, 0, qJDD(1), t204, t195, 0, 0, t104, t184, t105, t103, t251, 0, t122 * t168 + t222, -pkin(1) * t130 - t164 * t122 - t258, t203 + t221, pkin(1) * t122 + pkin(6) * t203, t104, t105, -t184, 0, -t251, t103, t168 * (pkin(2) * t133 + t173) + (t168 * t187 + t227) * qJ(3) + t222, t164 * (qJ(3) * t135 + t174) + t168 * (pkin(2) * t135 + t197 - t244) + t221, t164 * t173 + t258 + (pkin(1) + t247) * t130 + (t130 + t187) * t235, pkin(6) * (t164 * t81 + t168 * t79) + (pkin(1) - t188) * (qJ(3) * t187 + t173), t164 * (t167 * t83 + t191) + t168 * (-t163 * t83 + t189), t164 * (-t163 * t71 - t167 * t67) + t168 * (t163 * t67 - t167 * t71), t164 * (-t110 * t163 + t260) + t168 * (-t110 * t167 - t261), t164 * (-t163 * t82 - t190) + t168 * (-t167 * t82 + t192), t164 * (t109 * t167 + t240) + t168 * (-t109 * t163 + t236), t164 * (t190 - t191) + t168 * (-t192 - t189), t164 * (-pkin(7) * t64 + t241) + t168 * (-pkin(7) * t65 + t237) + pkin(6) * (t164 * t64 + t168 * t65) + t177 * t67, t164 * (-pkin(7) * t73 + t237) + t168 * (-pkin(7) * t74 - t241) + pkin(6) * (t164 * t73 + t168 * t74) + t177 * t71, t164 * (-pkin(7) * t39 - t15) + t168 * (-pkin(7) * t40 - t16) + pkin(6) * (t164 * t39 + t168 * t40) + t177 * (-t117 - t118), t177 * t54 + t254 * (t164 * t15 + t168 * t16), t164 * (t167 * t44 + t209) + t168 * (-t163 * t44 + t208), t164 * (t163 * t77 + t167 * t23) + t168 * (-t163 * t23 + t167 * t77), t164 * (t163 * t49 + t167 * t37) + t168 * (-t163 * t37 + t167 * t49), t164 * (t167 * t42 - t209) + t168 * (-t163 * t42 - t208), t164 * (-t163 * t45 + t167 * t38) + t168 * (-t163 * t38 - t167 * t45), t164 * (t163 * t80 + t167 * t56) + t168 * (-t163 * t56 + t167 * t80), t164 * (-pkin(7) * t17 + t11 * t167 - t163 * t6) + t168 * (-pkin(7) * t18 - t163 * t11 - t167 * t6) + pkin(6) * (t164 * t17 + t168 * t18) + t177 * t29, t164 * (-pkin(7) * t19 + t12 * t167 - t163 * t7) + t168 * (-pkin(7) * t20 - t163 * t12 - t167 * t7) + pkin(6) * (t164 * t19 + t168 * t20) + t177 * t31, t164 * (-pkin(7) * t13 + t167 * t1) + t168 * (-pkin(7) * t14 - t163 * t1) + pkin(6) * (t13 * t164 + t14 * t168) + (t164 * t205 + t168 * t200 + pkin(1)) * t22, (t164 * (-pkin(8) * t167 + t205) + t168 * (pkin(8) * t163 + t200) + pkin(1)) * t193 + t254 * (t164 * t2 + t168 * t3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t144, t136, t214, t144, t213, qJDD(2), -t107, -t108, 0, 0, -t144, t214, -t136, qJDD(2), -t213, t144, pkin(2) * t137 + qJ(3) * t142 - t81, (-pkin(2) * t164 + qJ(3) * t168) * qJDD(1), qJ(3) * t138 + (t140 - t170) * pkin(2) + t181, -pkin(2) * t81 + qJ(3) * t79, -t231, -t94, -t72, t231, t68, t156, qJ(3) * t65 - t248 * t64 + t33, qJ(3) * t74 - t248 * t73 + t34, qJ(3) * t40 - t248 * t39, qJ(3) * t16 - t15 * t248, t43, -t21, -t35, -t41, -t36, -t55, qJ(3) * t18 - t17 * t248 - t178, qJ(3) * t20 - t19 * t248 - t179, qJ(3) * t14 - t13 * t248 - t176, qJ(3) * t3 - t2 * t248 - t196; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t137, t214, -t140, t81, 0, 0, 0, 0, 0, 0, t64, t73, t39, t15, 0, 0, 0, 0, 0, 0, t17, t19, t13, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t231, t94, t72, -t231, -t68, -t156, -t33, -t34, 0, 0, -t43, t21, t35, t41, t36, t55, t178, t179, t176, t196; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, t77, t49, -t78, -t45, t80, -t9, -t10, 0, 0;];
tauJ_reg = t4;
