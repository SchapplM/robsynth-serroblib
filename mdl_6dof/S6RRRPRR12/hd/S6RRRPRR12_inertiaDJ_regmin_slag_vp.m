% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x35]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRRPRR12_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR12_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_inertiaDJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:46:29
% EndTime: 2019-03-09 19:46:44
% DurationCPUTime: 4.80s
% Computational Cost: add. (8636->437), mult. (23353->796), div. (0->0), fcn. (23601->12), ass. (0->193)
t191 = sin(qJ(5));
t189 = cos(pkin(12));
t252 = pkin(10) + qJ(4);
t167 = t252 * t189;
t187 = sin(pkin(12));
t203 = qJD(4) * t187 + qJD(5) * t167;
t261 = t191 * t203;
t192 = sin(qJ(3));
t196 = cos(qJ(3));
t206 = -t196 * pkin(3) - t192 * qJ(4);
t164 = -pkin(2) + t206;
t156 = t189 * t164;
t242 = t189 * t192;
t110 = -pkin(10) * t242 + t156 + (-pkin(9) * t187 - pkin(4)) * t196;
t241 = t189 * t196;
t177 = pkin(9) * t241;
t130 = t187 * t164 + t177;
t246 = t187 * t192;
t120 = -pkin(10) * t246 + t130;
t195 = cos(qJ(5));
t74 = t191 * t110 + t195 * t120;
t166 = t252 * t187;
t122 = -t191 * t166 + t195 * t167;
t180 = -t189 * pkin(4) - pkin(3);
t260 = 0.2e1 * t180;
t188 = sin(pkin(6));
t259 = 0.2e1 * t188;
t161 = t195 * t187 + t191 * t189;
t150 = t161 * qJD(5);
t258 = pkin(5) * t150;
t257 = pkin(9) * t188;
t193 = sin(qJ(2));
t244 = t188 * t193;
t247 = cos(pkin(6));
t152 = t247 * t192 + t196 * t244;
t197 = cos(qJ(2));
t236 = qJD(2) * t197;
t215 = t188 * t236;
t118 = qJD(3) * t152 + t192 * t215;
t256 = t118 * pkin(5);
t240 = t191 * t187;
t160 = -t195 * t189 + t240;
t149 = t160 * qJD(5);
t255 = t149 * pkin(11);
t151 = t192 * t244 - t247 * t196;
t254 = t151 * pkin(5);
t253 = t196 * pkin(5);
t243 = t188 * t197;
t117 = t152 * t189 - t187 * t243;
t219 = pkin(1) * t247;
t138 = -t247 * pkin(2) + pkin(8) * t244 - t197 * t219;
t83 = t151 * pkin(3) - t152 * qJ(4) + t138;
t199 = pkin(8) * t243 + t193 * t219;
t139 = t247 * pkin(9) + t199;
t140 = (-pkin(2) * t197 - pkin(9) * t193 - pkin(1)) * t188;
t239 = t196 * t139 + t192 * t140;
t84 = -qJ(4) * t243 + t239;
t54 = -t187 * t84 + t189 * t83;
t41 = t151 * pkin(4) - t117 * pkin(10) + t54;
t116 = t152 * t187 + t189 * t243;
t55 = t187 * t83 + t189 * t84;
t48 = -t116 * pkin(10) + t55;
t20 = t191 * t41 + t195 * t48;
t237 = qJD(2) * t193;
t144 = (pkin(2) * t193 - pkin(9) * t197) * t188 * qJD(2);
t212 = qJD(2) * t247;
t207 = t197 * t212;
t216 = t188 * t237;
t145 = -pkin(1) * t207 + pkin(8) * t216;
t182 = qJD(3) * t192;
t235 = qJD(3) * t196;
t66 = t139 * t182 - t140 * t235 - t192 * t144 + t196 * t145;
t57 = (qJ(4) * t237 - qJD(4) * t197) * t188 - t66;
t119 = -qJD(3) * t151 + t196 * t215;
t146 = t199 * qJD(2);
t62 = t118 * pkin(3) - t119 * qJ(4) - t152 * qJD(4) + t146;
t31 = t187 * t62 + t189 * t57;
t77 = t195 * t116 + t191 * t117;
t18 = -t77 * pkin(11) + t20;
t194 = cos(qJ(6));
t251 = t194 * t18;
t142 = t161 * t192;
t68 = -t142 * pkin(11) + t74;
t250 = t194 * t68;
t67 = -t139 * t235 - t140 * t182 + t196 * t144 + t192 * t145;
t60 = -pkin(3) * t216 - t67;
t249 = t60 * t187;
t248 = t60 * t189;
t245 = t187 * t196;
t229 = qJD(5) * t195;
t233 = qJD(4) * t195;
t238 = t166 * t229 - t189 * t233;
t147 = -t192 * qJD(4) + (pkin(3) * t192 - qJ(4) * t196) * qJD(3);
t223 = pkin(9) * t182;
t114 = t189 * t147 + t187 * t223;
t181 = pkin(9) * t235;
t217 = t187 * t235;
t153 = pkin(4) * t217 + t181;
t162 = pkin(4) * t246 + t192 * pkin(9);
t234 = qJD(3) * t197;
t232 = qJD(4) * t196;
t231 = qJD(5) * t191;
t230 = qJD(5) * t192;
t190 = sin(qJ(6));
t228 = qJD(6) * t190;
t227 = qJD(6) * t194;
t226 = -0.2e1 * pkin(2) * qJD(3);
t225 = pkin(9) * t245;
t224 = pkin(5) * t182;
t222 = pkin(5) * t228;
t221 = pkin(5) * t227;
t97 = t119 * t187 - t189 * t216;
t98 = t119 * t189 + t187 * t216;
t42 = -qJD(5) * t77 - t191 * t97 + t195 * t98;
t30 = -t187 * t57 + t189 * t62;
t200 = t118 * pkin(4) - t98 * pkin(10) + t30;
t25 = -t97 * pkin(10) + t31;
t9 = -t20 * qJD(5) - t191 * t25 + t195 * t200;
t5 = -t42 * pkin(11) + t256 + t9;
t78 = -t191 * t116 + t195 * t117;
t43 = qJD(5) * t78 + t191 * t98 + t195 * t97;
t8 = -t191 * t200 - t195 * t25 - t41 * t229 + t231 * t48;
t7 = -t43 * pkin(11) - t8;
t220 = -t190 * t7 + t194 * t5;
t185 = t188 ^ 2;
t218 = t185 * t236;
t214 = t192 * t235;
t101 = -t160 * t235 - t161 * t230;
t136 = t187 * t147;
t104 = t136 + (-pkin(9) * t242 - pkin(10) * t245) * qJD(3);
t198 = (pkin(4) * t192 - pkin(10) * t241) * qJD(3) + t114;
t45 = -qJD(5) * t74 - t191 * t104 + t195 * t198;
t32 = -t101 * pkin(11) + t224 + t45;
t102 = t161 * t235 + t229 * t242 - t230 * t240;
t44 = -t195 * t104 - t110 * t229 + t120 * t231 - t191 * t198;
t34 = -t102 * pkin(11) - t44;
t213 = -t190 * t34 + t194 * t32;
t19 = -t191 * t48 + t195 * t41;
t73 = t195 * t110 - t191 * t120;
t211 = -t192 * t139 + t196 * t140;
t121 = -t195 * t166 - t191 * t167;
t210 = t193 * t218;
t209 = 0.2e1 * (t187 ^ 2 + t189 ^ 2) * qJD(4);
t208 = -t150 * pkin(11) - t238;
t85 = pkin(3) * t243 - t211;
t17 = -t78 * pkin(11) + t19 + t254;
t11 = t190 * t17 + t251;
t143 = t160 * t192;
t63 = t143 * pkin(11) - t253 + t73;
t36 = t190 * t63 + t250;
t50 = t190 * t78 + t194 * t77;
t51 = -t190 * t77 + t194 * t78;
t100 = -t160 * pkin(11) + t122;
t99 = -t161 * pkin(11) + t121;
t65 = t194 * t100 + t190 * t99;
t115 = -t189 * t223 + t136;
t205 = -t114 * t187 + t115 * t189;
t95 = t194 * t142 - t190 * t143;
t96 = -t190 * t142 - t194 * t143;
t108 = t194 * t160 + t190 * t161;
t109 = -t190 * t160 + t194 * t161;
t204 = -qJ(4) * t118 - qJD(4) * t151;
t1 = -t17 * t227 + t18 * t228 - t190 * t5 - t194 * t7;
t12 = -t190 * t32 - t194 * t34 - t63 * t227 + t228 * t68;
t69 = t116 * pkin(4) + t85;
t202 = t192 * t234 + t196 * t237;
t201 = t192 * t237 - t196 * t234;
t49 = t97 * pkin(4) + t60;
t93 = -t161 * qJD(4) - qJD(5) * t122;
t172 = -0.2e1 * t214;
t135 = t160 * pkin(5) + t180;
t129 = t156 - t225;
t113 = t142 * pkin(5) + t162;
t92 = t238 + t261;
t91 = 0.2e1 * t151 * t118;
t86 = -t118 * t196 + t151 * t182;
t81 = t102 * pkin(5) + t153;
t71 = qJD(6) * t109 - t190 * t149 + t194 * t150;
t70 = -qJD(6) * t108 - t194 * t149 - t190 * t150;
t64 = -t190 * t100 + t194 * t99;
t53 = qJD(6) * t96 + t190 * t101 + t194 * t102;
t52 = -qJD(6) * t95 + t194 * t101 - t190 * t102;
t46 = t77 * pkin(5) + t69;
t35 = -t190 * t68 + t194 * t63;
t27 = -t190 * t208 + t194 * (-t167 * t229 - t187 * t233 + t255) + (t190 * t203 + t194 * (-qJD(4) * t189 + qJD(5) * t166)) * t191 - t65 * qJD(6);
t26 = t100 * t228 - t190 * (t93 + t255) - t194 * (t208 - t261) - t99 * t227;
t21 = t43 * pkin(5) + t49;
t15 = qJD(6) * t51 + t190 * t42 + t194 * t43;
t14 = -qJD(6) * t50 - t190 * t43 + t194 * t42;
t13 = -qJD(6) * t36 + t213;
t10 = t194 * t17 - t190 * t18;
t2 = -qJD(6) * t11 + t220;
t3 = [0, 0, 0, 0.2e1 * t210, 0.2e1 * (-t193 ^ 2 + t197 ^ 2) * t185 * qJD(2), t207 * t259, -0.2e1 * t212 * t244, 0, -0.2e1 * t185 * pkin(1) * t237 - 0.2e1 * t146 * t247, -0.2e1 * pkin(1) * t218 + 0.2e1 * t145 * t247, 0.2e1 * t152 * t119, -0.2e1 * t152 * t118 - 0.2e1 * t119 * t151 (-t119 * t197 + t152 * t237) * t259 (t118 * t197 - t151 * t237) * t259, -0.2e1 * t210, 0.2e1 * t138 * t118 + 0.2e1 * t146 * t151 + 0.2e1 * (-t67 * t197 + t211 * t237) * t188, 0.2e1 * t138 * t119 + 0.2e1 * t146 * t152 + 0.2e1 * (-t66 * t197 - t239 * t237) * t188, 0.2e1 * t60 * t116 + 0.2e1 * t54 * t118 + 0.2e1 * t30 * t151 + 0.2e1 * t85 * t97, 0.2e1 * t60 * t117 - 0.2e1 * t55 * t118 - 0.2e1 * t31 * t151 + 0.2e1 * t85 * t98, -0.2e1 * t31 * t116 - 0.2e1 * t30 * t117 - 0.2e1 * t54 * t98 - 0.2e1 * t55 * t97, 0.2e1 * t54 * t30 + 0.2e1 * t55 * t31 + 0.2e1 * t85 * t60, 0.2e1 * t78 * t42, -0.2e1 * t42 * t77 - 0.2e1 * t78 * t43, 0.2e1 * t78 * t118 + 0.2e1 * t42 * t151, -0.2e1 * t77 * t118 - 0.2e1 * t43 * t151, t91, 0.2e1 * t19 * t118 + 0.2e1 * t9 * t151 + 0.2e1 * t69 * t43 + 0.2e1 * t49 * t77, -0.2e1 * t20 * t118 + 0.2e1 * t8 * t151 + 0.2e1 * t69 * t42 + 0.2e1 * t49 * t78, 0.2e1 * t51 * t14, -0.2e1 * t14 * t50 - 0.2e1 * t15 * t51, 0.2e1 * t51 * t118 + 0.2e1 * t14 * t151, -0.2e1 * t50 * t118 - 0.2e1 * t15 * t151, t91, 0.2e1 * t10 * t118 + 0.2e1 * t46 * t15 + 0.2e1 * t2 * t151 + 0.2e1 * t21 * t50, 0.2e1 * t1 * t151 - 0.2e1 * t11 * t118 + 0.2e1 * t46 * t14 + 0.2e1 * t21 * t51; 0, 0, 0, 0, 0, t215, -t216, 0, -t146, t145, t119 * t192 + t152 * t235, -t192 * t118 + t119 * t196 + (-t151 * t196 - t152 * t192) * qJD(3), t201 * t188, t202 * t188, 0, -pkin(2) * t118 + t138 * t182 - t146 * t196 - t201 * t257, -pkin(2) * t119 + t138 * t235 + t146 * t192 - t202 * t257, t114 * t151 + t129 * t118 - t30 * t196 + (pkin(9) * t97 + t249) * t192 + (t192 * t54 + (pkin(9) * t116 + t187 * t85) * t196) * qJD(3), -t115 * t151 - t130 * t118 + t31 * t196 + (pkin(9) * t98 + t248) * t192 + (-t192 * t55 + (pkin(9) * t117 + t189 * t85) * t196) * qJD(3), -t114 * t117 - t115 * t116 - t129 * t98 - t130 * t97 + (-t187 * t31 - t189 * t30) * t192 + (-t187 * t55 - t189 * t54) * t235, t54 * t114 + t55 * t115 + t30 * t129 + t31 * t130 + (t192 * t60 + t235 * t85) * pkin(9), t78 * t101 - t42 * t143, -t101 * t77 - t78 * t102 - t42 * t142 + t143 * t43, t101 * t151 - t143 * t118 + t182 * t78 - t42 * t196, -t102 * t151 - t142 * t118 - t182 * t77 + t43 * t196, t86, t69 * t102 + t73 * t118 + t49 * t142 + t45 * t151 + t153 * t77 + t162 * t43 + t182 * t19 - t9 * t196, t69 * t101 - t74 * t118 - t49 * t143 + t44 * t151 + t153 * t78 + t162 * t42 - t182 * t20 - t8 * t196, t14 * t96 + t51 * t52, -t14 * t95 - t96 * t15 - t52 * t50 - t51 * t53, t96 * t118 - t14 * t196 + t52 * t151 + t182 * t51, -t95 * t118 + t15 * t196 - t53 * t151 - t182 * t50, t86, t10 * t182 + t113 * t15 + t35 * t118 + t13 * t151 - t2 * t196 + t21 * t95 + t46 * t53 + t81 * t50, -t1 * t196 - t11 * t182 + t113 * t14 - t36 * t118 + t12 * t151 + t21 * t96 + t46 * t52 + t81 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t214, 0.2e1 * (-t192 ^ 2 + t196 ^ 2) * qJD(3), 0, 0, 0, t192 * t226, t196 * t226, -0.2e1 * t114 * t196 + 0.2e1 * (t129 + 0.2e1 * t225) * t182, 0.2e1 * t115 * t196 + 0.2e1 * (-t130 + 0.2e1 * t177) * t182, 0.2e1 * (-t114 * t189 - t115 * t187) * t192 + 0.2e1 * (-t129 * t189 - t130 * t187) * t235, 0.2e1 * pkin(9) ^ 2 * t214 + 0.2e1 * t129 * t114 + 0.2e1 * t130 * t115, -0.2e1 * t143 * t101, -0.2e1 * t101 * t142 + 0.2e1 * t143 * t102, -0.2e1 * t101 * t196 - 0.2e1 * t143 * t182, 0.2e1 * t102 * t196 - 0.2e1 * t142 * t182, t172, 0.2e1 * t162 * t102 + 0.2e1 * t153 * t142 + 0.2e1 * t182 * t73 - 0.2e1 * t45 * t196, 0.2e1 * t162 * t101 - 0.2e1 * t153 * t143 - 0.2e1 * t182 * t74 - 0.2e1 * t44 * t196, 0.2e1 * t96 * t52, -0.2e1 * t52 * t95 - 0.2e1 * t96 * t53, 0.2e1 * t182 * t96 - 0.2e1 * t52 * t196, -0.2e1 * t182 * t95 + 0.2e1 * t53 * t196, t172, 0.2e1 * t113 * t53 - 0.2e1 * t13 * t196 + 0.2e1 * t182 * t35 + 0.2e1 * t81 * t95, 0.2e1 * t113 * t52 - 0.2e1 * t12 * t196 - 0.2e1 * t182 * t36 + 0.2e1 * t81 * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119, -t118, t216, t67, t66, -pkin(3) * t97 + t187 * t204 - t248, -pkin(3) * t98 + t189 * t204 + t249 (-qJ(4) * t97 - qJD(4) * t116 + t31) * t189 + (qJ(4) * t98 + qJD(4) * t117 - t30) * t187, -t60 * pkin(3) + (-t187 * t54 + t189 * t55) * qJD(4) + (-t30 * t187 + t31 * t189) * qJ(4), -t78 * t149 + t42 * t161, t149 * t77 - t78 * t150 - t42 * t160 - t161 * t43, t161 * t118 - t149 * t151, -t160 * t118 - t150 * t151, 0, t121 * t118 + t69 * t150 + t93 * t151 + t49 * t160 + t180 * t43, -t122 * t118 - t69 * t149 + t92 * t151 + t49 * t161 + t180 * t42, t14 * t109 + t51 * t70, -t14 * t108 - t109 * t15 - t70 * t50 - t51 * t71, t109 * t118 + t70 * t151, -t108 * t118 - t71 * t151, 0, t21 * t108 + t64 * t118 + t135 * t15 + t27 * t151 + t50 * t258 + t46 * t71, t21 * t109 - t65 * t118 + t135 * t14 + t26 * t151 + t51 * t258 + t46 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t235, -t182, 0, -t181, t223, t187 * t232 + (t187 * t206 - t177) * qJD(3), t189 * t232 + (t189 * t206 + t225) * qJD(3), t205, -pkin(3) * t181 + (-t129 * t187 + t130 * t189) * qJD(4) + t205 * qJ(4), t101 * t161 + t143 * t149, -t101 * t160 - t161 * t102 + t149 * t142 + t143 * t150, t149 * t196 + t161 * t182, t150 * t196 - t160 * t182, 0, t180 * t102 + t121 * t182 + t162 * t150 + t153 * t160 - t93 * t196, t180 * t101 - t122 * t182 - t162 * t149 + t153 * t161 - t92 * t196, t52 * t109 + t96 * t70, -t52 * t108 - t109 * t53 - t70 * t95 - t96 * t71, t109 * t182 - t70 * t196, -t108 * t182 + t71 * t196, 0, t81 * t108 + t113 * t71 + t135 * t53 + t182 * t64 - t27 * t196 + t95 * t258, t81 * t109 + t113 * t70 + t135 * t52 - t182 * t65 - t26 * t196 + t96 * t258; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t209, qJ(4) * t209, -0.2e1 * t161 * t149, 0.2e1 * t149 * t160 - 0.2e1 * t161 * t150, 0, 0, 0, t150 * t260, -t149 * t260, 0.2e1 * t109 * t70, -0.2e1 * t70 * t108 - 0.2e1 * t109 * t71, 0, 0, 0, 0.2e1 * t108 * t258 + 0.2e1 * t135 * t71, 0.2e1 * t109 * t258 + 0.2e1 * t135 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, t98, 0, t60, 0, 0, 0, 0, 0, t43, t42, 0, 0, 0, 0, 0, t15, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t217, t189 * t235, 0, t181, 0, 0, 0, 0, 0, t102, t101, 0, 0, 0, 0, 0, t53, t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t150, -t149, 0, 0, 0, 0, 0, t71, t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, -t43, t118, t9, t8, 0, 0, t14, -t15, t118, t194 * t256 + (-t251 + (-t17 - t254) * t190) * qJD(6) + t220 (-t118 * t190 - t151 * t227) * pkin(5) + t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, -t102, t182, t45, t44, 0, 0, t52, -t53, t182, t194 * t224 + (-t250 + (-t63 + t253) * t190) * qJD(6) + t213 (-t182 * t190 + t196 * t227) * pkin(5) + t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t149, -t150, 0, t93, t92, 0, 0, t70, -t71, 0, t27, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t222, -0.2e1 * t221; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t15, t118, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t53, t182, t13, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, -t71, 0, t27, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t222, -t221; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;