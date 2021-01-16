% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% 
% Output:
% tauc_reg [5x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:37
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRP10_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP10_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:35:52
% EndTime: 2021-01-16 00:36:13
% DurationCPUTime: 5.59s
% Computational Cost: add. (5778->436), mult. (15443->614), div. (0->0), fcn. (11664->8), ass. (0->196)
t159 = sin(pkin(5));
t165 = cos(qJ(2));
t229 = qJD(1) * t165;
t150 = t159 * t229;
t186 = t150 - qJD(3);
t162 = sin(qJ(2));
t247 = cos(pkin(5));
t203 = t247 * qJD(1);
t192 = pkin(1) * t203;
t148 = t162 * t192;
t121 = pkin(7) * t150 + t148;
t161 = sin(qJ(3));
t164 = cos(qJ(3));
t278 = -t121 - t186 * (pkin(3) * t161 - pkin(9) * t164);
t182 = t203 + qJD(2);
t230 = qJD(1) * t162;
t214 = t159 * t230;
t103 = t161 * t214 - t164 * t182;
t220 = qJD(1) * qJD(2);
t210 = t159 * t220;
t190 = t165 * t210;
t183 = t164 * t190;
t167 = -qJD(3) * t103 + t183;
t277 = -qJD(4) * t186 + t167;
t105 = t161 * t182 + t164 * t214;
t160 = sin(qJ(4));
t163 = cos(qJ(4));
t69 = t105 * t160 + t163 * t186;
t98 = qJD(4) + t103;
t262 = t69 * t98;
t191 = t162 * t210;
t223 = qJD(4) * t160;
t30 = t105 * t223 - t160 * t191 - t277 * t163;
t276 = -t30 - t262;
t71 = t163 * t105 - t160 * t186;
t261 = t71 * t98;
t222 = qJD(4) * t163;
t31 = t105 * t222 + t277 * t160 - t163 * t191;
t275 = t31 + t261;
t226 = qJD(3) * t161;
t219 = pkin(8) * t226;
t274 = t160 * t219 + t278 * t163;
t273 = t162 * t165;
t194 = t164 * t150;
t224 = qJD(3) * t164;
t272 = t194 - t224;
t215 = pkin(1) * t247;
t244 = t159 * t162;
t271 = -pkin(7) * t244 + t165 * t215;
t143 = -pkin(3) * t164 - pkin(9) * t161 - pkin(2);
t118 = -pkin(7) * t214 + t165 * t192;
t177 = (pkin(2) * t162 - pkin(8) * t165) * t159;
t119 = qJD(1) * t177;
t236 = t164 * t118 + t161 * t119;
t60 = pkin(9) * t214 + t236;
t270 = -t143 * t222 - t278 * t160 + t163 * t60;
t269 = t71 ^ 2;
t166 = qJD(1) ^ 2;
t86 = -t182 * pkin(2) - t118;
t44 = t103 * pkin(3) - t105 * pkin(9) + t86;
t87 = t182 * pkin(8) + t121;
t117 = (-pkin(2) * t165 - pkin(8) * t162 - pkin(1)) * t159;
t97 = qJD(1) * t117;
t53 = t161 * t97 + t164 * t87;
t47 = -t186 * pkin(9) + t53;
t14 = -t160 * t47 + t163 * t44;
t9 = -qJ(5) * t71 + t14;
t7 = pkin(4) * t98 + t9;
t268 = t7 - t9;
t184 = t161 * t190;
t73 = t105 * qJD(3) + t184;
t267 = pkin(4) * t73;
t266 = t69 * pkin(4);
t15 = t160 * t44 + t163 * t47;
t10 = -qJ(5) * t69 + t15;
t265 = t10 * t98;
t120 = qJD(2) * t177;
t112 = qJD(1) * t120;
t122 = t271 * qJD(2);
t113 = qJD(1) * t122;
t205 = t164 * t112 - t161 * t113 - t87 * t224 - t97 * t226;
t22 = -pkin(3) * t191 - t205;
t8 = pkin(4) * t31 + t22;
t264 = t160 * t8;
t263 = t163 * t8;
t260 = qJ(5) + pkin(9);
t239 = t163 * t164;
t153 = pkin(8) * t239;
t195 = t161 * t150;
t221 = qJD(5) * t163;
t246 = qJ(5) * t161;
t231 = qJD(1) * t159;
t238 = t163 * t165;
t92 = (t160 * t162 + t164 * t238) * t231;
t259 = pkin(4) * t195 - qJ(5) * t92 - t160 * t60 + t161 * t221 - (pkin(4) * t161 - qJ(5) * t239) * qJD(3) - (-t153 + (-t143 + t246) * t160) * qJD(4) - t274;
t241 = t161 * t163;
t91 = t160 * t194 - t163 * t214;
t258 = -qJ(5) * t91 - (-pkin(8) * qJD(3) - qJ(5) * qJD(4)) * t241 - (-qJD(5) * t161 + (-pkin(8) * qJD(4) - qJ(5) * qJD(3)) * t164) * t160 + t270;
t52 = -t161 * t87 + t164 * t97;
t63 = pkin(3) * t105 + pkin(9) * t103;
t257 = t160 * t63 + t163 * t52;
t115 = -t247 * pkin(2) - t271;
t126 = t161 * t244 - t247 * t164;
t127 = t247 * t161 + t164 * t244;
t56 = t126 * pkin(3) - t127 * pkin(9) + t115;
t243 = t159 * t165;
t172 = pkin(7) * t243 + t162 * t215;
t116 = t247 * pkin(8) + t172;
t237 = t164 * t116 + t161 * t117;
t58 = -pkin(9) * t243 + t237;
t256 = t160 * t56 + t163 * t58;
t254 = t160 * t30;
t253 = t160 * t71;
t252 = t160 * t73;
t251 = t163 * t73;
t106 = t161 * t118;
t59 = -pkin(3) * t214 - t119 * t164 + t106;
t250 = pkin(8) * t224 - t59 + (t160 * t224 + t161 * t222 - t91) * pkin(4);
t207 = qJD(4) * t260;
t242 = t160 * t103;
t249 = -qJ(5) * t242 - t160 * t207 + t221 - t257;
t240 = t163 * t103;
t62 = t163 * t63;
t248 = -pkin(4) * t105 - qJ(5) * t240 - t163 * t207 - t62 + (-qJD(5) + t52) * t160;
t156 = t159 ^ 2;
t245 = t156 * t166;
t233 = t160 * t143 + t153;
t232 = t162 ^ 2 - t165 ^ 2;
t228 = qJD(2) * t162;
t227 = qJD(3) * t160;
t225 = qJD(3) * t163;
t218 = t98 * t225;
t217 = t98 * t223;
t216 = t160 * t243;
t213 = t159 * t228;
t212 = qJD(2) * t243;
t211 = t156 * t220;
t209 = qJD(5) + t266;
t208 = -t160 * t58 + t163 * t56;
t176 = -t161 * t112 - t164 * t113 - t97 * t224 + t87 * t226;
t21 = pkin(9) * t191 - t176;
t185 = pkin(7) * t190;
t29 = t73 * pkin(3) - t167 * pkin(9) + qJD(2) * t148 + t185;
t206 = -t160 * t29 - t163 * t21 - t44 * t222 + t47 * t223;
t204 = -t161 * t116 + t117 * t164;
t202 = t163 * t98;
t201 = t165 * t186;
t200 = t186 * t159;
t199 = qJD(3) * t186;
t197 = 0.2e1 * t211;
t193 = qJD(4) * pkin(9) * t98 + t22;
t57 = pkin(3) * t243 - t204;
t189 = t159 * t166 * t247;
t188 = -0.2e1 * pkin(1) * t211;
t181 = 0.2e1 * t203 + qJD(2);
t180 = -t116 * t224 - t117 * t226 + t120 * t164 - t161 * t122;
t179 = qJ(5) * t31 + t206;
t178 = -t98 * t222 - t252;
t80 = t127 * t160 + t159 * t238;
t173 = -t116 * t226 + t117 * t224 + t161 * t120 + t164 * t122;
t25 = pkin(9) * t213 + t173;
t123 = t172 * qJD(2);
t78 = t127 * qJD(3) + t161 * t212;
t79 = -t126 * qJD(3) + t164 * t212;
t36 = t78 * pkin(3) - t79 * pkin(9) + t123;
t175 = t160 * t36 + t163 * t25 + t56 * t222 - t58 * t223;
t46 = t186 * pkin(3) - t52;
t174 = -pkin(9) * t73 + t98 * t46;
t171 = pkin(1) * (-qJD(2) * t203 + t245);
t26 = -pkin(3) * t213 - t180;
t6 = -t15 * qJD(4) - t160 * t21 + t163 * t29;
t170 = -t256 * qJD(4) - t160 * t25 + t163 * t36;
t168 = qJ(5) * t30 + t6;
t155 = -pkin(4) * t163 - pkin(3);
t145 = t260 * t163;
t144 = t260 * t160;
t136 = (pkin(4) * t160 + pkin(8)) * t161;
t133 = t163 * t143;
t114 = qJD(1) * t123;
t84 = -t160 * t246 + t233;
t81 = t127 * t163 - t216;
t74 = -qJ(5) * t241 + t133 + (-pkin(8) * t160 - pkin(4)) * t164;
t68 = t69 ^ 2;
t42 = -t80 * qJD(4) + t160 * t213 + t163 * t79;
t41 = -qJD(4) * t216 + t127 * t222 + t160 * t79 - t163 * t213;
t33 = -pkin(4) * t242 + t53;
t32 = pkin(4) * t80 + t57;
t23 = t209 + t46;
t16 = -qJ(5) * t80 + t256;
t12 = pkin(4) * t126 - qJ(5) * t81 + t208;
t11 = pkin(4) * t41 + t26;
t4 = -qJ(5) * t41 - qJD(5) * t80 + t175;
t3 = pkin(4) * t78 - qJ(5) * t42 - qJD(5) * t81 + t170;
t2 = -qJD(5) * t69 - t179;
t1 = -qJD(5) * t71 + t168 + t267;
t5 = [0, 0, 0, t197 * t273, -t232 * t197, t181 * t212, -t181 * t213, 0, -t114 * t247 - t123 * t182 + t162 * t188, -t113 * t247 - t122 * t182 + t165 * t188, t105 * t79 + t127 * t167, -t79 * t103 - t105 * t78 - t126 * t167 - t127 * t73, t105 * t213 + t127 * t191 - t167 * t243 - t186 * t79, t78 * t186 + (t73 * t165 + (-qJD(1) * t126 - t103) * t228) * t159, (-t156 * t229 - t200) * t228, -t180 * t186 + t123 * t103 + t115 * t73 + t114 * t126 + t86 * t78 + (-t205 * t165 + (qJD(1) * t204 + t52) * t228) * t159, t123 * t105 + t114 * t127 + t115 * t167 + t173 * t186 - t176 * t243 - t191 * t237 - t213 * t53 + t86 * t79, -t30 * t81 + t42 * t71, t30 * t80 - t31 * t81 - t41 * t71 - t42 * t69, -t126 * t30 + t42 * t98 + t71 * t78 + t73 * t81, -t126 * t31 - t41 * t98 - t69 * t78 - t73 * t80, t126 * t73 + t78 * t98, t6 * t126 + t14 * t78 + t170 * t98 + t208 * t73 + t22 * t80 + t26 * t69 + t57 * t31 + t46 * t41, t126 * t206 - t15 * t78 - t175 * t98 + t22 * t81 - t256 * t73 + t26 * t71 - t57 * t30 + t46 * t42, t1 * t126 + t11 * t69 + t12 * t73 + t23 * t41 + t3 * t98 + t31 * t32 + t7 * t78 + t8 * t80, -t10 * t78 + t11 * t71 - t126 * t2 - t16 * t73 + t23 * t42 - t30 * t32 - t4 * t98 + t8 * t81, -t1 * t81 - t10 * t41 + t12 * t30 - t16 * t31 - t2 * t80 - t3 * t71 - t4 * t69 - t42 * t7, t1 * t12 + t10 * t4 + t11 * t23 + t16 * t2 + t3 * t7 + t32 * t8; 0, 0, 0, -t245 * t273, t232 * t245, -t165 * t189, t162 * t189, 0, t121 * t182 + t162 * t171 - t185, pkin(7) * t191 + t118 * t182 + t165 * t171, -qJD(3) * t161 ^ 2 * t214 + ((qJD(3) * t182 + t190) * t161 - t186 * t105) * t164, -t161 * t73 + t164 * t167 + (t195 - t226) * t105 + t272 * t103, -t164 * t199 + (t164 * t201 + (t161 * qJD(2) - t105) * t162) * t231, t161 * t199 + (-t161 * t201 + (t164 * qJD(2) + t103) * t162) * t231, t200 * t230, -pkin(2) * t73 + t86 * t226 - t106 * t186 - t121 * t103 + (pkin(8) * t199 + t119 * t186 - t114) * t164 + (-t52 * t162 + (-pkin(8) * t228 - t165 * t86) * t161) * t231, -pkin(8) * t164 * t191 - pkin(2) * t167 - t121 * t105 + t114 * t161 + t214 * t53 - t272 * t86 + (-t219 - t236) * t186, -t30 * t241 + (-t161 * t223 + t163 * t224 - t92) * t71, t69 * t92 + t71 * t91 + (-t163 * t69 - t253) * t224 + (t254 - t163 * t31 + (t160 * t69 - t163 * t71) * qJD(4)) * t161, -t92 * t98 + (t30 + t218) * t164 + (-t186 * t71 - t217 + t251) * t161, t91 * t98 + (-t227 * t98 + t31) * t164 + (t186 * t69 + t178) * t161, -t161 * t186 * t98 - t164 * t73, t133 * t73 - t46 * t91 - t59 * t69 + ((-qJD(4) * t143 + t60) * t160 + t274) * t98 + (t46 * t227 - t6 + (qJD(3) * t69 + t178) * pkin(8)) * t164 + (pkin(8) * t31 - t14 * t186 + t22 * t160 + t222 * t46) * t161, -t233 * t73 - t59 * t71 - t46 * t92 + t270 * t98 + (t46 * t225 - t206 + (qJD(3) * t71 + t217) * pkin(8)) * t164 + (-t46 * t223 + t22 * t163 + t186 * t15 + (-t30 + t218) * pkin(8)) * t161, t136 * t31 - t23 * t91 + t73 * t74 - t259 * t98 + t250 * t69 + (t227 * t23 - t1) * t164 + (-t186 * t7 + t222 * t23 + t264) * t161, -t136 * t30 - t23 * t92 - t73 * t84 + t258 * t98 + t250 * t71 + (t225 * t23 + t2) * t164 + (t10 * t186 - t223 * t23 + t263) * t161, t10 * t91 + t30 * t74 - t31 * t84 + t7 * t92 + t259 * t71 + t258 * t69 + (-t10 * t160 - t163 * t7) * t224 + (-t1 * t163 - t160 * t2 + (-t10 * t163 + t160 * t7) * qJD(4)) * t161, t1 * t74 - t10 * t258 + t136 * t8 + t2 * t84 + t23 * t250 - t259 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105 * t103, -t103 ^ 2 + t105 ^ 2, -t103 * t150 + t183, -t105 * t150 - t184, t191, -t86 * t105 - t186 * t53 + t205, t86 * t103 - t186 * t52 + t176, t202 * t71 - t254, -t275 * t160 + t276 * t163, -t105 * t71 + t202 * t98 + t252, -t160 * t98 ^ 2 + t105 * t69 + t251, -t98 * t105, -pkin(3) * t31 - t14 * t105 - t53 * t69 - t62 * t98 - t193 * t163 + (t52 * t98 + t174) * t160, pkin(3) * t30 + t15 * t105 + t160 * t193 + t163 * t174 + t257 * t98 - t53 * t71, -t105 * t7 - t144 * t73 + t155 * t31 - t263 - t33 * t69 + t248 * t98 + (t103 * t23 + (t23 + t266) * qJD(4)) * t160, t23 * t240 + t10 * t105 - t145 * t73 - t155 * t30 + t264 - t33 * t71 - t249 * t98 + (pkin(4) * t253 + t163 * t23) * qJD(4), -t144 * t30 - t145 * t31 - t248 * t71 - t249 * t69 + (-t7 * t98 + t2) * t163 + (-t1 - t265) * t160, -t1 * t144 + t145 * t2 + t155 * t8 + t248 * t7 + (pkin(4) * t223 - t33) * t23 + t249 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71 * t69, -t68 + t269, -t30 + t262, -t31 + t261, t73, t15 * t98 - t46 * t71 + t6, t14 * t98 + t46 * t69 + t206, 0.2e1 * t267 + t265 + (-t209 - t23) * t71 + t168, -pkin(4) * t269 + t9 * t98 + (qJD(5) + t23) * t69 + t179, pkin(4) * t30 - t268 * t69, t268 * t10 + (-t23 * t71 + t1) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t275, t276, -t68 - t269, t10 * t69 + t7 * t71 + t8;];
tauc_reg = t5;
