% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:44
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPR3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:43:16
% EndTime: 2022-01-20 11:43:21
% DurationCPUTime: 2.11s
% Computational Cost: add. (4854->286), mult. (8387->381), div. (0->0), fcn. (5797->8), ass. (0->197)
t194 = cos(qJ(3));
t180 = t194 * qJD(4);
t191 = sin(qJ(3));
t260 = -qJ(4) - pkin(7);
t214 = qJD(3) * t260;
t142 = t191 * t214 + t180;
t143 = -t191 * qJD(4) + t194 * t214;
t188 = sin(pkin(9));
t189 = cos(pkin(9));
t154 = t188 * t194 + t189 * t191;
t195 = cos(qJ(2));
t253 = pkin(1) * qJD(1);
t223 = t195 * t253;
t249 = t188 * t142 - t189 * t143 - t154 * t223;
t241 = t189 * t194;
t243 = t188 * t191;
t153 = -t241 + t243;
t248 = -t189 * t142 - t188 * t143 - t153 * t223;
t146 = t154 * qJD(3);
t185 = qJD(1) + qJD(2);
t115 = t185 * t146;
t233 = qJD(3) * t191;
t219 = t185 * t233;
t163 = t188 * t219;
t232 = qJD(3) * t194;
t218 = t185 * t232;
t116 = t189 * t218 - t163;
t190 = sin(qJ(5));
t193 = cos(qJ(5));
t221 = t185 * t241;
t123 = t185 * t243 - t221;
t125 = t154 * t185;
t204 = t190 * t123 - t193 * t125;
t197 = t204 * qJD(5) - t193 * t115 - t190 * t116;
t184 = qJD(3) + qJD(5);
t251 = t204 * t184;
t278 = t197 - t251;
t229 = qJD(5) * t193;
t230 = qJD(5) * t190;
t202 = -t190 * t115 + t193 * t116 - t123 * t229 - t125 * t230;
t66 = -t193 * t123 - t190 * t125;
t250 = t66 * t184;
t277 = t202 - t250;
t262 = t204 ^ 2;
t263 = t66 ^ 2;
t276 = t262 - t263;
t261 = t66 * t204;
t147 = t153 * qJD(3);
t266 = t147 * pkin(8);
t275 = t266 - t249;
t139 = t146 * pkin(8);
t274 = -t139 - t248;
t186 = t191 ^ 2;
t187 = t194 ^ 2;
t234 = t186 + t187;
t267 = t125 * pkin(8);
t192 = sin(qJ(2));
t224 = t192 * t253;
t161 = t185 * pkin(7) + t224;
t211 = qJ(4) * t185 + t161;
t111 = t211 * t191;
t102 = qJD(3) * pkin(3) - t111;
t112 = t211 * t194;
t96 = t188 * t112;
t55 = t189 * t102 - t96;
t37 = qJD(3) * pkin(4) - t267 + t55;
t268 = t123 * pkin(8);
t242 = t189 * t112;
t56 = t188 * t102 + t242;
t38 = t56 - t268;
t15 = t190 * t37 + t193 * t38;
t252 = pkin(1) * qJD(2);
t222 = qJD(1) * t252;
t207 = t195 * t222;
t200 = qJD(4) * t185 + t207;
t203 = qJD(3) * t211;
t71 = -t191 * t203 + t200 * t194;
t72 = -t200 * t191 - t194 * t203;
t30 = -t188 * t71 + t189 * t72;
t20 = -t116 * pkin(8) + t30;
t31 = t188 * t72 + t189 * t71;
t21 = -t115 * pkin(8) + t31;
t5 = -t15 * qJD(5) - t190 * t21 + t193 * t20;
t177 = -t194 * pkin(3) - pkin(2);
t122 = t177 * t185 + qJD(4) - t223;
t73 = t123 * pkin(4) + t122;
t273 = t204 * t73 + t5;
t4 = (qJD(5) * t37 + t21) * t193 + t190 * t20 - t38 * t230;
t272 = -t73 * t66 - t4;
t271 = t125 ^ 2;
t270 = pkin(3) * t188;
t269 = pkin(3) * t191;
t265 = t154 * pkin(8);
t264 = t195 * pkin(1);
t169 = t260 * t191;
t182 = t194 * qJ(4);
t170 = t194 * pkin(7) + t182;
t100 = t189 * t169 - t188 * t170;
t76 = t100 - t265;
t101 = t188 * t169 + t189 * t170;
t150 = t153 * pkin(8);
t77 = -t150 + t101;
t32 = -t190 * t77 + t193 * t76;
t259 = t32 * qJD(5) + t275 * t190 + t274 * t193;
t33 = t190 * t76 + t193 * t77;
t258 = -t33 * qJD(5) - t274 * t190 + t275 * t193;
t84 = -t190 * t153 + t193 * t154;
t46 = t84 * qJD(5) + t193 * t146 - t190 * t147;
t172 = t192 * t222;
t144 = pkin(3) * t219 + t172;
t78 = t115 * pkin(4) + t144;
t83 = t193 * t153 + t190 * t154;
t257 = t73 * t46 + t78 * t83;
t45 = t190 * t146 + t193 * t147 + t153 * t229 + t154 * t230;
t256 = -t73 * t45 + t78 * t84;
t255 = t122 * t146 + t144 * t153;
t254 = -t122 * t147 + t144 * t154;
t175 = t192 * pkin(1) + pkin(7);
t238 = -qJ(4) - t175;
t209 = qJD(3) * t238;
t225 = t195 * t252;
t91 = t191 * t209 + t194 * t225 + t180;
t92 = (-qJD(4) - t225) * t191 + t194 * t209;
t52 = t188 * t92 + t189 * t91;
t58 = -t189 * t111 - t96;
t174 = t189 * pkin(3) + pkin(4);
t140 = t193 * t174 - t190 * t270;
t57 = t188 * t111 - t242;
t41 = t57 + t268;
t42 = t58 - t267;
t247 = t140 * qJD(5) - t190 * t41 - t193 * t42;
t141 = t190 * t174 + t193 * t270;
t246 = -t141 * qJD(5) + t190 * t42 - t193 * t41;
t245 = t125 * t123;
t244 = t185 * t191;
t240 = t192 * t194;
t196 = qJD(3) ^ 2;
t239 = t196 * t191;
t181 = t196 * t194;
t151 = t238 * t191;
t152 = t194 * t175 + t182;
t82 = t188 * t151 + t189 * t152;
t162 = -t185 * pkin(2) - t223;
t237 = t162 * t232 + t191 * t172;
t236 = t234 * t207;
t235 = t186 - t187;
t231 = qJD(3) * t195;
t228 = -qJD(1) - t185;
t227 = -qJD(2) + t185;
t226 = pkin(3) * t244;
t179 = t192 * t252;
t178 = pkin(3) * t233;
t183 = t185 ^ 2;
t220 = t191 * t183 * t194;
t217 = t191 * t231;
t14 = -t190 * t38 + t193 * t37;
t216 = t14 * t45 - t15 * t46 - t4 * t83 - t5 * t84;
t108 = t146 * pkin(4) + t178;
t51 = -t188 * t91 + t189 * t92;
t213 = -t56 * t146 + t55 * t147 - t31 * t153 - t30 * t154;
t81 = t189 * t151 - t188 * t152;
t208 = t234 * qJD(2);
t206 = t191 * t218;
t205 = t108 - t224;
t61 = t81 - t265;
t62 = -t150 + t82;
t24 = -t190 * t62 + t193 * t61;
t25 = t190 * t61 + t193 * t62;
t118 = t153 * pkin(4) + t177;
t201 = -t162 * t185 - t207;
t199 = -t192 * t244 + t194 * t231;
t176 = -pkin(2) - t264;
t164 = t177 - t264;
t160 = -0.2e1 * t206;
t159 = 0.2e1 * t206;
t158 = t179 + t178;
t148 = t162 * t233;
t134 = t147 * qJD(3);
t133 = t146 * qJD(3);
t129 = -0.2e1 * t235 * t185 * qJD(3);
t117 = t123 ^ 2;
t107 = t118 - t264;
t95 = t108 + t179;
t90 = t125 * pkin(4) + t226;
t48 = t116 * t154 - t125 * t147;
t47 = t115 * t153 + t123 * t146;
t40 = t46 * t184;
t39 = t45 * t184;
t35 = -t139 + t52;
t34 = t51 + t266;
t18 = -t154 * t115 - t116 * t153 + t147 * t123 - t125 * t146;
t9 = -t197 * t83 - t46 * t66;
t8 = t202 * t84 + t204 * t45;
t7 = -t25 * qJD(5) - t190 * t35 + t193 * t34;
t6 = t24 * qJD(5) + t190 * t34 + t193 * t35;
t1 = t197 * t84 - t202 * t83 + t204 * t46 - t45 * t66;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t185 * t179 - t172, t228 * t225, 0, 0, t159, t129, t181, t160, -t239, 0, t176 * t219 - t175 * t181 + t148 + (t228 * t240 - t217) * t252, t175 * t239 + t176 * t218 - t199 * t252 + t237, t185 * t208 * t264 + t236, ((qJD(1) * t176 + t162) * t192 + (qJD(1) * t175 + t161) * t195 * t234) * t252, t48, t18, -t134, t47, -t133, 0, t51 * qJD(3) + t164 * t115 + t158 * t123 + t255, -t52 * qJD(3) + t164 * t116 + t158 * t125 + t254, -t82 * t115 - t81 * t116 - t52 * t123 - t51 * t125 + t213, t122 * t158 + t144 * t164 + t30 * t81 + t31 * t82 + t55 * t51 + t56 * t52, t8, t1, -t39, t9, -t40, 0, -t107 * t197 + t7 * t184 - t66 * t95 + t257, t107 * t202 - t6 * t184 - t204 * t95 + t256, t197 * t25 - t202 * t24 + t204 * t7 + t6 * t66 + t216, t78 * t107 + t14 * t7 + t15 * t6 + t5 * t24 + t4 * t25 + t73 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t185 * t224 - t172, t227 * t223, 0, 0, t159, t129, t181, t160, -t239, 0, -pkin(2) * t219 - pkin(7) * t181 + t148 + (t227 * t240 + t217) * t253, -pkin(2) * t218 + pkin(7) * t239 + t199 * t253 + t237, -t185 * t223 * t234 + t236, ((-pkin(2) * qJD(2) - t162) * t192 + (pkin(7) * t208 - t161 * t234) * t195) * t253, t48, t18, -t134, t47, -t133, 0, -t123 * t224 + t177 * t115 + (t123 * t269 - t249) * qJD(3) + t255, -t125 * t224 + t177 * t116 + (t125 * t269 + t248) * qJD(3) + t254, -t100 * t116 - t101 * t115 + t123 * t248 + t125 * t249 + t213, t30 * t100 + t31 * t101 + t144 * t177 - t248 * t56 - t249 * t55 + (-t224 + t178) * t122, t8, t1, -t39, t9, -t40, 0, -t118 * t197 + t184 * t258 - t205 * t66 + t257, t118 * t202 - t184 * t259 - t204 * t205 + t256, t197 * t33 - t202 * t32 + t204 * t258 + t259 * t66 + t216, t78 * t118 + t14 * t258 + t15 * t259 + t205 * t73 + t5 * t32 + t4 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t220, t235 * t183, 0, t220, 0, 0, t201 * t191, t201 * t194, 0, 0, t245, -t117 + t271, -t163 + (t123 + t221) * qJD(3), -t245, 0, 0, -t57 * qJD(3) - t122 * t125 - t123 * t226 + t30, t58 * qJD(3) + t122 * t123 - t125 * t226 - t31, (t56 + t57) * t125 + (-t55 + t58) * t123 + (-t115 * t188 - t116 * t189) * pkin(3), -t55 * t57 - t56 * t58 + (-t122 * t244 + t188 * t31 + t189 * t30) * pkin(3), t261, t276, t277, -t261, t278, 0, t184 * t246 + t66 * t90 + t273, -t184 * t247 + t204 * t90 + t272, -t140 * t202 + t141 * t197 + (t14 + t247) * t66 + (-t15 + t246) * t204, t14 * t246 + t5 * t140 + t4 * t141 + t15 * t247 - t73 * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t125 * qJD(3), -t163 + (-t123 + t221) * qJD(3), -t117 - t271, t56 * t123 + t55 * t125 + t144, 0, 0, 0, 0, 0, 0, -t197 - t251, t202 + t250, -t262 - t263, -t14 * t204 - t15 * t66 + t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t261, t276, t277, -t261, t278, 0, t15 * t184 + t273, t14 * t184 + t272, 0, 0;];
tauc_reg = t2;
