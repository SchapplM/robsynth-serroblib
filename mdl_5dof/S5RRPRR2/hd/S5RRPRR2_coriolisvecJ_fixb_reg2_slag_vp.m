% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRR2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:28:27
% EndTime: 2019-12-05 18:28:37
% DurationCPUTime: 3.32s
% Computational Cost: add. (7243->312), mult. (19116->419), div. (0->0), fcn. (14526->8), ass. (0->184)
t192 = cos(qJ(5));
t189 = sin(qJ(5));
t221 = qJD(5) * t189;
t187 = sin(pkin(9));
t188 = cos(pkin(9));
t194 = cos(qJ(2));
t231 = t188 * t194;
t215 = qJD(1) * t231;
t191 = sin(qJ(2));
t224 = qJD(1) * t191;
t153 = t187 * t224 - t215;
t190 = sin(qJ(4));
t193 = cos(qJ(4));
t165 = t187 * t194 + t188 * t191;
t226 = qJD(1) * t165;
t105 = -t193 * t153 - t190 * t226;
t268 = t105 * pkin(8);
t248 = -qJ(3) - pkin(6);
t173 = t248 * t194;
t170 = qJD(1) * t173;
t159 = t187 * t170;
t172 = t248 * t191;
t169 = qJD(1) * t172;
t245 = qJD(2) * pkin(2);
t163 = t169 + t245;
t109 = t188 * t163 + t159;
t252 = t226 * pkin(7);
t78 = qJD(2) * pkin(3) + t109 - t252;
t232 = t188 * t170;
t110 = t187 * t163 - t232;
t253 = t153 * pkin(7);
t83 = t110 - t253;
t35 = t190 * t78 + t193 * t83;
t29 = t35 + t268;
t219 = qJD(1) * qJD(2);
t212 = t191 * t219;
t175 = t187 * t212;
t211 = t194 * t219;
t142 = t188 * t211 - t175;
t208 = qJD(2) * t248;
t149 = t194 * qJD(3) + t191 * t208;
t127 = t149 * qJD(1);
t150 = -t191 * qJD(3) + t194 * t208;
t128 = t150 * qJD(1);
t84 = -t187 * t127 + t188 * t128;
t66 = -t142 * pkin(7) + t84;
t155 = t165 * qJD(2);
t141 = qJD(1) * t155;
t85 = t188 * t127 + t187 * t128;
t67 = -t141 * pkin(7) + t85;
t19 = -t35 * qJD(4) - t190 * t67 + t193 * t66;
t222 = qJD(4) * t193;
t223 = qJD(4) * t190;
t200 = -t190 * t141 + t193 * t142 - t153 * t222 - t223 * t226;
t7 = -pkin(8) * t200 + t19;
t213 = t189 * t7 - t29 * t221;
t184 = qJD(2) + qJD(4);
t201 = t190 * t153 - t193 * t226;
t269 = pkin(8) * t201;
t34 = -t190 * t83 + t193 * t78;
t28 = t34 + t269;
t26 = t184 * pkin(4) + t28;
t197 = t201 * qJD(4) - t193 * t141 - t190 * t142;
t207 = -t190 * t66 - t193 * t67 - t78 * t222 + t83 * t223;
t6 = pkin(8) * t197 - t207;
t1 = (qJD(5) * t26 + t6) * t192 + t213;
t54 = t192 * t105 + t189 * t201;
t180 = -t194 * pkin(2) - pkin(1);
t225 = qJD(1) * t180;
t171 = qJD(3) + t225;
t115 = t153 * pkin(3) + t171;
t68 = -pkin(4) * t105 + t115;
t267 = t68 * t54;
t275 = -t1 - t267;
t216 = -t189 * t6 + t192 * t7;
t243 = t192 * t29;
t9 = t189 * t26 + t243;
t2 = -t9 * qJD(5) + t216;
t50 = t189 * t105 - t192 * t201;
t266 = t68 * t50;
t274 = t2 - t266;
t220 = qJD(5) * t192;
t16 = -t105 * t220 - t189 * t197 - t192 * t200 - t201 * t221;
t183 = qJD(5) + t184;
t241 = t54 * t183;
t262 = -t16 - t241;
t198 = -qJD(5) * t50 - t189 * t200 + t192 * t197;
t242 = t50 * t183;
t259 = t198 + t242;
t250 = t50 ^ 2;
t251 = t54 ^ 2;
t263 = t250 - t251;
t249 = t54 * t50;
t244 = t189 * t29;
t8 = t192 * t26 - t244;
t273 = t9 * t50 + t8 * t54;
t179 = t188 * pkin(2) + pkin(3);
t256 = pkin(2) * t187;
t147 = t193 * t179 - t190 * t256;
t113 = -t187 * t169 + t232;
t86 = t113 + t253;
t114 = t188 * t169 + t159;
t87 = t114 - t252;
t240 = t147 * qJD(4) - t190 * t86 - t193 * t87;
t148 = t190 * t179 + t193 * t256;
t239 = -t148 * qJD(4) + t190 * t87 - t193 * t86;
t238 = t201 * t184;
t272 = t197 - t238;
t236 = t105 * t184;
t271 = t200 - t236;
t235 = t201 ^ 2;
t237 = t105 ^ 2;
t270 = t235 - t237;
t254 = pkin(4) * t201;
t265 = -t269 + t240;
t264 = -t268 - t239;
t234 = t105 * t201;
t261 = t115 * t201 + t19;
t260 = -t115 * t105 + t207;
t258 = -0.2e1 * t219;
t257 = t226 ^ 2;
t255 = pkin(2) * t191;
t145 = pkin(4) + t147;
t94 = t189 * t145 + t192 * t148;
t247 = t94 * qJD(5) + t265 * t189 + t264 * t192;
t93 = t192 * t145 - t189 * t148;
t246 = -t93 * qJD(5) + t264 * t189 - t265 * t192;
t117 = t188 * t172 + t187 * t173;
t100 = -t165 * pkin(7) + t117;
t118 = t187 * t172 - t188 * t173;
t164 = t187 * t191 - t231;
t101 = -t164 * pkin(7) + t118;
t43 = t190 * t100 + t193 * t101;
t233 = t226 * t153;
t196 = qJD(1) ^ 2;
t230 = t194 * t196;
t195 = qJD(2) ^ 2;
t229 = t195 * t191;
t228 = t195 * t194;
t99 = t188 * t149 + t187 * t150;
t227 = t191 ^ 2 - t194 ^ 2;
t182 = t191 * t245;
t181 = pkin(2) * t224;
t217 = t191 * t230;
t214 = -pkin(4) * t183 - t26;
t177 = pkin(2) * t212;
t116 = t141 * pkin(3) + t177;
t123 = t155 * pkin(3) + t182;
t122 = pkin(3) * t226 + t181;
t42 = t193 * t100 - t190 * t101;
t206 = pkin(1) * t258;
t98 = -t187 * t149 + t188 * t150;
t204 = 0.2e1 * t226;
t203 = t191 * t211;
t112 = -t190 * t164 + t193 * t165;
t32 = -t112 * pkin(8) + t42;
t111 = t193 * t164 + t190 * t165;
t33 = -t111 * pkin(8) + t43;
t20 = -t189 * t33 + t192 * t32;
t21 = t189 * t32 + t192 * t33;
t62 = -t189 * t111 + t192 * t112;
t130 = t164 * pkin(3) + t180;
t36 = -pkin(4) * t197 + t116;
t158 = t164 * qJD(2);
t74 = t158 * pkin(7) + t98;
t75 = -t155 * pkin(7) + t99;
t24 = t100 * t222 - t101 * t223 + t190 * t74 + t193 * t75;
t25 = -t43 * qJD(4) - t190 * t75 + t193 * t74;
t151 = t153 ^ 2;
t79 = t111 * pkin(4) + t130;
t69 = t122 - t254;
t61 = t192 * t111 + t189 * t112;
t60 = t112 * qJD(4) + t193 * t155 - t190 * t158;
t59 = t190 * t155 + t193 * t158 + t164 * t222 + t165 * t223;
t39 = t60 * pkin(4) + t123;
t23 = t62 * qJD(5) - t189 * t59 + t192 * t60;
t22 = t111 * t220 + t112 * t221 + t189 * t60 + t192 * t59;
t15 = t59 * pkin(8) + t25;
t14 = -t60 * pkin(8) + t24;
t11 = t192 * t28 - t244;
t10 = -t189 * t28 - t243;
t4 = -t21 * qJD(5) - t189 * t14 + t192 * t15;
t3 = t20 * qJD(5) + t192 * t14 + t189 * t15;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t203, t227 * t258, t228, -0.2e1 * t203, -t229, 0, -pkin(6) * t228 + t191 * t206, pkin(6) * t229 + t194 * t206, 0, 0, t142 * t165 - t158 * t226, -t165 * t141 - t142 * t164 + t158 * t153 - t155 * t226, -t158 * qJD(2), t141 * t164 + t153 * t155, -t155 * qJD(2), 0, t180 * t141 + t171 * t155 + (t98 + (qJD(1) * t164 + t153) * t255) * qJD(2), t180 * t142 - t171 * t158 + (t204 * t255 - t99) * qJD(2), t109 * t158 - t110 * t155 - t117 * t142 - t118 * t141 - t99 * t153 - t85 * t164 - t84 * t165 - t226 * t98, t109 * t98 + t110 * t99 + t84 * t117 + t85 * t118 + (t171 + t225) * t182, t112 * t200 + t201 * t59, -t105 * t59 - t111 * t200 + t112 * t197 + t201 * t60, -t59 * t184, -t105 * t60 - t111 * t197, -t60 * t184, 0, -t105 * t123 + t116 * t111 + t115 * t60 - t130 * t197 + t25 * t184, t116 * t112 - t115 * t59 - t123 * t201 + t130 * t200 - t24 * t184, t105 * t24 + t111 * t207 - t19 * t112 + t197 * t43 - t200 * t42 + t201 * t25 + t34 * t59 - t35 * t60, t115 * t123 + t116 * t130 + t19 * t42 - t207 * t43 + t35 * t24 + t34 * t25, -t16 * t62 - t22 * t50, t16 * t61 + t198 * t62 - t22 * t54 - t23 * t50, -t22 * t183, -t198 * t61 - t23 * t54, -t23 * t183, 0, t4 * t183 - t198 * t79 + t68 * t23 + t36 * t61 - t39 * t54, -t79 * t16 - t3 * t183 - t68 * t22 + t36 * t62 + t39 * t50, -t1 * t61 + t20 * t16 + t198 * t21 - t2 * t62 + t8 * t22 - t9 * t23 + t3 * t54 - t4 * t50, t1 * t21 + t2 * t20 + t9 * t3 + t36 * t79 + t68 * t39 + t8 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t217, t227 * t196, 0, t217, 0, 0, t196 * pkin(1) * t191, pkin(1) * t230, 0, 0, t233, -t151 + t257, -t175 + (t153 + t215) * qJD(2), -t233, 0, 0, -t113 * qJD(2) - t153 * t181 - t171 * t226 + t84, t114 * qJD(2) + t171 * t153 - t181 * t226 - t85, (t110 + t113) * t226 + (-t109 + t114) * t153 + (-t141 * t187 - t142 * t188) * pkin(2), -t109 * t113 - t110 * t114 + (-t171 * t224 + t187 * t85 + t188 * t84) * pkin(2), t234, t270, t271, -t234, t272, 0, t105 * t122 + t239 * t184 + t261, t122 * t201 - t240 * t184 + t260, -t147 * t200 + t148 * t197 + (t239 - t35) * t201 + (t240 + t34) * t105, -t115 * t122 + t19 * t147 - t148 * t207 + t239 * t34 + t240 * t35, -t249, t263, t262, t249, t259, 0, -t247 * t183 + t54 * t69 + t274, t246 * t183 - t50 * t69 + t275, t93 * t16 + t198 * t94 - t246 * t54 + t247 * t50 + t273, t1 * t94 + t2 * t93 - t246 * t9 - t247 * t8 - t68 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t204 * qJD(2), -t175 + (-t153 + t215) * qJD(2), -t151 - t257, t109 * t226 + t110 * t153 + t177, 0, 0, 0, 0, 0, 0, -t197 - t238, t200 + t236, -t235 - t237, -t35 * t105 - t201 * t34 + t116, 0, 0, 0, 0, 0, 0, -t198 + t242, -t16 + t241, -t250 - t251, t8 * t50 - t9 * t54 + t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t234, t270, t271, -t234, t272, 0, t35 * t184 + t261, t34 * t184 + t260, 0, 0, -t249, t263, t262, t249, t259, 0, -t54 * t254 - t10 * t183 - t266 + (t189 * t214 - t243) * qJD(5) + t216, t50 * t254 + t11 * t183 - t267 + (qJD(5) * t214 - t6) * t192 - t213, t10 * t50 - t11 * t54 + (t16 * t192 + t198 * t189 + (t189 * t50 + t192 * t54) * qJD(5)) * pkin(4) + t273, -t8 * t10 - t9 * t11 + (t1 * t189 + t201 * t68 + t192 * t2 + (-t189 * t8 + t192 * t9) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t249, t263, t262, t249, t259, 0, t9 * t183 + t274, t8 * t183 + t275, 0, 0;];
tauc_reg = t5;
