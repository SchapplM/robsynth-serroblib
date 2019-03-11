% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% 
% Output:
% tauc_reg [6x35]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPRRR6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:53:25
% EndTime: 2019-03-09 13:53:38
% DurationCPUTime: 4.19s
% Computational Cost: add. (5554->358), mult. (13040->477), div. (0->0), fcn. (9496->8), ass. (0->205)
t183 = cos(qJ(6));
t241 = qJD(6) * t183;
t181 = sin(qJ(4));
t185 = cos(qJ(4));
t186 = cos(qJ(2));
t249 = qJD(1) * t186;
t182 = sin(qJ(2));
t250 = qJD(1) * t182;
t114 = -t181 * t250 - t185 * t249;
t116 = -t181 * t249 + t185 * t250;
t180 = sin(qJ(5));
t184 = cos(qJ(5));
t297 = -t184 * t114 + t116 * t180;
t311 = t183 * t297;
t332 = t241 + t311;
t164 = pkin(7) * t250;
t331 = -pkin(8) * t250 + qJD(3) + t164;
t173 = qJD(2) - qJD(4);
t168 = -qJD(5) + t173;
t243 = qJD(5) * t184;
t244 = qJD(5) * t180;
t129 = t181 * t182 + t185 * t186;
t198 = t129 * qJD(4);
t240 = qJD(1) * qJD(2);
t232 = t186 * t240;
t233 = t182 * t240;
t81 = -qJD(1) * t198 + t181 * t233 + t185 * t232;
t245 = qJD(4) * t185;
t246 = qJD(4) * t181;
t247 = qJD(2) * t186;
t324 = t181 * t247 + t182 * t245 - t186 * t246;
t82 = t324 * qJD(1) - t185 * t233;
t201 = -t114 * t243 + t116 * t244 + t180 * t82 - t184 * t81;
t208 = t114 * t180 + t184 * t116;
t179 = sin(qJ(6));
t242 = qJD(6) * t179;
t14 = -t168 * t241 - t183 * t201 - t208 * t242;
t211 = t168 * t179 - t183 * t208;
t15 = -qJD(6) * t211 - t179 * t201;
t307 = qJD(6) + t297;
t330 = t179 * t307;
t61 = t183 * t168 + t179 * t208;
t329 = t14 * t183 - t179 * t15 + t211 * t330 - t332 * t61;
t34 = qJD(5) * t208 + t180 * t81 + t184 * t82;
t30 = t183 * t34;
t328 = -t208 * t61 + t307 * t330 - t30;
t165 = pkin(7) * t249;
t136 = -pkin(8) * t249 + t165;
t187 = -pkin(2) - pkin(3);
t210 = -qJ(3) * t181 + t185 * t187;
t326 = qJD(4) * t210 - t181 * t136 + t185 * t331;
t139 = qJ(3) * t185 + t181 * t187;
t325 = qJD(4) * t139 + t185 * t136 + t181 * t331;
t47 = pkin(5) * t208 + pkin(10) * t297;
t323 = -t208 ^ 2 + t297 ^ 2;
t322 = qJD(6) - t307;
t12 = t14 * t179;
t320 = t211 * t332 - t12;
t28 = t179 * t34;
t66 = t307 * t241;
t278 = t28 + t66;
t319 = t208 * t211 + t307 * t311 + t278;
t236 = t187 * qJD(2);
t102 = t236 + t331;
t175 = qJD(2) * qJ(3);
t117 = t136 + t175;
t209 = t102 * t181 + t117 * t185;
t285 = pkin(9) * t114;
t57 = t209 + t285;
t272 = t180 * t57;
t226 = t185 * t102 - t117 * t181;
t284 = pkin(9) * t116;
t56 = t226 - t284;
t53 = -pkin(4) * t173 + t56;
t20 = t184 * t53 - t272;
t18 = pkin(5) * t168 - t20;
t316 = t18 * t297;
t315 = -t284 + t326;
t314 = t285 + t325;
t313 = t114 * t173 + t81;
t310 = t208 * t297;
t299 = t307 * t208;
t118 = -qJD(1) * pkin(1) - pkin(2) * t249 - qJ(3) * t250;
t100 = pkin(3) * t249 - t118;
t248 = qJD(2) * t182;
t287 = pkin(7) - pkin(8);
t135 = t287 * t248;
t174 = qJD(2) * qJD(3);
t107 = -qJD(1) * t135 + t174;
t157 = pkin(7) * t232;
t125 = -pkin(8) * t232 + t157;
t193 = -t209 * qJD(4) - t181 * t107 + t185 * t125;
t309 = t100 * t116 - t193;
t202 = -t102 * t245 - t185 * t107 + t117 * t246 - t181 * t125;
t31 = -pkin(9) * t82 - t202;
t32 = -t81 * pkin(9) + t193;
t3 = t180 * t31 - t184 * t32 + t57 * t243 + t53 * t244;
t80 = -pkin(4) * t114 + t100;
t293 = -t208 * t80 - t3;
t268 = t184 * t57;
t21 = t180 * t53 + t268;
t19 = -pkin(10) * t168 + t21;
t35 = pkin(5) * t297 - pkin(10) * t208 + t80;
t10 = t179 * t35 + t183 * t19;
t291 = t10 * t208 + t3 * t179 + t18 * t241;
t212 = t179 * t19 - t183 * t35;
t294 = t18 * t242 + t212 * t208;
t305 = t168 * t208 + t34;
t228 = -t180 * t32 - t184 * t31 - t53 * t243 + t57 * t244;
t304 = t297 * t80 + t228;
t303 = t168 * t297 + t201;
t301 = -0.2e1 * t240;
t298 = t116 * t173 + t82;
t130 = -t181 * t186 + t182 * t185;
t85 = t184 * t129 + t130 * t180;
t91 = -t185 * t248 + t324;
t92 = qJD(2) * t129 - t198;
t39 = -qJD(5) * t85 - t180 * t91 + t184 * t92;
t86 = -t129 * t180 + t130 * t184;
t141 = -t186 * pkin(2) - t182 * qJ(3) - pkin(1);
t126 = t186 * pkin(3) - t141;
t93 = pkin(4) * t129 + t126;
t42 = pkin(5) * t85 - pkin(10) * t86 + t93;
t144 = t287 * t182;
t145 = t287 * t186;
t67 = -pkin(9) * t130 + t144 * t185 - t145 * t181;
t205 = -t144 * t181 - t145 * t185;
t68 = -pkin(9) * t129 - t205;
t46 = t180 * t67 + t184 * t68;
t137 = qJD(2) * t145;
t199 = -t185 * t135 + t181 * t137 + t144 * t245 - t145 * t246;
t43 = -pkin(9) * t91 + t199;
t191 = qJD(4) * t205 + t181 * t135 + t185 * t137;
t44 = -t92 * pkin(9) + t191;
t45 = t180 * t68 - t184 * t67;
t6 = -qJD(5) * t45 + t180 * t44 + t184 * t43;
t288 = t18 * t39 - (qJD(6) * t42 + t6) * t307 - (qJD(6) * t35 - t228) * t85 + t3 * t86 - t46 * t34;
t283 = t116 * pkin(4);
t282 = t18 * t86;
t281 = t34 * t86;
t280 = t42 * t34;
t133 = -pkin(4) + t210;
t207 = t133 * t184 - t139 * t180;
t277 = -qJD(5) * t207 + t314 * t180 - t315 * t184;
t206 = t133 * t180 + t139 * t184;
t276 = qJD(5) * t206 + t315 * t180 + t314 * t184;
t275 = qJD(2) * pkin(2);
t128 = t180 * t181 - t184 * t185;
t267 = t168 * t128;
t131 = t180 * t185 + t181 * t184;
t266 = t168 * t131;
t263 = t116 * t114;
t189 = qJD(1) ^ 2;
t260 = t186 * t189;
t188 = qJD(2) ^ 2;
t259 = t188 * t182;
t258 = t188 * t186;
t169 = t182 * qJD(3);
t253 = qJ(3) * t232 + qJD(1) * t169;
t252 = qJ(3) * t247 + t169;
t176 = t182 ^ 2;
t251 = -t186 ^ 2 + t176;
t238 = t86 * t242;
t237 = t182 * t260;
t234 = t114 ^ 2 - t116 ^ 2;
t160 = qJ(3) * t249;
t110 = t187 * t250 + t160;
t83 = t110 - t283;
t88 = -pkin(10) + t206;
t227 = qJD(6) * t88 - t47 + t83;
t222 = pkin(1) * t301;
t221 = qJD(3) - t275;
t162 = pkin(4) * t180 + pkin(10);
t220 = qJD(6) * t162 + t283 + t47;
t217 = qJD(1) * t141 + t118;
t216 = t173 ^ 2;
t215 = t182 * t236;
t22 = t180 * t56 + t268;
t214 = pkin(4) * t244 - t22;
t213 = t307 * t39 + t281;
t204 = qJD(6) * t131 + t250;
t111 = pkin(2) * t248 - t252;
t99 = pkin(2) * t233 - t253;
t203 = -pkin(7) * t188 - qJD(1) * t111 - t99;
t96 = t215 + t252;
t94 = qJD(1) * t215 + t253;
t196 = t277 * t307 - t88 * t34 - t316;
t60 = t91 * pkin(4) + t96;
t195 = -t100 * t114 + t202;
t52 = t82 * pkin(4) + t94;
t23 = t184 * t56 - t272;
t192 = -t162 * t34 + t316 + (-pkin(4) * t243 + t23) * t307;
t138 = -pkin(7) * t233 + t174;
t140 = t164 + t221;
t142 = t165 + t175;
t190 = t138 * t186 + (t140 * t186 + (-t142 + t165) * t182) * qJD(2);
t163 = -pkin(4) * t184 - pkin(5);
t132 = pkin(2) * t250 - t160;
t87 = pkin(5) - t207;
t40 = qJD(5) * t86 + t180 * t92 + t184 * t91;
t8 = pkin(5) * t40 - pkin(10) * t39 + t60;
t7 = qJD(5) * t46 + t180 * t43 - t184 * t44;
t5 = pkin(5) * t34 + pkin(10) * t201 + t52;
t4 = t183 * t5;
t1 = [0, 0, 0, 0.2e1 * t182 * t232, t251 * t301, t258, -t259, 0, -pkin(7) * t258 + t182 * t222, pkin(7) * t259 + t186 * t222, t186 * t203 + t217 * t248, t190, t182 * t203 - t217 * t247, pkin(7) * t190 + t111 * t118 + t141 * t99, t116 * t92 + t130 * t81, t114 * t92 - t116 * t91 - t129 * t81 - t130 * t82, -t92 * t173, t91 * t173, 0, t100 * t91 - t96 * t114 + t126 * t82 + t94 * t129 - t173 * t191, t100 * t92 + t96 * t116 + t126 * t81 + t94 * t130 + t173 * t199, -t201 * t86 + t208 * t39, t201 * t85 - t208 * t40 - t297 * t39 - t281, -t39 * t168, t40 * t168, 0, t168 * t7 + t297 * t60 + t34 * t93 + t40 * t80 + t52 * t85, t168 * t6 - t201 * t93 + t208 * t60 + t39 * t80 + t52 * t86, t211 * t238 + (t14 * t86 - t211 * t39) * t183 (t179 * t211 - t183 * t61) * t39 + (-t12 - t15 * t183 + (t179 * t61 + t183 * t211) * qJD(6)) * t86, t14 * t85 + t183 * t213 - t211 * t40 - t238 * t307, -t15 * t85 - t179 * t213 - t61 * t40 - t66 * t86, t307 * t40 + t34 * t85, t45 * t15 + t4 * t85 - t212 * t40 + t7 * t61 + (t280 + t8 * t307 + (-t19 * t85 - t307 * t46 + t282) * qJD(6)) * t183 + t288 * t179, -t10 * t40 + t45 * t14 - t7 * t211 + (-(-qJD(6) * t46 + t8) * t307 - t280 - (-qJD(6) * t19 + t5) * t85 - qJD(6) * t282) * t179 + t288 * t183; 0, 0, 0, -t237, t251 * t189, 0, 0, 0, t189 * pkin(1) * t182, pkin(1) * t260 (-t118 * t182 + t132 * t186) * qJD(1) ((t142 - t175) * t182 + (-t140 + t221) * t186) * qJD(1), 0.2e1 * t174 + (t118 * t186 + t132 * t182) * qJD(1), qJ(3) * t138 + qJD(3) * t142 - t118 * t132 + (t142 * t182 + (-t140 - t275) * t186) * qJD(1) * pkin(7), t263, t234, -t313, t298, 0, t110 * t114 + t325 * t173 + t309, -t110 * t116 + t326 * t173 - t195, -t310, t323, t303, t305, 0, t168 * t276 - t297 * t83 - t293, -t168 * t277 - t208 * t83 - t304, t320, -t329, -t319, t328, t299, t87 * t15 + t276 * t61 + (-t227 * t307 + t3) * t183 + t196 * t179 - t294, t87 * t14 + t183 * t196 - t211 * t276 + t227 * t330 - t291; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t237, 0, -t176 * t189 - t188, -qJD(2) * t142 + t118 * t250 + t157, 0, 0, 0, 0, 0, t114 * t250 - t181 * t216, -t116 * t250 - t185 * t216, 0, 0, 0, 0, 0, -t168 * t266 - t250 * t297, t168 * t267 - t208 * t250, 0, 0, 0, 0, 0, -t131 * t28 + t128 * t15 - t266 * t61 + (-t179 * t267 - t183 * t204) * t307, -t131 * t30 + t128 * t14 + t266 * t211 + (t179 * t204 - t183 * t267) * t307; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t263, -t234, t313, -t298, 0, -t173 * t209 - t309, -t173 * t226 + t195, t310, -t323, -t303, -t305, 0, -t168 * t22 + (-t116 * t297 + t168 * t244) * pkin(4) + t293, -t168 * t23 + (-t116 * t208 + t168 * t243) * pkin(4) + t304, -t320, t329, t319, -t328, -t299, t163 * t15 + t214 * t61 + (-t220 * t307 - t3) * t183 + t192 * t179 + t294, t163 * t14 + t183 * t192 - t211 * t214 + t220 * t330 + t291; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t310, -t323, -t303, -t305, 0, -t168 * t21 + t293, -t168 * t20 + t304, -t320, t329, t319, -t328, -t299, -pkin(5) * t15 - t3 * t183 - (-t179 * t20 + t183 * t47) * t307 - t21 * t61 + t179 * t316 - t278 * pkin(10) + t294, -pkin(5) * t14 + (t179 * t47 + t183 * t20) * t307 + t21 * t211 + t18 * t311 + (t242 * t307 - t30) * pkin(10) + t291; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t211 * t61, t211 ^ 2 - t61 ^ 2, t307 * t61 + t14, -t211 * t307 - t15, t34, -t322 * t10 + t179 * t228 + t18 * t211 + t4, -t179 * t5 + t18 * t61 + t183 * t228 + t322 * t212;];
tauc_reg  = t1;
