% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:19
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPPR1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR1_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR1_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR1_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:18:10
% EndTime: 2019-12-05 18:18:18
% DurationCPUTime: 4.01s
% Computational Cost: add. (8624->350), mult. (5274->409), div. (0->0), fcn. (3876->10), ass. (0->214)
t161 = pkin(9) + qJ(5);
t154 = sin(t161);
t155 = cos(t161);
t123 = rSges(6,1) * t154 + rSges(6,2) * t155;
t322 = qJD(5) * t123;
t163 = qJ(1) + qJ(2);
t156 = pkin(8) + t163;
t150 = sin(t156);
t166 = -pkin(7) - qJ(4);
t142 = t150 * t166;
t151 = cos(t156);
t264 = qJ(4) * t150;
t165 = cos(pkin(9));
t152 = pkin(4) * t165 + pkin(3);
t293 = pkin(3) - t152;
t256 = t151 * t154;
t132 = rSges(6,2) * t256;
t255 = t151 * t155;
t198 = rSges(6,1) * t255 + rSges(6,3) * t150;
t79 = -t132 + t198;
t324 = t151 * t293 + t142 + t264 - t79;
t162 = qJD(1) + qJD(2);
t148 = Icges(6,4) * t155;
t205 = -Icges(6,2) * t154 + t148;
t314 = Icges(6,1) * t154 + t148;
t248 = t314 + t205;
t269 = Icges(6,4) * t154;
t113 = Icges(6,2) * t155 + t269;
t116 = Icges(6,1) * t155 - t269;
t249 = t113 - t116;
t323 = (t154 * t248 + t155 * t249) * t162;
t321 = 0.2e1 * qJD(5);
t282 = rSges(5,2) * sin(pkin(9));
t139 = t151 * t282;
t118 = t162 * t139;
t280 = rSges(5,3) * t150;
t285 = rSges(5,1) * t165;
t199 = -t151 * t285 - t280;
t83 = -t139 - t199;
t320 = -t162 * t83 - t118;
t145 = qJD(4) * t150;
t319 = -t145 - t162 * (-pkin(3) * t150 + qJ(4) * t151);
t147 = t150 * rSges(4,2);
t286 = rSges(4,1) * t151;
t107 = -t147 + t286;
t258 = t150 * t162;
t135 = rSges(4,2) * t258;
t158 = cos(t163);
t295 = pkin(2) * t158;
t216 = -t286 - t295;
t318 = -t135 + (-t107 - t216) * t162;
t106 = pkin(3) * t151 + t264;
t146 = qJD(4) * t151;
t317 = t162 * t106 - t146;
t167 = sin(qJ(1));
t278 = pkin(1) * qJD(1);
t238 = t167 * t278;
t157 = sin(t163);
t213 = rSges(3,1) * t157 + rSges(3,2) * t158;
t99 = t213 * t162;
t92 = t99 + t238;
t74 = Icges(6,4) * t255 - Icges(6,2) * t256 + Icges(6,6) * t150;
t131 = Icges(6,4) * t256;
t76 = Icges(6,1) * t255 + Icges(6,5) * t150 - t131;
t208 = t154 * t74 - t155 * t76;
t316 = t208 * t151;
t315 = t132 + t142;
t253 = t157 * t162;
t242 = pkin(2) * t253;
t202 = t242 + t319;
t243 = qJD(5) * t151;
t283 = rSges(6,1) * t155;
t210 = -rSges(6,2) * t154 + t283;
t276 = t151 * rSges(6,3);
t309 = t150 * t210 - t276;
t177 = t123 * t243 + t202 + (-(-qJ(4) - t166) * t151 - t293 * t150 + t309) * t162;
t25 = t177 + t238;
t168 = cos(qJ(1));
t237 = t168 * t278;
t219 = t146 - t237;
t224 = -t106 - t295;
t244 = qJD(5) * t150;
t94 = t123 * t244;
t26 = t94 + (t224 + t324) * t162 + t219;
t313 = t150 * t26 + t151 * t25;
t251 = t162 * t166;
t126 = t150 * t251;
t233 = t162 * t132 + t322 * t150;
t312 = t324 * t162 - t126 - t146 - t233 + t94;
t111 = Icges(6,5) * t154 + Icges(6,6) * t155;
t308 = -Icges(6,3) * t162 + qJD(5) * t111;
t307 = -Icges(6,6) * t162 + qJD(5) * t113;
t103 = t205 * qJD(5);
t104 = t116 * qJD(5);
t306 = t103 * t154 - t104 * t155 - t111 * t162 + (t113 * t155 + t154 * t314) * qJD(5);
t305 = -Icges(6,5) * t162 + qJD(5) * t314;
t287 = -Icges(6,2) * t255 - t131 + t76;
t289 = t151 * t314 + t74;
t303 = t154 * t287 + t155 * t289;
t260 = t150 * t154;
t130 = Icges(6,4) * t260;
t259 = t150 * t155;
t75 = -Icges(6,1) * t259 + Icges(6,5) * t151 + t130;
t288 = Icges(6,2) * t259 + t130 + t75;
t73 = Icges(6,6) * t151 - t150 * t205;
t290 = -t150 * t314 + t73;
t302 = -t154 * t288 - t155 * t290;
t299 = pkin(1) * t167;
t298 = pkin(1) * t168;
t297 = pkin(1) * qJD(1) ^ 2;
t296 = pkin(2) * t157;
t294 = pkin(2) * t162 ^ 2;
t112 = Icges(6,5) * t155 - Icges(6,6) * t154;
t71 = Icges(6,3) * t151 - t112 * t150;
t292 = t151 * t71 + t73 * t260;
t291 = t150 * t71 + t75 * t255;
t279 = rSges(5,3) * t151;
t274 = t154 * t73;
t273 = t155 * t75;
t203 = t113 * t154 - t155 * t314;
t84 = t111 * t150;
t46 = -t151 * t203 + t84;
t272 = t46 * t162;
t271 = rSges(5,3) + qJ(4);
t263 = t111 * t151;
t262 = t112 * t162;
t90 = t123 * t150;
t261 = t123 * t151;
t257 = t151 * t152;
t254 = t151 * t162;
t252 = t158 * t162;
t250 = t151 * t251 + t152 * t258;
t247 = -rSges(4,1) * t258 - rSges(4,2) * t254;
t153 = t167 * t297;
t246 = t157 * t294 + t153;
t138 = pkin(3) * t258;
t245 = t145 - t138;
t241 = pkin(2) * t252;
t240 = t168 * t297;
t239 = t150 * t282;
t234 = -t322 * t151 - t258 * t283;
t117 = t258 * t285;
t232 = t117 - t245;
t229 = -pkin(3) - t285;
t228 = -t244 / 0.2e1;
t226 = -t243 / 0.2e1;
t225 = t243 / 0.2e1;
t72 = Icges(6,5) * t255 - Icges(6,6) * t256 + Icges(6,3) * t150;
t223 = -t72 - t273;
t125 = rSges(3,1) * t158 - t157 * rSges(3,2);
t222 = -qJ(4) * t254 - t145 - t245;
t100 = -rSges(3,1) * t252 + rSges(3,2) * t253;
t218 = -t296 - t299;
t217 = -t295 - t298;
t212 = -rSges(4,1) * t150 - rSges(4,2) * t151;
t211 = t282 - t285;
t40 = t154 * t75 + t155 * t73;
t209 = -t273 + t274;
t41 = t154 * t76 + t155 * t74;
t201 = -t241 - t317;
t28 = t151 * t72 - t259 * t76 + t74 * t260;
t200 = t147 + t216;
t197 = rSges(6,2) * t260 + t276;
t196 = -t145 - t234 + t250;
t93 = -t125 * t162 - t237;
t27 = -t259 * t75 + t292;
t195 = (t150 * t28 + t151 * t27) * qJD(5);
t29 = -t256 * t73 + t291;
t30 = t150 * t72 - t316;
t194 = (t150 * t30 + t151 * t29) * qJD(5);
t193 = -t158 * t294 - t240;
t192 = t116 * t162;
t191 = t205 * t162;
t189 = t238 + t242;
t188 = t237 + t241;
t187 = t212 - t296;
t182 = -t262 * t150 - t151 * t308 + t208 * t162;
t181 = t150 * t308 - t262 * t151 + t209 * t162;
t180 = t112 * qJD(5) + t162 * t203;
t179 = t193 + (t146 - t317) * t162;
t178 = t188 + t317;
t10 = t194 + t272;
t16 = -qJD(5) * t209 + t154 * (t150 * t305 - t151 * t192) + t155 * (t150 * t307 - t151 * t191);
t17 = -qJD(5) * t208 + t154 * (-t150 * t192 - t151 * t305) + t155 * (-t150 * t191 - t151 * t307);
t20 = t180 * t150 - t151 * t306;
t21 = t150 * t306 + t180 * t151;
t45 = t150 * t203 + t263;
t44 = t45 * t162;
t9 = t44 + t195;
t176 = (t44 + ((t292 + t30 + t316) * t151 + (-t29 + (t223 - t274) * t151 + t28 + t291) * t150) * qJD(5)) * t228 + (t10 - t272 + ((t28 + (-t72 + t274) * t151 - t291) * t151 + (t150 * t223 - t27 + t292) * t150) * qJD(5)) * t226 + (t16 + t21) * t225 + (t17 + t20 + t9) * t244 / 0.2e1 + (-qJD(5) * t203 + t103 * t155 + t104 * t154 + (t40 + t45) * t228 + (t41 + t46) * t225) * t162;
t173 = t150 * t309 + t151 * t79;
t53 = t162 * t197 + t234;
t54 = -t162 * t198 + t233;
t172 = (-rSges(6,3) * t254 + t53) * t151 + (-t54 + (t151 * t210 - t79) * t162) * t150;
t105 = t210 * qJD(5);
t14 = t105 * t244 + (-t138 - t53 + (qJ(4) * t162 + t322) * t151 + t222 + t250) * t162 + t246;
t15 = -t105 * t243 + (t94 + t126 + t54 + (t106 - t257) * t162) * t162 + t179;
t171 = (-t14 * rSges(6,3) + t15 * (-t152 - t210)) * t150 + (t14 * (-t152 - t283) + t15 * (rSges(6,3) - t166)) * t151 + (t26 * (-t197 + t296) - t25 * (-t198 - t257 - t295)) * t162;
t31 = (t117 + (-t239 - t279) * t162 + t222) * t162 + t246;
t32 = (t162 * t199 + t118) * t162 + t179;
t77 = t162 * (t150 * t211 + t279);
t42 = t189 - t77 + t319;
t43 = (t224 - t83) * t162 + t219;
t170 = (-t31 * t271 + t32 * (-pkin(3) + t211)) * t150 + (t31 * t229 + t32 * t271) * t151 + (t43 * (-t239 + t296) - t42 * (-t264 - t280 - t295) + (-t229 * t42 - t271 * t43) * t151) * t162;
t96 = t162 * t212;
t82 = t100 * t162 - t240;
t81 = t162 * t99 + t153;
t69 = -t237 + (-t107 - t295) * t162;
t68 = t189 - t96;
t59 = t162 * (-rSges(4,1) * t254 + t135) + t193;
t58 = -t162 * t247 + t246;
t39 = qJD(5) * t173 + qJD(3);
t11 = t172 * qJD(5);
t1 = [t176 + m(3) * (t81 * (-t125 - t298) + t82 * (-t213 - t299) + (t93 - t100 + t237) * t92) + (t14 * (t217 + t315) + t26 * (t196 + t238) + t15 * t218 + t171 + (t237 - t178 - t26 + t312) * t25) * m(6) + (t31 * (t139 + t217) + t43 * (t232 + t238) + t32 * t218 + t170 + (-t219 - t178 - t43 + t320) * t42) * m(5) + (t58 * (t200 - t298) + t69 * (t189 - t247) + t59 * (t187 - t299) + (t237 - t188 - t69 + t318) * t68) * m(4); t176 + (t14 * (-t295 + t315) - t15 * t296 + t171 + (-t177 + t196) * t26 + (t201 + t312) * t25) * m(6) + (t31 * (t139 - t295) - t32 * t296 + t170 + (-t202 + t77 + t232) * t43 + (t201 - t146 + t320) * t42) * m(5) + (t59 * t187 + t58 * t200 + (-t241 + t318) * t68 + (-t247 + t96) * t69) * m(4) + (-(t92 * t125 + t213 * t93) * t162 - t92 * t100 - t81 * t125 - t82 * t213 + t93 * t99) * m(3); m(6) * t11; m(5) * (t150 * t32 + t151 * t31) + m(6) * (t14 * t151 + t15 * t150); t162 * ((t162 * t41 + t16) * t151 + (-t162 * t40 + t17) * t150) / 0.2e1 - t162 * ((-t249 * t154 + t248 * t155) * t162 + ((t150 * t287 + t151 * t288) * t155 + (-t150 * t289 - t151 * t290) * t154) * qJD(5)) / 0.2e1 + ((t84 * t243 + t262) * t151 + (t323 + (t303 * t150 + (-t302 - t263) * t151) * qJD(5)) * t150) * t226 + ((-t244 * t263 + t262) * t150 + (-t323 + (t302 * t151 + (-t303 + t84) * t150) * qJD(5)) * t151) * t228 + (t162 * t20 + ((t181 * t150 + t162 * t30) * t151 + (t182 * t150 - t162 * t29) * t150) * t321) * t150 / 0.2e1 + (t162 * t21 + ((t181 * t151 + t162 * t28) * t151 + (t182 * t151 - t162 * t27) * t150) * t321) * t151 / 0.2e1 - (t9 + t195) * t258 / 0.2e1 + (t10 + t194) * t254 / 0.2e1 + (t11 * t173 + t14 * t90 - t15 * t261 + t39 * t172 - (-t25 * t90 + t26 * t261) * t162 - (t39 * (-t150 * t90 - t151 * t261) + t313 * t210) * qJD(5) + (-t25 * t258 + t254 * t26) * t123 + t313 * t105) * m(6);];
tauc = t1(:);
