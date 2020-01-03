% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRPR6_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR6_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR6_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR6_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR6_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:20
% EndTime: 2019-12-31 16:24:33
% DurationCPUTime: 11.31s
% Computational Cost: add. (7275->480), mult. (12172->762), div. (0->0), fcn. (11569->8), ass. (0->249)
t207 = sin(qJ(2));
t208 = cos(qJ(2));
t201 = pkin(7) + qJ(4);
t196 = sin(t201);
t197 = cos(t201);
t273 = rSges(5,1) * t197 - rSges(5,2) * t196;
t132 = -rSges(5,3) * t208 + t207 * t273;
t205 = cos(pkin(6));
t298 = qJD(4) * t207;
t203 = sin(pkin(6));
t303 = qJD(2) * t203;
t180 = t205 * t298 + t303;
t371 = t132 * t180;
t302 = qJD(2) * t205;
t181 = t203 * t298 - t302;
t370 = t132 * t181;
t199 = t203 ^ 2;
t200 = t205 ^ 2;
t304 = t199 + t200;
t202 = sin(pkin(7));
t204 = cos(pkin(7));
t300 = qJD(2) * t208;
t301 = qJD(2) * t207;
t331 = Icges(3,4) * t208;
t332 = Icges(3,4) * t207;
t369 = -(-Icges(3,1) * t207 - t331) * t300 + (-(-Icges(4,5) * t204 + Icges(4,6) * t202) * t207 - t332 + (-Icges(4,3) - Icges(3,2)) * t208) * t301;
t356 = -rSges(4,1) * t204 + rSges(4,2) * t202;
t354 = rSges(4,3) * t208 + t207 * t356;
t317 = t205 * t208;
t169 = -t202 * t317 + t203 * t204;
t322 = t202 * t203;
t170 = t204 * t317 + t322;
t260 = -Icges(3,2) * t207 + t331;
t263 = Icges(3,1) * t208 - t332;
t318 = t205 * t207;
t368 = t369 * t205 + ((-Icges(4,5) * t170 + Icges(3,6) * t203 - Icges(4,6) * t169 - Icges(4,3) * t318 + t205 * t260) * t208 + (-t202 * (Icges(4,4) * t170 + Icges(4,2) * t169 + Icges(4,6) * t318) + t204 * (Icges(4,1) * t170 + Icges(4,4) * t169 + Icges(4,5) * t318) + Icges(3,5) * t203 + t205 * t263) * t207) * qJD(2);
t319 = t203 * t208;
t167 = -t202 * t319 - t204 * t205;
t321 = t202 * t205;
t168 = t204 * t319 - t321;
t320 = t203 * t207;
t367 = t369 * t203 + ((-Icges(4,5) * t168 - Icges(3,6) * t205 - Icges(4,6) * t167 - Icges(4,3) * t320 + t203 * t260) * t208 + (-t202 * (Icges(4,4) * t168 + Icges(4,2) * t167 + Icges(4,6) * t320) + t204 * (Icges(4,1) * t168 + Icges(4,4) * t167 + Icges(4,5) * t320) - Icges(3,5) * t205 + t203 * t263) * t207) * qJD(2);
t198 = qJD(3) * t207;
t190 = t203 * t198;
t183 = pkin(2) * t207 - qJ(3) * t208;
t241 = qJD(2) * t183;
t121 = -t203 * t241 + t190;
t192 = t205 * t198;
t122 = -t205 * t241 + t192;
t366 = t121 * t203 + t122 * t205 + t241 * t304 - t198;
t140 = -t196 * t319 - t197 * t205;
t141 = -t196 * t205 + t197 * t319;
t54 = Icges(5,5) * t141 + Icges(5,6) * t140 + Icges(5,3) * t320;
t330 = Icges(5,4) * t141;
t56 = Icges(5,2) * t140 + Icges(5,6) * t320 + t330;
t134 = Icges(5,4) * t140;
t58 = Icges(5,1) * t141 + Icges(5,5) * t320 + t134;
t23 = t140 * t56 + t141 * t58 + t320 * t54;
t142 = -t196 * t317 + t197 * t203;
t143 = t196 * t203 + t197 * t317;
t55 = Icges(5,5) * t143 + Icges(5,6) * t142 + Icges(5,3) * t318;
t329 = Icges(5,4) * t143;
t57 = Icges(5,2) * t142 + Icges(5,6) * t318 + t329;
t135 = Icges(5,4) * t142;
t59 = Icges(5,1) * t143 + Icges(5,5) * t318 + t135;
t24 = t140 * t57 + t141 * t59 + t320 * t55;
t25 = t142 * t56 + t143 * t58 + t318 * t54;
t26 = t142 * t57 + t143 * t59 + t318 * t55;
t297 = qJD(4) * t208;
t256 = Icges(5,5) * t197 - Icges(5,6) * t196;
t125 = -Icges(5,3) * t208 + t207 * t256;
t327 = Icges(5,4) * t197;
t258 = -Icges(5,2) * t196 + t327;
t127 = -Icges(5,6) * t208 + t207 * t258;
t328 = Icges(5,4) * t196;
t261 = Icges(5,1) * t197 - t328;
t129 = -Icges(5,5) * t208 + t207 * t261;
t39 = t125 * t320 + t127 * t140 + t129 * t141;
t40 = t125 * t318 + t127 * t142 + t129 * t143;
t364 = t205 * (t180 * t26 + t181 * t25 - t297 * t40) + t203 * (t180 * t24 + t181 * t23 - t297 * t39);
t363 = t304 * qJD(2) * t354;
t347 = t180 / 0.2e1;
t345 = t181 / 0.2e1;
t360 = qJD(4) / 0.2e1;
t184 = rSges(3,1) * t207 + rSges(3,2) * t208;
t242 = qJD(2) * t184;
t270 = -t196 * t57 + t197 * t59;
t271 = -t196 * t56 + t197 * t58;
t353 = -(-t125 * t205 - t270) * t180 - (-t125 * t203 - t271) * t181;
t206 = -pkin(5) - qJ(3);
t316 = qJ(3) + t206;
t194 = pkin(3) * t204 + pkin(2);
t343 = pkin(2) - t194;
t351 = t207 * t343 - t208 * t316;
t153 = (-Icges(5,2) * t197 - t328) * t207;
t217 = t180 * (-Icges(5,2) * t143 + t135 + t59) + t181 * (-Icges(5,2) * t141 + t134 + t58) - t297 * (t129 + t153);
t287 = t205 * t300;
t288 = t205 * t301;
t92 = -qJD(4) * t143 + t196 * t288;
t93 = qJD(4) * t142 - t197 * t288;
t43 = Icges(5,5) * t93 + Icges(5,6) * t92 + Icges(5,3) * t287;
t246 = t207 * t43 + t300 * t55;
t45 = Icges(5,4) * t93 + Icges(5,2) * t92 + Icges(5,6) * t287;
t47 = Icges(5,1) * t93 + Icges(5,4) * t92 + Icges(5,5) * t287;
t10 = t142 * t45 + t143 * t47 + t205 * t246 + t57 * t92 + t59 * t93;
t126 = Icges(5,3) * t207 + t208 * t256;
t152 = (-Icges(5,5) * t196 - Icges(5,6) * t197) * t207;
t63 = qJD(2) * t126 + qJD(4) * t152;
t245 = t125 * t300 + t207 * t63;
t336 = t203 * t25;
t266 = t205 * t26 + t336;
t333 = t40 * t207;
t128 = Icges(5,6) * t207 + t208 * t258;
t64 = qJD(2) * t128 + qJD(4) * t153;
t130 = Icges(5,5) * t207 + t208 * t261;
t154 = (-Icges(5,1) * t196 - t327) * t207;
t65 = qJD(2) * t130 + qJD(4) * t154;
t215 = -(t127 * t92 + t129 * t93 + t142 * t64 + t143 * t65 + t205 * t245) * t208 + (t208 * t266 + t333) * qJD(2);
t289 = t203 * t300;
t290 = t203 * t301;
t90 = -qJD(4) * t141 + t196 * t290;
t91 = qJD(4) * t140 - t197 * t290;
t42 = Icges(5,5) * t91 + Icges(5,6) * t90 + Icges(5,3) * t289;
t247 = t207 * t42 + t300 * t54;
t44 = Icges(5,4) * t91 + Icges(5,2) * t90 + Icges(5,6) * t289;
t46 = Icges(5,1) * t91 + Icges(5,4) * t90 + Icges(5,5) * t289;
t9 = t142 * t44 + t143 * t46 + t205 * t247 + t56 * t92 + t58 * t93;
t350 = t10 * t347 + t215 * t360 + t345 * t9;
t349 = -t304 * t301 / 0.2e1;
t348 = -t180 / 0.2e1;
t346 = -t181 / 0.2e1;
t344 = -t208 / 0.2e1;
t335 = t205 * t24;
t334 = t39 * t207;
t323 = t125 * t208;
t124 = -t207 * t316 - t208 * t343;
t185 = pkin(2) * t208 + qJ(3) * t207;
t299 = qJD(3) * t208;
t166 = qJD(2) * t185 - t299;
t315 = -qJD(2) * t124 - t166;
t313 = -t183 + t351;
t312 = -t124 - t185;
t151 = rSges(4,3) * t207 - t208 * t356;
t310 = -qJD(2) * t151 - t166;
t307 = -t183 + t354;
t306 = -t151 - t185;
t178 = t185 * t203;
t179 = t185 * t205;
t305 = t178 * t203 + t179 * t205;
t296 = qJD(2) * qJD(3);
t133 = rSges(5,3) * t207 + t208 * t273;
t155 = (-rSges(5,1) * t196 - rSges(5,2) * t197) * t207;
t66 = qJD(2) * t133 + qJD(4) * t155;
t295 = -t66 + t315;
t293 = t121 * t303 + t122 * t302 + t207 * t296;
t292 = -t132 + t313;
t286 = t208 * t296;
t282 = -t297 / 0.2e1;
t281 = t297 / 0.2e1;
t280 = qJD(2) * t360;
t279 = qJD(2) * t313;
t278 = qJD(2) * t310;
t277 = qJD(2) * t307;
t276 = t207 * t280;
t275 = t208 * t280;
t186 = rSges(3,1) * t208 - rSges(3,2) * t207;
t272 = t180 * t55 + t181 * t54;
t267 = t203 * t23 + t335;
t27 = t207 * t271 - t208 * t54;
t28 = t207 * t270 - t208 * t55;
t265 = t203 * t27 + t28 * t205;
t60 = rSges(5,1) * t141 + rSges(5,2) * t140 + rSges(5,3) * t320;
t61 = rSges(5,1) * t143 + rSges(5,2) * t142 + rSges(5,3) * t318;
t264 = -t203 * t61 + t205 * t60;
t257 = -Icges(3,5) * t207 - Icges(3,6) * t208;
t255 = -t127 * t196 + t129 * t197;
t254 = t304 * t186;
t253 = t304 * t242;
t252 = t203 * t275;
t251 = t205 * t275;
t248 = t178 * t303 + t179 * t302 + qJD(1) - t299;
t244 = t126 - t255;
t86 = -pkin(3) * t321 + t124 * t203;
t87 = pkin(3) * t322 + t124 * t205;
t22 = t180 * t60 - t181 * t61 + (t203 * t86 + t205 * t87) * qJD(2) + t248;
t243 = t22 * t264;
t238 = qJD(2) * t257;
t236 = t152 * t297 - t180 * (Icges(5,5) * t142 - Icges(5,6) * t143) - t181 * (Icges(5,5) * t140 - Icges(5,6) * t141);
t234 = qJD(2) * t351;
t233 = Icges(4,5) * t208 + (-Icges(4,1) * t204 + Icges(4,4) * t202) * t207;
t231 = Icges(4,6) * t208 + (-Icges(4,4) * t204 + Icges(4,2) * t202) * t207;
t227 = t207 * t236;
t84 = rSges(4,1) * t168 + rSges(4,2) * t167 + rSges(4,3) * t320;
t85 = rSges(4,1) * t170 + rSges(4,2) * t169 + rSges(4,3) * t318;
t33 = (t203 * t84 + t205 * t85) * qJD(2) + t248;
t226 = t33 * t354;
t224 = t22 * (-t194 * t207 - t206 * t208 + t183);
t223 = qJD(2) * t233;
t222 = qJD(2) * t231;
t218 = (Icges(5,1) * t142 - t329 - t57) * t180 + (Icges(5,1) * t140 - t330 - t56) * t181 - (-t127 + t154) * t297;
t216 = -(t127 * t90 + t129 * t91 + t140 * t64 + t141 * t65 + t203 * t245) * t208 + (t208 * t267 + t334) * qJD(2);
t41 = t207 * t255 - t323;
t214 = -((qJD(2) * t255 - t63) * t208 + (qJD(2) * t125 - t196 * t64 + t197 * t65 + (-t127 * t197 - t129 * t196) * qJD(4)) * t207) * t208 + (t41 * t207 + t208 * t265) * qJD(2);
t210 = (-t244 * t297 - t353) * t207;
t193 = t205 * t299;
t191 = t203 * t299;
t188 = t205 * t286;
t187 = t203 * t286;
t173 = t257 * t205;
t172 = t257 * t203;
t159 = t205 * t238;
t158 = t203 * t238;
t120 = t233 * t205;
t119 = t233 * t203;
t118 = t231 * t205;
t117 = t231 * t203;
t105 = t129 * t205;
t104 = t129 * t203;
t103 = t127 * t205;
t102 = t127 * t203;
t99 = t205 * t223;
t98 = t203 * t223;
t97 = t205 * t222;
t96 = t203 * t222;
t89 = t205 * t234;
t88 = t203 * t234;
t83 = rSges(5,1) * t142 - rSges(5,2) * t143;
t82 = rSges(5,1) * t140 - rSges(5,2) * t141;
t69 = t253 * qJD(2);
t68 = t205 * t277 + t192;
t67 = t203 * t277 + t190;
t62 = qJD(2) * t254 + qJD(1);
t51 = t205 * t278 + t188;
t50 = t203 * t278 + t187;
t49 = rSges(5,1) * t93 + rSges(5,2) * t92 + rSges(5,3) * t287;
t48 = rSges(5,1) * t91 + rSges(5,2) * t90 + rSges(5,3) * t289;
t34 = qJD(2) * t363 + t293;
t32 = t205 * t279 + t297 * t60 + t192 + t370;
t31 = t203 * t279 - t297 * t61 + t190 - t371;
t17 = t48 * t297 + t181 * t66 + t188 + (t315 * t205 + (t132 * t319 - t207 * t60) * qJD(4)) * qJD(2);
t16 = -t49 * t297 - t180 * t66 + t187 + (t315 * t203 + (-t132 * t317 + t207 * t61) * qJD(4)) * qJD(2);
t12 = t180 * t48 - t181 * t49 + (t203 * t88 + t205 * t89 + t264 * t297) * qJD(2) + t293;
t11 = t180 * t28 + t181 * t27 - t297 * t41;
t8 = t140 * t45 + t141 * t47 + t203 * t246 + t57 * t90 + t59 * t91;
t7 = t140 * t44 + t141 * t46 + t203 * t247 + t56 * t90 + t58 * t91;
t4 = (qJD(2) * t270 - t43) * t208 + (qJD(2) * t55 - t196 * t45 + t197 * t47 + (-t196 * t59 - t197 * t57) * qJD(4)) * t207;
t3 = (qJD(2) * t271 - t42) * t208 + (qJD(2) * t54 - t196 * t44 + t197 * t46 + (-t196 * t58 - t197 * t56) * qJD(4)) * t207;
t1 = qJD(4) * t216 + t180 * t8 + t181 * t7;
t2 = [-m(3) * t69 + m(4) * t34 + m(5) * t12; (((t103 * t196 - t105 * t197 + t55) * t180 + (t102 * t196 - t104 * t197 + t54) * t181 + t41 * qJD(4)) * t207 + ((t244 * t208 + (t128 * t196 - t130 * t197 - t125) * t207 + t265) * qJD(4) + t353) * t208) * t281 + (t203 * t8 - t205 * t7) * t345 + ((-t103 * t140 - t105 * t141) * t180 + (-t102 * t140 - t104 * t141) * t181 + (t334 + (-t128 * t140 - t130 * t141 + t335) * t208) * qJD(4) + (((t23 - t323) * qJD(4) + t272) * t208 + t210) * t203) * t346 + (t10 * t203 - t205 * t9) * t347 + ((-t103 * t142 - t105 * t143) * t180 + (-t102 * t142 - t104 * t143) * t181 + (t333 + (-t128 * t142 - t130 * t143 + t336) * t208) * qJD(4) + (((t26 - t323) * qJD(4) + t272) * t208 + t210) * t205) * t348 + t203 * t350 - t11 * t298 / 0.2e1 + (t203 * t28 - t205 * t27) * t276 - t205 * t1 / 0.2e1 + (t203 * t26 - t205 * t25) * t251 + (t203 * t24 - t205 * t23) * t252 + ((t159 * t203 + t169 * t97 + t170 * t99) * t203 + (-t169 * t96 - t170 * t98 + t367 * t205 + (-t158 - t368) * t203) * t205) * t303 + ((-t158 * t205 + t167 * t96 + t168 * t98) * t205 + (-t167 * t97 - t168 * t99 + t368 * t203 + (t159 - t367) * t205) * t203) * t302 - ((t118 * t169 + t120 * t170) * t303 + t173 * qJD(2) * t199 + (-t169 * t117 - t170 * t119 - t172 * t203) * t302) * t303 / 0.2e1 + (-(t117 * t167 + t119 * t168) * t302 + t172 * qJD(2) * t200 + (t167 * t118 + t168 * t120 - t173 * t205) * t303) * t302 / 0.2e1 + (t12 * t305 + (t17 * t292 + t32 * t295 + t12 * (t61 + t87)) * t205 + (t16 * t292 + t31 * t295 + t12 * (t60 + t86)) * t203 - t32 * (t133 * t181 + t193) - t31 * (-t133 * t180 + t191) - (t205 * t224 + t312 * t32) * t302 - (t203 * t224 + t31 * t312) * t303 - ((t31 * t61 - t32 * t60) * t207 + t243 * t208) * qJD(4) + (t366 + (t49 + t89 - t370) * t205 + (t48 + t88 + t371) * t203) * t22) * m(5) + (t34 * t305 + (t307 * t51 + t310 * t68 + t34 * t85) * t205 + (t307 * t50 + t310 * t67 + t34 * t84) * t203 - t68 * t193 - t67 * t191 - ((t205 * t226 + t306 * t68) * t205 + (t203 * t226 + t306 * t67) * t203) * qJD(2) + (t363 + t366) * t33) * m(4) + (t203 * t4 - t205 * t3 + t364) * t282 + (-t62 * t253 - t69 * t254 + (qJD(2) ^ 2 * t184 * t186 + t242 * t62) * t304) * m(3); 0.2e1 * (t12 * t344 + t22 * t349) * m(5) + 0.2e1 * (t33 * t349 + t34 * t344) * m(4) + 0.2e1 * (m(4) * (qJD(2) * t33 + t203 * t50 + t205 * t51) / 0.2e1 + m(5) * (qJD(2) * t22 + t16 * t203 + t17 * t205) / 0.2e1) * t207; t318 * t350 + (t207 * t266 - t208 * t40) * t251 + ((t10 * t205 + t203 * t9) * t207 + t215) * t347 + t1 * t320 / 0.2e1 + (t207 * t267 - t208 * t39) * t252 + ((t203 * t7 + t205 * t8) * t207 + t216) * t345 + t11 * t301 / 0.2e1 + (qJD(4) * t214 + t180 * t4 + t181 * t3) * t344 + (t207 * t265 - t208 * t41) * t276 + ((t203 * t3 + t205 * t4) * t207 + t214) * t282 + (t142 * t217 + t143 * t218 - t205 * t227) * t348 + (t140 * t217 + t141 * t218 - t203 * t227) * t346 + (t236 * t208 + (-t196 * t217 + t218 * t197) * t207) * t281 + t364 * t300 / 0.2e1 + ((-t16 * t61 + t17 * t60 - t31 * t49 + t32 * t48 + (t243 + (t203 * t32 - t205 * t31) * t132) * qJD(2)) * t208 + (t32 * (-qJD(2) * t60 + t203 * t66) + t31 * (qJD(2) * t61 - t205 * t66) + t12 * t264 + t22 * (-t203 * t49 + t205 * t48) + (-t16 * t205 + t17 * t203) * t132) * t207 - t32 * (t155 * t181 + t297 * t82) - t31 * (-t155 * t180 - t297 * t83) - t22 * (t180 * t82 - t181 * t83)) * m(5);];
tauc = t2(:);
