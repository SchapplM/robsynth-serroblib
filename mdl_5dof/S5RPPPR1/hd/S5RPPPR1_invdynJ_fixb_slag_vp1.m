% Calculate vector of inverse dynamics joint torques for
% S5RPPPR1
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% m [6x1]
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:13
% Revision: 008671b0a00594318b890887636eaaff83cd5e2f (2021-12-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPPR1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR1_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR1_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:12:13
% EndTime: 2022-01-20 09:12:29
% DurationCPUTime: 9.56s
% Computational Cost: add. (9352->488), mult. (9298->649), div. (0->0), fcn. (8638->10), ass. (0->231)
t207 = cos(qJ(1));
t198 = t207 * pkin(1);
t200 = qJ(1) + pkin(7);
t195 = sin(t200);
t197 = cos(t200);
t153 = rSges(3,1) * t195 + rSges(3,2) * t197;
t206 = sin(qJ(1));
t327 = pkin(1) * t206;
t146 = -t153 - t327;
t199 = pkin(9) + qJ(5);
t194 = sin(t199);
t196 = cos(t199);
t202 = sin(pkin(8));
t204 = cos(pkin(8));
t115 = -rSges(6,3) * t204 + (rSges(6,1) * t196 - rSges(6,2) * t194) * t202;
t182 = -qJD(5) * t204 + qJD(1);
t279 = qJD(5) * t202;
t261 = t195 * t279;
t301 = t195 * t204;
t116 = t194 * t301 + t196 * t197;
t117 = -t197 * t194 + t196 * t301;
t294 = -t117 * rSges(6,1) + t116 * rSges(6,2);
t302 = t195 * t202;
t72 = rSges(6,3) * t302 - t294;
t349 = t115 * t261 - t182 * t72;
t208 = qJD(1) ^ 2;
t348 = t208 * t198;
t63 = Icges(6,5) * t117 - Icges(6,6) * t116 + Icges(6,3) * t302;
t103 = Icges(6,4) * t117;
t66 = -Icges(6,2) * t116 + Icges(6,6) * t302 + t103;
t102 = Icges(6,4) * t116;
t70 = -Icges(6,1) * t117 - Icges(6,5) * t302 + t102;
t25 = -(t194 * t66 + t196 * t70) * t202 - t204 * t63;
t346 = -t116 * t66 - t117 * t70;
t298 = t197 * t204;
t118 = -t194 * t298 + t195 * t196;
t119 = t194 * t195 + t196 * t298;
t345 = t118 * t66 - t119 * t70;
t299 = t197 * t202;
t342 = t63 * t299;
t145 = pkin(3) * t298 + qJ(4) * t299;
t154 = t197 * pkin(2) + t195 * qJ(3);
t336 = t198 + t154;
t241 = t145 + t336;
t203 = cos(pkin(9));
t296 = t203 * t204;
t201 = sin(pkin(9));
t297 = t201 * t204;
t303 = t195 * t201;
t341 = -(t197 * t296 + t303) * rSges(5,1) - (t195 * t203 - t197 * t297) * rSges(5,2);
t82 = rSges(5,3) * t299 - t341;
t51 = t241 + t82;
t101 = rSges(4,1) * t298 - rSges(4,2) * t299 + t195 * rSges(4,3);
t96 = t101 + t336;
t300 = t197 * t201;
t233 = (t195 * t296 - t300) * rSges(5,1) - (t195 * t297 + t197 * t203) * rSges(5,2);
t18 = t342 + t345;
t337 = t18 - t342;
t155 = t197 * rSges(3,1) - rSges(3,2) * t195;
t147 = t155 + t198;
t65 = Icges(6,5) * t119 + Icges(6,6) * t118 + Icges(6,3) * t299;
t308 = Icges(6,4) * t119;
t68 = Icges(6,2) * t118 + Icges(6,6) * t299 + t308;
t104 = Icges(6,4) * t118;
t71 = Icges(6,1) * t119 + Icges(6,5) * t299 + t104;
t17 = -t116 * t68 + t117 * t71 + t65 * t302;
t334 = t195 * (-Icges(6,2) * t117 - t102 - t70) + t197 * (-Icges(6,2) * t119 + t104 + t71);
t332 = -m(5) - m(6);
t276 = qJD(1) * qJD(5);
t134 = (qJDD(5) * t195 + t197 * t276) * t202;
t331 = t134 / 0.2e1;
t135 = (qJDD(5) * t197 - t195 * t276) * t202;
t330 = t135 / 0.2e1;
t329 = t195 / 0.2e1;
t328 = -t197 / 0.2e1;
t326 = pkin(3) * t204;
t325 = pkin(4) * t201;
t324 = g(1) * t195;
t144 = (-rSges(6,1) * t194 - rSges(6,2) * t196) * t202;
t127 = qJD(5) * t144;
t181 = -qJDD(5) * t204 + qJDD(1);
t277 = qJD(1) * qJD(3);
t224 = qJDD(3) * t195 + t197 * t277 - t348;
t274 = qJDD(4) * t202;
t218 = t197 * t274 + t224;
t305 = qJ(4) * t202;
t231 = t305 + t326;
t143 = t231 * t195;
t187 = t197 * qJ(3);
t152 = pkin(2) * t195 - t187;
t255 = -t152 - t327;
t242 = -t143 + t255;
t192 = pkin(4) * t203 + pkin(3);
t245 = pkin(4) * t300 - t192 * t301;
t205 = -pkin(6) - qJ(4);
t252 = (qJ(4) + t205) * t202;
t79 = (t252 + t326) * t195 + t245;
t226 = t242 + t79;
t258 = (pkin(3) - t192) * t204;
t269 = qJD(1) * t325;
t280 = qJD(4) * t202;
t283 = qJD(1) * t197;
t185 = qJD(3) * t197;
t121 = qJD(1) * t154 - t185;
t262 = t195 * t280;
t309 = -t231 * t283 - t121 - t262;
t85 = qJD(1) * t118 - qJD(5) * t117;
t86 = qJD(1) * t119 - qJD(5) * t116;
t240 = -rSges(6,1) * t86 - rSges(6,2) * t85;
t282 = qJD(1) * t202;
t263 = t197 * t282;
t45 = rSges(6,3) * t263 - t240;
t7 = t127 * t261 + t134 * t115 - t181 * t72 - t182 * t45 + t226 * qJDD(1) + ((-t269 - t280) * t195 + (t258 + t252) * t283 + t309) * qJD(1) + t218;
t323 = t7 * t195;
t166 = t197 * t280;
t246 = qJDD(1) * t198 - t208 * t327;
t284 = qJD(1) * t195;
t184 = qJD(3) * t195;
t286 = qJ(3) * t283 + t184;
t213 = qJD(1) * (-pkin(2) * t284 + t286) + qJDD(1) * t154 + t195 * t277 + t246;
t209 = qJDD(1) * t145 + t195 * t274 + t213 + (-t231 * t284 + 0.2e1 * t166) * qJD(1);
t264 = t195 * t282;
t290 = t197 * t269 + t205 * t264;
t83 = qJD(1) * t116 - qJD(5) * t119;
t84 = -qJD(1) * t117 + qJD(5) * t118;
t318 = t84 * rSges(6,1) + t83 * rSges(6,2);
t44 = -rSges(6,3) * t264 + t318;
t74 = t119 * rSges(6,1) + t118 * rSges(6,2) + rSges(6,3) * t299;
t223 = pkin(4) * t303 + t192 * t298 - t205 * t299;
t80 = t223 - t145;
t8 = ((t258 + t305) * t284 + t290) * qJD(1) + (-t127 * t279 - qJDD(3)) * t197 + t209 + qJDD(1) * t80 - t135 * t115 + t181 * t74 + t182 * t44;
t322 = t8 * t197;
t306 = Icges(6,4) * t196;
t113 = -Icges(6,6) * t204 + (-Icges(6,2) * t194 + t306) * t202;
t307 = Icges(6,4) * t194;
t114 = -Icges(6,5) * t204 + (Icges(6,1) * t196 - t307) * t202;
t140 = (-Icges(6,5) * t194 - Icges(6,6) * t196) * t202;
t124 = qJD(5) * t140;
t141 = (-Icges(6,2) * t196 - t307) * t202;
t125 = qJD(5) * t141;
t142 = (-Icges(6,1) * t194 - t306) * t202;
t126 = qJD(5) * t142;
t28 = -t124 * t204 + (-t125 * t194 + t126 * t196 + (-t113 * t196 - t114 * t194) * qJD(5)) * t202;
t112 = -Icges(6,3) * t204 + (Icges(6,5) * t196 - Icges(6,6) * t194) * t202;
t50 = -t112 * t204 + (-t113 * t194 + t114 * t196) * t202;
t321 = t50 * t181 + t28 * t182;
t16 = t302 * t63 + t346;
t316 = t16 * t195;
t31 = t112 * t302 - t113 * t116 + t114 * t117;
t315 = t182 * t31;
t288 = t166 + t184;
t23 = qJD(1) * t226 + t288 + t349;
t314 = t197 * t23;
t313 = t25 * t134;
t26 = -t204 * t65 + (-t194 * t68 + t196 * t71) * t202;
t312 = t26 * t135;
t311 = -rSges(5,3) - qJ(4);
t310 = -rSges(6,3) + t205;
t304 = t192 * t204;
t293 = t233 * qJD(1);
t292 = -t113 + t142;
t291 = t114 + t141;
t289 = rSges(4,2) * t264 + rSges(4,3) * t283;
t287 = rSges(4,2) * t302 + t197 * rSges(4,3);
t285 = -qJD(1) * t152 + t184;
t281 = qJD(4) * t195;
t278 = -m(4) + t332;
t275 = qJDD(3) * t197;
t19 = t118 * t68 + t119 * t71 + t65 * t299;
t273 = rSges(4,1) * t301;
t271 = rSges(5,3) * t302;
t267 = t166 + t286;
t265 = -pkin(2) - t326;
t259 = -rSges(4,1) * t204 - pkin(2);
t257 = -t279 / 0.2e1;
t256 = t279 / 0.2e1;
t253 = t187 - t327;
t251 = -qJD(1) * t143 + t166 + t285;
t180 = -qJDD(4) * t204 + qJDD(2);
t250 = t195 * t257;
t249 = t195 * t256;
t248 = t197 * t257;
t247 = t197 * t256;
t100 = t273 - t287;
t244 = -t100 + t255;
t243 = -t327 - t233;
t38 = Icges(6,5) * t84 + Icges(6,6) * t83 - Icges(6,3) * t264;
t39 = Icges(6,5) * t86 + Icges(6,6) * t85 + Icges(6,3) * t263;
t40 = Icges(6,4) * t84 + Icges(6,2) * t83 - Icges(6,6) * t264;
t41 = Icges(6,4) * t86 + Icges(6,2) * t85 + Icges(6,6) * t263;
t42 = Icges(6,1) * t84 + Icges(6,4) * t83 - Icges(6,5) * t264;
t43 = Icges(6,1) * t86 + Icges(6,4) * t85 + Icges(6,5) * t263;
t239 = (t118 * t41 + t119 * t43 + t66 * t83 - t70 * t84 + (t197 * t39 - t284 * t63) * t202) * t195 + t197 * (t118 * t40 + t119 * t42 + t68 * t83 + t71 * t84 + (t197 * t38 - t284 * t65) * t202);
t238 = t195 * (-t116 * t41 + t117 * t43 + t66 * t85 - t70 * t86 + (t195 * t39 + t283 * t63) * t202) + t197 * (-t116 * t40 + t117 * t42 + t68 * t85 + t71 * t86 + (t195 * t38 + t283 * t65) * t202);
t237 = t202 * t310 - pkin(2);
t10 = -t204 * t38 + (-t194 * t40 + t196 * t42 + (-t194 * t71 - t196 * t68) * qJD(5)) * t202;
t9 = -t204 * t39 + (-t194 * t41 + t196 * t43 + (t194 * t70 - t196 * t66) * qJD(5)) * t202;
t236 = t10 * t197 + t195 * t9;
t173 = rSges(2,1) * t207 - rSges(2,2) * t206;
t172 = rSges(2,1) * t206 + rSges(2,2) * t207;
t234 = t341 * qJD(1);
t24 = t182 * t74 - t185 + (-qJD(5) * t115 * t197 + t281) * t202 + (t241 + t80) * qJD(1);
t230 = t195 * t23 - t197 * t24;
t229 = -t195 * t44 + t197 * t45;
t228 = -t195 * t74 + t197 * t72;
t227 = t195 * (-Icges(6,5) * t116 - Icges(6,6) * t117) + t197 * (Icges(6,5) * t118 - Icges(6,6) * t119);
t225 = t242 - t233 - t271;
t217 = t202 * t227;
t216 = (t17 * t197 + t316) * t202;
t215 = (t18 * t195 + t19 * t197) * t202;
t214 = t311 * t202 + t265;
t211 = (Icges(6,1) * t118 - t308 - t68) * t197 + (-Icges(6,1) * t116 - t103 - t66) * t195;
t94 = rSges(6,1) * t118 - rSges(6,2) * t119;
t93 = -rSges(6,1) * t116 - rSges(6,2) * t117;
t76 = qJD(1) * t96 - t185;
t75 = qJD(1) * t244 + t184;
t47 = qJD(1) * t51 - t185 + t262;
t46 = qJD(1) * t225 + t288;
t34 = -t275 + qJDD(1) * t101 + qJD(1) * (-qJD(1) * t273 + t289) + t213;
t33 = t244 * qJDD(1) + (-qJD(1) * t101 - t121) * qJD(1) + t224;
t32 = t112 * t299 + t113 * t118 + t114 * t119;
t30 = t32 * t182;
t29 = -qJD(4) * t204 + t228 * t279 + qJD(2);
t15 = -t275 + qJDD(1) * t82 + qJD(1) * (-rSges(5,3) * t264 - t293) + t209;
t14 = t225 * qJDD(1) + ((-rSges(5,3) * t283 - t281) * t202 + t234 + t309) * qJD(1) + t218;
t13 = t113 * t85 + t114 * t86 - t116 * t125 + t117 * t126 + (t112 * t283 + t124 * t195) * t202;
t12 = t113 * t83 + t114 * t84 + t118 * t125 + t119 * t126 + (-t112 * t284 + t124 * t197) * t202;
t11 = -t134 * t74 + t135 * t72 + t229 * t279 + t180;
t6 = qJD(5) * t215 + t30;
t5 = qJD(5) * t216 + t315;
t1 = [(t30 + ((t16 + t19 - t346) * t197 + t337 * t195) * t279) * t250 + t32 * t330 + t31 * t331 + t321 + t313 / 0.2e1 + t312 / 0.2e1 - m(2) * (-g(1) * t172 + g(2) * t173) + (t10 + t12) * t247 + ((-t153 * t208 - g(2) + t246) * t147 + (-t348 + (-0.2e1 * t155 - t198 + t147) * t208 - g(1)) * t146) * m(3) + (t13 + t9 + t6) * t249 + (-t315 + ((-t17 + t337 - t345) * t197 - t316) * t279 + t5) * t248 + (-(-t23 + (t79 - t327) * qJD(1) + t251 + t349) * t24 - t237 * t324 + t23 * (t185 + t240) + t24 * (t267 + t290 + t318) + (-t7 * pkin(2) + (-t23 * qJD(4) + t310 * t7) * t202) * t195 + ((-t206 * t24 - t207 * t23) * pkin(1) + (t237 - t304) * t314 + (t23 * (-qJ(3) - t325) + t24 * (-rSges(6,3) * t202 - pkin(2) - t304)) * t195) * qJD(1) + (-g(2) + t8) * (t223 + t336 + t74) + (-g(1) + t7) * (t245 + t253 + t294)) * m(6) + (-g(1) * (t187 + t243) - t214 * t324 + t14 * (-t233 + t253) + t46 * (t185 + t234) + t47 * (t267 - t293) + (t14 * t265 + (-t46 * qJD(4) + t14 * t311) * t202) * t195 + ((-t206 * t47 - t207 * t46) * pkin(1) + t46 * t214 * t197 + (-t46 * qJ(3) + t47 * (-rSges(5,3) * t202 - pkin(2) - t231)) * t195) * qJD(1) - (-t46 + (t243 - t271) * qJD(1) + t251) * t47 + (-g(2) + t15) * t51) * m(5) + (t75 * t185 + t76 * (t286 + t289) + ((-t206 * t76 - t207 * t75) * pkin(1) + t75 * (rSges(4,2) * t202 + t259) * t197 + (t75 * (-rSges(4,3) - qJ(3)) + t76 * t259) * t195) * qJD(1) - (-t75 + (-t100 - t327) * qJD(1) + t285) * t76 + (t34 - g(2)) * t96 + (t33 - g(1)) * (t259 * t195 + t253 + t287)) * m(4) + (m(3) * (t146 ^ 2 + t155 * t147) + Icges(2,3) + Icges(3,3) + m(2) * (t172 ^ 2 + t173 ^ 2) + (Icges(5,3) + Icges(4,2)) * t204 ^ 2 + ((Icges(5,1) * t203 ^ 2 + (-0.2e1 * Icges(5,4) * t203 + Icges(5,2) * t201) * t201 + Icges(4,1)) * t202 + 0.2e1 * (-Icges(5,5) * t203 + Icges(5,6) * t201 + Icges(4,4)) * t204) * t202) * qJDD(1); (m(3) + m(4)) * qJDD(2) + m(5) * t180 + m(6) * t11 + (-m(3) + t278) * g(3); t278 * (-g(2) * t197 + t324) + 0.2e1 * (t323 / 0.2e1 - t322 / 0.2e1) * m(6) + 0.2e1 * (t14 * t329 + t15 * t328) * m(5) + 0.2e1 * (t328 * t34 + t329 * t33) * m(4); t332 * (-g(3) * t204 + (g(1) * t197 + g(2) * t195) * t202) + m(5) * (t14 * t299 + t15 * t302 - t180 * t204) + m(6) * (-t11 * t204 + t299 * t7 + t302 * t8); -t6 * t264 / 0.2e1 + (-t204 * t32 + t215) * t330 + (-t12 * t204 + ((t18 * t197 - t19 * t195) * qJD(1) + t239) * t202) * t247 + (t13 * t182 + t134 * t16 + t135 * t17 + t181 * t31 + t238 * t279) * t302 / 0.2e1 + (-t204 * t31 + t216) * t331 + (-t13 * t204 + ((t16 * t197 - t17 * t195) * qJD(1) + t238) * t202) * t249 - t204 * (t236 * t279 + t312 + t313 + t321) / 0.2e1 + t181 * (-t204 * t50 + (t195 * t25 + t197 * t26) * t202) / 0.2e1 + t182 * (-t204 * t28 + ((-t195 * t26 + t25 * t197) * qJD(1) + t236) * t202) / 0.2e1 + ((t118 * t291 + t119 * t292 + t140 * t299) * t182 + (t118 * t334 + t211 * t119 + t197 * t217) * t279) * t248 + ((-t116 * t291 + t117 * t292 + t140 * t302) * t182 + (-t116 * t334 + t117 * t211 + t195 * t217) * t279) * t250 - t182 * (-t204 * t140 * t182 + ((-t194 * t291 + t196 * t292) * t182 + ((-t194 * t334 + t196 * t211) * t202 - t227 * t204) * qJD(5)) * t202) / 0.2e1 + (qJD(1) * t5 + t12 * t182 + t134 * t18 + t135 * t19 + t181 * t32 + t239 * t279) * t299 / 0.2e1 + ((t23 * t45 - t24 * t44 + t7 * t72 - t74 * t8) * t204 + (t11 * t228 + t29 * (-t283 * t74 - t284 * t72 + t229) + t230 * t127 + (t323 - t322 + (t195 * t24 + t314) * qJD(1)) * t115) * t202 - (-t23 * t93 + t24 * t94) * t182 - (t29 * (-t195 * t94 + t197 * t93) + t230 * t144) * t279 - g(1) * t94 - g(2) * t93 - g(3) * t144) * m(6);];
tau = t1;
