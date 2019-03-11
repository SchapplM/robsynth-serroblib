% Calculate kinetic energy for
% S6RRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRP4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP4_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP4_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRRP4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP4_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP4_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP4_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP4_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:52:29
% EndTime: 2019-03-09 11:52:32
% DurationCPUTime: 3.39s
% Computational Cost: add. (2254->330), mult. (2320->472), div. (0->0), fcn. (2224->10), ass. (0->166)
t352 = Icges(6,1) + Icges(7,1);
t351 = -Icges(6,4) + Icges(7,5);
t350 = Icges(7,4) + Icges(6,5);
t349 = Icges(6,2) + Icges(7,3);
t348 = -Icges(7,6) + Icges(6,6);
t347 = Icges(3,3) + Icges(4,3);
t346 = -Icges(6,3) - Icges(7,2);
t262 = qJ(2) + pkin(10);
t253 = sin(t262);
t254 = cos(t262);
t266 = sin(qJ(2));
t269 = cos(qJ(2));
t345 = Icges(3,5) * t269 + Icges(4,5) * t254 - Icges(3,6) * t266 - Icges(4,6) * t253;
t344 = rSges(7,1) + pkin(5);
t343 = rSges(7,3) + qJ(6);
t263 = qJ(4) + qJ(5);
t259 = cos(t263);
t270 = cos(qJ(1));
t311 = t259 * t270;
t258 = sin(t263);
t267 = sin(qJ(1));
t313 = t258 * t267;
t199 = t254 * t313 + t311;
t308 = t267 * t259;
t312 = t258 * t270;
t200 = t254 * t308 - t312;
t315 = t253 * t267;
t342 = t349 * t199 + t351 * t200 - t348 * t315;
t201 = t254 * t312 - t308;
t202 = t254 * t311 + t313;
t314 = t253 * t270;
t341 = t349 * t201 + t351 * t202 - t348 * t314;
t340 = -t348 * t199 + t350 * t200 - t346 * t315;
t339 = -t348 * t201 + t350 * t202 - t346 * t314;
t338 = t351 * t199 + t352 * t200 + t350 * t315;
t337 = t351 * t201 + t352 * t202 + t350 * t314;
t336 = t348 * t254 + (t349 * t258 + t351 * t259) * t253;
t335 = t346 * t254 + (-t348 * t258 + t350 * t259) * t253;
t334 = -t350 * t254 + (t351 * t258 + t352 * t259) * t253;
t316 = Icges(4,4) * t254;
t289 = -Icges(4,2) * t253 + t316;
t191 = -Icges(4,6) * t270 + t267 * t289;
t192 = Icges(4,6) * t267 + t270 * t289;
t317 = Icges(4,4) * t253;
t291 = Icges(4,1) * t254 - t317;
t193 = -Icges(4,5) * t270 + t267 * t291;
t194 = Icges(4,5) * t267 + t270 * t291;
t318 = Icges(3,4) * t269;
t290 = -Icges(3,2) * t266 + t318;
t205 = -Icges(3,6) * t270 + t267 * t290;
t206 = Icges(3,6) * t267 + t270 * t290;
t319 = Icges(3,4) * t266;
t292 = Icges(3,1) * t269 - t319;
t207 = -Icges(3,5) * t270 + t267 * t292;
t208 = Icges(3,5) * t267 + t270 * t292;
t222 = Icges(4,2) * t254 + t317;
t223 = Icges(4,1) * t253 + t316;
t235 = Icges(3,2) * t269 + t319;
t238 = Icges(3,1) * t266 + t318;
t248 = -qJD(2) * t270 + V_base(5);
t249 = qJD(2) * t267 + V_base(4);
t255 = V_base(6) + qJD(1);
t333 = (-t222 * t253 + t223 * t254 - t235 * t266 + t238 * t269) * t255 + (-t192 * t253 + t194 * t254 - t206 * t266 + t208 * t269) * t249 + (-t191 * t253 + t193 * t254 - t205 * t266 + t207 * t269) * t248;
t332 = (Icges(3,5) * t266 + Icges(4,5) * t253 + Icges(3,6) * t269 + Icges(4,6) * t254) * t255 + (t347 * t267 + t345 * t270) * t249 + (t345 * t267 - t347 * t270) * t248;
t325 = pkin(2) * t266;
t324 = pkin(2) * t269;
t268 = cos(qJ(4));
t323 = pkin(4) * t268;
t320 = Icges(2,4) * t267;
t265 = sin(qJ(4));
t310 = t265 * t267;
t309 = t265 * t270;
t307 = t267 * t268;
t306 = t268 * t270;
t305 = rSges(7,2) * t315 + t343 * t199 + t200 * t344;
t304 = rSges(7,2) * t314 + t343 * t201 + t202 * t344;
t303 = -rSges(7,2) * t254 + (t343 * t258 + t259 * t344) * t253;
t186 = -qJ(3) * t270 + t267 * t324;
t246 = pkin(1) * t267 - pkin(7) * t270;
t302 = -t186 - t246;
t301 = qJD(4) * t253;
t300 = qJD(5) * t253;
t299 = V_base(5) * pkin(6) + V_base(1);
t215 = t270 * t301 + t249;
t296 = qJD(3) * t267 + t248 * t325 + t299;
t295 = pkin(3) * t254 + pkin(8) * t253;
t294 = rSges(3,1) * t269 - rSges(3,2) * t266;
t293 = rSges(4,1) * t254 - rSges(4,2) * t253;
t214 = t267 * t301 + t248;
t247 = pkin(1) * t270 + pkin(7) * t267;
t286 = -V_base(4) * pkin(6) + t255 * t247 + V_base(2);
t285 = V_base(4) * t246 - t247 * V_base(5) + V_base(3);
t284 = t249 * t186 + t285;
t283 = pkin(9) * t253 + t254 * t323;
t187 = qJ(3) * t267 + t270 * t324;
t280 = -qJD(3) * t270 + t255 * t187 + t286;
t211 = t295 * t267;
t225 = t253 * pkin(3) - t254 * pkin(8);
t279 = t248 * t225 + (-t211 + t302) * t255 + t296;
t212 = t295 * t270;
t278 = t249 * t211 + (-t187 - t212) * t248 + t284;
t153 = -pkin(4) * t309 + t267 * t283;
t165 = -pkin(9) * t254 + t253 * t323;
t228 = -qJD(4) * t254 + t255;
t277 = -t153 * t228 + t214 * t165 + t279;
t154 = pkin(4) * t310 + t270 * t283;
t276 = t215 * t153 - t154 * t214 + t278;
t275 = t255 * t212 + (-t225 - t325) * t249 + t280;
t274 = t228 * t154 - t165 * t215 + t275;
t260 = Icges(2,4) * t270;
t245 = rSges(2,1) * t270 - rSges(2,2) * t267;
t244 = rSges(2,1) * t267 + rSges(2,2) * t270;
t243 = rSges(3,1) * t266 + rSges(3,2) * t269;
t240 = Icges(2,1) * t270 - t320;
t239 = Icges(2,1) * t267 + t260;
t237 = -Icges(2,2) * t267 + t260;
t236 = Icges(2,2) * t270 + t320;
t231 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t230 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t229 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t224 = rSges(4,1) * t253 + rSges(4,2) * t254;
t219 = t254 * t306 + t310;
t218 = -t254 * t309 + t307;
t217 = t254 * t307 - t309;
t216 = -t254 * t310 - t306;
t213 = (-qJD(4) - qJD(5)) * t254 + t255;
t210 = rSges(3,3) * t267 + t270 * t294;
t209 = -rSges(3,3) * t270 + t267 * t294;
t196 = rSges(4,3) * t267 + t270 * t293;
t195 = -rSges(4,3) * t270 + t267 * t293;
t185 = -rSges(5,3) * t254 + (rSges(5,1) * t268 - rSges(5,2) * t265) * t253;
t184 = V_base(5) * rSges(2,3) - t244 * t255 + t299;
t183 = t245 * t255 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t182 = -Icges(5,5) * t254 + (Icges(5,1) * t268 - Icges(5,4) * t265) * t253;
t181 = -Icges(5,6) * t254 + (Icges(5,4) * t268 - Icges(5,2) * t265) * t253;
t180 = -Icges(5,3) * t254 + (Icges(5,5) * t268 - Icges(5,6) * t265) * t253;
t179 = t270 * t300 + t215;
t178 = t267 * t300 + t214;
t176 = t244 * V_base(4) - t245 * V_base(5) + V_base(3);
t174 = -rSges(6,3) * t254 + (rSges(6,1) * t259 - rSges(6,2) * t258) * t253;
t162 = rSges(5,1) * t219 + rSges(5,2) * t218 + rSges(5,3) * t314;
t161 = rSges(5,1) * t217 + rSges(5,2) * t216 + rSges(5,3) * t315;
t160 = Icges(5,1) * t219 + Icges(5,4) * t218 + Icges(5,5) * t314;
t159 = Icges(5,1) * t217 + Icges(5,4) * t216 + Icges(5,5) * t315;
t158 = Icges(5,4) * t219 + Icges(5,2) * t218 + Icges(5,6) * t314;
t157 = Icges(5,4) * t217 + Icges(5,2) * t216 + Icges(5,6) * t315;
t156 = Icges(5,5) * t219 + Icges(5,6) * t218 + Icges(5,3) * t314;
t155 = Icges(5,5) * t217 + Icges(5,6) * t216 + Icges(5,3) * t315;
t151 = rSges(6,1) * t202 - rSges(6,2) * t201 + rSges(6,3) * t314;
t149 = rSges(6,1) * t200 - rSges(6,2) * t199 + rSges(6,3) * t315;
t134 = t243 * t248 + (-t209 - t246) * t255 + t299;
t133 = t210 * t255 - t243 * t249 + t286;
t131 = t209 * t249 - t210 * t248 + t285;
t130 = t224 * t248 + (-t195 + t302) * t255 + t296;
t129 = t196 * t255 + (-t224 - t325) * t249 + t280;
t128 = t195 * t249 + (-t187 - t196) * t248 + t284;
t127 = -t161 * t228 + t185 * t214 + t279;
t126 = t162 * t228 - t185 * t215 + t275;
t125 = t161 * t215 - t162 * t214 + t278;
t124 = -t149 * t213 + t174 * t178 + t277;
t123 = t151 * t213 - t174 * t179 + t274;
t122 = t149 * t179 - t151 * t178 + t276;
t121 = qJD(6) * t201 + t178 * t303 - t213 * t305 + t277;
t120 = qJD(6) * t199 - t179 * t303 + t213 * t304 + t274;
t119 = qJD(6) * t253 * t258 - t178 * t304 + t179 * t305 + t276;
t1 = t215 * ((t156 * t314 + t218 * t158 + t219 * t160) * t215 + (t155 * t314 + t218 * t157 + t219 * t159) * t214 + (t180 * t314 + t218 * t181 + t219 * t182) * t228) / 0.2e1 + t214 * ((t156 * t315 + t216 * t158 + t217 * t160) * t215 + (t155 * t315 + t216 * t157 + t217 * t159) * t214 + (t180 * t315 + t216 * t181 + t217 * t182) * t228) / 0.2e1 + t228 * ((-t155 * t214 - t156 * t215 - t180 * t228) * t254 + ((-t158 * t265 + t160 * t268) * t215 + (-t157 * t265 + t159 * t268) * t214 + (-t181 * t265 + t182 * t268) * t228) * t253) / 0.2e1 + m(1) * (t229 ^ 2 + t230 ^ 2 + t231 ^ 2) / 0.2e1 + m(2) * (t176 ^ 2 + t183 ^ 2 + t184 ^ 2) / 0.2e1 + m(3) * (t131 ^ 2 + t133 ^ 2 + t134 ^ 2) / 0.2e1 + m(5) * (t125 ^ 2 + t126 ^ 2 + t127 ^ 2) / 0.2e1 + m(4) * (t128 ^ 2 + t129 ^ 2 + t130 ^ 2) / 0.2e1 + m(7) * (t119 ^ 2 + t120 ^ 2 + t121 ^ 2) / 0.2e1 + m(6) * (t122 ^ 2 + t123 ^ 2 + t124 ^ 2) / 0.2e1 + ((t199 * t336 + t200 * t334 + t315 * t335) * t213 + (t199 * t341 + t200 * t337 + t315 * t339) * t179 + (t342 * t199 + t338 * t200 + t340 * t315) * t178) * t178 / 0.2e1 + ((t201 * t336 + t202 * t334 + t314 * t335) * t213 + (t341 * t201 + t337 * t202 + t339 * t314) * t179 + (t201 * t342 + t338 * t202 + t340 * t314) * t178) * t179 / 0.2e1 + ((-t178 * t340 - t179 * t339 - t213 * t335) * t254 + ((t258 * t336 + t259 * t334) * t213 + (t258 * t341 + t259 * t337) * t179 + (t258 * t342 + t338 * t259) * t178) * t253) * t213 / 0.2e1 + (t333 * t267 - t332 * t270) * t248 / 0.2e1 + (t332 * t267 + t333 * t270) * t249 / 0.2e1 + ((-t236 * t267 + t239 * t270 + Icges(1,4)) * V_base(5) + (-t267 * t237 + t270 * t240 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t270 * t236 + t267 * t239 + Icges(1,2)) * V_base(5) + (t237 * t270 + t240 * t267 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t254 * t192 + t253 * t194 + t269 * t206 + t266 * t208) * t249 + (t254 * t191 + t253 * t193 + t269 * t205 + t266 * t207) * t248 + (t254 * t222 + t253 * t223 + t269 * t235 + t266 * t238 + Icges(2,3)) * t255) * t255 / 0.2e1 + t255 * V_base(4) * (Icges(2,5) * t270 - Icges(2,6) * t267) + t255 * V_base(5) * (Icges(2,5) * t267 + Icges(2,6) * t270) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
