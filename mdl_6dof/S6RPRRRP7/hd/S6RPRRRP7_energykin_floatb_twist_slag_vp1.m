% Calculate kinetic energy for
% S6RPRRRP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRP7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP7_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP7_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRRP7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP7_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP7_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP7_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRP7_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:18:33
% EndTime: 2019-03-09 06:18:37
% DurationCPUTime: 3.52s
% Computational Cost: add. (2228->336), mult. (2294->486), div. (0->0), fcn. (2198->10), ass. (0->171)
t347 = Icges(6,1) + Icges(7,1);
t346 = -Icges(6,4) + Icges(7,5);
t345 = Icges(7,4) + Icges(6,5);
t344 = Icges(6,2) + Icges(7,3);
t343 = -Icges(7,6) + Icges(6,6);
t342 = -Icges(6,3) - Icges(7,2);
t341 = rSges(7,1) + pkin(5);
t340 = rSges(7,3) + qJ(6);
t262 = pkin(10) + qJ(3);
t254 = cos(t262);
t263 = qJ(4) + qJ(5);
t259 = cos(t263);
t270 = cos(qJ(1));
t313 = t259 * t270;
t258 = sin(t263);
t268 = sin(qJ(1));
t315 = t258 * t268;
t199 = t254 * t315 + t313;
t310 = t268 * t259;
t314 = t258 * t270;
t200 = t254 * t310 - t314;
t253 = sin(t262);
t317 = t253 * t268;
t339 = t344 * t199 + t346 * t200 - t343 * t317;
t201 = t254 * t314 - t310;
t202 = t254 * t313 + t315;
t316 = t253 * t270;
t338 = t344 * t201 + t346 * t202 - t343 * t316;
t337 = -t343 * t199 + t345 * t200 - t342 * t317;
t336 = -t343 * t201 + t345 * t202 - t342 * t316;
t335 = t346 * t199 + t347 * t200 + t345 * t317;
t334 = t346 * t201 + t347 * t202 + t345 * t316;
t333 = t343 * t254 + (t344 * t258 + t346 * t259) * t253;
t332 = t342 * t254 + (-t343 * t258 + t345 * t259) * t253;
t331 = -t345 * t254 + (t346 * t258 + t347 * t259) * t253;
t264 = sin(pkin(10));
t326 = pkin(2) * t264;
t265 = cos(pkin(10));
t325 = pkin(2) * t265;
t269 = cos(qJ(4));
t324 = pkin(4) * t269;
t322 = Icges(2,4) * t268;
t321 = Icges(3,4) * t264;
t320 = Icges(3,4) * t265;
t319 = Icges(4,4) * t253;
t318 = Icges(4,4) * t254;
t267 = sin(qJ(4));
t312 = t267 * t268;
t311 = t267 * t270;
t309 = t268 * t269;
t308 = t269 * t270;
t306 = rSges(7,2) * t317 + t340 * t199 + t341 * t200;
t305 = rSges(7,2) * t316 + t340 * t201 + t341 * t202;
t304 = -rSges(7,2) * t254 + (t340 * t258 + t341 * t259) * t253;
t186 = -pkin(7) * t270 + t268 * t325;
t243 = pkin(1) * t268 - qJ(2) * t270;
t303 = -t186 - t243;
t302 = qJD(4) * t253;
t301 = qJD(5) * t253;
t300 = V_base(4) * t243 + V_base(3);
t299 = V_base(5) * pkin(6) + V_base(1);
t248 = qJD(3) * t268 + V_base(4);
t255 = V_base(6) + qJD(1);
t296 = qJD(2) * t268 + t299;
t215 = t270 * t302 + t248;
t295 = V_base(5) * t326 + t296;
t294 = pkin(3) * t254 + pkin(8) * t253;
t247 = -qJD(3) * t270 + V_base(5);
t293 = rSges(3,1) * t265 - rSges(3,2) * t264;
t292 = rSges(4,1) * t254 - rSges(4,2) * t253;
t291 = Icges(3,1) * t265 - t321;
t290 = Icges(4,1) * t254 - t319;
t289 = -Icges(3,2) * t264 + t320;
t288 = -Icges(4,2) * t253 + t318;
t287 = Icges(3,5) * t265 - Icges(3,6) * t264;
t286 = Icges(4,5) * t254 - Icges(4,6) * t253;
t245 = pkin(1) * t270 + qJ(2) * t268;
t285 = -qJD(2) * t270 + t255 * t245 + V_base(2);
t214 = t268 * t302 + t247;
t284 = pkin(9) * t253 + t254 * t324;
t283 = (-Icges(4,3) * t270 + t268 * t286) * t247 + (Icges(4,3) * t268 + t270 * t286) * t248 + (Icges(4,5) * t253 + Icges(4,6) * t254) * t255;
t187 = pkin(7) * t268 + t270 * t325;
t282 = V_base(4) * t186 + (-t187 - t245) * V_base(5) + t300;
t281 = (-Icges(3,3) * t270 + t268 * t287) * V_base(5) + (Icges(3,3) * t268 + t270 * t287) * V_base(4) + (Icges(3,5) * t264 + Icges(3,6) * t265) * t255;
t211 = t294 * t268;
t225 = t253 * pkin(3) - t254 * pkin(8);
t280 = t247 * t225 + (-t211 + t303) * t255 + t295;
t212 = t294 * t270;
t279 = t248 * t211 - t212 * t247 + t282;
t278 = t255 * t187 + (-pkin(6) - t326) * V_base(4) + t285;
t153 = -pkin(4) * t311 + t268 * t284;
t165 = -pkin(9) * t254 + t253 * t324;
t227 = -qJD(4) * t254 + t255;
t277 = -t153 * t227 + t214 * t165 + t280;
t154 = pkin(4) * t312 + t270 * t284;
t276 = t215 * t153 - t154 * t214 + t279;
t275 = t255 * t212 - t225 * t248 + t278;
t274 = t227 * t154 - t165 * t215 + t275;
t191 = -Icges(4,6) * t270 + t268 * t288;
t192 = Icges(4,6) * t268 + t270 * t288;
t193 = -Icges(4,5) * t270 + t268 * t290;
t194 = Icges(4,5) * t268 + t270 * t290;
t222 = Icges(4,2) * t254 + t319;
t223 = Icges(4,1) * t253 + t318;
t273 = (-t192 * t253 + t194 * t254) * t248 + (-t191 * t253 + t193 * t254) * t247 + (-t222 * t253 + t223 * t254) * t255;
t205 = -Icges(3,6) * t270 + t268 * t289;
t206 = Icges(3,6) * t268 + t270 * t289;
t207 = -Icges(3,5) * t270 + t268 * t291;
t208 = Icges(3,5) * t268 + t270 * t291;
t232 = Icges(3,2) * t265 + t321;
t233 = Icges(3,1) * t264 + t320;
t272 = (-t206 * t264 + t208 * t265) * V_base(4) + (-t205 * t264 + t207 * t265) * V_base(5) + (-t232 * t264 + t233 * t265) * t255;
t260 = Icges(2,4) * t270;
t246 = rSges(2,1) * t270 - rSges(2,2) * t268;
t244 = rSges(2,1) * t268 + rSges(2,2) * t270;
t240 = Icges(2,1) * t270 - t322;
t239 = Icges(2,1) * t268 + t260;
t238 = -Icges(2,2) * t268 + t260;
t237 = Icges(2,2) * t270 + t322;
t236 = Icges(2,5) * t270 - Icges(2,6) * t268;
t235 = Icges(2,5) * t268 + Icges(2,6) * t270;
t234 = rSges(3,1) * t264 + rSges(3,2) * t265;
t230 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t229 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t228 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t224 = rSges(4,1) * t253 + rSges(4,2) * t254;
t219 = t254 * t308 + t312;
t218 = -t254 * t311 + t309;
t217 = t254 * t309 - t311;
t216 = -t254 * t312 - t308;
t213 = (-qJD(4) - qJD(5)) * t254 + t255;
t210 = rSges(3,3) * t268 + t270 * t293;
t209 = -rSges(3,3) * t270 + t268 * t293;
t196 = rSges(4,3) * t268 + t270 * t292;
t195 = -rSges(4,3) * t270 + t268 * t292;
t185 = -rSges(5,3) * t254 + (rSges(5,1) * t269 - rSges(5,2) * t267) * t253;
t184 = V_base(5) * rSges(2,3) - t244 * t255 + t299;
t183 = t246 * t255 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t182 = -Icges(5,5) * t254 + (Icges(5,1) * t269 - Icges(5,4) * t267) * t253;
t181 = -Icges(5,6) * t254 + (Icges(5,4) * t269 - Icges(5,2) * t267) * t253;
t180 = -Icges(5,3) * t254 + (Icges(5,5) * t269 - Icges(5,6) * t267) * t253;
t179 = t270 * t301 + t215;
t178 = t268 * t301 + t214;
t176 = t244 * V_base(4) - t246 * V_base(5) + V_base(3);
t174 = -rSges(6,3) * t254 + (rSges(6,1) * t259 - rSges(6,2) * t258) * t253;
t162 = rSges(5,1) * t219 + rSges(5,2) * t218 + rSges(5,3) * t316;
t161 = rSges(5,1) * t217 + rSges(5,2) * t216 + rSges(5,3) * t317;
t160 = Icges(5,1) * t219 + Icges(5,4) * t218 + Icges(5,5) * t316;
t159 = Icges(5,1) * t217 + Icges(5,4) * t216 + Icges(5,5) * t317;
t158 = Icges(5,4) * t219 + Icges(5,2) * t218 + Icges(5,6) * t316;
t157 = Icges(5,4) * t217 + Icges(5,2) * t216 + Icges(5,6) * t317;
t156 = Icges(5,5) * t219 + Icges(5,6) * t218 + Icges(5,3) * t316;
t155 = Icges(5,5) * t217 + Icges(5,6) * t216 + Icges(5,3) * t317;
t151 = rSges(6,1) * t202 - rSges(6,2) * t201 + rSges(6,3) * t316;
t149 = rSges(6,1) * t200 - rSges(6,2) * t199 + rSges(6,3) * t317;
t134 = t234 * V_base(5) + (-t209 - t243) * t255 + t296;
t133 = t210 * t255 + (-pkin(6) - t234) * V_base(4) + t285;
t131 = t209 * V_base(4) + (-t210 - t245) * V_base(5) + t300;
t130 = t224 * t247 + (-t195 + t303) * t255 + t295;
t129 = t196 * t255 - t224 * t248 + t278;
t128 = t195 * t248 - t196 * t247 + t282;
t127 = -t161 * t227 + t185 * t214 + t280;
t126 = t162 * t227 - t185 * t215 + t275;
t125 = t161 * t215 - t162 * t214 + t279;
t124 = -t149 * t213 + t174 * t178 + t277;
t123 = t151 * t213 - t174 * t179 + t274;
t122 = t149 * t179 - t151 * t178 + t276;
t121 = qJD(6) * t201 + t178 * t304 - t213 * t306 + t277;
t120 = qJD(6) * t199 - t179 * t304 + t213 * t305 + t274;
t119 = qJD(6) * t253 * t258 - t178 * t305 + t179 * t306 + t276;
t1 = t248 * (t268 * t283 + t270 * t273) / 0.2e1 + t247 * (t268 * t273 - t283 * t270) / 0.2e1 + t227 * ((-t155 * t214 - t156 * t215 - t180 * t227) * t254 + ((-t158 * t267 + t160 * t269) * t215 + (-t157 * t267 + t159 * t269) * t214 + (-t181 * t267 + t182 * t269) * t227) * t253) / 0.2e1 + m(1) * (t228 ^ 2 + t229 ^ 2 + t230 ^ 2) / 0.2e1 + m(2) * (t176 ^ 2 + t183 ^ 2 + t184 ^ 2) / 0.2e1 + m(3) * (t131 ^ 2 + t133 ^ 2 + t134 ^ 2) / 0.2e1 + m(6) * (t122 ^ 2 + t123 ^ 2 + t124 ^ 2) / 0.2e1 + m(5) * (t125 ^ 2 + t126 ^ 2 + t127 ^ 2) / 0.2e1 + m(4) * (t128 ^ 2 + t129 ^ 2 + t130 ^ 2) / 0.2e1 + m(7) * (t119 ^ 2 + t120 ^ 2 + t121 ^ 2) / 0.2e1 + t215 * ((t156 * t316 + t218 * t158 + t219 * t160) * t215 + (t155 * t316 + t157 * t218 + t159 * t219) * t214 + (t180 * t316 + t181 * t218 + t182 * t219) * t227) / 0.2e1 + t214 * ((t156 * t317 + t158 * t216 + t160 * t217) * t215 + (t155 * t317 + t216 * t157 + t217 * t159) * t214 + (t180 * t317 + t181 * t216 + t182 * t217) * t227) / 0.2e1 + ((t333 * t199 + t331 * t200 + t332 * t317) * t213 + (t338 * t199 + t334 * t200 + t336 * t317) * t179 + (t339 * t199 + t335 * t200 + t337 * t317) * t178) * t178 / 0.2e1 + ((t333 * t201 + t331 * t202 + t332 * t316) * t213 + (t338 * t201 + t334 * t202 + t336 * t316) * t179 + (t339 * t201 + t335 * t202 + t337 * t316) * t178) * t179 / 0.2e1 + ((-t337 * t178 - t336 * t179 - t332 * t213) * t254 + ((t333 * t258 + t331 * t259) * t213 + (t338 * t258 + t334 * t259) * t179 + (t339 * t258 + t335 * t259) * t178) * t253) * t213 / 0.2e1 + ((t192 * t254 + t194 * t253) * t248 + (t191 * t254 + t193 * t253) * t247 + (t205 * t265 + t207 * t264 + t235) * V_base(5) + (t206 * t265 + t208 * t264 + t236) * V_base(4) + (t254 * t222 + t253 * t223 + t265 * t232 + t264 * t233 + Icges(2,3)) * t255) * t255 / 0.2e1 + (t236 * t255 + t268 * t281 + t270 * t272 + (-t237 * t268 + t239 * t270 + Icges(1,4)) * V_base(5) + (-t268 * t238 + t270 * t240 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t235 * t255 + t268 * t272 - t270 * t281 + (t270 * t237 + t268 * t239 + Icges(1,2)) * V_base(5) + (t238 * t270 + t240 * t268 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
