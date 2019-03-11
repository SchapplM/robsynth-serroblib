% Calculate kinetic energy for
% S6RRPRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-03-09 11:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPR12_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR12_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRPR12_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR12_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR12_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR12_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:17:10
% EndTime: 2019-03-09 11:17:14
% DurationCPUTime: 4.32s
% Computational Cost: add. (2587->391), mult. (4666->544), div. (0->0), fcn. (5358->12), ass. (0->172)
t370 = Icges(3,1) + Icges(4,2);
t369 = Icges(3,4) + Icges(4,6);
t368 = Icges(3,5) - Icges(4,4);
t367 = Icges(3,2) + Icges(4,3);
t366 = Icges(3,6) - Icges(4,5);
t365 = Icges(3,3) + Icges(4,1);
t364 = Icges(5,3) + Icges(6,3);
t302 = cos(pkin(6));
t307 = sin(qJ(1));
t310 = cos(qJ(2));
t332 = t307 * t310;
t306 = sin(qJ(2));
t311 = cos(qJ(1));
t333 = t306 * t311;
t268 = t302 * t332 + t333;
t329 = qJ(4) + pkin(11);
t297 = sin(t329);
t325 = cos(t329);
t301 = sin(pkin(6));
t338 = t301 * t307;
t228 = -t268 * t325 + t297 * t338;
t324 = t301 * t325;
t229 = t268 * t297 + t307 * t324;
t305 = sin(qJ(4));
t309 = cos(qJ(4));
t234 = t268 * t309 - t305 * t338;
t340 = t268 * t305;
t235 = t309 * t338 + t340;
t331 = t310 * t311;
t334 = t306 * t307;
t269 = -t302 * t334 + t331;
t363 = Icges(5,5) * t235 + Icges(6,5) * t229 + Icges(5,6) * t234 - Icges(6,6) * t228 + t364 * t269;
t266 = -t302 * t331 + t334;
t336 = t301 * t311;
t230 = t266 * t325 + t297 * t336;
t231 = t266 * t297 - t311 * t324;
t236 = t266 * t309 + t305 * t336;
t341 = t266 * t305;
t237 = -t309 * t336 + t341;
t267 = t302 * t333 + t332;
t362 = Icges(5,5) * t237 + Icges(6,5) * t231 + Icges(5,6) * t236 + Icges(6,6) * t230 + t364 * t267;
t252 = t302 * t297 + t310 * t324;
t337 = t301 * t310;
t253 = -t297 * t337 + t302 * t325;
t264 = -t302 * t305 - t309 * t337;
t335 = t305 * t310;
t265 = -t301 * t335 + t302 * t309;
t339 = t301 * t306;
t361 = Icges(5,5) * t265 + Icges(6,5) * t253 + Icges(5,6) * t264 - Icges(6,6) * t252 + t364 * t339;
t360 = t367 * t268 - t369 * t269 - t366 * t338;
t359 = t367 * t266 - t369 * t267 + t366 * t336;
t358 = -t369 * t268 + t370 * t269 + t368 * t338;
t357 = -t369 * t266 + t370 * t267 - t368 * t336;
t356 = -t366 * t268 + t368 * t269 + t365 * t338;
t355 = -t366 * t266 + t368 * t267 - t365 * t336;
t354 = t365 * t302 + (t368 * t306 + t366 * t310) * t301;
t353 = t366 * t302 + (t369 * t306 + t367 * t310) * t301;
t352 = t368 * t302 + (t370 * t306 + t369 * t310) * t301;
t345 = pkin(8) * t302;
t344 = pkin(4) * t309;
t342 = Icges(2,4) * t307;
t330 = qJD(2) * t301;
t328 = V_base(5) * pkin(7) + V_base(1);
t279 = t307 * t330 + V_base(4);
t298 = V_base(6) + qJD(1);
t233 = qJD(4) * t269 + t279;
t280 = qJD(2) * t302 + t298;
t262 = qJD(4) * t339 + t280;
t278 = -t311 * t330 + V_base(5);
t273 = t307 * pkin(1) - pkin(8) * t336;
t323 = -t273 * t298 + V_base(5) * t345 + t328;
t274 = pkin(1) * t311 + pkin(8) * t338;
t322 = V_base(4) * t273 - t274 * V_base(5) + V_base(3);
t232 = qJD(4) * t267 + t278;
t270 = (pkin(2) * t306 - qJ(3) * t310) * t301;
t321 = qJD(3) * t268 + t278 * t270 + t323;
t320 = t298 * t274 + V_base(2) + (-pkin(7) - t345) * V_base(4);
t225 = pkin(2) * t269 + qJ(3) * t268;
t319 = qJD(3) * t266 + t280 * t225 + t320;
t224 = pkin(2) * t267 + qJ(3) * t266;
t318 = -qJD(3) * t337 + t279 * t224 + t322;
t243 = -pkin(3) * t336 + t267 * pkin(9);
t272 = pkin(3) * t302 + pkin(9) * t339;
t317 = t278 * t272 + (-t224 - t243) * t280 + t321;
t220 = t344 * t302 + (-pkin(4) * t335 + qJ(5) * t306) * t301;
t316 = qJD(5) * t269 + t232 * t220 + t317;
t242 = pkin(3) * t338 + pkin(9) * t269;
t315 = t280 * t242 + (-t270 - t272) * t279 + t319;
t314 = t279 * t243 + (-t225 - t242) * t278 + t318;
t182 = pkin(4) * t340 + qJ(5) * t269 + t338 * t344;
t313 = qJD(5) * t267 + t262 * t182 + t315;
t183 = pkin(4) * t341 + qJ(5) * t267 - t336 * t344;
t312 = qJD(5) * t339 + t233 * t183 + t314;
t308 = cos(qJ(6));
t304 = sin(qJ(6));
t299 = Icges(2,4) * t311;
t288 = rSges(2,1) * t311 - t307 * rSges(2,2);
t287 = t307 * rSges(2,1) + rSges(2,2) * t311;
t286 = Icges(2,1) * t311 - t342;
t285 = Icges(2,1) * t307 + t299;
t284 = -Icges(2,2) * t307 + t299;
t283 = Icges(2,2) * t311 + t342;
t277 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t276 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t275 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t255 = rSges(4,1) * t302 + (-rSges(4,2) * t306 - rSges(4,3) * t310) * t301;
t254 = rSges(3,3) * t302 + (rSges(3,1) * t306 + rSges(3,2) * t310) * t301;
t241 = V_base(5) * rSges(2,3) - t287 * t298 + t328;
t240 = t288 * t298 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t238 = t287 * V_base(4) - t288 * V_base(5) + V_base(3);
t227 = t253 * t308 + t304 * t339;
t226 = -t253 * t304 + t308 * t339;
t221 = qJD(6) * t252 + t262;
t219 = rSges(3,1) * t269 - rSges(3,2) * t268 + rSges(3,3) * t338;
t218 = t267 * rSges(3,1) - t266 * rSges(3,2) - rSges(3,3) * t336;
t217 = -rSges(4,1) * t336 - t267 * rSges(4,2) + t266 * rSges(4,3);
t216 = rSges(4,1) * t338 - rSges(4,2) * t269 + rSges(4,3) * t268;
t203 = rSges(5,1) * t265 + rSges(5,2) * t264 + rSges(5,3) * t339;
t202 = Icges(5,1) * t265 + Icges(5,4) * t264 + Icges(5,5) * t339;
t201 = Icges(5,4) * t265 + Icges(5,2) * t264 + Icges(5,6) * t339;
t199 = pkin(5) * t253 + pkin(10) * t252;
t196 = rSges(6,1) * t253 - rSges(6,2) * t252 + rSges(6,3) * t339;
t195 = Icges(6,1) * t253 - Icges(6,4) * t252 + Icges(6,5) * t339;
t194 = Icges(6,4) * t253 - Icges(6,2) * t252 + Icges(6,6) * t339;
t192 = t231 * t308 + t267 * t304;
t191 = -t231 * t304 + t267 * t308;
t190 = t229 * t308 + t269 * t304;
t189 = -t229 * t304 + t269 * t308;
t188 = qJD(6) * t228 + t233;
t187 = -qJD(6) * t230 + t232;
t186 = pkin(5) * t231 - pkin(10) * t230;
t185 = pkin(5) * t229 + pkin(10) * t228;
t181 = rSges(5,1) * t237 + rSges(5,2) * t236 + rSges(5,3) * t267;
t180 = rSges(5,1) * t235 + rSges(5,2) * t234 + rSges(5,3) * t269;
t179 = Icges(5,1) * t237 + Icges(5,4) * t236 + Icges(5,5) * t267;
t178 = Icges(5,1) * t235 + Icges(5,4) * t234 + Icges(5,5) * t269;
t177 = Icges(5,4) * t237 + Icges(5,2) * t236 + Icges(5,6) * t267;
t176 = Icges(5,4) * t235 + Icges(5,2) * t234 + Icges(5,6) * t269;
t173 = rSges(6,1) * t231 + rSges(6,2) * t230 + rSges(6,3) * t267;
t172 = rSges(6,1) * t229 - rSges(6,2) * t228 + rSges(6,3) * t269;
t171 = Icges(6,1) * t231 + Icges(6,4) * t230 + Icges(6,5) * t267;
t170 = Icges(6,1) * t229 - Icges(6,4) * t228 + Icges(6,5) * t269;
t169 = Icges(6,4) * t231 + Icges(6,2) * t230 + Icges(6,6) * t267;
t168 = Icges(6,4) * t229 - Icges(6,2) * t228 + Icges(6,6) * t269;
t165 = rSges(7,1) * t227 + rSges(7,2) * t226 + rSges(7,3) * t252;
t164 = Icges(7,1) * t227 + Icges(7,4) * t226 + Icges(7,5) * t252;
t163 = Icges(7,4) * t227 + Icges(7,2) * t226 + Icges(7,6) * t252;
t162 = Icges(7,5) * t227 + Icges(7,6) * t226 + Icges(7,3) * t252;
t159 = -t218 * t280 + t254 * t278 + t323;
t158 = t219 * t280 - t254 * t279 + t320;
t157 = t218 * t279 - t219 * t278 + t322;
t156 = rSges(7,1) * t192 + rSges(7,2) * t191 - rSges(7,3) * t230;
t155 = rSges(7,1) * t190 + rSges(7,2) * t189 + rSges(7,3) * t228;
t154 = Icges(7,1) * t192 + Icges(7,4) * t191 - Icges(7,5) * t230;
t153 = Icges(7,1) * t190 + Icges(7,4) * t189 + Icges(7,5) * t228;
t152 = Icges(7,4) * t192 + Icges(7,2) * t191 - Icges(7,6) * t230;
t151 = Icges(7,4) * t190 + Icges(7,2) * t189 + Icges(7,6) * t228;
t150 = Icges(7,5) * t192 + Icges(7,6) * t191 - Icges(7,3) * t230;
t149 = Icges(7,5) * t190 + Icges(7,6) * t189 + Icges(7,3) * t228;
t148 = t255 * t278 + (-t217 - t224) * t280 + t321;
t147 = t216 * t280 + (-t255 - t270) * t279 + t319;
t146 = t217 * t279 + (-t216 - t225) * t278 + t318;
t145 = -t181 * t262 + t203 * t232 + t317;
t144 = t180 * t262 - t203 * t233 + t315;
t143 = -t180 * t232 + t181 * t233 + t314;
t142 = t196 * t232 + (-t173 - t183) * t262 + t316;
t141 = t172 * t262 + (-t196 - t220) * t233 + t313;
t140 = t173 * t233 + (-t172 - t182) * t232 + t312;
t139 = -t156 * t221 + t165 * t187 + t199 * t232 + (-t183 - t186) * t262 + t316;
t138 = t155 * t221 - t165 * t188 + t185 * t262 + (-t199 - t220) * t233 + t313;
t137 = -t155 * t187 + t156 * t188 + t186 * t233 + (-t182 - t185) * t232 + t312;
t1 = t187 * ((-t149 * t230 + t151 * t191 + t153 * t192) * t188 + (-t230 * t150 + t191 * t152 + t192 * t154) * t187 + (-t162 * t230 + t163 * t191 + t164 * t192) * t221) / 0.2e1 + t188 * ((t228 * t149 + t189 * t151 + t190 * t153) * t188 + (t150 * t228 + t152 * t189 + t154 * t190) * t187 + (t162 * t228 + t163 * t189 + t164 * t190) * t221) / 0.2e1 + m(3) * (t157 ^ 2 + t158 ^ 2 + t159 ^ 2) / 0.2e1 + m(4) * (t146 ^ 2 + t147 ^ 2 + t148 ^ 2) / 0.2e1 + m(6) * (t140 ^ 2 + t141 ^ 2 + t142 ^ 2) / 0.2e1 + m(5) * (t143 ^ 2 + t144 ^ 2 + t145 ^ 2) / 0.2e1 + m(7) * (t137 ^ 2 + t138 ^ 2 + t139 ^ 2) / 0.2e1 + m(1) * (t275 ^ 2 + t276 ^ 2 + t277 ^ 2) / 0.2e1 + m(2) * (t238 ^ 2 + t240 ^ 2 + t241 ^ 2) / 0.2e1 + t221 * ((t149 * t252 + t151 * t226 + t153 * t227) * t188 + (t150 * t252 + t152 * t226 + t154 * t227) * t187 + (t252 * t162 + t226 * t163 + t227 * t164) * t221) / 0.2e1 + ((t194 * t230 + t195 * t231 + t201 * t236 + t202 * t237 + t267 * t361) * t262 + (t168 * t230 + t170 * t231 + t176 * t236 + t178 * t237 + t267 * t363) * t233 + (t230 * t169 + t231 * t171 + t236 * t177 + t237 * t179 + t362 * t267) * t232) * t232 / 0.2e1 + ((-t194 * t228 + t195 * t229 + t201 * t234 + t202 * t235 + t269 * t361) * t262 + (-t228 * t168 + t229 * t170 + t234 * t176 + t235 * t178 + t363 * t269) * t233 + (-t169 * t228 + t171 * t229 + t177 * t234 + t179 * t235 + t269 * t362) * t232) * t233 / 0.2e1 + ((-t252 * t194 + t253 * t195 + t264 * t201 + t265 * t202 + t361 * t339) * t262 + (-t168 * t252 + t170 * t253 + t176 * t264 + t178 * t265 + t339 * t363) * t233 + (-t169 * t252 + t171 * t253 + t177 * t264 + t179 * t265 + t339 * t362) * t232) * t262 / 0.2e1 + ((-t266 * t353 + t267 * t352 - t336 * t354) * t280 + (t266 * t360 + t267 * t358 - t336 * t356) * t279 + (t359 * t266 + t357 * t267 - t355 * t336) * t278) * t278 / 0.2e1 + ((-t268 * t353 + t269 * t352 + t338 * t354) * t280 + (t360 * t268 + t358 * t269 + t356 * t338) * t279 + (t268 * t359 + t269 * t357 + t338 * t355) * t278) * t279 / 0.2e1 + ((t278 * t355 + t279 * t356 + t280 * t354) * t302 + ((t306 * t352 + t310 * t353) * t280 + (t306 * t358 - t310 * t360) * t279 + (t306 * t357 - t310 * t359) * t278) * t301) * t280 / 0.2e1 + ((-t307 * t283 + t285 * t311 + Icges(1,4)) * V_base(5) + (-t307 * t284 + t311 * t286 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t283 * t311 + t307 * t285 + Icges(1,2)) * V_base(5) + (t284 * t311 + t307 * t286 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t307 + Icges(2,6) * t311) * V_base(5) + (Icges(2,5) * t311 - Icges(2,6) * t307) * V_base(4) + Icges(2,3) * t298 / 0.2e1) * t298;
T  = t1;
