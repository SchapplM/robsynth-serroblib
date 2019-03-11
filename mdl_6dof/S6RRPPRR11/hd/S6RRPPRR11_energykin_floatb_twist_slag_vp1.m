% Calculate kinetic energy for
% S6RRPPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-03-09 09:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRR11_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR11_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPPRR11_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR11_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR11_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRR11_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:38:57
% EndTime: 2019-03-09 09:39:01
% DurationCPUTime: 3.90s
% Computational Cost: add. (2542->396), mult. (4576->553), div. (0->0), fcn. (5250->12), ass. (0->179)
t368 = Icges(3,1) + Icges(4,2);
t367 = Icges(3,4) + Icges(4,6);
t366 = Icges(3,5) - Icges(4,4);
t365 = Icges(3,2) + Icges(4,3);
t364 = Icges(3,6) - Icges(4,5);
t363 = Icges(3,3) + Icges(4,1);
t304 = cos(pkin(6));
t308 = sin(qJ(1));
t310 = cos(qJ(2));
t336 = t308 * t310;
t307 = sin(qJ(2));
t311 = cos(qJ(1));
t337 = t307 * t311;
t268 = t304 * t336 + t337;
t335 = t310 * t311;
t338 = t307 * t308;
t269 = -t304 * t338 + t335;
t302 = sin(pkin(6));
t341 = t302 * t308;
t362 = t365 * t268 - t367 * t269 - t364 * t341;
t266 = -t304 * t335 + t338;
t267 = t304 * t337 + t336;
t339 = t302 * t311;
t361 = t365 * t266 - t367 * t267 + t364 * t339;
t360 = -t367 * t268 + t368 * t269 + t366 * t341;
t359 = -t367 * t266 + t368 * t267 - t366 * t339;
t358 = -t364 * t268 + t366 * t269 + t363 * t341;
t357 = -t364 * t266 + t366 * t267 - t363 * t339;
t356 = t363 * t304 + (t366 * t307 + t364 * t310) * t302;
t355 = t364 * t304 + (t367 * t307 + t365 * t310) * t302;
t354 = t366 * t304 + (t368 * t307 + t367 * t310) * t302;
t301 = sin(pkin(11));
t303 = cos(pkin(11));
t232 = t268 * t303 - t301 * t341;
t343 = t268 * t301;
t233 = t303 * t341 + t343;
t174 = Icges(5,5) * t233 + Icges(5,6) * t232 + Icges(5,3) * t269;
t353 = t174 + t360;
t234 = t266 * t303 + t301 * t339;
t344 = t266 * t301;
t235 = -t303 * t339 + t344;
t175 = Icges(5,5) * t235 + Icges(5,6) * t234 + Icges(5,3) * t267;
t352 = t175 + t359;
t340 = t302 * t310;
t264 = -t301 * t304 - t303 * t340;
t265 = -t301 * t340 + t303 * t304;
t342 = t302 * t307;
t199 = Icges(5,5) * t265 + Icges(5,6) * t264 + Icges(5,3) * t342;
t351 = t199 + t354;
t347 = pkin(8) * t304;
t346 = pkin(4) * t303;
t345 = Icges(2,4) * t308;
t224 = pkin(2) * t267 + qJ(3) * t266;
t243 = -pkin(3) * t339 + t267 * qJ(4);
t333 = -t224 - t243;
t225 = pkin(2) * t269 + qJ(3) * t268;
t242 = pkin(3) * t341 + qJ(4) * t269;
t332 = -t225 - t242;
t270 = (pkin(2) * t307 - qJ(3) * t310) * t302;
t272 = pkin(3) * t304 + qJ(4) * t342;
t331 = -t270 - t272;
t330 = qJD(2) * t302;
t329 = pkin(11) + qJ(5);
t328 = V_base(5) * pkin(7) + V_base(1);
t279 = t308 * t330 + V_base(4);
t298 = V_base(6) + qJD(1);
t325 = cos(t329);
t237 = qJD(5) * t269 + t279;
t280 = qJD(2) * t304 + t298;
t324 = t302 * t325;
t262 = qJD(5) * t342 + t280;
t278 = -t311 * t330 + V_base(5);
t273 = t308 * pkin(1) - pkin(8) * t339;
t323 = -t273 * t298 + V_base(5) * t347 + t328;
t274 = pkin(1) * t311 + pkin(8) * t341;
t322 = V_base(4) * t273 - t274 * V_base(5) + V_base(3);
t236 = qJD(5) * t267 + t278;
t321 = qJD(3) * t268 + t278 * t270 + t323;
t320 = t298 * t274 + V_base(2) + (-pkin(7) - t347) * V_base(4);
t319 = qJD(4) * t269 + t278 * t272 + t321;
t318 = qJD(3) * t266 + t280 * t225 + t320;
t317 = -qJD(3) * t340 + t279 * t224 + t322;
t316 = qJD(4) * t267 + t280 * t242 + t318;
t315 = qJD(4) * t342 + t279 * t243 + t317;
t183 = pkin(4) * t344 + pkin(9) * t267 - t339 * t346;
t220 = t346 * t304 + (-pkin(4) * t301 * t310 + pkin(9) * t307) * t302;
t314 = t278 * t220 + (-t183 + t333) * t280 + t319;
t182 = pkin(4) * t343 + pkin(9) * t269 + t341 * t346;
t313 = t280 * t182 + (-t220 + t331) * t279 + t316;
t312 = t279 * t183 + (-t182 + t332) * t278 + t315;
t309 = cos(qJ(6));
t306 = sin(qJ(6));
t299 = Icges(2,4) * t311;
t297 = sin(t329);
t288 = rSges(2,1) * t311 - t308 * rSges(2,2);
t287 = t308 * rSges(2,1) + rSges(2,2) * t311;
t286 = Icges(2,1) * t311 - t345;
t285 = Icges(2,1) * t308 + t299;
t284 = -Icges(2,2) * t308 + t299;
t283 = Icges(2,2) * t311 + t345;
t277 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t276 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t275 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t255 = rSges(4,1) * t304 + (-rSges(4,2) * t307 - rSges(4,3) * t310) * t302;
t254 = rSges(3,3) * t304 + (rSges(3,1) * t307 + rSges(3,2) * t310) * t302;
t253 = -t297 * t340 + t304 * t325;
t252 = t304 * t297 + t310 * t324;
t241 = V_base(5) * rSges(2,3) - t287 * t298 + t328;
t240 = t288 * t298 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t238 = t287 * V_base(4) - t288 * V_base(5) + V_base(3);
t231 = t266 * t297 - t311 * t324;
t230 = t266 * t325 + t297 * t339;
t229 = t268 * t297 + t308 * t324;
t228 = -t268 * t325 + t297 * t341;
t227 = t253 * t309 + t306 * t342;
t226 = -t253 * t306 + t309 * t342;
t221 = qJD(6) * t252 + t262;
t219 = rSges(3,1) * t269 - rSges(3,2) * t268 + rSges(3,3) * t341;
t218 = t267 * rSges(3,1) - t266 * rSges(3,2) - rSges(3,3) * t339;
t217 = -rSges(4,1) * t339 - t267 * rSges(4,2) + t266 * rSges(4,3);
t216 = rSges(4,1) * t341 - rSges(4,2) * t269 + rSges(4,3) * t268;
t203 = pkin(5) * t253 + pkin(10) * t252;
t202 = rSges(5,1) * t265 + rSges(5,2) * t264 + rSges(5,3) * t342;
t201 = Icges(5,1) * t265 + Icges(5,4) * t264 + Icges(5,5) * t342;
t200 = Icges(5,4) * t265 + Icges(5,2) * t264 + Icges(5,6) * t342;
t196 = rSges(6,1) * t253 - rSges(6,2) * t252 + rSges(6,3) * t342;
t195 = Icges(6,1) * t253 - Icges(6,4) * t252 + Icges(6,5) * t342;
t194 = Icges(6,4) * t253 - Icges(6,2) * t252 + Icges(6,6) * t342;
t193 = Icges(6,5) * t253 - Icges(6,6) * t252 + Icges(6,3) * t342;
t191 = t231 * t309 + t267 * t306;
t190 = -t231 * t306 + t267 * t309;
t189 = t229 * t309 + t269 * t306;
t188 = -t229 * t306 + t269 * t309;
t187 = qJD(6) * t228 + t237;
t186 = -qJD(6) * t230 + t236;
t185 = pkin(5) * t231 - pkin(10) * t230;
t184 = pkin(5) * t229 + pkin(10) * t228;
t181 = rSges(5,1) * t235 + rSges(5,2) * t234 + rSges(5,3) * t267;
t180 = rSges(5,1) * t233 + rSges(5,2) * t232 + rSges(5,3) * t269;
t179 = Icges(5,1) * t235 + Icges(5,4) * t234 + Icges(5,5) * t267;
t178 = Icges(5,1) * t233 + Icges(5,4) * t232 + Icges(5,5) * t269;
t177 = Icges(5,4) * t235 + Icges(5,2) * t234 + Icges(5,6) * t267;
t176 = Icges(5,4) * t233 + Icges(5,2) * t232 + Icges(5,6) * t269;
t172 = rSges(6,1) * t231 + rSges(6,2) * t230 + rSges(6,3) * t267;
t171 = rSges(6,1) * t229 - rSges(6,2) * t228 + rSges(6,3) * t269;
t170 = Icges(6,1) * t231 + Icges(6,4) * t230 + Icges(6,5) * t267;
t169 = Icges(6,1) * t229 - Icges(6,4) * t228 + Icges(6,5) * t269;
t168 = Icges(6,4) * t231 + Icges(6,2) * t230 + Icges(6,6) * t267;
t167 = Icges(6,4) * t229 - Icges(6,2) * t228 + Icges(6,6) * t269;
t166 = Icges(6,5) * t231 + Icges(6,6) * t230 + Icges(6,3) * t267;
t165 = Icges(6,5) * t229 - Icges(6,6) * t228 + Icges(6,3) * t269;
t163 = rSges(7,1) * t227 + rSges(7,2) * t226 + rSges(7,3) * t252;
t162 = Icges(7,1) * t227 + Icges(7,4) * t226 + Icges(7,5) * t252;
t161 = Icges(7,4) * t227 + Icges(7,2) * t226 + Icges(7,6) * t252;
t160 = Icges(7,5) * t227 + Icges(7,6) * t226 + Icges(7,3) * t252;
t159 = -t218 * t280 + t254 * t278 + t323;
t158 = t219 * t280 - t254 * t279 + t320;
t157 = t218 * t279 - t219 * t278 + t322;
t156 = rSges(7,1) * t191 + rSges(7,2) * t190 - rSges(7,3) * t230;
t155 = rSges(7,1) * t189 + rSges(7,2) * t188 + rSges(7,3) * t228;
t154 = Icges(7,1) * t191 + Icges(7,4) * t190 - Icges(7,5) * t230;
t153 = Icges(7,1) * t189 + Icges(7,4) * t188 + Icges(7,5) * t228;
t152 = Icges(7,4) * t191 + Icges(7,2) * t190 - Icges(7,6) * t230;
t151 = Icges(7,4) * t189 + Icges(7,2) * t188 + Icges(7,6) * t228;
t150 = Icges(7,5) * t191 + Icges(7,6) * t190 - Icges(7,3) * t230;
t149 = Icges(7,5) * t189 + Icges(7,6) * t188 + Icges(7,3) * t228;
t148 = t255 * t278 + (-t217 - t224) * t280 + t321;
t147 = t216 * t280 + (-t255 - t270) * t279 + t318;
t146 = t217 * t279 + (-t216 - t225) * t278 + t317;
t145 = t202 * t278 + (-t181 + t333) * t280 + t319;
t144 = t180 * t280 + (-t202 + t331) * t279 + t316;
t143 = t181 * t279 + (-t180 + t332) * t278 + t315;
t142 = -t172 * t262 + t196 * t236 + t314;
t141 = t171 * t262 - t196 * t237 + t313;
t140 = -t171 * t236 + t172 * t237 + t312;
t139 = -t156 * t221 + t163 * t186 - t185 * t262 + t203 * t236 + t314;
t138 = t155 * t221 - t163 * t187 + t184 * t262 - t203 * t237 + t313;
t137 = -t155 * t186 + t156 * t187 - t184 * t236 + t185 * t237 + t312;
t1 = m(1) * (t275 ^ 2 + t276 ^ 2 + t277 ^ 2) / 0.2e1 + t237 * ((t269 * t165 - t228 * t167 + t229 * t169) * t237 + (t166 * t269 - t168 * t228 + t170 * t229) * t236 + (t193 * t269 - t194 * t228 + t195 * t229) * t262) / 0.2e1 + t236 * ((t165 * t267 + t167 * t230 + t169 * t231) * t237 + (t267 * t166 + t230 * t168 + t231 * t170) * t236 + (t193 * t267 + t194 * t230 + t195 * t231) * t262) / 0.2e1 + t221 * ((t149 * t252 + t151 * t226 + t153 * t227) * t187 + (t150 * t252 + t152 * t226 + t154 * t227) * t186 + (t252 * t160 + t226 * t161 + t227 * t162) * t221) / 0.2e1 + m(2) * (t238 ^ 2 + t240 ^ 2 + t241 ^ 2) / 0.2e1 + t186 * ((-t149 * t230 + t151 * t190 + t153 * t191) * t187 + (-t230 * t150 + t190 * t152 + t191 * t154) * t186 + (-t160 * t230 + t161 * t190 + t162 * t191) * t221) / 0.2e1 + t187 * ((t228 * t149 + t188 * t151 + t189 * t153) * t187 + (t150 * t228 + t152 * t188 + t154 * t189) * t186 + (t160 * t228 + t161 * t188 + t162 * t189) * t221) / 0.2e1 + m(3) * (t157 ^ 2 + t158 ^ 2 + t159 ^ 2) / 0.2e1 + m(5) * (t143 ^ 2 + t144 ^ 2 + t145 ^ 2) / 0.2e1 + m(4) * (t146 ^ 2 + t147 ^ 2 + t148 ^ 2) / 0.2e1 + m(6) * (t140 ^ 2 + t141 ^ 2 + t142 ^ 2) / 0.2e1 + m(7) * (t137 ^ 2 + t138 ^ 2 + t139 ^ 2) / 0.2e1 + t262 * ((t165 * t342 - t167 * t252 + t169 * t253) * t237 + (t166 * t342 - t168 * t252 + t170 * t253) * t236 + (t193 * t342 - t252 * t194 + t253 * t195) * t262) / 0.2e1 + ((-t308 * t283 + t285 * t311 + Icges(1,4)) * V_base(5) + (-t308 * t284 + t286 * t311 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t283 * t311 + t308 * t285 + Icges(1,2)) * V_base(5) + (t284 * t311 + t308 * t286 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t200 * t234 + t201 * t235 - t266 * t355 + t267 * t351 - t339 * t356) * t280 + (t176 * t234 + t178 * t235 + t266 * t362 + t353 * t267 - t358 * t339) * t279 + (t234 * t177 + t235 * t179 + t266 * t361 + t267 * t352 - t339 * t357) * t278) * t278 / 0.2e1 + ((t200 * t232 + t201 * t233 - t268 * t355 + t269 * t351 + t341 * t356) * t280 + (t232 * t176 + t233 * t178 + t268 * t362 + t353 * t269 + t358 * t341) * t279 + (t177 * t232 + t179 * t233 + t268 * t361 + t269 * t352 + t341 * t357) * t278) * t279 / 0.2e1 + ((t174 * t342 + t176 * t264 + t178 * t265) * t279 + (t175 * t342 + t177 * t264 + t179 * t265) * t278 + (t199 * t342 + t264 * t200 + t265 * t201) * t280 + (t278 * t357 + t279 * t358 + t280 * t356) * t304 + ((t307 * t354 + t310 * t355) * t280 + (t360 * t307 - t310 * t362) * t279 + (t307 * t359 - t310 * t361) * t278) * t302) * t280 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t308 + Icges(2,6) * t311) * V_base(5) + (Icges(2,5) * t311 - Icges(2,6) * t308) * V_base(4) + Icges(2,3) * t298 / 0.2e1) * t298;
T  = t1;
