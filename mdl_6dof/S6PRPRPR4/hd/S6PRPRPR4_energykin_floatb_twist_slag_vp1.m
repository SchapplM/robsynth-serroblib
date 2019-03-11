% Calculate kinetic energy for
% S6PRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2019-03-08 19:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRPR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR4_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRPRPR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR4_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR4_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRPR4_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:38:56
% EndTime: 2019-03-08 19:39:01
% DurationCPUTime: 4.60s
% Computational Cost: add. (3467->437), mult. (5732->607), div. (0->0), fcn. (6800->14), ass. (0->189)
t368 = Icges(5,2) + Icges(6,3);
t305 = sin(pkin(10));
t309 = cos(pkin(10));
t314 = cos(qJ(2));
t310 = cos(pkin(6));
t313 = sin(qJ(2));
t342 = t310 * t313;
t266 = t305 * t314 + t309 * t342;
t337 = pkin(11) + qJ(4);
t298 = sin(t337);
t328 = cos(t337);
t306 = sin(pkin(6));
t345 = t306 * t309;
t239 = t266 * t328 - t298 * t345;
t341 = t310 * t314;
t265 = t305 * t313 - t309 * t341;
t303 = sin(pkin(12));
t307 = cos(pkin(12));
t204 = -t239 * t303 + t265 * t307;
t349 = t265 * t303;
t205 = t239 * t307 + t349;
t327 = t306 * t328;
t238 = t266 * t298 + t309 * t327;
t367 = -Icges(5,4) * t239 + Icges(6,5) * t205 - Icges(5,6) * t265 + Icges(6,6) * t204 + t368 * t238;
t268 = -t305 * t342 + t309 * t314;
t346 = t305 * t306;
t241 = t268 * t328 + t298 * t346;
t267 = t305 * t341 + t309 * t313;
t206 = -t241 * t303 + t267 * t307;
t348 = t267 * t303;
t207 = t241 * t307 + t348;
t240 = t268 * t298 - t305 * t327;
t366 = -Icges(5,4) * t241 + Icges(6,5) * t207 - Icges(5,6) * t267 + Icges(6,6) * t206 + t368 * t240;
t257 = t310 * t298 + t313 * t327;
t343 = t306 * t314;
t236 = -t257 * t303 - t307 * t343;
t331 = t303 * t343;
t237 = t257 * t307 - t331;
t344 = t306 * t313;
t256 = t298 * t344 - t310 * t328;
t365 = -Icges(5,4) * t257 + Icges(6,5) * t237 + Icges(5,6) * t343 + Icges(6,6) * t236 + t368 * t256;
t304 = sin(pkin(11));
t308 = cos(pkin(11));
t243 = -t266 * t304 - t308 * t345;
t329 = t304 * t345;
t244 = t266 * t308 - t329;
t187 = Icges(4,5) * t244 + Icges(4,6) * t243 + Icges(4,3) * t265;
t221 = Icges(3,4) * t266 - Icges(3,2) * t265 - Icges(3,6) * t345;
t364 = t187 - t221;
t245 = -t268 * t304 + t308 * t346;
t330 = t304 * t346;
t246 = t268 * t308 + t330;
t188 = Icges(4,5) * t246 + Icges(4,6) * t245 + Icges(4,3) * t267;
t222 = Icges(3,4) * t268 - Icges(3,2) * t267 + Icges(3,6) * t346;
t363 = t188 - t222;
t263 = -t304 * t344 + t308 * t310;
t347 = t304 * t310;
t264 = t308 * t344 + t347;
t215 = Icges(4,5) * t264 + Icges(4,6) * t263 - Icges(4,3) * t343;
t254 = Icges(3,6) * t310 + (Icges(3,4) * t313 + Icges(3,2) * t314) * t306;
t362 = t215 - t254;
t353 = pkin(7) * t310;
t352 = pkin(3) * t308;
t351 = pkin(5) * t307;
t350 = Icges(2,4) * t305;
t338 = qJD(2) * t306;
t336 = V_base(5) * qJ(1) + V_base(1);
t332 = qJD(1) + V_base(3);
t280 = t305 * t338 + V_base(4);
t291 = qJD(2) * t310 + V_base(6);
t248 = qJD(4) * t267 + t280;
t279 = -t309 * t338 + V_base(5);
t247 = qJD(4) * t265 + t279;
t269 = -qJD(4) * t343 + t291;
t273 = pkin(1) * t305 - pkin(7) * t345;
t326 = -t273 * V_base(6) + V_base(5) * t353 + t336;
t274 = pkin(1) * t309 + pkin(7) * t346;
t325 = V_base(4) * t273 - V_base(5) * t274 + t332;
t324 = V_base(6) * t274 + V_base(2) + (-qJ(1) - t353) * V_base(4);
t272 = (pkin(2) * t313 - qJ(3) * t314) * t306;
t323 = qJD(3) * t267 + t279 * t272 + t326;
t233 = pkin(2) * t268 + qJ(3) * t267;
t322 = qJD(3) * t265 + t291 * t233 + t324;
t232 = pkin(2) * t266 + qJ(3) * t265;
t321 = -qJD(3) * t343 + t280 * t232 + t325;
t184 = -pkin(3) * t329 + pkin(8) * t265 + t266 * t352;
t226 = pkin(3) * t347 + (-pkin(8) * t314 + t313 * t352) * t306;
t320 = t279 * t226 + (-t184 - t232) * t291 + t323;
t185 = pkin(3) * t330 + pkin(8) * t267 + t268 * t352;
t319 = t291 * t185 + (-t226 - t272) * t280 + t322;
t218 = t257 * pkin(4) + t256 * qJ(5);
t318 = qJD(5) * t240 + t247 * t218 + t320;
t317 = t280 * t184 + (-t185 - t233) * t279 + t321;
t197 = pkin(4) * t241 + qJ(5) * t240;
t316 = qJD(5) * t238 + t269 * t197 + t319;
t196 = pkin(4) * t239 + qJ(5) * t238;
t315 = qJD(5) * t256 + t248 * t196 + t317;
t302 = pkin(12) + qJ(6);
t300 = Icges(2,4) * t309;
t299 = cos(t302);
t297 = sin(t302);
t288 = rSges(2,1) * t309 - rSges(2,2) * t305;
t287 = rSges(2,1) * t305 + rSges(2,2) * t309;
t286 = Icges(2,1) * t309 - t350;
t285 = Icges(2,1) * t305 + t300;
t284 = -Icges(2,2) * t305 + t300;
t283 = Icges(2,2) * t309 + t350;
t278 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t277 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t276 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t258 = t310 * rSges(3,3) + (rSges(3,1) * t313 + rSges(3,2) * t314) * t306;
t255 = Icges(3,5) * t310 + (Icges(3,1) * t313 + Icges(3,4) * t314) * t306;
t253 = Icges(3,3) * t310 + (Icges(3,5) * t313 + Icges(3,6) * t314) * t306;
t251 = V_base(5) * rSges(2,3) - t287 * V_base(6) + t336;
t250 = t288 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t242 = t287 * V_base(4) - t288 * V_base(5) + t332;
t231 = t257 * t299 - t297 * t343;
t230 = -t257 * t297 - t299 * t343;
t229 = qJD(6) * t256 + t269;
t228 = rSges(3,1) * t268 - rSges(3,2) * t267 + rSges(3,3) * t346;
t227 = rSges(3,1) * t266 - rSges(3,2) * t265 - rSges(3,3) * t345;
t225 = t264 * rSges(4,1) + t263 * rSges(4,2) - rSges(4,3) * t343;
t224 = Icges(3,1) * t268 - Icges(3,4) * t267 + Icges(3,5) * t346;
t223 = Icges(3,1) * t266 - Icges(3,4) * t265 - Icges(3,5) * t345;
t220 = Icges(3,5) * t268 - Icges(3,6) * t267 + Icges(3,3) * t346;
t219 = Icges(3,5) * t266 - Icges(3,6) * t265 - Icges(3,3) * t345;
t217 = Icges(4,1) * t264 + Icges(4,4) * t263 - Icges(4,5) * t343;
t216 = Icges(4,4) * t264 + Icges(4,2) * t263 - Icges(4,6) * t343;
t212 = t257 * rSges(5,1) - t256 * rSges(5,2) - rSges(5,3) * t343;
t211 = Icges(5,1) * t257 - Icges(5,4) * t256 - Icges(5,5) * t343;
t209 = Icges(5,5) * t257 - Icges(5,6) * t256 - Icges(5,3) * t343;
t203 = t241 * t299 + t267 * t297;
t202 = -t241 * t297 + t267 * t299;
t201 = t239 * t299 + t265 * t297;
t200 = -t239 * t297 + t265 * t299;
t199 = qJD(6) * t240 + t248;
t198 = qJD(6) * t238 + t247;
t194 = rSges(4,1) * t246 + rSges(4,2) * t245 + rSges(4,3) * t267;
t193 = rSges(4,1) * t244 + rSges(4,2) * t243 + rSges(4,3) * t265;
t192 = Icges(4,1) * t246 + Icges(4,4) * t245 + Icges(4,5) * t267;
t191 = Icges(4,1) * t244 + Icges(4,4) * t243 + Icges(4,5) * t265;
t190 = Icges(4,4) * t246 + Icges(4,2) * t245 + Icges(4,6) * t267;
t189 = Icges(4,4) * t244 + Icges(4,2) * t243 + Icges(4,6) * t265;
t183 = rSges(5,1) * t241 - rSges(5,2) * t240 + rSges(5,3) * t267;
t182 = rSges(5,1) * t239 - rSges(5,2) * t238 + rSges(5,3) * t265;
t181 = Icges(5,1) * t241 - Icges(5,4) * t240 + Icges(5,5) * t267;
t180 = Icges(5,1) * t239 - Icges(5,4) * t238 + Icges(5,5) * t265;
t177 = Icges(5,5) * t241 - Icges(5,6) * t240 + Icges(5,3) * t267;
t176 = Icges(5,5) * t239 - Icges(5,6) * t238 + Icges(5,3) * t265;
t175 = rSges(6,1) * t237 + rSges(6,2) * t236 + rSges(6,3) * t256;
t174 = Icges(6,1) * t237 + Icges(6,4) * t236 + Icges(6,5) * t256;
t173 = Icges(6,4) * t237 + Icges(6,2) * t236 + Icges(6,6) * t256;
t169 = rSges(7,1) * t231 + rSges(7,2) * t230 + rSges(7,3) * t256;
t167 = Icges(7,1) * t231 + Icges(7,4) * t230 + Icges(7,5) * t256;
t166 = Icges(7,4) * t231 + Icges(7,2) * t230 + Icges(7,6) * t256;
t165 = Icges(7,5) * t231 + Icges(7,6) * t230 + Icges(7,3) * t256;
t164 = -pkin(5) * t331 + pkin(9) * t256 + t257 * t351;
t163 = -t227 * t291 + t258 * t279 + t326;
t162 = t228 * t291 - t258 * t280 + t324;
t161 = t227 * t280 - t228 * t279 + t325;
t160 = rSges(6,1) * t207 + rSges(6,2) * t206 + rSges(6,3) * t240;
t159 = rSges(6,1) * t205 + rSges(6,2) * t204 + rSges(6,3) * t238;
t158 = Icges(6,1) * t207 + Icges(6,4) * t206 + Icges(6,5) * t240;
t157 = Icges(6,1) * t205 + Icges(6,4) * t204 + Icges(6,5) * t238;
t156 = Icges(6,4) * t207 + Icges(6,2) * t206 + Icges(6,6) * t240;
t155 = Icges(6,4) * t205 + Icges(6,2) * t204 + Icges(6,6) * t238;
t152 = rSges(7,1) * t203 + rSges(7,2) * t202 + rSges(7,3) * t240;
t151 = rSges(7,1) * t201 + rSges(7,2) * t200 + rSges(7,3) * t238;
t150 = Icges(7,1) * t203 + Icges(7,4) * t202 + Icges(7,5) * t240;
t149 = Icges(7,1) * t201 + Icges(7,4) * t200 + Icges(7,5) * t238;
t148 = Icges(7,4) * t203 + Icges(7,2) * t202 + Icges(7,6) * t240;
t147 = Icges(7,4) * t201 + Icges(7,2) * t200 + Icges(7,6) * t238;
t146 = Icges(7,5) * t203 + Icges(7,6) * t202 + Icges(7,3) * t240;
t145 = Icges(7,5) * t201 + Icges(7,6) * t200 + Icges(7,3) * t238;
t144 = pkin(5) * t348 + pkin(9) * t240 + t241 * t351;
t143 = pkin(5) * t349 + pkin(9) * t238 + t239 * t351;
t142 = t225 * t279 + (-t193 - t232) * t291 + t323;
t141 = t194 * t291 + (-t225 - t272) * t280 + t322;
t140 = t280 * t193 + (-t194 - t233) * t279 + t321;
t139 = -t182 * t269 + t212 * t247 + t320;
t138 = t183 * t269 - t212 * t248 + t319;
t137 = t248 * t182 - t247 * t183 + t317;
t136 = t175 * t247 + (-t159 - t196) * t269 + t318;
t135 = t160 * t269 + (-t175 - t218) * t248 + t316;
t134 = t248 * t159 + (-t160 - t197) * t247 + t315;
t133 = -t151 * t229 + t164 * t247 + t169 * t198 + (-t143 - t196) * t269 + t318;
t132 = t144 * t269 + t152 * t229 - t169 * t199 + (-t164 - t218) * t248 + t316;
t131 = (-t144 - t197) * t247 + t248 * t143 - t198 * t152 + t199 * t151 + t315;
t1 = m(1) * (t276 ^ 2 + t277 ^ 2 + t278 ^ 2) / 0.2e1 + t229 * ((t146 * t256 + t148 * t230 + t150 * t231) * t199 + (t145 * t256 + t147 * t230 + t149 * t231) * t198 + (t256 * t165 + t230 * t166 + t231 * t167) * t229) / 0.2e1 + m(2) * (t242 ^ 2 + t250 ^ 2 + t251 ^ 2) / 0.2e1 + t199 * ((t240 * t146 + t148 * t202 + t150 * t203) * t199 + (t145 * t240 + t147 * t202 + t149 * t203) * t198 + (t165 * t240 + t166 * t202 + t167 * t203) * t229) / 0.2e1 + t198 * ((t146 * t238 + t148 * t200 + t150 * t201) * t199 + (t145 * t238 + t200 * t147 + t149 * t201) * t198 + (t165 * t238 + t166 * t200 + t167 * t201) * t229) / 0.2e1 + m(3) * (t161 ^ 2 + t162 ^ 2 + t163 ^ 2) / 0.2e1 + m(5) * (t137 ^ 2 + t138 ^ 2 + t139 ^ 2) / 0.2e1 + m(4) * (t140 ^ 2 + t141 ^ 2 + t142 ^ 2) / 0.2e1 + m(6) * (t134 ^ 2 + t135 ^ 2 + t136 ^ 2) / 0.2e1 + m(7) * (t131 ^ 2 + t132 ^ 2 + t133 ^ 2) / 0.2e1 + ((t173 * t204 + t174 * t205 + t209 * t265 + t211 * t239 + t238 * t365) * t269 + (t156 * t204 + t158 * t205 + t177 * t265 + t181 * t239 + t238 * t366) * t248 + (t155 * t204 + t157 * t205 + t265 * t176 + t239 * t180 + t367 * t238) * t247) * t247 / 0.2e1 + ((t173 * t206 + t174 * t207 + t209 * t267 + t211 * t241 + t240 * t365) * t269 + (t156 * t206 + t158 * t207 + t177 * t267 + t181 * t241 + t366 * t240) * t248 + (t155 * t206 + t157 * t207 + t176 * t267 + t180 * t241 + t240 * t367) * t247) * t248 / 0.2e1 + ((t173 * t236 + t174 * t237 - t209 * t343 + t257 * t211 + t365 * t256) * t269 + (t156 * t236 + t158 * t237 - t177 * t343 + t257 * t181 + t256 * t366) * t248 + (t155 * t236 + t157 * t237 - t176 * t343 + t257 * t180 + t256 * t367) * t247) * t269 / 0.2e1 + ((t216 * t243 + t217 * t244 - t253 * t345 + t255 * t266 + t265 * t362) * t291 + (t190 * t243 + t192 * t244 - t220 * t345 + t224 * t266 + t265 * t363) * t280 + (t189 * t243 + t191 * t244 - t219 * t345 + t223 * t266 + t364 * t265) * t279) * t279 / 0.2e1 + ((t216 * t245 + t217 * t246 + t253 * t346 + t255 * t268 + t267 * t362) * t291 + (t190 * t245 + t192 * t246 + t220 * t346 + t224 * t268 + t363 * t267) * t280 + (t189 * t245 + t191 * t246 + t219 * t346 + t223 * t268 + t267 * t364) * t279) * t280 / 0.2e1 + ((t219 * t279 + t220 * t280 + t253 * t291) * t310 + ((t222 * t314 + t224 * t313) * t280 + (t221 * t314 + t223 * t313) * t279 + (t254 * t314 + t255 * t313) * t291) * t306 + (-t188 * t343 + t263 * t190 + t264 * t192) * t280 + (-t187 * t343 + t263 * t189 + t264 * t191) * t279 + (-t215 * t343 + t263 * t216 + t264 * t217) * t291) * t291 / 0.2e1 + ((-t283 * t305 + t285 * t309 + Icges(1,4)) * V_base(5) + (-t284 * t305 + t286 * t309 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t283 * t309 + t285 * t305 + Icges(1,2)) * V_base(5) + (t284 * t309 + t286 * t305 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t309 - Icges(2,6) * t305 + Icges(1,5)) * V_base(4) + (Icges(2,5) * t305 + Icges(2,6) * t309 + Icges(1,6)) * V_base(5) + (Icges(2,3) / 0.2e1 + Icges(1,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
