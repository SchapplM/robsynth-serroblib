% Calculate kinetic energy for
% S5RRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d2,d3,d4,d5]';
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR12_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR12_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR12_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRRR12_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR12_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR12_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR12_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR12_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:46:30
% EndTime: 2019-12-31 22:46:33
% DurationCPUTime: 3.32s
% Computational Cost: add. (3607->365), mult. (9334->555), div. (0->0), fcn. (11822->14), ass. (0->162)
t349 = cos(qJ(3));
t348 = cos(qJ(4));
t310 = cos(pkin(5));
t347 = pkin(8) * t310;
t346 = cos(pkin(6));
t345 = sin(pkin(6));
t315 = sin(qJ(1));
t344 = Icges(2,4) * t315;
t309 = sin(pkin(5));
t314 = sin(qJ(2));
t343 = t309 * t314;
t342 = t309 * t315;
t318 = cos(qJ(1));
t341 = t309 * t318;
t340 = t314 * t318;
t339 = t315 * t314;
t317 = cos(qJ(2));
t338 = t315 * t317;
t337 = t317 * t318;
t336 = qJD(2) * t309;
t335 = V_base(5) * pkin(7) + V_base(1);
t292 = t315 * t336 + V_base(4);
t306 = V_base(6) + qJD(1);
t332 = t309 * t346;
t331 = t309 * t345;
t281 = -t310 * t338 - t340;
t263 = -t281 * t345 + t315 * t332;
t253 = qJD(3) * t263 + t292;
t293 = qJD(2) * t310 + t306;
t330 = t346 * t349;
t329 = t349 * t345;
t282 = -t310 * t339 + t337;
t313 = sin(qJ(3));
t328 = t309 * t329;
t240 = -t281 * t330 + t282 * t313 - t315 * t328;
t219 = qJD(4) * t240 + t253;
t278 = t310 * t346 - t317 * t331;
t264 = qJD(3) * t278 + t293;
t291 = -t318 * t336 + V_base(5);
t284 = pkin(1) * t315 - pkin(8) * t341;
t327 = -t284 * t306 + t347 * V_base(5) + t335;
t260 = -t309 * t317 * t330 - t310 * t329 + t313 * t343;
t231 = qJD(4) * t260 + t264;
t285 = pkin(1) * t318 + pkin(8) * t342;
t326 = t284 * V_base(4) - t285 * V_base(5) + V_base(3);
t279 = t310 * t337 - t339;
t262 = -t279 * t345 - t318 * t332;
t252 = qJD(3) * t262 + t291;
t280 = t310 * t340 + t338;
t238 = -t279 * t330 + t280 * t313 + t318 * t328;
t218 = qJD(4) * t238 + t252;
t325 = t306 * t285 + V_base(2) + (-pkin(7) - t347) * V_base(4);
t242 = t280 * pkin(2) + pkin(9) * t262;
t268 = pkin(2) * t343 + pkin(9) * t278;
t324 = -t242 * t293 + t268 * t291 + t327;
t243 = pkin(2) * t282 + pkin(9) * t263;
t323 = t242 * t292 - t243 * t291 + t326;
t322 = t243 * t293 - t268 * t292 + t325;
t239 = t280 * t349 + (t279 * t346 - t318 * t331) * t313;
t215 = pkin(3) * t239 + pkin(10) * t238;
t261 = t310 * t345 * t313 + (t313 * t317 * t346 + t314 * t349) * t309;
t230 = pkin(3) * t261 + pkin(10) * t260;
t321 = -t215 * t264 + t230 * t252 + t324;
t241 = t282 * t349 + (t281 * t346 + t315 * t331) * t313;
t216 = pkin(3) * t241 + pkin(10) * t240;
t320 = t215 * t253 - t216 * t252 + t323;
t319 = t216 * t264 - t230 * t253 + t322;
t316 = cos(qJ(5));
t312 = sin(qJ(4));
t311 = sin(qJ(5));
t307 = Icges(2,4) * t318;
t301 = rSges(2,1) * t318 - rSges(2,2) * t315;
t300 = rSges(2,1) * t315 + rSges(2,2) * t318;
t299 = Icges(2,1) * t318 - t344;
t298 = Icges(2,1) * t315 + t307;
t297 = -Icges(2,2) * t315 + t307;
t296 = Icges(2,2) * t318 + t344;
t290 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t289 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t288 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t274 = rSges(3,3) * t310 + (rSges(3,1) * t314 + rSges(3,2) * t317) * t309;
t273 = Icges(3,5) * t310 + (Icges(3,1) * t314 + Icges(3,4) * t317) * t309;
t272 = Icges(3,6) * t310 + (Icges(3,4) * t314 + Icges(3,2) * t317) * t309;
t271 = Icges(3,3) * t310 + (Icges(3,5) * t314 + Icges(3,6) * t317) * t309;
t267 = V_base(5) * rSges(2,3) - t300 * t306 + t335;
t266 = t301 * t306 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t265 = t300 * V_base(4) - t301 * V_base(5) + V_base(3);
t251 = rSges(3,1) * t282 + rSges(3,2) * t281 + rSges(3,3) * t342;
t250 = rSges(3,1) * t280 + rSges(3,2) * t279 - rSges(3,3) * t341;
t249 = Icges(3,1) * t282 + Icges(3,4) * t281 + Icges(3,5) * t342;
t248 = Icges(3,1) * t280 + Icges(3,4) * t279 - Icges(3,5) * t341;
t247 = Icges(3,4) * t282 + Icges(3,2) * t281 + Icges(3,6) * t342;
t246 = Icges(3,4) * t280 + Icges(3,2) * t279 - Icges(3,6) * t341;
t245 = Icges(3,5) * t282 + Icges(3,6) * t281 + Icges(3,3) * t342;
t244 = Icges(3,5) * t280 + Icges(3,6) * t279 - Icges(3,3) * t341;
t237 = t261 * t348 + t278 * t312;
t236 = t261 * t312 - t278 * t348;
t229 = rSges(4,1) * t261 - rSges(4,2) * t260 + rSges(4,3) * t278;
t228 = Icges(4,1) * t261 - Icges(4,4) * t260 + Icges(4,5) * t278;
t227 = Icges(4,4) * t261 - Icges(4,2) * t260 + Icges(4,6) * t278;
t226 = Icges(4,5) * t261 - Icges(4,6) * t260 + Icges(4,3) * t278;
t225 = t241 * t348 + t263 * t312;
t224 = t241 * t312 - t263 * t348;
t223 = t239 * t348 + t262 * t312;
t222 = t239 * t312 - t262 * t348;
t221 = t237 * t316 + t260 * t311;
t220 = -t237 * t311 + t260 * t316;
t214 = pkin(4) * t237 + pkin(11) * t236;
t213 = -t250 * t293 + t274 * t291 + t327;
t212 = t251 * t293 - t274 * t292 + t325;
t211 = qJD(5) * t236 + t231;
t209 = rSges(4,1) * t241 - rSges(4,2) * t240 + rSges(4,3) * t263;
t208 = rSges(4,1) * t239 - rSges(4,2) * t238 + rSges(4,3) * t262;
t207 = t250 * t292 - t251 * t291 + t326;
t206 = Icges(4,1) * t241 - Icges(4,4) * t240 + Icges(4,5) * t263;
t205 = Icges(4,1) * t239 - Icges(4,4) * t238 + Icges(4,5) * t262;
t204 = Icges(4,4) * t241 - Icges(4,2) * t240 + Icges(4,6) * t263;
t203 = Icges(4,4) * t239 - Icges(4,2) * t238 + Icges(4,6) * t262;
t202 = Icges(4,5) * t241 - Icges(4,6) * t240 + Icges(4,3) * t263;
t201 = Icges(4,5) * t239 - Icges(4,6) * t238 + Icges(4,3) * t262;
t200 = rSges(5,1) * t237 - rSges(5,2) * t236 + rSges(5,3) * t260;
t199 = Icges(5,1) * t237 - Icges(5,4) * t236 + Icges(5,5) * t260;
t198 = Icges(5,4) * t237 - Icges(5,2) * t236 + Icges(5,6) * t260;
t197 = Icges(5,5) * t237 - Icges(5,6) * t236 + Icges(5,3) * t260;
t196 = t225 * t316 + t240 * t311;
t195 = -t225 * t311 + t240 * t316;
t194 = t223 * t316 + t238 * t311;
t193 = -t223 * t311 + t238 * t316;
t191 = pkin(4) * t225 + pkin(11) * t224;
t190 = pkin(4) * t223 + pkin(11) * t222;
t189 = qJD(5) * t224 + t219;
t188 = qJD(5) * t222 + t218;
t187 = rSges(5,1) * t225 - rSges(5,2) * t224 + rSges(5,3) * t240;
t186 = rSges(5,1) * t223 - rSges(5,2) * t222 + rSges(5,3) * t238;
t185 = Icges(5,1) * t225 - Icges(5,4) * t224 + Icges(5,5) * t240;
t184 = Icges(5,1) * t223 - Icges(5,4) * t222 + Icges(5,5) * t238;
t183 = Icges(5,4) * t225 - Icges(5,2) * t224 + Icges(5,6) * t240;
t182 = Icges(5,4) * t223 - Icges(5,2) * t222 + Icges(5,6) * t238;
t181 = Icges(5,5) * t225 - Icges(5,6) * t224 + Icges(5,3) * t240;
t180 = Icges(5,5) * t223 - Icges(5,6) * t222 + Icges(5,3) * t238;
t179 = rSges(6,1) * t221 + rSges(6,2) * t220 + rSges(6,3) * t236;
t178 = Icges(6,1) * t221 + Icges(6,4) * t220 + Icges(6,5) * t236;
t177 = Icges(6,4) * t221 + Icges(6,2) * t220 + Icges(6,6) * t236;
t176 = Icges(6,5) * t221 + Icges(6,6) * t220 + Icges(6,3) * t236;
t175 = rSges(6,1) * t196 + rSges(6,2) * t195 + rSges(6,3) * t224;
t174 = rSges(6,1) * t194 + rSges(6,2) * t193 + rSges(6,3) * t222;
t173 = Icges(6,1) * t196 + Icges(6,4) * t195 + Icges(6,5) * t224;
t172 = Icges(6,1) * t194 + Icges(6,4) * t193 + Icges(6,5) * t222;
t171 = Icges(6,4) * t196 + Icges(6,2) * t195 + Icges(6,6) * t224;
t170 = Icges(6,4) * t194 + Icges(6,2) * t193 + Icges(6,6) * t222;
t169 = Icges(6,5) * t196 + Icges(6,6) * t195 + Icges(6,3) * t224;
t168 = Icges(6,5) * t194 + Icges(6,6) * t193 + Icges(6,3) * t222;
t167 = -t208 * t264 + t229 * t252 + t324;
t166 = t209 * t264 - t229 * t253 + t322;
t165 = t208 * t253 - t209 * t252 + t323;
t164 = -t186 * t231 + t200 * t218 + t321;
t163 = t187 * t231 - t200 * t219 + t319;
t162 = t186 * t219 - t187 * t218 + t320;
t161 = -t174 * t211 + t179 * t188 - t190 * t231 + t214 * t218 + t321;
t160 = t175 * t211 - t179 * t189 + t191 * t231 - t214 * t219 + t319;
t159 = t174 * t189 - t175 * t188 + t190 * t219 - t191 * t218 + t320;
t1 = m(1) * (t288 ^ 2 + t289 ^ 2 + t290 ^ 2) / 0.2e1 + m(2) * (t265 ^ 2 + t266 ^ 2 + t267 ^ 2) / 0.2e1 + m(3) * (t207 ^ 2 + t212 ^ 2 + t213 ^ 2) / 0.2e1 + t292 * ((t245 * t342 + t247 * t281 + t249 * t282) * t292 + (t244 * t342 + t246 * t281 + t248 * t282) * t291 + (t271 * t342 + t272 * t281 + t273 * t282) * t293) / 0.2e1 + t291 * ((-t245 * t341 + t247 * t279 + t249 * t280) * t292 + (-t244 * t341 + t279 * t246 + t280 * t248) * t291 + (-t271 * t341 + t272 * t279 + t273 * t280) * t293) / 0.2e1 + t293 * ((t244 * t291 + t245 * t292 + t271 * t293) * t310 + ((t247 * t317 + t249 * t314) * t292 + (t246 * t317 + t248 * t314) * t291 + (t272 * t317 + t273 * t314) * t293) * t309) / 0.2e1 + m(4) * (t165 ^ 2 + t166 ^ 2 + t167 ^ 2) / 0.2e1 + t253 * ((t263 * t202 - t240 * t204 + t241 * t206) * t253 + (t201 * t263 - t203 * t240 + t205 * t241) * t252 + (t226 * t263 - t227 * t240 + t228 * t241) * t264) / 0.2e1 + t252 * ((t202 * t262 - t204 * t238 + t206 * t239) * t253 + (t262 * t201 - t238 * t203 + t239 * t205) * t252 + (t226 * t262 - t227 * t238 + t228 * t239) * t264) / 0.2e1 + t264 * ((t202 * t278 - t204 * t260 + t206 * t261) * t253 + (t201 * t278 - t203 * t260 + t205 * t261) * t252 + (t226 * t278 - t227 * t260 + t228 * t261) * t264) / 0.2e1 + m(5) * (t162 ^ 2 + t163 ^ 2 + t164 ^ 2) / 0.2e1 + t219 * ((t240 * t181 - t224 * t183 + t225 * t185) * t219 + (t180 * t240 - t182 * t224 + t184 * t225) * t218 + (t197 * t240 - t198 * t224 + t199 * t225) * t231) / 0.2e1 + t218 * ((t181 * t238 - t183 * t222 + t185 * t223) * t219 + (t238 * t180 - t222 * t182 + t223 * t184) * t218 + (t197 * t238 - t198 * t222 + t199 * t223) * t231) / 0.2e1 + t231 * ((t181 * t260 - t183 * t236 + t185 * t237) * t219 + (t180 * t260 - t182 * t236 + t184 * t237) * t218 + (t260 * t197 - t236 * t198 + t237 * t199) * t231) / 0.2e1 + m(6) * (t159 ^ 2 + t160 ^ 2 + t161 ^ 2) / 0.2e1 + t189 * ((t224 * t169 + t195 * t171 + t196 * t173) * t189 + (t168 * t224 + t170 * t195 + t172 * t196) * t188 + (t176 * t224 + t177 * t195 + t178 * t196) * t211) / 0.2e1 + t188 * ((t169 * t222 + t171 * t193 + t173 * t194) * t189 + (t222 * t168 + t193 * t170 + t194 * t172) * t188 + (t176 * t222 + t177 * t193 + t178 * t194) * t211) / 0.2e1 + t211 * ((t169 * t236 + t171 * t220 + t173 * t221) * t189 + (t168 * t236 + t170 * t220 + t172 * t221) * t188 + (t236 * t176 + t220 * t177 + t221 * t178) * t211) / 0.2e1 + ((-t296 * t315 + t298 * t318 + Icges(1,4)) * V_base(5) + (-t315 * t297 + t299 * t318 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t296 * t318 + t315 * t298 + Icges(1,2)) * V_base(5) + (t297 * t318 + t299 * t315 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t315 + Icges(2,6) * t318) * V_base(5) + (Icges(2,5) * t318 - Icges(2,6) * t315) * V_base(4) + Icges(2,3) * t306 / 0.2e1) * t306;
T = t1;
