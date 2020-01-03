% Calculate kinetic energy for
% S5RPRRR14
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
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR14_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR14_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR14_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRRR14_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR14_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR14_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR14_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:16:48
% EndTime: 2019-12-31 19:16:52
% DurationCPUTime: 3.51s
% Computational Cost: add. (3511->370), mult. (9174->543), div. (0->0), fcn. (11662->14), ass. (0->161)
t346 = cos(qJ(3));
t345 = cos(qJ(4));
t344 = cos(pkin(6));
t343 = sin(pkin(6));
t313 = sin(qJ(1));
t342 = Icges(2,4) * t313;
t309 = cos(pkin(5));
t341 = qJ(2) * t309;
t306 = sin(pkin(11));
t307 = sin(pkin(5));
t340 = t306 * t307;
t339 = t307 * t313;
t315 = cos(qJ(1));
t338 = t307 * t315;
t337 = t309 * t315;
t336 = t313 * t306;
t308 = cos(pkin(11));
t335 = t313 * t308;
t334 = qJD(2) * t307;
t333 = V_base(5) * pkin(7) + V_base(1);
t281 = -t306 * t315 - t309 * t335;
t329 = t307 * t344;
t262 = -t281 * t343 + t313 * t329;
t253 = qJD(3) * t262 + V_base(4);
t279 = t308 * t337 - t336;
t261 = -t279 * t343 - t315 * t329;
t252 = qJD(3) * t261 + V_base(5);
t303 = V_base(6) + qJD(1);
t330 = -pkin(7) - t341;
t328 = t307 * t343;
t284 = t313 * pkin(1) - qJ(2) * t338;
t327 = qJD(2) * t309 + V_base(4) * t284 + V_base(3);
t282 = t308 * t315 - t309 * t336;
t312 = sin(qJ(3));
t325 = t346 * t343;
t323 = t307 * t325;
t326 = t344 * t346;
t239 = -t281 * t326 + t282 * t312 - t313 * t323;
t221 = qJD(4) * t239 + t253;
t280 = t306 * t337 + t335;
t237 = -t279 * t326 + t280 * t312 + t315 * t323;
t220 = qJD(4) * t237 + t252;
t278 = -t308 * t328 + t309 * t344;
t268 = qJD(3) * t278 + t303;
t324 = t313 * t334 + V_base(5) * t341 + t333;
t259 = -t307 * t308 * t326 - t309 * t325 + t312 * t340;
t231 = qJD(4) * t259 + t268;
t285 = pkin(1) * t315 + qJ(2) * t339;
t322 = t303 * t285 - t315 * t334 + V_base(2);
t242 = t280 * pkin(2) + pkin(8) * t261;
t265 = pkin(2) * t340 + pkin(8) * t278;
t321 = V_base(5) * t265 + (-t242 - t284) * t303 + t324;
t243 = t282 * pkin(2) + pkin(8) * t262;
t320 = V_base(4) * t242 + (-t243 - t285) * V_base(5) + t327;
t238 = t280 * t346 + (t279 * t344 - t315 * t328) * t312;
t215 = pkin(3) * t238 + pkin(9) * t237;
t260 = t309 * t343 * t312 + (t308 * t312 * t344 + t306 * t346) * t307;
t230 = pkin(3) * t260 + pkin(9) * t259;
t319 = -t215 * t268 + t252 * t230 + t321;
t240 = t282 * t346 + (t281 * t344 + t313 * t328) * t312;
t216 = pkin(3) * t240 + pkin(9) * t239;
t318 = t253 * t215 - t216 * t252 + t320;
t317 = t303 * t243 + (-t265 + t330) * V_base(4) + t322;
t316 = t268 * t216 - t253 * t230 + t317;
t314 = cos(qJ(5));
t311 = sin(qJ(4));
t310 = sin(qJ(5));
t304 = Icges(2,4) * t315;
t298 = rSges(2,1) * t315 - t313 * rSges(2,2);
t297 = t313 * rSges(2,1) + rSges(2,2) * t315;
t296 = Icges(2,1) * t315 - t342;
t295 = Icges(2,1) * t313 + t304;
t294 = -Icges(2,2) * t313 + t304;
t293 = Icges(2,2) * t315 + t342;
t292 = Icges(2,5) * t315 - Icges(2,6) * t313;
t291 = Icges(2,5) * t313 + Icges(2,6) * t315;
t290 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t289 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t288 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t274 = rSges(3,3) * t309 + (rSges(3,1) * t306 + rSges(3,2) * t308) * t307;
t273 = Icges(3,5) * t309 + (Icges(3,1) * t306 + Icges(3,4) * t308) * t307;
t272 = Icges(3,6) * t309 + (Icges(3,4) * t306 + Icges(3,2) * t308) * t307;
t271 = Icges(3,3) * t309 + (Icges(3,5) * t306 + Icges(3,6) * t308) * t307;
t267 = V_base(5) * rSges(2,3) - t297 * t303 + t333;
t266 = t298 * t303 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t264 = t297 * V_base(4) - t298 * V_base(5) + V_base(3);
t251 = rSges(3,1) * t282 + rSges(3,2) * t281 + rSges(3,3) * t339;
t250 = t280 * rSges(3,1) + t279 * rSges(3,2) - rSges(3,3) * t338;
t249 = Icges(3,1) * t282 + Icges(3,4) * t281 + Icges(3,5) * t339;
t248 = Icges(3,1) * t280 + Icges(3,4) * t279 - Icges(3,5) * t338;
t247 = Icges(3,4) * t282 + Icges(3,2) * t281 + Icges(3,6) * t339;
t246 = Icges(3,4) * t280 + Icges(3,2) * t279 - Icges(3,6) * t338;
t245 = Icges(3,5) * t282 + Icges(3,6) * t281 + Icges(3,3) * t339;
t244 = Icges(3,5) * t280 + Icges(3,6) * t279 - Icges(3,3) * t338;
t235 = t260 * t345 + t278 * t311;
t234 = t260 * t311 - t278 * t345;
t229 = rSges(4,1) * t260 - rSges(4,2) * t259 + rSges(4,3) * t278;
t228 = Icges(4,1) * t260 - Icges(4,4) * t259 + Icges(4,5) * t278;
t227 = Icges(4,4) * t260 - Icges(4,2) * t259 + Icges(4,6) * t278;
t226 = Icges(4,5) * t260 - Icges(4,6) * t259 + Icges(4,3) * t278;
t225 = t240 * t345 + t262 * t311;
t224 = t240 * t311 - t262 * t345;
t223 = t238 * t345 + t261 * t311;
t222 = t238 * t311 - t261 * t345;
t219 = t235 * t314 + t259 * t310;
t218 = -t235 * t310 + t259 * t314;
t214 = t274 * V_base(5) + (-t250 - t284) * t303 + t324;
t213 = t303 * t251 + (-t274 + t330) * V_base(4) + t322;
t212 = qJD(5) * t234 + t231;
t211 = pkin(4) * t235 + pkin(10) * t234;
t210 = t250 * V_base(4) + (-t251 - t285) * V_base(5) + t327;
t208 = rSges(4,1) * t240 - rSges(4,2) * t239 + rSges(4,3) * t262;
t207 = rSges(4,1) * t238 - rSges(4,2) * t237 + rSges(4,3) * t261;
t206 = Icges(4,1) * t240 - Icges(4,4) * t239 + Icges(4,5) * t262;
t205 = Icges(4,1) * t238 - Icges(4,4) * t237 + Icges(4,5) * t261;
t204 = Icges(4,4) * t240 - Icges(4,2) * t239 + Icges(4,6) * t262;
t203 = Icges(4,4) * t238 - Icges(4,2) * t237 + Icges(4,6) * t261;
t202 = Icges(4,5) * t240 - Icges(4,6) * t239 + Icges(4,3) * t262;
t201 = Icges(4,5) * t238 - Icges(4,6) * t237 + Icges(4,3) * t261;
t200 = rSges(5,1) * t235 - rSges(5,2) * t234 + rSges(5,3) * t259;
t198 = Icges(5,1) * t235 - Icges(5,4) * t234 + Icges(5,5) * t259;
t197 = Icges(5,4) * t235 - Icges(5,2) * t234 + Icges(5,6) * t259;
t196 = Icges(5,5) * t235 - Icges(5,6) * t234 + Icges(5,3) * t259;
t195 = t225 * t314 + t239 * t310;
t194 = -t225 * t310 + t239 * t314;
t193 = t223 * t314 + t237 * t310;
t192 = -t223 * t310 + t237 * t314;
t191 = qJD(5) * t224 + t221;
t190 = qJD(5) * t222 + t220;
t189 = pkin(4) * t225 + pkin(10) * t224;
t188 = pkin(4) * t223 + pkin(10) * t222;
t187 = rSges(5,1) * t225 - rSges(5,2) * t224 + rSges(5,3) * t239;
t186 = rSges(5,1) * t223 - rSges(5,2) * t222 + rSges(5,3) * t237;
t185 = Icges(5,1) * t225 - Icges(5,4) * t224 + Icges(5,5) * t239;
t184 = Icges(5,1) * t223 - Icges(5,4) * t222 + Icges(5,5) * t237;
t183 = Icges(5,4) * t225 - Icges(5,2) * t224 + Icges(5,6) * t239;
t182 = Icges(5,4) * t223 - Icges(5,2) * t222 + Icges(5,6) * t237;
t181 = Icges(5,5) * t225 - Icges(5,6) * t224 + Icges(5,3) * t239;
t180 = Icges(5,5) * t223 - Icges(5,6) * t222 + Icges(5,3) * t237;
t179 = rSges(6,1) * t219 + rSges(6,2) * t218 + rSges(6,3) * t234;
t178 = Icges(6,1) * t219 + Icges(6,4) * t218 + Icges(6,5) * t234;
t177 = Icges(6,4) * t219 + Icges(6,2) * t218 + Icges(6,6) * t234;
t176 = Icges(6,5) * t219 + Icges(6,6) * t218 + Icges(6,3) * t234;
t175 = rSges(6,1) * t195 + rSges(6,2) * t194 + rSges(6,3) * t224;
t174 = rSges(6,1) * t193 + rSges(6,2) * t192 + rSges(6,3) * t222;
t173 = Icges(6,1) * t195 + Icges(6,4) * t194 + Icges(6,5) * t224;
t172 = Icges(6,1) * t193 + Icges(6,4) * t192 + Icges(6,5) * t222;
t171 = Icges(6,4) * t195 + Icges(6,2) * t194 + Icges(6,6) * t224;
t170 = Icges(6,4) * t193 + Icges(6,2) * t192 + Icges(6,6) * t222;
t169 = Icges(6,5) * t195 + Icges(6,6) * t194 + Icges(6,3) * t224;
t168 = Icges(6,5) * t193 + Icges(6,6) * t192 + Icges(6,3) * t222;
t167 = -t207 * t268 + t229 * t252 + t321;
t166 = t268 * t208 - t253 * t229 + t317;
t165 = t207 * t253 - t208 * t252 + t320;
t164 = -t186 * t231 + t200 * t220 + t319;
t163 = t231 * t187 - t221 * t200 + t316;
t162 = t186 * t221 - t187 * t220 + t318;
t161 = -t174 * t212 + t179 * t190 - t188 * t231 + t211 * t220 + t319;
t160 = t212 * t175 - t191 * t179 + t231 * t189 - t221 * t211 + t316;
t159 = t174 * t191 - t175 * t190 + t188 * t221 - t189 * t220 + t318;
t1 = m(1) * (t288 ^ 2 + t289 ^ 2 + t290 ^ 2) / 0.2e1 + m(2) * (t264 ^ 2 + t266 ^ 2 + t267 ^ 2) / 0.2e1 + m(3) * (t210 ^ 2 + t213 ^ 2 + t214 ^ 2) / 0.2e1 + m(4) * (t165 ^ 2 + t166 ^ 2 + t167 ^ 2) / 0.2e1 + t253 * ((t262 * t202 - t239 * t204 + t240 * t206) * t253 + (t201 * t262 - t203 * t239 + t205 * t240) * t252 + (t226 * t262 - t227 * t239 + t228 * t240) * t268) / 0.2e1 + t252 * ((t202 * t261 - t204 * t237 + t206 * t238) * t253 + (t261 * t201 - t237 * t203 + t238 * t205) * t252 + (t226 * t261 - t227 * t237 + t228 * t238) * t268) / 0.2e1 + t268 * ((t202 * t278 - t204 * t259 + t206 * t260) * t253 + (t201 * t278 - t203 * t259 + t205 * t260) * t252 + (t278 * t226 - t259 * t227 + t260 * t228) * t268) / 0.2e1 + m(5) * (t162 ^ 2 + t163 ^ 2 + t164 ^ 2) / 0.2e1 + t221 * ((t239 * t181 - t224 * t183 + t225 * t185) * t221 + (t180 * t239 - t182 * t224 + t184 * t225) * t220 + (t196 * t239 - t197 * t224 + t198 * t225) * t231) / 0.2e1 + t220 * ((t181 * t237 - t183 * t222 + t185 * t223) * t221 + (t237 * t180 - t222 * t182 + t223 * t184) * t220 + (t196 * t237 - t197 * t222 + t198 * t223) * t231) / 0.2e1 + t231 * ((t181 * t259 - t183 * t234 + t185 * t235) * t221 + (t180 * t259 - t182 * t234 + t184 * t235) * t220 + (t259 * t196 - t234 * t197 + t235 * t198) * t231) / 0.2e1 + m(6) * (t159 ^ 2 + t160 ^ 2 + t161 ^ 2) / 0.2e1 + t191 * ((t224 * t169 + t194 * t171 + t195 * t173) * t191 + (t168 * t224 + t170 * t194 + t172 * t195) * t190 + (t176 * t224 + t177 * t194 + t178 * t195) * t212) / 0.2e1 + t190 * ((t169 * t222 + t171 * t192 + t173 * t193) * t191 + (t222 * t168 + t192 * t170 + t193 * t172) * t190 + (t176 * t222 + t177 * t192 + t178 * t193) * t212) / 0.2e1 + t212 * ((t169 * t234 + t171 * t218 + t173 * t219) * t191 + (t168 * t234 + t170 * t218 + t172 * t219) * t190 + (t234 * t176 + t218 * t177 + t219 * t178) * t212) / 0.2e1 + (Icges(2,3) * t303 + t291 * V_base(5) + t292 * V_base(4) + (t244 * V_base(5) + t245 * V_base(4) + t271 * t303) * t309 + ((t247 * t308 + t249 * t306) * V_base(4) + (t246 * t308 + t248 * t306) * V_base(5) + (t272 * t308 + t273 * t306) * t303) * t307) * t303 / 0.2e1 + ((t271 * t339 + t272 * t281 + t273 * t282 + t292) * t303 + (t244 * t339 + t246 * t281 + t248 * t282 - t313 * t293 + t295 * t315 + Icges(1,4)) * V_base(5) + (t245 * t339 + t247 * t281 + t249 * t282 - t313 * t294 + t296 * t315 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((-t271 * t338 + t279 * t272 + t280 * t273 + t291) * t303 + (-t244 * t338 + t279 * t246 + t280 * t248 + t293 * t315 + t313 * t295 + Icges(1,2)) * V_base(5) + (-t245 * t338 + t279 * t247 + t280 * t249 + t294 * t315 + t313 * t296 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
