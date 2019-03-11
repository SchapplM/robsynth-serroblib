% Calculate kinetic energy for
% S6RRPPRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2019-03-09 09:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR5_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPPRR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR5_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR5_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRR5_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:07:08
% EndTime: 2019-03-09 09:07:11
% DurationCPUTime: 3.81s
% Computational Cost: add. (1876->344), mult. (4163->487), div. (0->0), fcn. (4711->10), ass. (0->156)
t347 = Icges(3,1) + Icges(4,1) + Icges(5,1);
t346 = Icges(3,4) - Icges(5,4) - Icges(4,5);
t345 = -Icges(5,5) + Icges(4,4) + Icges(3,5);
t344 = Icges(3,2) + Icges(5,2) + Icges(4,3);
t343 = Icges(4,6) - Icges(5,6) - Icges(3,6);
t342 = -Icges(5,3) - Icges(3,3) - Icges(4,2);
t290 = cos(pkin(6));
t296 = cos(qJ(2));
t297 = cos(qJ(1));
t319 = t296 * t297;
t293 = sin(qJ(2));
t294 = sin(qJ(1));
t322 = t293 * t294;
t255 = -t290 * t319 + t322;
t320 = t294 * t296;
t321 = t293 * t297;
t256 = t290 * t321 + t320;
t289 = sin(pkin(6));
t323 = t289 * t297;
t341 = -t255 * t343 - t256 * t345 - t323 * t342;
t257 = t290 * t320 + t321;
t258 = -t290 * t322 + t319;
t325 = t289 * t294;
t340 = -t257 * t343 - t258 * t345 + t325 * t342;
t339 = t255 * t344 - t256 * t346 - t323 * t343;
t338 = t257 * t344 - t258 * t346 + t325 * t343;
t337 = -t346 * t255 + t256 * t347 - t345 * t323;
t336 = -t346 * t257 + t258 * t347 + t345 * t325;
t335 = t342 * t290 + (-t293 * t345 + t296 * t343) * t289;
t334 = t343 * t290 + (-t293 * t346 - t296 * t344) * t289;
t333 = t345 * t290 + (t293 * t347 + t346 * t296) * t289;
t329 = cos(qJ(5));
t328 = pkin(8) * t290;
t327 = Icges(2,4) * t294;
t326 = t289 * t293;
t324 = t289 * t296;
t215 = pkin(2) * t256 + qJ(3) * t255;
t232 = t256 * pkin(3) + qJ(4) * t323;
t318 = -t215 - t232;
t216 = pkin(2) * t258 + qJ(3) * t257;
t233 = pkin(3) * t258 - qJ(4) * t325;
t317 = -t216 - t233;
t259 = (pkin(2) * t293 - qJ(3) * t296) * t289;
t262 = pkin(3) * t326 - qJ(4) * t290;
t316 = -t259 - t262;
t315 = qJD(2) * t289;
t314 = qJD(4) * t289;
t313 = V_base(5) * pkin(7) + V_base(1);
t310 = t289 * t329;
t269 = t294 * t315 + V_base(4);
t286 = V_base(6) + qJD(1);
t222 = -qJD(5) * t257 + t269;
t270 = qJD(2) * t290 + t286;
t251 = qJD(5) * t324 + t270;
t268 = -t297 * t315 + V_base(5);
t263 = t294 * pkin(1) - pkin(8) * t323;
t309 = -t263 * t286 + V_base(5) * t328 + t313;
t264 = pkin(1) * t297 + pkin(8) * t325;
t308 = V_base(4) * t263 - t264 * V_base(5) + V_base(3);
t221 = -qJD(5) * t255 + t268;
t307 = qJD(3) * t257 + t268 * t259 + t309;
t306 = t286 * t264 + V_base(2) + (-pkin(7) - t328) * V_base(4);
t305 = qJD(3) * t255 + t270 * t216 + t306;
t304 = -qJD(3) * t324 + t269 * t215 + t308;
t303 = t270 * t233 + t297 * t314 + t305;
t302 = t268 * t262 - t294 * t314 + t307;
t301 = -qJD(4) * t290 + t269 * t232 + t304;
t218 = pkin(4) * t258 - pkin(9) * t257;
t261 = (pkin(4) * t293 + pkin(9) * t296) * t289;
t300 = t270 * t218 + (-t261 + t316) * t269 + t303;
t217 = pkin(4) * t256 - pkin(9) * t255;
t299 = t268 * t261 + (-t217 + t318) * t270 + t302;
t298 = t269 * t217 + (-t218 + t317) * t268 + t301;
t295 = cos(qJ(6));
t292 = sin(qJ(5));
t291 = sin(qJ(6));
t287 = Icges(2,4) * t297;
t278 = rSges(2,1) * t297 - t294 * rSges(2,2);
t277 = t294 * rSges(2,1) + rSges(2,2) * t297;
t276 = Icges(2,1) * t297 - t327;
t275 = Icges(2,1) * t294 + t287;
t274 = -Icges(2,2) * t294 + t287;
t273 = Icges(2,2) * t297 + t327;
t267 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t266 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t265 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t254 = -t290 * t292 + t293 * t310;
t253 = t290 * t329 + t292 * t326;
t246 = rSges(3,3) * t290 + (rSges(3,1) * t293 + rSges(3,2) * t296) * t289;
t245 = rSges(4,2) * t290 + (rSges(4,1) * t293 - rSges(4,3) * t296) * t289;
t244 = -rSges(5,3) * t290 + (rSges(5,1) * t293 - rSges(5,2) * t296) * t289;
t231 = V_base(5) * rSges(2,3) - t277 * t286 + t313;
t230 = t278 * t286 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t227 = t277 * V_base(4) - t278 * V_base(5) + V_base(3);
t226 = t258 * t329 - t292 * t325;
t225 = t258 * t292 + t294 * t310;
t224 = t256 * t329 + t292 * t323;
t223 = t256 * t292 - t297 * t310;
t220 = t254 * t295 + t291 * t324;
t219 = -t254 * t291 + t295 * t324;
t214 = qJD(6) * t253 + t251;
t213 = pkin(5) * t254 + pkin(10) * t253;
t210 = rSges(3,1) * t258 - rSges(3,2) * t257 + rSges(3,3) * t325;
t209 = rSges(4,1) * t258 + rSges(4,2) * t325 + rSges(4,3) * t257;
t208 = rSges(5,1) * t258 + rSges(5,2) * t257 - rSges(5,3) * t325;
t207 = t256 * rSges(3,1) - t255 * rSges(3,2) - rSges(3,3) * t323;
t206 = t256 * rSges(4,1) - rSges(4,2) * t323 + t255 * rSges(4,3);
t205 = t256 * rSges(5,1) + t255 * rSges(5,2) + rSges(5,3) * t323;
t186 = rSges(6,1) * t254 - rSges(6,2) * t253 + rSges(6,3) * t324;
t185 = Icges(6,1) * t254 - Icges(6,4) * t253 + Icges(6,5) * t324;
t184 = Icges(6,4) * t254 - Icges(6,2) * t253 + Icges(6,6) * t324;
t183 = Icges(6,5) * t254 - Icges(6,6) * t253 + Icges(6,3) * t324;
t178 = t226 * t295 - t257 * t291;
t177 = -t226 * t291 - t257 * t295;
t176 = t224 * t295 - t255 * t291;
t175 = -t224 * t291 - t255 * t295;
t174 = qJD(6) * t225 + t222;
t173 = qJD(6) * t223 + t221;
t172 = pkin(5) * t226 + pkin(10) * t225;
t171 = pkin(5) * t224 + pkin(10) * t223;
t170 = rSges(6,1) * t226 - rSges(6,2) * t225 - rSges(6,3) * t257;
t169 = rSges(6,1) * t224 - rSges(6,2) * t223 - rSges(6,3) * t255;
t168 = Icges(6,1) * t226 - Icges(6,4) * t225 - Icges(6,5) * t257;
t167 = Icges(6,1) * t224 - Icges(6,4) * t223 - Icges(6,5) * t255;
t166 = Icges(6,4) * t226 - Icges(6,2) * t225 - Icges(6,6) * t257;
t165 = Icges(6,4) * t224 - Icges(6,2) * t223 - Icges(6,6) * t255;
t164 = Icges(6,5) * t226 - Icges(6,6) * t225 - Icges(6,3) * t257;
t163 = Icges(6,5) * t224 - Icges(6,6) * t223 - Icges(6,3) * t255;
t162 = rSges(7,1) * t220 + rSges(7,2) * t219 + rSges(7,3) * t253;
t161 = Icges(7,1) * t220 + Icges(7,4) * t219 + Icges(7,5) * t253;
t160 = Icges(7,4) * t220 + Icges(7,2) * t219 + Icges(7,6) * t253;
t159 = Icges(7,5) * t220 + Icges(7,6) * t219 + Icges(7,3) * t253;
t158 = -t207 * t270 + t246 * t268 + t309;
t157 = t210 * t270 - t246 * t269 + t306;
t156 = rSges(7,1) * t178 + rSges(7,2) * t177 + rSges(7,3) * t225;
t155 = rSges(7,1) * t176 + rSges(7,2) * t175 + rSges(7,3) * t223;
t154 = Icges(7,1) * t178 + Icges(7,4) * t177 + Icges(7,5) * t225;
t153 = Icges(7,1) * t176 + Icges(7,4) * t175 + Icges(7,5) * t223;
t152 = Icges(7,4) * t178 + Icges(7,2) * t177 + Icges(7,6) * t225;
t151 = Icges(7,4) * t176 + Icges(7,2) * t175 + Icges(7,6) * t223;
t150 = Icges(7,5) * t178 + Icges(7,6) * t177 + Icges(7,3) * t225;
t149 = Icges(7,5) * t176 + Icges(7,6) * t175 + Icges(7,3) * t223;
t148 = t207 * t269 - t210 * t268 + t308;
t147 = t245 * t268 + (-t206 - t215) * t270 + t307;
t146 = t209 * t270 + (-t245 - t259) * t269 + t305;
t145 = t206 * t269 + (-t209 - t216) * t268 + t304;
t144 = t244 * t268 + (-t205 + t318) * t270 + t302;
t143 = t208 * t270 + (-t244 + t316) * t269 + t303;
t142 = t205 * t269 + (-t208 + t317) * t268 + t301;
t141 = -t169 * t251 + t186 * t221 + t299;
t140 = t170 * t251 - t186 * t222 + t300;
t139 = t169 * t222 - t170 * t221 + t298;
t138 = -t155 * t214 + t162 * t173 - t171 * t251 + t213 * t221 + t299;
t137 = t156 * t214 - t162 * t174 + t172 * t251 - t213 * t222 + t300;
t136 = t155 * t174 - t156 * t173 + t171 * t222 - t172 * t221 + t298;
t1 = m(7) * (t136 ^ 2 + t137 ^ 2 + t138 ^ 2) / 0.2e1 + m(1) * (t265 ^ 2 + t266 ^ 2 + t267 ^ 2) / 0.2e1 + t222 * ((-t257 * t164 - t225 * t166 + t226 * t168) * t222 + (-t163 * t257 - t165 * t225 + t167 * t226) * t221 + (-t183 * t257 - t184 * t225 + t185 * t226) * t251) / 0.2e1 + t251 * ((t164 * t324 - t166 * t253 + t168 * t254) * t222 + (t163 * t324 - t165 * t253 + t167 * t254) * t221 + (t183 * t324 - t253 * t184 + t254 * t185) * t251) / 0.2e1 + t221 * ((-t164 * t255 - t166 * t223 + t168 * t224) * t222 + (-t255 * t163 - t223 * t165 + t224 * t167) * t221 + (-t183 * t255 - t184 * t223 + t185 * t224) * t251) / 0.2e1 + t214 * ((t150 * t253 + t152 * t219 + t154 * t220) * t174 + (t149 * t253 + t151 * t219 + t153 * t220) * t173 + (t253 * t159 + t219 * t160 + t220 * t161) * t214) / 0.2e1 + m(2) * (t227 ^ 2 + t230 ^ 2 + t231 ^ 2) / 0.2e1 + t173 * ((t150 * t223 + t152 * t175 + t154 * t176) * t174 + (t223 * t149 + t175 * t151 + t176 * t153) * t173 + (t159 * t223 + t160 * t175 + t161 * t176) * t214) / 0.2e1 + t174 * ((t225 * t150 + t177 * t152 + t178 * t154) * t174 + (t149 * t225 + t151 * t177 + t153 * t178) * t173 + (t159 * t225 + t160 * t177 + t161 * t178) * t214) / 0.2e1 + m(6) * (t139 ^ 2 + t140 ^ 2 + t141 ^ 2) / 0.2e1 + m(5) * (t142 ^ 2 + t143 ^ 2 + t144 ^ 2) / 0.2e1 + m(4) * (t145 ^ 2 + t146 ^ 2 + t147 ^ 2) / 0.2e1 + m(3) * (t148 ^ 2 + t157 ^ 2 + t158 ^ 2) / 0.2e1 + ((-t294 * t273 + t275 * t297 + Icges(1,4)) * V_base(5) + (-t294 * t274 + t297 * t276 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t297 * t273 + t294 * t275 + Icges(1,2)) * V_base(5) + (t274 * t297 + t294 * t276 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t255 * t334 + t256 * t333 + t323 * t335) * t270 + (t255 * t338 + t256 * t336 + t323 * t340) * t269 + (t339 * t255 + t337 * t256 + t341 * t323) * t268) * t268 / 0.2e1 + ((t257 * t334 + t258 * t333 - t325 * t335) * t270 + (t338 * t257 + t336 * t258 - t340 * t325) * t269 + (t339 * t257 + t337 * t258 - t325 * t341) * t268) * t269 / 0.2e1 + ((-t268 * t341 - t340 * t269 - t335 * t270) * t290 + ((t293 * t333 - t296 * t334) * t270 + (t293 * t336 - t296 * t338) * t269 + (t293 * t337 - t296 * t339) * t268) * t289) * t270 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t294 + Icges(2,6) * t297) * V_base(5) + (Icges(2,5) * t297 - Icges(2,6) * t294) * V_base(4) + Icges(2,3) * t286 / 0.2e1) * t286;
T  = t1;
