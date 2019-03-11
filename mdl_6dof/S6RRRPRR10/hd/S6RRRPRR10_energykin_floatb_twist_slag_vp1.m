% Calculate kinetic energy for
% S6RRRPRR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRR10_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR10_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR10_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPRR10_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR10_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR10_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR10_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR10_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:15:32
% EndTime: 2019-03-09 19:15:36
% DurationCPUTime: 4.10s
% Computational Cost: add. (1813->372), mult. (3400->547), div. (0->0), fcn. (3680->10), ass. (0->171)
t331 = Icges(4,1) + Icges(5,1);
t330 = -Icges(4,4) + Icges(5,5);
t329 = Icges(5,4) + Icges(4,5);
t328 = Icges(4,2) + Icges(5,3);
t327 = -Icges(5,6) + Icges(4,6);
t326 = -Icges(4,3) - Icges(5,2);
t268 = sin(qJ(3));
t272 = cos(qJ(3));
t274 = cos(qJ(1));
t270 = sin(qJ(1));
t273 = cos(qJ(2));
t302 = t270 * t273;
t222 = t268 * t302 + t272 * t274;
t223 = -t268 * t274 + t272 * t302;
t269 = sin(qJ(2));
t304 = t269 * t270;
t325 = t222 * t328 + t223 * t330 - t304 * t327;
t301 = t273 * t274;
t224 = t268 * t301 - t270 * t272;
t225 = t268 * t270 + t272 * t301;
t303 = t269 * t274;
t324 = t224 * t328 + t225 * t330 - t303 * t327;
t323 = -t222 * t327 + t223 * t329 - t304 * t326;
t322 = -t224 * t327 + t225 * t329 - t303 * t326;
t321 = t330 * t222 + t223 * t331 + t329 * t304;
t320 = t330 * t224 + t225 * t331 + t329 * t303;
t319 = t327 * t273 + (t268 * t328 + t272 * t330) * t269;
t318 = t326 * t273 + (-t268 * t327 + t272 * t329) * t269;
t317 = -t329 * t273 + (t330 * t268 + t272 * t331) * t269;
t271 = cos(qJ(5));
t312 = pkin(5) * t271;
t310 = Icges(2,4) * t270;
t309 = Icges(3,4) * t269;
t308 = Icges(3,4) * t273;
t267 = sin(qJ(5));
t307 = t222 * t267;
t306 = t224 * t267;
t305 = t267 * t268;
t300 = qJD(3) * t269;
t299 = qJD(5) * t269;
t298 = V_base(5) * pkin(6) + V_base(1);
t253 = qJD(2) * t270 + V_base(4);
t295 = t269 * pkin(10);
t259 = V_base(6) + qJD(1);
t294 = t269 * (-qJD(5) - qJD(6));
t221 = t274 * t300 + t253;
t293 = pkin(2) * t273 + pkin(8) * t269;
t252 = -qJD(2) * t274 + V_base(5);
t292 = rSges(3,1) * t273 - rSges(3,2) * t269;
t291 = Icges(3,1) * t273 - t309;
t290 = -Icges(3,2) * t269 + t308;
t289 = Icges(3,5) * t273 - Icges(3,6) * t269;
t220 = t270 * t300 + t252;
t251 = pkin(1) * t274 + pkin(7) * t270;
t288 = -V_base(4) * pkin(6) + t259 * t251 + V_base(2);
t245 = -qJD(3) * t273 + t259;
t250 = pkin(1) * t270 - pkin(7) * t274;
t287 = V_base(4) * t250 - t251 * V_base(5) + V_base(3);
t228 = t293 * t270;
t249 = pkin(2) * t269 - pkin(8) * t273;
t286 = t252 * t249 + (-t228 - t250) * t259 + t298;
t285 = (-Icges(3,3) * t274 + t270 * t289) * t252 + (Icges(3,3) * t270 + t274 * t289) * t253 + (Icges(3,5) * t269 + Icges(3,6) * t273) * t259;
t229 = t293 * t274;
t284 = t259 * t229 - t249 * t253 + t288;
t226 = (pkin(3) * t272 + qJ(4) * t268) * t269;
t283 = qJD(4) * t224 + t220 * t226 + t286;
t282 = t253 * t228 - t229 * t252 + t287;
t181 = pkin(3) * t225 + qJ(4) * t224;
t281 = qJD(4) * t222 + t245 * t181 + t284;
t180 = pkin(3) * t223 + qJ(4) * t222;
t280 = qJD(4) * t269 * t268 + t221 * t180 + t282;
t190 = t223 * pkin(4) - pkin(9) * t304;
t232 = t269 * t272 * pkin(4) + t273 * pkin(9);
t279 = t220 * t232 + (-t180 - t190) * t245 + t283;
t191 = t225 * pkin(4) - pkin(9) * t303;
t278 = t245 * t191 + (-t226 - t232) * t221 + t281;
t277 = t221 * t190 + (-t181 - t191) * t220 + t280;
t204 = -Icges(3,6) * t274 + t270 * t290;
t205 = Icges(3,6) * t270 + t274 * t290;
t208 = -Icges(3,5) * t274 + t270 * t291;
t209 = Icges(3,5) * t270 + t274 * t291;
t239 = Icges(3,2) * t273 + t309;
t242 = Icges(3,1) * t269 + t308;
t276 = (-t205 * t269 + t209 * t273) * t253 + (-t204 * t269 + t208 * t273) * t252 + (-t239 * t269 + t242 * t273) * t259;
t266 = qJ(5) + qJ(6);
t264 = Icges(2,4) * t274;
t263 = cos(t266);
t262 = sin(t266);
t260 = qJD(5) * t273;
t248 = rSges(2,1) * t274 - rSges(2,2) * t270;
t247 = rSges(2,1) * t270 + rSges(2,2) * t274;
t246 = rSges(3,1) * t269 + rSges(3,2) * t273;
t244 = Icges(2,1) * t274 - t310;
t243 = Icges(2,1) * t270 + t264;
t241 = -Icges(2,2) * t270 + t264;
t240 = Icges(2,2) * t274 + t310;
t235 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t234 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t233 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t230 = t245 + t260;
t216 = (t271 * t272 + t305) * t269;
t215 = (-t267 * t272 + t268 * t271) * t269;
t214 = t260 + (-qJD(3) + qJD(6)) * t273 + t259;
t213 = rSges(3,3) * t270 + t274 * t292;
t212 = -rSges(3,3) * t274 + t270 * t292;
t211 = -rSges(4,3) * t273 + (rSges(4,1) * t272 - rSges(4,2) * t268) * t269;
t210 = -rSges(5,2) * t273 + (rSges(5,1) * t272 + rSges(5,3) * t268) * t269;
t196 = (t262 * t268 + t263 * t272) * t269;
t195 = (-t262 * t272 + t263 * t268) * t269;
t194 = -t274 * t299 + t221;
t193 = -t270 * t299 + t220;
t189 = V_base(5) * rSges(2,3) - t247 * t259 + t298;
t188 = t248 * t259 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t186 = t247 * V_base(4) - t248 * V_base(5) + V_base(3);
t185 = t274 * t294 + t221;
t184 = t270 * t294 + t220;
t179 = t225 * t271 + t306;
t178 = t224 * t271 - t225 * t267;
t177 = t223 * t271 + t307;
t176 = t222 * t271 - t223 * t267;
t175 = pkin(10) * t273 + (pkin(5) * t305 + t272 * t312) * t269;
t173 = t224 * t262 + t225 * t263;
t172 = t224 * t263 - t225 * t262;
t171 = t222 * t262 + t223 * t263;
t170 = t222 * t263 - t223 * t262;
t169 = rSges(4,1) * t225 - rSges(4,2) * t224 + rSges(4,3) * t303;
t168 = rSges(5,1) * t225 + rSges(5,2) * t303 + rSges(5,3) * t224;
t167 = rSges(4,1) * t223 - rSges(4,2) * t222 + rSges(4,3) * t304;
t166 = rSges(5,1) * t223 + rSges(5,2) * t304 + rSges(5,3) * t222;
t153 = rSges(6,1) * t216 + rSges(6,2) * t215 + rSges(6,3) * t273;
t151 = Icges(6,1) * t216 + Icges(6,4) * t215 + Icges(6,5) * t273;
t150 = Icges(6,4) * t216 + Icges(6,2) * t215 + Icges(6,6) * t273;
t149 = Icges(6,5) * t216 + Icges(6,6) * t215 + Icges(6,3) * t273;
t147 = rSges(7,1) * t196 + rSges(7,2) * t195 + rSges(7,3) * t273;
t146 = Icges(7,1) * t196 + Icges(7,4) * t195 + Icges(7,5) * t273;
t145 = Icges(7,4) * t196 + Icges(7,2) * t195 + Icges(7,6) * t273;
t144 = Icges(7,5) * t196 + Icges(7,6) * t195 + Icges(7,3) * t273;
t143 = t246 * t252 + (-t212 - t250) * t259 + t298;
t142 = t213 * t259 - t246 * t253 + t288;
t141 = pkin(5) * t306 + t225 * t312 - t274 * t295;
t140 = pkin(5) * t307 + t223 * t312 - t270 * t295;
t139 = t212 * t253 - t213 * t252 + t287;
t138 = rSges(6,1) * t179 + rSges(6,2) * t178 - rSges(6,3) * t303;
t137 = rSges(6,1) * t177 + rSges(6,2) * t176 - rSges(6,3) * t304;
t136 = Icges(6,1) * t179 + Icges(6,4) * t178 - Icges(6,5) * t303;
t135 = Icges(6,1) * t177 + Icges(6,4) * t176 - Icges(6,5) * t304;
t134 = Icges(6,4) * t179 + Icges(6,2) * t178 - Icges(6,6) * t303;
t133 = Icges(6,4) * t177 + Icges(6,2) * t176 - Icges(6,6) * t304;
t132 = Icges(6,5) * t179 + Icges(6,6) * t178 - Icges(6,3) * t303;
t131 = Icges(6,5) * t177 + Icges(6,6) * t176 - Icges(6,3) * t304;
t130 = rSges(7,1) * t173 + rSges(7,2) * t172 - rSges(7,3) * t303;
t129 = rSges(7,1) * t171 + rSges(7,2) * t170 - rSges(7,3) * t304;
t128 = Icges(7,1) * t173 + Icges(7,4) * t172 - Icges(7,5) * t303;
t127 = Icges(7,1) * t171 + Icges(7,4) * t170 - Icges(7,5) * t304;
t126 = Icges(7,4) * t173 + Icges(7,2) * t172 - Icges(7,6) * t303;
t125 = Icges(7,4) * t171 + Icges(7,2) * t170 - Icges(7,6) * t304;
t124 = Icges(7,5) * t173 + Icges(7,6) * t172 - Icges(7,3) * t303;
t123 = Icges(7,5) * t171 + Icges(7,6) * t170 - Icges(7,3) * t304;
t122 = -t167 * t245 + t211 * t220 + t286;
t121 = t169 * t245 - t211 * t221 + t284;
t120 = t167 * t221 - t169 * t220 + t282;
t119 = t210 * t220 + (-t166 - t180) * t245 + t283;
t118 = t168 * t245 + (-t210 - t226) * t221 + t281;
t117 = t166 * t221 + (-t168 - t181) * t220 + t280;
t116 = -t137 * t230 + t153 * t193 + t279;
t115 = t138 * t230 - t153 * t194 + t278;
t114 = t137 * t194 - t138 * t193 + t277;
t113 = -t129 * t214 - t140 * t230 + t147 * t184 + t175 * t193 + t279;
t112 = t130 * t214 + t141 * t230 - t147 * t185 - t175 * t194 + t278;
t111 = t129 * t185 - t130 * t184 + t140 * t194 - t141 * t193 + t277;
t1 = t193 * ((-t132 * t304 + t134 * t176 + t136 * t177) * t194 + (-t131 * t304 + t176 * t133 + t177 * t135) * t193 + (-t149 * t304 + t150 * t176 + t151 * t177) * t230) / 0.2e1 + t184 * ((-t124 * t304 + t126 * t170 + t128 * t171) * t185 + (-t123 * t304 + t170 * t125 + t171 * t127) * t184 + (-t144 * t304 + t145 * t170 + t146 * t171) * t214) / 0.2e1 + t194 * ((-t132 * t303 + t178 * t134 + t179 * t136) * t194 + (-t131 * t303 + t133 * t178 + t135 * t179) * t193 + (-t149 * t303 + t150 * t178 + t151 * t179) * t230) / 0.2e1 + t185 * ((-t124 * t303 + t172 * t126 + t173 * t128) * t185 + (-t123 * t303 + t125 * t172 + t127 * t173) * t184 + (-t144 * t303 + t145 * t172 + t146 * t173) * t214) / 0.2e1 + t230 * ((t132 * t273 + t134 * t215 + t136 * t216) * t194 + (t131 * t273 + t133 * t215 + t135 * t216) * t193 + (t273 * t149 + t215 * t150 + t216 * t151) * t230) / 0.2e1 + t214 * ((t124 * t273 + t126 * t195 + t128 * t196) * t185 + (t123 * t273 + t125 * t195 + t127 * t196) * t184 + (t273 * t144 + t195 * t145 + t196 * t146) * t214) / 0.2e1 + m(1) * (t233 ^ 2 + t234 ^ 2 + t235 ^ 2) / 0.2e1 + m(2) * (t186 ^ 2 + t188 ^ 2 + t189 ^ 2) / 0.2e1 + m(3) * (t139 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + t253 * (t285 * t270 + t276 * t274) / 0.2e1 + t252 * (t276 * t270 - t285 * t274) / 0.2e1 + m(4) * (t120 ^ 2 + t121 ^ 2 + t122 ^ 2) / 0.2e1 + m(7) * (t111 ^ 2 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + m(6) * (t114 ^ 2 + t115 ^ 2 + t116 ^ 2) / 0.2e1 + m(5) * (t117 ^ 2 + t118 ^ 2 + t119 ^ 2) / 0.2e1 + ((t222 * t319 + t223 * t317 + t304 * t318) * t245 + (t222 * t324 + t223 * t320 + t304 * t322) * t221 + (t222 * t325 + t321 * t223 + t323 * t304) * t220) * t220 / 0.2e1 + ((t224 * t319 + t225 * t317 + t303 * t318) * t245 + (t224 * t324 + t225 * t320 + t303 * t322) * t221 + (t224 * t325 + t321 * t225 + t323 * t303) * t220) * t221 / 0.2e1 + ((-t220 * t323 - t221 * t322 - t245 * t318) * t273 + ((t268 * t319 + t272 * t317) * t245 + (t268 * t324 + t272 * t320) * t221 + (t268 * t325 + t321 * t272) * t220) * t269) * t245 / 0.2e1 + ((t205 * t273 + t209 * t269) * t253 + (t204 * t273 + t208 * t269) * t252 + (t273 * t239 + t269 * t242 + Icges(2,3)) * t259) * t259 / 0.2e1 + ((-t240 * t270 + t243 * t274 + Icges(1,4)) * V_base(5) + (-t270 * t241 + t274 * t244 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t274 * t240 + t270 * t243 + Icges(1,2)) * V_base(5) + (t241 * t274 + t244 * t270 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t259 * (Icges(2,5) * t274 - Icges(2,6) * t270) + V_base(5) * t259 * (Icges(2,5) * t270 + Icges(2,6) * t274) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
