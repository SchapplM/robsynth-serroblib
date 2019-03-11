% Calculate kinetic energy for
% S6RRPPRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 09:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRR8_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR8_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPPRR8_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR8_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR8_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRR8_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:22:47
% EndTime: 2019-03-09 09:22:50
% DurationCPUTime: 3.61s
% Computational Cost: add. (1717->373), mult. (3240->530), div. (0->0), fcn. (3520->10), ass. (0->171)
t334 = Icges(4,1) + Icges(5,1);
t333 = -Icges(4,4) + Icges(5,5);
t332 = Icges(5,4) + Icges(4,5);
t331 = Icges(4,2) + Icges(5,3);
t330 = -Icges(5,6) + Icges(4,6);
t329 = -Icges(4,3) - Icges(5,2);
t267 = sin(pkin(10));
t268 = cos(pkin(10));
t274 = cos(qJ(1));
t271 = sin(qJ(1));
t273 = cos(qJ(2));
t306 = t271 * t273;
t220 = t267 * t306 + t268 * t274;
t221 = -t267 * t274 + t268 * t306;
t270 = sin(qJ(2));
t308 = t270 * t271;
t328 = t331 * t220 + t333 * t221 - t330 * t308;
t305 = t273 * t274;
t222 = t267 * t305 - t271 * t268;
t223 = t267 * t271 + t268 * t305;
t307 = t270 * t274;
t327 = t331 * t222 + t333 * t223 - t330 * t307;
t326 = -t330 * t220 + t332 * t221 - t329 * t308;
t325 = -t330 * t222 + t332 * t223 - t329 * t307;
t324 = t333 * t220 + t334 * t221 + t332 * t308;
t323 = t333 * t222 + t334 * t223 + t332 * t307;
t322 = t330 * t273 + (t331 * t267 + t333 * t268) * t270;
t321 = t329 * t273 + (-t330 * t267 + t332 * t268) * t270;
t320 = -t332 * t273 + (t333 * t267 + t334 * t268) * t270;
t272 = cos(qJ(5));
t316 = pkin(5) * t272;
t314 = Icges(2,4) * t271;
t313 = Icges(3,4) * t270;
t312 = Icges(3,4) * t273;
t269 = sin(qJ(5));
t311 = t220 * t269;
t310 = t222 * t269;
t309 = t267 * t269;
t183 = pkin(3) * t223 + qJ(4) * t222;
t291 = pkin(2) * t273 + qJ(3) * t270;
t228 = t291 * t274;
t304 = -t183 - t228;
t226 = (pkin(3) * t268 + qJ(4) * t267) * t270;
t246 = pkin(2) * t270 - qJ(3) * t273;
t303 = -t226 - t246;
t227 = t291 * t271;
t250 = pkin(1) * t271 - pkin(7) * t274;
t302 = -t227 - t250;
t301 = qJD(3) * t270;
t300 = qJD(5) * t270;
t299 = V_base(5) * pkin(6) + V_base(1);
t182 = pkin(3) * t221 + qJ(4) * t220;
t296 = -t182 + t302;
t253 = qJD(2) * t271 + V_base(4);
t295 = t270 * pkin(9);
t259 = V_base(6) + qJD(1);
t294 = t270 * (-qJD(5) - qJD(6));
t245 = qJD(5) * t273 + t259;
t252 = -qJD(2) * t274 + V_base(5);
t293 = t252 * t246 + t274 * t301 + t299;
t292 = rSges(3,1) * t273 - rSges(3,2) * t270;
t290 = Icges(3,1) * t273 - t313;
t289 = -Icges(3,2) * t270 + t312;
t288 = Icges(3,5) * t273 - Icges(3,6) * t270;
t251 = pkin(1) * t274 + pkin(7) * t271;
t287 = -V_base(4) * pkin(6) + t259 * t251 + V_base(2);
t286 = V_base(4) * t250 - t251 * V_base(5) + V_base(3);
t285 = qJD(4) * t222 + t252 * t226 + t293;
t284 = (-Icges(3,3) * t274 + t271 * t288) * t252 + (Icges(3,3) * t271 + t274 * t288) * t253 + (Icges(3,5) * t270 + Icges(3,6) * t273) * t259;
t283 = t259 * t228 + t271 * t301 + t287;
t282 = -qJD(3) * t273 + t253 * t227 + t286;
t281 = qJD(4) * t220 + t259 * t183 + t283;
t280 = qJD(4) * t270 * t267 + t253 * t182 + t282;
t189 = t221 * pkin(4) - pkin(8) * t308;
t232 = t270 * t268 * pkin(4) + t273 * pkin(8);
t279 = t252 * t232 + (-t189 + t296) * t259 + t285;
t190 = t223 * pkin(4) - pkin(8) * t307;
t278 = t259 * t190 + (-t232 + t303) * t253 + t281;
t277 = t253 * t189 + (-t190 + t304) * t252 + t280;
t209 = -Icges(3,6) * t274 + t271 * t289;
t210 = Icges(3,6) * t271 + t274 * t289;
t211 = -Icges(3,5) * t274 + t271 * t290;
t212 = Icges(3,5) * t271 + t274 * t290;
t239 = Icges(3,2) * t273 + t313;
t242 = Icges(3,1) * t270 + t312;
t276 = (-t210 * t270 + t212 * t273) * t253 + (-t209 * t270 + t211 * t273) * t252 + (-t239 * t270 + t242 * t273) * t259;
t266 = qJ(5) + qJ(6);
t264 = Icges(2,4) * t274;
t263 = cos(t266);
t262 = sin(t266);
t249 = rSges(2,1) * t274 - rSges(2,2) * t271;
t248 = rSges(2,1) * t271 + rSges(2,2) * t274;
t247 = rSges(3,1) * t270 + rSges(3,2) * t273;
t244 = Icges(2,1) * t274 - t314;
t243 = Icges(2,1) * t271 + t264;
t241 = -Icges(2,2) * t271 + t264;
t240 = Icges(2,2) * t274 + t314;
t235 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t234 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t233 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t230 = qJD(6) * t273 + t245;
t225 = -t274 * t300 + t253;
t224 = -t271 * t300 + t252;
t216 = (t268 * t272 + t309) * t270;
t215 = (t267 * t272 - t268 * t269) * t270;
t214 = rSges(3,3) * t271 + t274 * t292;
t213 = -rSges(3,3) * t274 + t271 * t292;
t206 = -rSges(4,3) * t273 + (rSges(4,1) * t268 - rSges(4,2) * t267) * t270;
t205 = -rSges(5,2) * t273 + (rSges(5,1) * t268 + rSges(5,3) * t267) * t270;
t196 = (t262 * t267 + t263 * t268) * t270;
t195 = (-t262 * t268 + t263 * t267) * t270;
t194 = t274 * t294 + t253;
t193 = t271 * t294 + t252;
t188 = V_base(5) * rSges(2,3) - t248 * t259 + t299;
t187 = t249 * t259 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t186 = t248 * V_base(4) - t249 * V_base(5) + V_base(3);
t181 = t223 * t272 + t310;
t180 = t222 * t272 - t223 * t269;
t179 = t221 * t272 + t311;
t178 = t220 * t272 - t221 * t269;
t177 = pkin(9) * t273 + (pkin(5) * t309 + t268 * t316) * t270;
t175 = t222 * t262 + t223 * t263;
t174 = t222 * t263 - t223 * t262;
t173 = t220 * t262 + t221 * t263;
t172 = t220 * t263 - t221 * t262;
t171 = rSges(4,1) * t223 - rSges(4,2) * t222 + rSges(4,3) * t307;
t170 = rSges(5,1) * t223 + rSges(5,2) * t307 + rSges(5,3) * t222;
t169 = rSges(4,1) * t221 - rSges(4,2) * t220 + rSges(4,3) * t308;
t168 = rSges(5,1) * t221 + rSges(5,2) * t308 + rSges(5,3) * t220;
t154 = rSges(6,1) * t216 + rSges(6,2) * t215 + rSges(6,3) * t273;
t153 = Icges(6,1) * t216 + Icges(6,4) * t215 + Icges(6,5) * t273;
t152 = Icges(6,4) * t216 + Icges(6,2) * t215 + Icges(6,6) * t273;
t151 = Icges(6,5) * t216 + Icges(6,6) * t215 + Icges(6,3) * t273;
t150 = rSges(7,1) * t196 + rSges(7,2) * t195 + rSges(7,3) * t273;
t149 = Icges(7,1) * t196 + Icges(7,4) * t195 + Icges(7,5) * t273;
t148 = Icges(7,4) * t196 + Icges(7,2) * t195 + Icges(7,6) * t273;
t147 = Icges(7,5) * t196 + Icges(7,6) * t195 + Icges(7,3) * t273;
t146 = t247 * t252 + (-t213 - t250) * t259 + t299;
t145 = t214 * t259 - t247 * t253 + t287;
t144 = pkin(5) * t310 + t223 * t316 - t274 * t295;
t143 = pkin(5) * t311 + t221 * t316 - t271 * t295;
t142 = t213 * t253 - t214 * t252 + t286;
t141 = rSges(6,1) * t181 + rSges(6,2) * t180 - rSges(6,3) * t307;
t140 = rSges(6,1) * t179 + rSges(6,2) * t178 - rSges(6,3) * t308;
t139 = Icges(6,1) * t181 + Icges(6,4) * t180 - Icges(6,5) * t307;
t138 = Icges(6,1) * t179 + Icges(6,4) * t178 - Icges(6,5) * t308;
t137 = Icges(6,4) * t181 + Icges(6,2) * t180 - Icges(6,6) * t307;
t136 = Icges(6,4) * t179 + Icges(6,2) * t178 - Icges(6,6) * t308;
t135 = Icges(6,5) * t181 + Icges(6,6) * t180 - Icges(6,3) * t307;
t134 = Icges(6,5) * t179 + Icges(6,6) * t178 - Icges(6,3) * t308;
t133 = rSges(7,1) * t175 + rSges(7,2) * t174 - rSges(7,3) * t307;
t132 = rSges(7,1) * t173 + rSges(7,2) * t172 - rSges(7,3) * t308;
t131 = Icges(7,1) * t175 + Icges(7,4) * t174 - Icges(7,5) * t307;
t130 = Icges(7,1) * t173 + Icges(7,4) * t172 - Icges(7,5) * t308;
t129 = Icges(7,4) * t175 + Icges(7,2) * t174 - Icges(7,6) * t307;
t128 = Icges(7,4) * t173 + Icges(7,2) * t172 - Icges(7,6) * t308;
t127 = Icges(7,5) * t175 + Icges(7,6) * t174 - Icges(7,3) * t307;
t126 = Icges(7,5) * t173 + Icges(7,6) * t172 - Icges(7,3) * t308;
t125 = t206 * t252 + (-t169 + t302) * t259 + t293;
t124 = t171 * t259 + (-t206 - t246) * t253 + t283;
t123 = t169 * t253 + (-t171 - t228) * t252 + t282;
t122 = t205 * t252 + (-t168 + t296) * t259 + t285;
t121 = t170 * t259 + (-t205 + t303) * t253 + t281;
t120 = t168 * t253 + (-t170 + t304) * t252 + t280;
t119 = -t140 * t245 + t154 * t224 + t279;
t118 = t141 * t245 - t154 * t225 + t278;
t117 = t140 * t225 - t141 * t224 + t277;
t116 = -t132 * t230 - t143 * t245 + t150 * t193 + t177 * t224 + t279;
t115 = t133 * t230 + t144 * t245 - t150 * t194 - t177 * t225 + t278;
t114 = t132 * t194 - t133 * t193 + t143 * t225 - t144 * t224 + t277;
t1 = t224 * ((-t135 * t308 + t137 * t178 + t139 * t179) * t225 + (-t134 * t308 + t178 * t136 + t179 * t138) * t224 + (-t151 * t308 + t152 * t178 + t153 * t179) * t245) / 0.2e1 + t193 * ((-t127 * t308 + t129 * t172 + t131 * t173) * t194 + (-t126 * t308 + t172 * t128 + t173 * t130) * t193 + (-t147 * t308 + t148 * t172 + t149 * t173) * t230) / 0.2e1 + t225 * ((-t135 * t307 + t180 * t137 + t181 * t139) * t225 + (-t134 * t307 + t136 * t180 + t138 * t181) * t224 + (-t151 * t307 + t152 * t180 + t153 * t181) * t245) / 0.2e1 + t194 * ((-t127 * t307 + t174 * t129 + t175 * t131) * t194 + (-t126 * t307 + t128 * t174 + t130 * t175) * t193 + (-t147 * t307 + t148 * t174 + t149 * t175) * t230) / 0.2e1 + t245 * ((t135 * t273 + t137 * t215 + t139 * t216) * t225 + (t134 * t273 + t136 * t215 + t138 * t216) * t224 + (t273 * t151 + t215 * t152 + t216 * t153) * t245) / 0.2e1 + t230 * ((t127 * t273 + t129 * t195 + t131 * t196) * t194 + (t126 * t273 + t128 * t195 + t130 * t196) * t193 + (t273 * t147 + t195 * t148 + t196 * t149) * t230) / 0.2e1 + m(1) * (t233 ^ 2 + t234 ^ 2 + t235 ^ 2) / 0.2e1 + m(2) * (t186 ^ 2 + t187 ^ 2 + t188 ^ 2) / 0.2e1 + m(3) * (t142 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(6) * (t117 ^ 2 + t118 ^ 2 + t119 ^ 2) / 0.2e1 + m(5) * (t120 ^ 2 + t121 ^ 2 + t122 ^ 2) / 0.2e1 + m(4) * (t123 ^ 2 + t124 ^ 2 + t125 ^ 2) / 0.2e1 + m(7) * (t114 ^ 2 + t115 ^ 2 + t116 ^ 2) / 0.2e1 + ((-t240 * t271 + t243 * t274 + Icges(1,4)) * V_base(5) + (-t241 * t271 + t244 * t274 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t240 * t274 + t243 * t271 + Icges(1,2)) * V_base(5) + (t241 * t274 + t244 * t271 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (t271 * t276 - t284 * t274 + (t220 * t322 + t221 * t320 + t308 * t321) * t259 + (t220 * t327 + t221 * t323 + t308 * t325) * t253 + (t220 * t328 + t324 * t221 + t326 * t308) * t252) * t252 / 0.2e1 + (t271 * t284 + t274 * t276 + (t222 * t322 + t223 * t320 + t307 * t321) * t259 + (t222 * t327 + t223 * t323 + t307 * t325) * t253 + (t222 * t328 + t324 * t223 + t326 * t307) * t252) * t253 / 0.2e1 + (((t210 - t325) * t253 + (t209 - t326) * t252) * t273 + ((t267 * t327 + t268 * t323 + t212) * t253 + (t267 * t328 + t324 * t268 + t211) * t252) * t270 + (Icges(2,3) + (t239 - t321) * t273 + (t267 * t322 + t268 * t320 + t242) * t270) * t259) * t259 / 0.2e1 + t259 * V_base(4) * (Icges(2,5) * t274 - Icges(2,6) * t271) + t259 * V_base(5) * (Icges(2,5) * t271 + Icges(2,6) * t274) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
