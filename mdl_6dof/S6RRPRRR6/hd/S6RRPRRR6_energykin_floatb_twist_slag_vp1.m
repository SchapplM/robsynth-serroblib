% Calculate kinetic energy for
% S6RRPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 13:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR6_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR6_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRRR6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR6_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR6_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR6_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR6_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:50:23
% EndTime: 2019-03-09 13:50:27
% DurationCPUTime: 4.67s
% Computational Cost: add. (1811->344), mult. (2628->499), div. (0->0), fcn. (2722->10), ass. (0->167)
t345 = Icges(3,4) - Icges(4,5);
t344 = Icges(3,1) + Icges(4,1);
t343 = Icges(3,2) + Icges(4,3);
t269 = sin(qJ(2));
t342 = t345 * t269;
t273 = cos(qJ(2));
t341 = t345 * t273;
t340 = Icges(4,4) + Icges(3,5);
t339 = Icges(3,6) - Icges(4,6);
t338 = t343 * t269 - t341;
t337 = t344 * t273 - t342;
t336 = Icges(4,2) + Icges(3,3);
t270 = sin(qJ(1));
t274 = cos(qJ(1));
t335 = t338 * t270 + t339 * t274;
t334 = -t339 * t270 + t338 * t274;
t333 = t337 * t270 - t340 * t274;
t332 = t340 * t270 + t337 * t274;
t331 = -t343 * t273 - t342;
t330 = t344 * t269 + t341;
t329 = -t339 * t269 + t340 * t273;
t255 = -qJD(2) * t274 + V_base(5);
t256 = qJD(2) * t270 + V_base(4);
t260 = V_base(6) + qJD(1);
t328 = (t269 * t331 + t273 * t330) * t260 + (t269 * t334 + t273 * t332) * t256 + (t269 * t335 + t273 * t333) * t255;
t327 = (t340 * t269 + t339 * t273) * t260 + (t270 * t336 + t329 * t274) * t256 + (t329 * t270 - t336 * t274) * t255;
t311 = qJ(4) + qJ(5);
t264 = sin(t311);
t304 = cos(t311);
t301 = t269 * t304;
t219 = -t273 * t264 + t301;
t323 = pkin(3) * t269;
t272 = cos(qJ(4));
t322 = pkin(4) * t272;
t320 = Icges(2,4) * t270;
t268 = sin(qJ(4));
t315 = t268 * t269;
t314 = t268 * t273;
t312 = t273 * t274;
t298 = pkin(2) * t273 + qJ(3) * t269;
t220 = t298 * t270;
t253 = pkin(1) * t270 - pkin(7) * t274;
t310 = -t220 - t253;
t309 = qJD(3) * t269;
t308 = V_base(5) * pkin(6) + V_base(1);
t229 = t270 * t273 * pkin(3) + t274 * pkin(8);
t305 = -t229 + t310;
t286 = pkin(4) * t315 + t273 * t322;
t166 = pkin(9) * t274 + t270 * t286;
t303 = -t166 + t305;
t248 = pkin(2) * t269 - qJ(3) * t273;
t302 = t255 * t248 + t274 * t309 + t308;
t300 = rSges(3,1) * t273 - rSges(3,2) * t269;
t299 = rSges(4,1) * t273 + rSges(4,3) * t269;
t226 = t269 * t272 - t314;
t291 = t272 * t273 + t315;
t290 = t255 * t323 + t302;
t223 = qJD(4) * t274 + t255;
t254 = pkin(1) * t274 + pkin(7) * t270;
t289 = -V_base(4) * pkin(6) + t260 * t254 + V_base(2);
t288 = V_base(4) * t253 - t254 * V_base(5) + V_base(3);
t191 = -pkin(4) * t314 + t269 * t322;
t287 = t223 * t191 + t290;
t215 = qJD(5) * t274 + t223;
t216 = (-qJD(4) - qJD(5)) * t270 + t256;
t218 = t269 * t264 + t273 * t304;
t221 = t298 * t274;
t283 = t260 * t221 + t270 * t309 + t289;
t282 = -qJD(3) * t273 + t256 * t220 + t288;
t230 = pkin(3) * t312 - t270 * pkin(8);
t281 = t260 * t230 + (-t248 - t323) * t256 + t283;
t280 = t256 * t229 + (-t221 - t230) * t255 + t282;
t167 = -pkin(9) * t270 + t274 * t286;
t224 = -qJD(4) * t270 + t256;
t279 = t260 * t167 - t191 * t224 + t281;
t278 = t224 * t166 - t167 * t223 + t280;
t271 = cos(qJ(6));
t267 = sin(qJ(6));
t265 = Icges(2,4) * t274;
t252 = rSges(2,1) * t274 - rSges(2,2) * t270;
t251 = rSges(2,1) * t270 + rSges(2,2) * t274;
t250 = rSges(3,1) * t269 + rSges(3,2) * t273;
t249 = rSges(4,1) * t269 - rSges(4,3) * t273;
t247 = Icges(2,1) * t274 - t320;
t246 = Icges(2,1) * t270 + t265;
t243 = -Icges(2,2) * t270 + t265;
t242 = Icges(2,2) * t274 + t320;
t233 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t232 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t231 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t214 = t291 * t274;
t213 = t226 * t274;
t212 = t291 * t270;
t211 = t226 * t270;
t209 = rSges(3,3) * t270 + t274 * t300;
t208 = rSges(4,2) * t270 + t274 * t299;
t207 = -rSges(3,3) * t274 + t270 * t300;
t206 = -rSges(4,2) * t274 + t270 * t299;
t190 = t218 * t274;
t189 = t264 * t312 - t274 * t301;
t188 = t218 * t270;
t187 = t219 * t270;
t186 = qJD(6) * t218 + t260;
t184 = V_base(5) * rSges(2,3) - t251 * t260 + t308;
t183 = t252 * t260 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t182 = t251 * V_base(4) - t252 * V_base(5) + V_base(3);
t181 = rSges(5,1) * t226 - rSges(5,2) * t291;
t180 = Icges(5,1) * t226 - Icges(5,4) * t291;
t179 = Icges(5,4) * t226 - Icges(5,2) * t291;
t178 = Icges(5,5) * t226 - Icges(5,6) * t291;
t177 = t190 * t271 - t267 * t270;
t176 = -t190 * t267 - t270 * t271;
t175 = t188 * t271 + t267 * t274;
t174 = -t188 * t267 + t271 * t274;
t172 = pkin(5) * t219 + pkin(10) * t218;
t171 = rSges(6,1) * t219 - rSges(6,2) * t218;
t170 = Icges(6,1) * t219 - Icges(6,4) * t218;
t169 = Icges(6,4) * t219 - Icges(6,2) * t218;
t168 = Icges(6,5) * t219 - Icges(6,6) * t218;
t165 = qJD(6) * t189 + t216;
t164 = -qJD(6) * t187 + t215;
t162 = rSges(5,1) * t214 + rSges(5,2) * t213 - rSges(5,3) * t270;
t161 = rSges(5,1) * t212 + rSges(5,2) * t211 + rSges(5,3) * t274;
t160 = Icges(5,1) * t214 + Icges(5,4) * t213 - Icges(5,5) * t270;
t159 = Icges(5,1) * t212 + Icges(5,4) * t211 + Icges(5,5) * t274;
t158 = Icges(5,4) * t214 + Icges(5,2) * t213 - Icges(5,6) * t270;
t157 = Icges(5,4) * t212 + Icges(5,2) * t211 + Icges(5,6) * t274;
t156 = Icges(5,5) * t214 + Icges(5,6) * t213 - Icges(5,3) * t270;
t155 = Icges(5,5) * t212 + Icges(5,6) * t211 + Icges(5,3) * t274;
t154 = pkin(5) * t190 + pkin(10) * t189;
t153 = pkin(5) * t188 - pkin(10) * t187;
t152 = rSges(6,1) * t190 - rSges(6,2) * t189 - rSges(6,3) * t270;
t151 = rSges(6,1) * t188 + rSges(6,2) * t187 + rSges(6,3) * t274;
t149 = Icges(6,1) * t190 - Icges(6,4) * t189 - Icges(6,5) * t270;
t148 = Icges(6,1) * t188 + Icges(6,4) * t187 + Icges(6,5) * t274;
t147 = Icges(6,4) * t190 - Icges(6,2) * t189 - Icges(6,6) * t270;
t146 = Icges(6,4) * t188 + Icges(6,2) * t187 + Icges(6,6) * t274;
t145 = Icges(6,5) * t190 - Icges(6,6) * t189 - Icges(6,3) * t270;
t144 = Icges(6,5) * t188 + Icges(6,6) * t187 + Icges(6,3) * t274;
t143 = rSges(7,3) * t218 + (rSges(7,1) * t271 - rSges(7,2) * t267) * t219;
t142 = Icges(7,5) * t218 + (Icges(7,1) * t271 - Icges(7,4) * t267) * t219;
t141 = Icges(7,6) * t218 + (Icges(7,4) * t271 - Icges(7,2) * t267) * t219;
t140 = Icges(7,3) * t218 + (Icges(7,5) * t271 - Icges(7,6) * t267) * t219;
t139 = t250 * t255 + (-t207 - t253) * t260 + t308;
t138 = t209 * t260 - t250 * t256 + t289;
t137 = t207 * t256 - t209 * t255 + t288;
t136 = rSges(7,1) * t177 + rSges(7,2) * t176 + rSges(7,3) * t189;
t135 = rSges(7,1) * t175 + rSges(7,2) * t174 - rSges(7,3) * t187;
t134 = Icges(7,1) * t177 + Icges(7,4) * t176 + Icges(7,5) * t189;
t133 = Icges(7,1) * t175 + Icges(7,4) * t174 - Icges(7,5) * t187;
t132 = Icges(7,4) * t177 + Icges(7,2) * t176 + Icges(7,6) * t189;
t131 = Icges(7,4) * t175 + Icges(7,2) * t174 - Icges(7,6) * t187;
t130 = Icges(7,5) * t177 + Icges(7,6) * t176 + Icges(7,3) * t189;
t129 = Icges(7,5) * t175 + Icges(7,6) * t174 - Icges(7,3) * t187;
t128 = t249 * t255 + (-t206 + t310) * t260 + t302;
t127 = t208 * t260 + (-t248 - t249) * t256 + t283;
t126 = t206 * t256 + (-t208 - t221) * t255 + t282;
t125 = t181 * t223 + (-t161 + t305) * t260 + t290;
t124 = t162 * t260 - t181 * t224 + t281;
t123 = t161 * t224 - t162 * t223 + t280;
t122 = t171 * t215 + (-t151 + t303) * t260 + t287;
t121 = t152 * t260 - t171 * t216 + t279;
t120 = t151 * t216 - t152 * t215 + t278;
t119 = -t135 * t186 + t143 * t164 + t172 * t215 + (-t153 + t303) * t260 + t287;
t118 = t136 * t186 - t143 * t165 + t154 * t260 - t172 * t216 + t279;
t117 = t135 * t165 - t136 * t164 + t153 * t216 - t154 * t215 + t278;
t1 = m(7) * (t117 ^ 2 + t118 ^ 2 + t119 ^ 2) / 0.2e1 + m(6) * (t120 ^ 2 + t121 ^ 2 + t122 ^ 2) / 0.2e1 + t186 * ((t129 * t164 + t130 * t165 + t140 * t186) * t218 + ((-t132 * t267 + t134 * t271) * t165 + (-t131 * t267 + t133 * t271) * t164 + (-t141 * t267 + t142 * t271) * t186) * t219) / 0.2e1 + t164 * ((-t130 * t187 + t132 * t174 + t134 * t175) * t165 + (-t187 * t129 + t174 * t131 + t175 * t133) * t164 + (-t140 * t187 + t141 * t174 + t142 * t175) * t186) / 0.2e1 + t165 * ((t130 * t189 + t176 * t132 + t177 * t134) * t165 + (t129 * t189 + t131 * t176 + t133 * t177) * t164 + (t140 * t189 + t141 * t176 + t142 * t177) * t186) / 0.2e1 + m(2) * (t182 ^ 2 + t183 ^ 2 + t184 ^ 2) / 0.2e1 + t215 * ((t145 * t274 + t147 * t187 + t149 * t188) * t216 + (t144 * t274 + t187 * t146 + t188 * t148) * t215 + (t168 * t274 + t169 * t187 + t170 * t188) * t260) / 0.2e1 + t223 * ((t156 * t274 + t158 * t211 + t160 * t212) * t224 + (t155 * t274 + t211 * t157 + t212 * t159) * t223 + (t178 * t274 + t179 * t211 + t180 * t212) * t260) / 0.2e1 + t216 * ((-t270 * t145 - t189 * t147 + t149 * t190) * t216 + (-t144 * t270 - t146 * t189 + t148 * t190) * t215 + (-t168 * t270 - t169 * t189 + t170 * t190) * t260) / 0.2e1 + t224 * ((-t270 * t156 + t213 * t158 + t214 * t160) * t224 + (-t155 * t270 + t157 * t213 + t159 * t214) * t223 + (-t178 * t270 + t179 * t213 + t180 * t214) * t260) / 0.2e1 + m(3) * (t137 ^ 2 + t138 ^ 2 + t139 ^ 2) / 0.2e1 + m(5) * (t123 ^ 2 + t124 ^ 2 + t125 ^ 2) / 0.2e1 + m(4) * (t126 ^ 2 + t127 ^ 2 + t128 ^ 2) / 0.2e1 + m(1) * (t231 ^ 2 + t232 ^ 2 + t233 ^ 2) / 0.2e1 + (t328 * t270 - t327 * t274) * t255 / 0.2e1 + (t327 * t270 + t328 * t274) * t256 / 0.2e1 + ((-t242 * t270 + t246 * t274 + Icges(1,4)) * V_base(5) + (-t270 * t243 + t274 * t247 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t274 * t242 + t270 * t246 + Icges(1,2)) * V_base(5) + (t243 * t274 + t247 * t270 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t158 * t291 + t160 * t226) * t224 + (-t157 * t291 + t159 * t226) * t223 + (-t147 * t218 + t149 * t219) * t216 + (-t146 * t218 + t148 * t219) * t215 + (t269 * t332 - t273 * t334) * t256 + (t269 * t333 - t273 * t335) * t255 + (-t218 * t169 + t219 * t170 - t291 * t179 + t226 * t180 + t269 * t330 - t273 * t331 + Icges(2,3)) * t260) * t260 / 0.2e1 + t260 * V_base(4) * (Icges(2,5) * t274 - Icges(2,6) * t270) + t260 * V_base(5) * (Icges(2,5) * t270 + Icges(2,6) * t274) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
