% Calculate kinetic energy for
% S6RRPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
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
% Datum: 2019-03-09 09:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPP3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP3_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRPP3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP3_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPP3_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPP3_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:54:09
% EndTime: 2019-03-09 09:54:12
% DurationCPUTime: 3.13s
% Computational Cost: add. (1796->311), mult. (2532->439), div. (0->0), fcn. (2508->8), ass. (0->144)
t319 = Icges(5,1) + Icges(6,2) + Icges(7,3);
t318 = -Icges(5,4) - Icges(6,6) + Icges(7,6);
t317 = -Icges(5,5) - Icges(7,5) + Icges(6,4);
t316 = Icges(5,2) + Icges(7,2) + Icges(6,3);
t315 = Icges(5,6) - Icges(6,5) - Icges(7,4);
t314 = -Icges(5,3) - Icges(7,1) - Icges(6,1);
t313 = rSges(7,1) + pkin(5);
t312 = rSges(7,3) + qJ(6);
t249 = pkin(9) + qJ(4);
t243 = sin(t249);
t244 = cos(t249);
t256 = cos(qJ(1));
t254 = sin(qJ(1));
t255 = cos(qJ(2));
t287 = t254 * t255;
t191 = t243 * t287 + t244 * t256;
t192 = -t243 * t256 + t244 * t287;
t253 = sin(qJ(2));
t290 = t253 * t254;
t311 = t318 * t191 + t319 * t192 - t317 * t290;
t286 = t255 * t256;
t193 = t243 * t286 - t254 * t244;
t194 = t254 * t243 + t244 * t286;
t289 = t253 * t256;
t310 = t318 * t193 + t319 * t194 - t317 * t289;
t309 = t316 * t191 + t318 * t192 - t315 * t290;
t308 = t316 * t193 + t318 * t194 - t315 * t289;
t307 = -t315 * t191 - t317 * t192 - t314 * t290;
t306 = -t315 * t193 - t317 * t194 - t314 * t289;
t305 = t314 * t255 + (-t315 * t243 - t317 * t244) * t253;
t304 = t315 * t255 + (t316 * t243 + t318 * t244) * t253;
t303 = t317 * t255 + (t318 * t243 + t319 * t244) * t253;
t251 = cos(pkin(9));
t296 = pkin(3) * t251;
t295 = Icges(2,4) * t254;
t294 = Icges(3,4) * t253;
t293 = Icges(3,4) * t255;
t292 = t244 * t253;
t250 = sin(pkin(9));
t291 = t250 * t256;
t288 = t254 * t250;
t284 = rSges(7,2) * t191 + t312 * t192 + t290 * t313;
t283 = t193 * rSges(7,2) + t312 * t194 + t289 * t313;
t282 = (rSges(7,2) * t243 + rSges(7,3) * t244) * t253 + qJ(6) * t292 - t313 * t255;
t273 = pkin(2) * t255 + qJ(3) * t253;
t211 = t273 * t254;
t234 = t254 * pkin(1) - pkin(7) * t256;
t281 = -t211 - t234;
t280 = qJD(3) * t253;
t279 = qJD(4) * t253;
t278 = V_base(5) * pkin(6) + V_base(1);
t239 = qJD(2) * t254 + V_base(4);
t245 = V_base(6) + qJD(1);
t230 = pkin(2) * t253 - qJ(3) * t255;
t238 = -qJD(2) * t256 + V_base(5);
t275 = t238 * t230 + t256 * t280 + t278;
t274 = rSges(3,1) * t255 - rSges(3,2) * t253;
t272 = Icges(3,1) * t255 - t294;
t271 = -Icges(3,2) * t253 + t293;
t270 = Icges(3,5) * t255 - Icges(3,6) * t253;
t235 = pkin(1) * t256 + t254 * pkin(7);
t269 = -V_base(4) * pkin(6) + t245 * t235 + V_base(2);
t268 = V_base(4) * t234 - t235 * V_base(5) + V_base(3);
t267 = (-Icges(3,3) * t256 + t254 * t270) * t238 + (Icges(3,3) * t254 + t256 * t270) * t239 + (Icges(3,5) * t253 + Icges(3,6) * t255) * t245;
t266 = pkin(8) * t253 + t255 * t296;
t212 = t273 * t256;
t265 = t245 * t212 + t254 * t280 + t269;
t264 = -qJD(3) * t255 + t239 * t211 + t268;
t161 = -pkin(3) * t291 + t254 * t266;
t170 = -pkin(8) * t255 + t253 * t296;
t263 = t238 * t170 + (-t161 + t281) * t245 + t275;
t201 = (pkin(4) * t244 + qJ(5) * t243) * t253;
t209 = t254 * t279 + t238;
t262 = qJD(5) * t193 + t209 * t201 + t263;
t162 = pkin(3) * t288 + t256 * t266;
t261 = t245 * t162 + (-t170 - t230) * t239 + t265;
t260 = t239 * t161 + (-t162 - t212) * t238 + t264;
t152 = pkin(4) * t194 + qJ(5) * t193;
t229 = -qJD(4) * t255 + t245;
t259 = qJD(5) * t191 + t229 * t152 + t261;
t151 = pkin(4) * t192 + qJ(5) * t191;
t210 = t256 * t279 + t239;
t258 = qJD(5) * t253 * t243 + t210 * t151 + t260;
t197 = -Icges(3,6) * t256 + t254 * t271;
t198 = Icges(3,6) * t254 + t256 * t271;
t199 = -Icges(3,5) * t256 + t254 * t272;
t200 = Icges(3,5) * t254 + t256 * t272;
t222 = Icges(3,2) * t255 + t294;
t225 = Icges(3,1) * t253 + t293;
t257 = (-t198 * t253 + t200 * t255) * t239 + (-t197 * t253 + t199 * t255) * t238 + (-t222 * t253 + t225 * t255) * t245;
t247 = Icges(2,4) * t256;
t233 = rSges(2,1) * t256 - t254 * rSges(2,2);
t232 = t254 * rSges(2,1) + rSges(2,2) * t256;
t231 = rSges(3,1) * t253 + rSges(3,2) * t255;
t227 = Icges(2,1) * t256 - t295;
t226 = Icges(2,1) * t254 + t247;
t224 = -Icges(2,2) * t254 + t247;
t223 = Icges(2,2) * t256 + t295;
t218 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t217 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t216 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t208 = t251 * t286 + t288;
t207 = -t250 * t286 + t254 * t251;
t206 = t251 * t287 - t291;
t205 = -t250 * t287 - t251 * t256;
t203 = t254 * rSges(3,3) + t256 * t274;
t202 = -rSges(3,3) * t256 + t254 * t274;
t190 = -rSges(4,3) * t255 + (rSges(4,1) * t251 - rSges(4,2) * t250) * t253;
t188 = -Icges(4,5) * t255 + (Icges(4,1) * t251 - Icges(4,4) * t250) * t253;
t187 = -Icges(4,6) * t255 + (Icges(4,4) * t251 - Icges(4,2) * t250) * t253;
t186 = -Icges(4,3) * t255 + (Icges(4,5) * t251 - Icges(4,6) * t250) * t253;
t182 = -rSges(6,1) * t255 + (-rSges(6,2) * t244 + rSges(6,3) * t243) * t253;
t180 = -rSges(5,3) * t255 + (rSges(5,1) * t244 - rSges(5,2) * t243) * t253;
t169 = V_base(5) * rSges(2,3) - t232 * t245 + t278;
t168 = t233 * t245 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t167 = t232 * V_base(4) - t233 * V_base(5) + V_base(3);
t160 = t208 * rSges(4,1) + t207 * rSges(4,2) + rSges(4,3) * t289;
t159 = rSges(4,1) * t206 + rSges(4,2) * t205 + rSges(4,3) * t290;
t158 = Icges(4,1) * t208 + Icges(4,4) * t207 + Icges(4,5) * t289;
t157 = Icges(4,1) * t206 + Icges(4,4) * t205 + Icges(4,5) * t290;
t156 = Icges(4,4) * t208 + Icges(4,2) * t207 + Icges(4,6) * t289;
t155 = Icges(4,4) * t206 + Icges(4,2) * t205 + Icges(4,6) * t290;
t154 = Icges(4,5) * t208 + Icges(4,6) * t207 + Icges(4,3) * t289;
t153 = Icges(4,5) * t206 + Icges(4,6) * t205 + Icges(4,3) * t290;
t149 = t194 * rSges(5,1) - t193 * rSges(5,2) + rSges(5,3) * t289;
t148 = rSges(5,1) * t192 - rSges(5,2) * t191 + rSges(5,3) * t290;
t147 = rSges(6,1) * t289 - t194 * rSges(6,2) + t193 * rSges(6,3);
t145 = rSges(6,1) * t290 - rSges(6,2) * t192 + rSges(6,3) * t191;
t122 = t231 * t238 + (-t202 - t234) * t245 + t278;
t121 = t203 * t245 - t231 * t239 + t269;
t120 = t202 * t239 - t203 * t238 + t268;
t119 = t190 * t238 + (-t159 + t281) * t245 + t275;
t118 = t160 * t245 + (-t190 - t230) * t239 + t265;
t117 = t159 * t239 + (-t160 - t212) * t238 + t264;
t116 = -t148 * t229 + t180 * t209 + t263;
t115 = t149 * t229 - t180 * t210 + t261;
t114 = t148 * t210 - t149 * t209 + t260;
t113 = t182 * t209 + (-t145 - t151) * t229 + t262;
t112 = t147 * t229 + (-t182 - t201) * t210 + t259;
t111 = t145 * t210 + (-t147 - t152) * t209 + t258;
t110 = qJD(6) * t194 + t282 * t209 + (-t151 - t284) * t229 + t262;
t109 = qJD(6) * t192 + t283 * t229 + (-t201 - t282) * t210 + t259;
t108 = qJD(6) * t292 + t284 * t210 + (-t152 - t283) * t209 + t258;
t1 = m(1) * (t216 ^ 2 + t217 ^ 2 + t218 ^ 2) / 0.2e1 + m(2) * (t167 ^ 2 + t168 ^ 2 + t169 ^ 2) / 0.2e1 + m(4) * (t117 ^ 2 + t118 ^ 2 + t119 ^ 2) / 0.2e1 + m(3) * (t120 ^ 2 + t121 ^ 2 + t122 ^ 2) / 0.2e1 + m(7) * (t108 ^ 2 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + m(6) * (t111 ^ 2 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + m(5) * (t114 ^ 2 + t115 ^ 2 + t116 ^ 2) / 0.2e1 + ((t154 * t290 + t156 * t205 + t158 * t206) * t239 + (t153 * t290 + t155 * t205 + t157 * t206) * t238 + (t186 * t290 + t187 * t205 + t188 * t206) * t245 + t254 * t257 - t256 * t267) * t238 / 0.2e1 + ((t154 * t289 + t207 * t156 + t208 * t158) * t239 + (t153 * t289 + t207 * t155 + t208 * t157) * t238 + (t186 * t289 + t207 * t187 + t208 * t188) * t245 + t254 * t267 + t256 * t257) * t239 / 0.2e1 + ((-t254 * t223 + t226 * t256 + Icges(1,4)) * V_base(5) + (-t254 * t224 + t227 * t256 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t223 * t256 + t254 * t226 + Icges(1,2)) * V_base(5) + (t224 * t256 + t254 * t227 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t191 * t304 + t192 * t303 + t290 * t305) * t229 + (t191 * t308 + t192 * t310 + t290 * t306) * t210 + (t309 * t191 + t311 * t192 + t307 * t290) * t209) * t209 / 0.2e1 + ((t193 * t304 + t194 * t303 + t289 * t305) * t229 + (t308 * t193 + t310 * t194 + t306 * t289) * t210 + (t309 * t193 + t194 * t311 + t307 * t289) * t209) * t210 / 0.2e1 + ((-t209 * t307 - t210 * t306 - t229 * t305) * t255 + ((t243 * t304 + t244 * t303) * t229 + (t243 * t308 + t244 * t310) * t210 + (t309 * t243 + t244 * t311) * t209) * t253) * t229 / 0.2e1 + ((-t153 * t238 - t154 * t239) * t255 + ((-t156 * t250 + t158 * t251) * t239 + (-t155 * t250 + t157 * t251) * t238) * t253 + (t198 * t255 + t200 * t253) * t239 + (t197 * t255 + t199 * t253) * t238 + (Icges(2,3) + (-t186 + t222) * t255 + (-t187 * t250 + t188 * t251 + t225) * t253) * t245) * t245 / 0.2e1 + t245 * V_base(4) * (Icges(2,5) * t256 - Icges(2,6) * t254) + t245 * V_base(5) * (Icges(2,5) * t254 + Icges(2,6) * t256) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
