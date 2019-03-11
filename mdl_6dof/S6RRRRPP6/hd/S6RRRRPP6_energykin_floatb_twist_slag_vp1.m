% Calculate kinetic energy for
% S6RRRRPP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 21:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPP6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP6_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP6_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRPP6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP6_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP6_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP6_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPP6_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:11:04
% EndTime: 2019-03-09 21:11:07
% DurationCPUTime: 3.26s
% Computational Cost: add. (1892->309), mult. (2692->450), div. (0->0), fcn. (2668->8), ass. (0->147)
t321 = Icges(5,1) + Icges(6,2) + Icges(7,3);
t320 = -Icges(5,4) - Icges(6,6) + Icges(7,6);
t319 = -Icges(5,5) - Icges(7,5) + Icges(6,4);
t318 = Icges(5,2) + Icges(7,2) + Icges(6,3);
t317 = Icges(5,6) - Icges(6,5) - Icges(7,4);
t316 = -Icges(5,3) - Icges(7,1) - Icges(6,1);
t315 = rSges(7,1) + pkin(5);
t314 = rSges(7,3) + qJ(6);
t253 = qJ(3) + qJ(4);
t249 = sin(t253);
t250 = cos(t253);
t259 = cos(qJ(1));
t256 = sin(qJ(1));
t258 = cos(qJ(2));
t290 = t256 * t258;
t202 = t249 * t290 + t250 * t259;
t288 = t259 * t249;
t203 = t250 * t290 - t288;
t255 = sin(qJ(2));
t292 = t255 * t256;
t313 = t320 * t202 + t203 * t321 - t319 * t292;
t204 = -t256 * t250 + t258 * t288;
t289 = t258 * t259;
t205 = t249 * t256 + t250 * t289;
t291 = t255 * t259;
t312 = t320 * t204 + t205 * t321 - t319 * t291;
t311 = t202 * t318 + t203 * t320 - t292 * t317;
t310 = t204 * t318 + t205 * t320 - t291 * t317;
t309 = -t202 * t317 - t203 * t319 - t292 * t316;
t308 = -t204 * t317 - t205 * t319 - t291 * t316;
t307 = t316 * t258 + (-t249 * t317 - t250 * t319) * t255;
t306 = t317 * t258 + (t249 * t318 + t250 * t320) * t255;
t305 = t319 * t258 + (t320 * t249 + t250 * t321) * t255;
t257 = cos(qJ(3));
t300 = pkin(3) * t257;
t298 = Icges(2,4) * t256;
t297 = Icges(3,4) * t255;
t296 = Icges(3,4) * t258;
t295 = t250 * t255;
t254 = sin(qJ(3));
t294 = t254 * t256;
t293 = t254 * t259;
t287 = rSges(7,2) * t202 + t314 * t203 + t292 * t315;
t286 = rSges(7,2) * t204 + t314 * t205 + t291 * t315;
t285 = (rSges(7,2) * t249 + rSges(7,3) * t250) * t255 + qJ(6) * t295 - t315 * t258;
t284 = qJD(3) * t255;
t283 = qJD(4) * t255;
t282 = V_base(5) * pkin(6) + V_base(1);
t241 = qJD(2) * t256 + V_base(4);
t247 = V_base(6) + qJD(1);
t209 = t259 * t284 + t241;
t279 = pkin(2) * t258 + pkin(8) * t255;
t240 = -qJD(2) * t259 + V_base(5);
t278 = rSges(3,1) * t258 - rSges(3,2) * t255;
t277 = Icges(3,1) * t258 - t297;
t276 = -Icges(3,2) * t255 + t296;
t275 = Icges(3,5) * t258 - Icges(3,6) * t255;
t208 = t256 * t284 + t240;
t239 = pkin(1) * t259 + pkin(7) * t256;
t274 = -V_base(4) * pkin(6) + t247 * t239 + V_base(2);
t238 = pkin(1) * t256 - pkin(7) * t259;
t273 = V_base(4) * t238 - t239 * V_base(5) + V_base(3);
t272 = pkin(9) * t255 + t258 * t300;
t215 = t279 * t256;
t237 = t255 * pkin(2) - t258 * pkin(8);
t271 = t240 * t237 + (-t215 - t238) * t247 + t282;
t270 = (-Icges(3,3) * t259 + t256 * t275) * t240 + (Icges(3,3) * t256 + t259 * t275) * t241 + (Icges(3,5) * t255 + Icges(3,6) * t258) * t247;
t216 = t279 * t259;
t269 = t247 * t216 - t237 * t241 + t274;
t268 = t241 * t215 - t216 * t240 + t273;
t163 = -pkin(3) * t293 + t256 * t272;
t171 = -pkin(9) * t258 + t255 * t300;
t232 = -qJD(3) * t258 + t247;
t267 = -t163 * t232 + t208 * t171 + t271;
t164 = pkin(3) * t294 + t259 * t272;
t266 = t232 * t164 - t171 * t209 + t269;
t185 = t256 * t283 + t208;
t206 = (pkin(4) * t250 + qJ(5) * t249) * t255;
t265 = qJD(5) * t204 + t185 * t206 + t267;
t264 = t209 * t163 - t164 * t208 + t268;
t153 = pkin(4) * t205 + qJ(5) * t204;
t217 = (-qJD(3) - qJD(4)) * t258 + t247;
t263 = qJD(5) * t202 + t217 * t153 + t266;
t152 = pkin(4) * t203 + qJ(5) * t202;
t186 = t259 * t283 + t209;
t262 = qJD(5) * t255 * t249 + t186 * t152 + t264;
t194 = -Icges(3,6) * t259 + t256 * t276;
t195 = Icges(3,6) * t256 + t259 * t276;
t197 = -Icges(3,5) * t259 + t256 * t277;
t198 = Icges(3,5) * t256 + t259 * t277;
t226 = Icges(3,2) * t258 + t297;
t229 = Icges(3,1) * t255 + t296;
t261 = (-t195 * t255 + t198 * t258) * t241 + (-t194 * t255 + t197 * t258) * t240 + (-t226 * t255 + t229 * t258) * t247;
t251 = Icges(2,4) * t259;
t235 = rSges(2,1) * t259 - rSges(2,2) * t256;
t234 = rSges(2,1) * t256 + rSges(2,2) * t259;
t233 = rSges(3,1) * t255 + rSges(3,2) * t258;
t231 = Icges(2,1) * t259 - t298;
t230 = Icges(2,1) * t256 + t251;
t228 = -Icges(2,2) * t256 + t251;
t227 = Icges(2,2) * t259 + t298;
t222 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t221 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t220 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t213 = t257 * t289 + t294;
t212 = -t254 * t289 + t256 * t257;
t211 = t257 * t290 - t293;
t210 = -t254 * t290 - t257 * t259;
t201 = rSges(3,3) * t256 + t259 * t278;
t200 = -rSges(3,3) * t259 + t256 * t278;
t199 = -rSges(4,3) * t258 + (rSges(4,1) * t257 - rSges(4,2) * t254) * t255;
t196 = -Icges(4,5) * t258 + (Icges(4,1) * t257 - Icges(4,4) * t254) * t255;
t193 = -Icges(4,6) * t258 + (Icges(4,4) * t257 - Icges(4,2) * t254) * t255;
t190 = -Icges(4,3) * t258 + (Icges(4,5) * t257 - Icges(4,6) * t254) * t255;
t184 = -rSges(6,1) * t258 + (-rSges(6,2) * t250 + rSges(6,3) * t249) * t255;
t182 = -rSges(5,3) * t258 + (rSges(5,1) * t250 - rSges(5,2) * t249) * t255;
t170 = V_base(5) * rSges(2,3) - t234 * t247 + t282;
t169 = t235 * t247 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t168 = t234 * V_base(4) - t235 * V_base(5) + V_base(3);
t162 = rSges(4,1) * t213 + rSges(4,2) * t212 + rSges(4,3) * t291;
t161 = rSges(4,1) * t211 + rSges(4,2) * t210 + rSges(4,3) * t292;
t159 = Icges(4,1) * t213 + Icges(4,4) * t212 + Icges(4,5) * t291;
t158 = Icges(4,1) * t211 + Icges(4,4) * t210 + Icges(4,5) * t292;
t157 = Icges(4,4) * t213 + Icges(4,2) * t212 + Icges(4,6) * t291;
t156 = Icges(4,4) * t211 + Icges(4,2) * t210 + Icges(4,6) * t292;
t155 = Icges(4,5) * t213 + Icges(4,6) * t212 + Icges(4,3) * t291;
t154 = Icges(4,5) * t211 + Icges(4,6) * t210 + Icges(4,3) * t292;
t151 = rSges(5,1) * t205 - rSges(5,2) * t204 + rSges(5,3) * t291;
t150 = rSges(5,1) * t203 - rSges(5,2) * t202 + rSges(5,3) * t292;
t149 = rSges(6,1) * t291 - rSges(6,2) * t205 + rSges(6,3) * t204;
t147 = rSges(6,1) * t292 - rSges(6,2) * t203 + rSges(6,3) * t202;
t124 = t233 * t240 + (-t200 - t238) * t247 + t282;
t123 = t201 * t247 - t233 * t241 + t274;
t121 = t200 * t241 - t201 * t240 + t273;
t120 = -t161 * t232 + t199 * t208 + t271;
t119 = t162 * t232 - t199 * t209 + t269;
t118 = t161 * t209 - t162 * t208 + t268;
t117 = -t150 * t217 + t182 * t185 + t267;
t116 = t151 * t217 - t182 * t186 + t266;
t115 = t150 * t186 - t151 * t185 + t264;
t114 = t184 * t185 + (-t147 - t152) * t217 + t265;
t113 = t149 * t217 + (-t184 - t206) * t186 + t263;
t112 = t147 * t186 + (-t149 - t153) * t185 + t262;
t111 = qJD(6) * t205 + t285 * t185 + (-t152 - t287) * t217 + t265;
t110 = qJD(6) * t203 + t286 * t217 + (-t206 - t285) * t186 + t263;
t109 = qJD(6) * t295 + t287 * t186 + (-t153 - t286) * t185 + t262;
t1 = t241 * (t270 * t256 + t261 * t259) / 0.2e1 + t232 * ((-t154 * t208 - t155 * t209 - t190 * t232) * t258 + ((-t157 * t254 + t159 * t257) * t209 + (-t156 * t254 + t158 * t257) * t208 + (-t193 * t254 + t196 * t257) * t232) * t255) / 0.2e1 + m(1) * (t220 ^ 2 + t221 ^ 2 + t222 ^ 2) / 0.2e1 + m(2) * (t168 ^ 2 + t169 ^ 2 + t170 ^ 2) / 0.2e1 + m(3) * (t121 ^ 2 + t123 ^ 2 + t124 ^ 2) / 0.2e1 + m(4) * (t118 ^ 2 + t119 ^ 2 + t120 ^ 2) / 0.2e1 + m(7) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + m(6) * (t112 ^ 2 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + m(5) * (t115 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + t240 * (t261 * t256 - t270 * t259) / 0.2e1 + t209 * ((t155 * t291 + t157 * t212 + t159 * t213) * t209 + (t154 * t291 + t156 * t212 + t158 * t213) * t208 + (t190 * t291 + t193 * t212 + t196 * t213) * t232) / 0.2e1 + t208 * ((t155 * t292 + t157 * t210 + t159 * t211) * t209 + (t154 * t292 + t156 * t210 + t158 * t211) * t208 + (t190 * t292 + t193 * t210 + t196 * t211) * t232) / 0.2e1 + ((t195 * t258 + t198 * t255) * t241 + (t194 * t258 + t197 * t255) * t240 + (t226 * t258 + t229 * t255 + Icges(2,3)) * t247) * t247 / 0.2e1 + ((-t227 * t256 + t230 * t259 + Icges(1,4)) * V_base(5) + (-t228 * t256 + t231 * t259 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t227 * t259 + t230 * t256 + Icges(1,2)) * V_base(5) + (t228 * t259 + t231 * t256 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t202 * t306 + t203 * t305 + t292 * t307) * t217 + (t202 * t310 + t203 * t312 + t292 * t308) * t186 + (t311 * t202 + t313 * t203 + t309 * t292) * t185) * t185 / 0.2e1 + ((t204 * t306 + t205 * t305 + t291 * t307) * t217 + (t310 * t204 + t312 * t205 + t308 * t291) * t186 + (t311 * t204 + t205 * t313 + t309 * t291) * t185) * t186 / 0.2e1 + ((-t185 * t309 - t186 * t308 - t217 * t307) * t258 + ((t249 * t306 + t250 * t305) * t217 + (t249 * t310 + t250 * t312) * t186 + (t311 * t249 + t250 * t313) * t185) * t255) * t217 / 0.2e1 + V_base(4) * t247 * (Icges(2,5) * t259 - Icges(2,6) * t256) + V_base(5) * t247 * (Icges(2,5) * t256 + Icges(2,6) * t259) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
