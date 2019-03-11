% Calculate kinetic energy for
% S6RPRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
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
% Datum: 2019-03-09 02:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPPR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR4_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPPR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR4_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR4_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPPR4_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:46:41
% EndTime: 2019-03-09 02:46:44
% DurationCPUTime: 3.21s
% Computational Cost: add. (1913->348), mult. (2485->477), div. (0->0), fcn. (2543->10), ass. (0->168)
t323 = Icges(5,1) + Icges(6,1);
t322 = -Icges(5,4) + Icges(6,5);
t321 = Icges(6,4) + Icges(5,5);
t320 = Icges(5,2) + Icges(6,3);
t319 = -Icges(6,6) + Icges(5,6);
t318 = -Icges(5,3) - Icges(6,2);
t247 = pkin(9) + qJ(3);
t241 = cos(t247);
t248 = sin(pkin(10));
t254 = sin(qJ(1));
t294 = t254 * t248;
t250 = cos(pkin(10));
t256 = cos(qJ(1));
t295 = t250 * t256;
t200 = t241 * t294 + t295;
t293 = t254 * t250;
t296 = t248 * t256;
t201 = t241 * t293 - t296;
t240 = sin(t247);
t298 = t240 * t254;
t317 = t320 * t200 + t322 * t201 - t319 * t298;
t202 = t241 * t296 - t293;
t203 = t241 * t295 + t294;
t297 = t240 * t256;
t316 = t320 * t202 + t322 * t203 - t319 * t297;
t315 = -t319 * t200 + t321 * t201 - t318 * t298;
t314 = -t319 * t202 + t321 * t203 - t318 * t297;
t313 = t322 * t200 + t323 * t201 + t321 * t298;
t312 = t322 * t202 + t323 * t203 + t321 * t297;
t311 = t319 * t241 + (t320 * t248 + t322 * t250) * t240;
t310 = t318 * t241 + (-t319 * t248 + t321 * t250) * t240;
t309 = -t321 * t241 + (t322 * t248 + t323 * t250) * t240;
t249 = sin(pkin(9));
t305 = pkin(2) * t249;
t251 = cos(pkin(9));
t304 = pkin(2) * t251;
t303 = Icges(2,4) * t254;
t302 = Icges(3,4) * t249;
t301 = Icges(3,4) * t251;
t300 = Icges(4,4) * t240;
t299 = Icges(4,4) * t241;
t155 = pkin(4) * t203 + qJ(5) * t202;
t276 = pkin(3) * t241 + qJ(4) * t240;
t199 = t276 * t256;
t291 = -t155 - t199;
t173 = -pkin(7) * t256 + t254 * t304;
t231 = t254 * pkin(1) - qJ(2) * t256;
t290 = -t173 - t231;
t195 = (pkin(4) * t250 + qJ(5) * t248) * t240;
t211 = pkin(3) * t240 - qJ(4) * t241;
t289 = -t195 - t211;
t288 = qJD(4) * t240;
t287 = qJD(6) * t240;
t286 = V_base(4) * t231 + V_base(3);
t285 = V_base(5) * pkin(6) + V_base(1);
t198 = t276 * t254;
t282 = -t198 + t290;
t236 = qJD(3) * t254 + V_base(4);
t242 = V_base(6) + qJD(1);
t281 = qJD(2) * t254 + t285;
t154 = pkin(4) * t201 + qJ(5) * t200;
t280 = -t154 + t282;
t279 = V_base(5) * t305 + t281;
t235 = -qJD(3) * t256 + V_base(5);
t278 = rSges(3,1) * t251 - rSges(3,2) * t249;
t277 = rSges(4,1) * t241 - rSges(4,2) * t240;
t275 = Icges(3,1) * t251 - t302;
t274 = Icges(4,1) * t241 - t300;
t273 = -Icges(3,2) * t249 + t301;
t272 = -Icges(4,2) * t240 + t299;
t271 = Icges(3,5) * t251 - Icges(3,6) * t249;
t270 = Icges(4,5) * t241 - Icges(4,6) * t240;
t233 = pkin(1) * t256 + t254 * qJ(2);
t269 = -qJD(2) * t256 + t242 * t233 + V_base(2);
t268 = t235 * t211 + t256 * t288 + t279;
t267 = (-Icges(4,3) * t256 + t254 * t270) * t235 + (Icges(4,3) * t254 + t256 * t270) * t236 + (Icges(4,5) * t240 + Icges(4,6) * t241) * t242;
t266 = qJD(5) * t202 + t235 * t195 + t268;
t174 = pkin(7) * t254 + t256 * t304;
t265 = V_base(4) * t173 + (-t174 - t233) * V_base(5) + t286;
t264 = (-Icges(3,3) * t256 + t254 * t271) * V_base(5) + (Icges(3,3) * t254 + t256 * t271) * V_base(4) + (Icges(3,5) * t249 + Icges(3,6) * t251) * t242;
t263 = t242 * t174 + (-pkin(6) - t305) * V_base(4) + t269;
t262 = -qJD(4) * t241 + t236 * t198 + t265;
t261 = t242 * t199 + t254 * t288 + t263;
t260 = qJD(5) * t240 * t248 + t236 * t154 + t262;
t259 = qJD(5) * t200 + t242 * t155 + t261;
t178 = -Icges(4,6) * t256 + t254 * t272;
t179 = Icges(4,6) * t254 + t256 * t272;
t180 = -Icges(4,5) * t256 + t254 * t274;
t181 = Icges(4,5) * t254 + t256 * t274;
t209 = Icges(4,2) * t241 + t300;
t210 = Icges(4,1) * t240 + t299;
t258 = (-t179 * t240 + t181 * t241) * t236 + (-t178 * t240 + t180 * t241) * t235 + (-t209 * t240 + t210 * t241) * t242;
t189 = -Icges(3,6) * t256 + t254 * t273;
t190 = Icges(3,6) * t254 + t256 * t273;
t191 = -Icges(3,5) * t256 + t254 * t275;
t192 = Icges(3,5) * t254 + t256 * t275;
t219 = Icges(3,2) * t251 + t302;
t220 = Icges(3,1) * t249 + t301;
t257 = (-t190 * t249 + t192 * t251) * V_base(4) + (-t189 * t249 + t191 * t251) * V_base(5) + (-t219 * t249 + t220 * t251) * t242;
t255 = cos(qJ(6));
t253 = sin(qJ(6));
t245 = Icges(2,4) * t256;
t234 = rSges(2,1) * t256 - t254 * rSges(2,2);
t232 = t254 * rSges(2,1) + rSges(2,2) * t256;
t228 = Icges(2,1) * t256 - t303;
t227 = Icges(2,1) * t254 + t245;
t226 = -Icges(2,2) * t254 + t245;
t225 = Icges(2,2) * t256 + t303;
t224 = Icges(2,5) * t256 - Icges(2,6) * t254;
t223 = Icges(2,5) * t254 + Icges(2,6) * t256;
t221 = rSges(3,1) * t249 + rSges(3,2) * t251;
t217 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t216 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t215 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t214 = qJD(6) * t241 + t242;
t212 = rSges(4,1) * t240 + rSges(4,2) * t241;
t206 = pkin(5) * t240 * t250 + pkin(8) * t241;
t205 = -t256 * t287 + t236;
t204 = -t254 * t287 + t235;
t197 = t254 * rSges(3,3) + t256 * t278;
t196 = -rSges(3,3) * t256 + t254 * t278;
t185 = (t248 * t253 + t250 * t255) * t240;
t184 = (t248 * t255 - t250 * t253) * t240;
t183 = t254 * rSges(4,3) + t256 * t277;
t182 = -rSges(4,3) * t256 + t254 * t277;
t172 = V_base(5) * rSges(2,3) - t232 * t242 + t285;
t171 = t234 * t242 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t170 = -rSges(5,3) * t241 + (rSges(5,1) * t250 - rSges(5,2) * t248) * t240;
t169 = -rSges(6,2) * t241 + (rSges(6,1) * t250 + rSges(6,3) * t248) * t240;
t161 = t232 * V_base(4) - t234 * V_base(5) + V_base(3);
t158 = t203 * pkin(5) - pkin(8) * t297;
t157 = pkin(5) * t201 - pkin(8) * t298;
t153 = t202 * t253 + t203 * t255;
t152 = t202 * t255 - t203 * t253;
t151 = t200 * t253 + t201 * t255;
t150 = t200 * t255 - t201 * t253;
t147 = t203 * rSges(5,1) - t202 * rSges(5,2) + rSges(5,3) * t297;
t146 = t203 * rSges(6,1) + rSges(6,2) * t297 + t202 * rSges(6,3);
t145 = rSges(5,1) * t201 - rSges(5,2) * t200 + rSges(5,3) * t298;
t144 = rSges(6,1) * t201 + rSges(6,2) * t298 + rSges(6,3) * t200;
t131 = rSges(7,1) * t185 + rSges(7,2) * t184 + rSges(7,3) * t241;
t130 = Icges(7,1) * t185 + Icges(7,4) * t184 + Icges(7,5) * t241;
t129 = Icges(7,4) * t185 + Icges(7,2) * t184 + Icges(7,6) * t241;
t128 = Icges(7,5) * t185 + Icges(7,6) * t184 + Icges(7,3) * t241;
t127 = t221 * V_base(5) + (-t196 - t231) * t242 + t281;
t126 = t242 * t197 + (-pkin(6) - t221) * V_base(4) + t269;
t125 = t196 * V_base(4) + (-t197 - t233) * V_base(5) + t286;
t124 = t153 * rSges(7,1) + t152 * rSges(7,2) - rSges(7,3) * t297;
t123 = rSges(7,1) * t151 + rSges(7,2) * t150 - rSges(7,3) * t298;
t122 = Icges(7,1) * t153 + Icges(7,4) * t152 - Icges(7,5) * t297;
t121 = Icges(7,1) * t151 + Icges(7,4) * t150 - Icges(7,5) * t298;
t120 = Icges(7,4) * t153 + Icges(7,2) * t152 - Icges(7,6) * t297;
t119 = Icges(7,4) * t151 + Icges(7,2) * t150 - Icges(7,6) * t298;
t118 = Icges(7,5) * t153 + Icges(7,6) * t152 - Icges(7,3) * t297;
t117 = Icges(7,5) * t151 + Icges(7,6) * t150 - Icges(7,3) * t298;
t116 = t212 * t235 + (-t182 + t290) * t242 + t279;
t115 = t242 * t183 - t236 * t212 + t263;
t114 = t182 * t236 - t183 * t235 + t265;
t113 = t170 * t235 + (-t145 + t282) * t242 + t268;
t112 = t242 * t147 + (-t170 - t211) * t236 + t261;
t111 = t145 * t236 + (-t147 - t199) * t235 + t262;
t110 = t169 * t235 + (-t144 + t280) * t242 + t266;
t109 = t242 * t146 + (-t169 + t289) * t236 + t259;
t108 = t144 * t236 + (-t146 + t291) * t235 + t260;
t107 = -t123 * t214 + t131 * t204 + t206 * t235 + (-t157 + t280) * t242 + t266;
t106 = t214 * t124 - t205 * t131 + t242 * t158 + (-t206 + t289) * t236 + t259;
t105 = t123 * t205 - t124 * t204 + t157 * t236 + (-t158 + t291) * t235 + t260;
t1 = t214 * ((t118 * t241 + t120 * t184 + t122 * t185) * t205 + (t117 * t241 + t119 * t184 + t121 * t185) * t204 + (t241 * t128 + t184 * t129 + t185 * t130) * t214) / 0.2e1 + m(1) * (t215 ^ 2 + t216 ^ 2 + t217 ^ 2) / 0.2e1 + t205 * ((-t118 * t297 + t152 * t120 + t153 * t122) * t205 + (-t117 * t297 + t152 * t119 + t153 * t121) * t204 + (-t128 * t297 + t152 * t129 + t153 * t130) * t214) / 0.2e1 + t204 * ((-t118 * t298 + t120 * t150 + t122 * t151) * t205 + (-t117 * t298 + t150 * t119 + t151 * t121) * t204 + (-t128 * t298 + t129 * t150 + t130 * t151) * t214) / 0.2e1 + m(2) * (t161 ^ 2 + t171 ^ 2 + t172 ^ 2) / 0.2e1 + m(3) * (t125 ^ 2 + t126 ^ 2 + t127 ^ 2) / 0.2e1 + m(6) * (t108 ^ 2 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + m(5) * (t111 ^ 2 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + m(4) * (t114 ^ 2 + t115 ^ 2 + t116 ^ 2) / 0.2e1 + m(7) * (t105 ^ 2 + t106 ^ 2 + t107 ^ 2) / 0.2e1 + (t254 * t258 - t256 * t267 + (t200 * t311 + t201 * t309 + t298 * t310) * t242 + (t200 * t316 + t201 * t312 + t298 * t314) * t236 + (t317 * t200 + t313 * t201 + t315 * t298) * t235) * t235 / 0.2e1 + (t254 * t267 + t256 * t258 + (t202 * t311 + t203 * t309 + t297 * t310) * t242 + (t316 * t202 + t312 * t203 + t314 * t297) * t236 + (t202 * t317 + t313 * t203 + t315 * t297) * t235) * t236 / 0.2e1 + (t224 * t242 + t254 * t264 + t256 * t257 + (-t254 * t225 + t227 * t256 + Icges(1,4)) * V_base(5) + (-t254 * t226 + t228 * t256 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t223 * t242 + t254 * t257 - t264 * t256 + (t225 * t256 + t254 * t227 + Icges(1,2)) * V_base(5) + (t226 * t256 + t254 * t228 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t189 * t251 + t191 * t249 + t223) * V_base(5) + (t190 * t251 + t192 * t249 + t224) * V_base(4) + ((t179 - t314) * t236 + (t178 - t315) * t235) * t241 + ((t248 * t316 + t250 * t312 + t181) * t236 + (t248 * t317 + t313 * t250 + t180) * t235) * t240 + (t219 * t251 + t220 * t249 + Icges(2,3) + (t209 - t310) * t241 + (t248 * t311 + t250 * t309 + t210) * t240) * t242) * t242 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
