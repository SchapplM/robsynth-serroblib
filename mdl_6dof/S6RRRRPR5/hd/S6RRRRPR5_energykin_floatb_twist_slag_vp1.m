% Calculate kinetic energy for
% S6RRRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-09 22:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR5_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRPR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR5_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR5_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR5_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:11:56
% EndTime: 2019-03-09 22:11:59
% DurationCPUTime: 3.29s
% Computational Cost: add. (2129->343), mult. (2680->502), div. (0->0), fcn. (2738->10), ass. (0->166)
t322 = Icges(5,1) + Icges(6,1);
t321 = -Icges(5,4) + Icges(6,5);
t320 = Icges(6,4) + Icges(5,5);
t319 = Icges(5,2) + Icges(6,3);
t318 = -Icges(6,6) + Icges(5,6);
t317 = -Icges(5,3) - Icges(6,2);
t250 = qJ(2) + qJ(3);
t247 = cos(t250);
t256 = cos(qJ(4));
t258 = cos(qJ(1));
t290 = t256 * t258;
t252 = sin(qJ(4));
t254 = sin(qJ(1));
t293 = t252 * t254;
t203 = t247 * t293 + t290;
t291 = t254 * t256;
t292 = t252 * t258;
t204 = t247 * t291 - t292;
t246 = sin(t250);
t295 = t246 * t254;
t316 = t319 * t203 + t321 * t204 - t318 * t295;
t205 = t247 * t292 - t291;
t206 = t247 * t290 + t293;
t294 = t246 * t258;
t315 = t319 * t205 + t321 * t206 - t318 * t294;
t314 = -t318 * t203 + t320 * t204 - t317 * t295;
t313 = -t318 * t205 + t320 * t206 - t317 * t294;
t312 = t321 * t203 + t322 * t204 + t320 * t295;
t311 = t321 * t205 + t322 * t206 + t320 * t294;
t310 = t318 * t247 + (t319 * t252 + t321 * t256) * t246;
t309 = t317 * t247 + (-t318 * t252 + t320 * t256) * t246;
t308 = -t320 * t247 + (t321 * t252 + t322 * t256) * t246;
t253 = sin(qJ(2));
t303 = pkin(2) * t253;
t257 = cos(qJ(2));
t302 = pkin(2) * t257;
t300 = Icges(2,4) * t254;
t299 = Icges(3,4) * t253;
t298 = Icges(3,4) * t257;
t297 = Icges(4,4) * t246;
t296 = Icges(4,4) * t247;
t174 = -pkin(8) * t258 + t254 * t302;
t237 = t254 * pkin(1) - t258 * pkin(7);
t289 = -t174 - t237;
t288 = qJD(4) * t246;
t287 = qJD(6) * t246;
t286 = V_base(5) * pkin(6) + V_base(1);
t240 = qJD(2) * t254 + V_base(4);
t243 = V_base(6) + qJD(1);
t239 = -qJD(2) * t258 + V_base(5);
t283 = t239 * t303 + t286;
t215 = qJD(3) * t254 + t240;
t282 = pkin(3) * t247 + pkin(9) * t246;
t281 = rSges(3,1) * t257 - rSges(3,2) * t253;
t280 = rSges(4,1) * t247 - rSges(4,2) * t246;
t186 = t258 * t288 + t215;
t279 = Icges(3,1) * t257 - t299;
t278 = Icges(4,1) * t247 - t297;
t277 = -Icges(3,2) * t253 + t298;
t276 = -Icges(4,2) * t246 + t296;
t275 = Icges(3,5) * t257 - Icges(3,6) * t253;
t274 = Icges(4,5) * t247 - Icges(4,6) * t246;
t238 = t258 * pkin(1) + t254 * pkin(7);
t273 = -V_base(4) * pkin(6) + t243 * t238 + V_base(2);
t272 = V_base(4) * t237 - t238 * V_base(5) + V_base(3);
t214 = V_base(5) + (-qJD(2) - qJD(3)) * t258;
t185 = t254 * t288 + t214;
t271 = (-Icges(4,3) * t258 + t254 * t274) * t214 + (Icges(4,3) * t254 + t258 * t274) * t215 + (Icges(4,5) * t246 + Icges(4,6) * t247) * t243;
t270 = (-Icges(3,3) * t258 + t254 * t275) * t239 + (Icges(3,3) * t254 + t258 * t275) * t240 + (Icges(3,5) * t253 + Icges(3,6) * t257) * t243;
t200 = t282 * t254;
t213 = pkin(3) * t246 - pkin(9) * t247;
t269 = t214 * t213 + (-t200 + t289) * t243 + t283;
t175 = pkin(8) * t254 + t258 * t302;
t268 = t240 * t174 - t175 * t239 + t272;
t267 = t243 * t175 - t240 * t303 + t273;
t199 = (pkin(4) * t256 + qJ(5) * t252) * t246;
t266 = qJD(5) * t205 + t185 * t199 + t269;
t201 = t282 * t258;
t265 = t215 * t200 - t201 * t214 + t268;
t264 = t243 * t201 - t213 * t215 + t267;
t153 = pkin(4) * t204 + qJ(5) * t203;
t263 = qJD(5) * t246 * t252 + t186 * t153 + t265;
t154 = pkin(4) * t206 + qJ(5) * t205;
t221 = -qJD(4) * t247 + t243;
t262 = qJD(5) * t203 + t221 * t154 + t264;
t179 = -Icges(4,6) * t258 + t254 * t276;
t180 = Icges(4,6) * t254 + t258 * t276;
t181 = -Icges(4,5) * t258 + t254 * t278;
t182 = Icges(4,5) * t254 + t258 * t278;
t210 = Icges(4,2) * t247 + t297;
t211 = Icges(4,1) * t246 + t296;
t261 = (-t180 * t246 + t182 * t247) * t215 + (-t179 * t246 + t181 * t247) * t214 + (-t210 * t246 + t211 * t247) * t243;
t191 = -Icges(3,6) * t258 + t254 * t277;
t192 = Icges(3,6) * t254 + t258 * t277;
t193 = -Icges(3,5) * t258 + t254 * t279;
t194 = Icges(3,5) * t254 + t258 * t279;
t225 = Icges(3,2) * t257 + t299;
t228 = Icges(3,1) * t253 + t298;
t260 = (-t192 * t253 + t194 * t257) * t240 + (-t191 * t253 + t193 * t257) * t239 + (-t225 * t253 + t228 * t257) * t243;
t255 = cos(qJ(6));
t251 = sin(qJ(6));
t248 = Icges(2,4) * t258;
t233 = rSges(2,1) * t258 - rSges(2,2) * t254;
t232 = rSges(2,1) * t254 + rSges(2,2) * t258;
t231 = rSges(3,1) * t253 + rSges(3,2) * t257;
t230 = Icges(2,1) * t258 - t300;
t229 = Icges(2,1) * t254 + t248;
t227 = -Icges(2,2) * t254 + t248;
t226 = Icges(2,2) * t258 + t300;
t220 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t219 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t218 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t212 = rSges(4,1) * t246 + rSges(4,2) * t247;
t207 = pkin(5) * t246 * t256 + pkin(10) * t247;
t202 = (-qJD(4) + qJD(6)) * t247 + t243;
t198 = rSges(3,3) * t254 + t258 * t281;
t197 = -rSges(3,3) * t258 + t254 * t281;
t188 = (t251 * t252 + t255 * t256) * t246;
t187 = (-t251 * t256 + t252 * t255) * t246;
t184 = rSges(4,3) * t254 + t258 * t280;
t183 = -rSges(4,3) * t258 + t254 * t280;
t173 = -rSges(5,3) * t247 + (rSges(5,1) * t256 - rSges(5,2) * t252) * t246;
t172 = -rSges(6,2) * t247 + (rSges(6,1) * t256 + rSges(6,3) * t252) * t246;
t165 = V_base(5) * rSges(2,3) - t232 * t243 + t286;
t164 = t233 * t243 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t162 = pkin(5) * t206 - pkin(10) * t294;
t161 = pkin(5) * t204 - pkin(10) * t295;
t160 = t232 * V_base(4) - t233 * V_base(5) + V_base(3);
t158 = -t258 * t287 + t186;
t157 = -t254 * t287 + t185;
t151 = t205 * t251 + t206 * t255;
t150 = t205 * t255 - t206 * t251;
t149 = t203 * t251 + t204 * t255;
t148 = t203 * t255 - t204 * t251;
t147 = rSges(5,1) * t206 - rSges(5,2) * t205 + rSges(5,3) * t294;
t146 = rSges(6,1) * t206 + rSges(6,2) * t294 + rSges(6,3) * t205;
t145 = rSges(5,1) * t204 - rSges(5,2) * t203 + rSges(5,3) * t295;
t144 = rSges(6,1) * t204 + rSges(6,2) * t295 + rSges(6,3) * t203;
t130 = rSges(7,1) * t188 + rSges(7,2) * t187 + rSges(7,3) * t247;
t129 = Icges(7,1) * t188 + Icges(7,4) * t187 + Icges(7,5) * t247;
t128 = Icges(7,4) * t188 + Icges(7,2) * t187 + Icges(7,6) * t247;
t127 = Icges(7,5) * t188 + Icges(7,6) * t187 + Icges(7,3) * t247;
t125 = t231 * t239 + (-t197 - t237) * t243 + t286;
t124 = t198 * t243 - t231 * t240 + t273;
t123 = t197 * t240 - t198 * t239 + t272;
t122 = rSges(7,1) * t151 + rSges(7,2) * t150 - rSges(7,3) * t294;
t121 = rSges(7,1) * t149 + rSges(7,2) * t148 - rSges(7,3) * t295;
t120 = Icges(7,1) * t151 + Icges(7,4) * t150 - Icges(7,5) * t294;
t119 = Icges(7,1) * t149 + Icges(7,4) * t148 - Icges(7,5) * t295;
t118 = Icges(7,4) * t151 + Icges(7,2) * t150 - Icges(7,6) * t294;
t117 = Icges(7,4) * t149 + Icges(7,2) * t148 - Icges(7,6) * t295;
t116 = Icges(7,5) * t151 + Icges(7,6) * t150 - Icges(7,3) * t294;
t115 = Icges(7,5) * t149 + Icges(7,6) * t148 - Icges(7,3) * t295;
t114 = t212 * t214 + (-t183 + t289) * t243 + t283;
t113 = t184 * t243 - t212 * t215 + t267;
t112 = t183 * t215 - t184 * t214 + t268;
t111 = -t145 * t221 + t173 * t185 + t269;
t110 = t147 * t221 - t173 * t186 + t264;
t109 = t145 * t186 - t147 * t185 + t265;
t108 = t172 * t185 + (-t144 - t153) * t221 + t266;
t107 = t146 * t221 + (-t172 - t199) * t186 + t262;
t106 = t144 * t186 + (-t146 - t154) * t185 + t263;
t105 = -t121 * t202 + t130 * t157 + t185 * t207 + (-t153 - t161) * t221 + t266;
t104 = t122 * t202 - t130 * t158 + t162 * t221 + (-t199 - t207) * t186 + t262;
t103 = t121 * t158 - t122 * t157 + t161 * t186 + (-t154 - t162) * t185 + t263;
t1 = t202 * ((t116 * t247 + t118 * t187 + t120 * t188) * t158 + (t115 * t247 + t117 * t187 + t119 * t188) * t157 + (t247 * t127 + t187 * t128 + t188 * t129) * t202) / 0.2e1 + m(1) * (t218 ^ 2 + t219 ^ 2 + t220 ^ 2) / 0.2e1 + t157 * ((-t116 * t295 + t118 * t148 + t120 * t149) * t158 + (-t115 * t295 + t148 * t117 + t149 * t119) * t157 + (-t127 * t295 + t128 * t148 + t129 * t149) * t202) / 0.2e1 + t158 * ((-t116 * t294 + t150 * t118 + t151 * t120) * t158 + (-t115 * t294 + t117 * t150 + t119 * t151) * t157 + (-t127 * t294 + t128 * t150 + t129 * t151) * t202) / 0.2e1 + m(2) * (t160 ^ 2 + t164 ^ 2 + t165 ^ 2) / 0.2e1 + m(3) * (t123 ^ 2 + t124 ^ 2 + t125 ^ 2) / 0.2e1 + m(4) * (t112 ^ 2 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + m(6) * (t106 ^ 2 + t107 ^ 2 + t108 ^ 2) / 0.2e1 + m(5) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + m(7) * (t103 ^ 2 + t104 ^ 2 + t105 ^ 2) / 0.2e1 + t240 * (t254 * t270 + t258 * t260) / 0.2e1 + t239 * (t254 * t260 - t258 * t270) / 0.2e1 + t215 * (t254 * t271 + t258 * t261) / 0.2e1 + t214 * (t254 * t261 - t271 * t258) / 0.2e1 + ((t203 * t310 + t204 * t308 + t295 * t309) * t221 + (t203 * t315 + t204 * t311 + t295 * t313) * t186 + (t316 * t203 + t312 * t204 + t314 * t295) * t185) * t185 / 0.2e1 + ((t205 * t310 + t206 * t308 + t294 * t309) * t221 + (t315 * t205 + t311 * t206 + t313 * t294) * t186 + (t205 * t316 + t312 * t206 + t314 * t294) * t185) * t186 / 0.2e1 + ((-t185 * t314 - t186 * t313 - t221 * t309) * t247 + ((t252 * t310 + t256 * t308) * t221 + (t252 * t315 + t256 * t311) * t186 + (t252 * t316 + t312 * t256) * t185) * t246) * t221 / 0.2e1 + ((-t226 * t254 + t229 * t258 + Icges(1,4)) * V_base(5) + (-t227 * t254 + t230 * t258 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t226 * t258 + t229 * t254 + Icges(1,2)) * V_base(5) + (t227 * t258 + t230 * t254 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t192 * t257 + t194 * t253) * t240 + (t191 * t257 + t193 * t253) * t239 + (t180 * t247 + t182 * t246) * t215 + (t179 * t247 + t181 * t246) * t214 + (t210 * t247 + t211 * t246 + t225 * t257 + t228 * t253 + Icges(2,3)) * t243) * t243 / 0.2e1 + t243 * V_base(4) * (Icges(2,5) * t258 - Icges(2,6) * t254) + t243 * V_base(5) * (Icges(2,5) * t254 + Icges(2,6) * t258) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
