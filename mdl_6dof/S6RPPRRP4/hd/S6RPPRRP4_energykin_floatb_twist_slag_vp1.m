% Calculate kinetic energy for
% S6RPPRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2019-03-09 02:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRP4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP4_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPPRRP4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP4_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP4_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRP4_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:05:05
% EndTime: 2019-03-09 02:05:08
% DurationCPUTime: 2.83s
% Computational Cost: add. (1383->268), mult. (2610->358), div. (0->0), fcn. (2998->8), ass. (0->135)
t308 = Icges(2,4) - Icges(3,5);
t307 = Icges(2,1) + Icges(3,1);
t306 = Icges(6,1) + Icges(7,1);
t305 = Icges(3,4) + Icges(2,5);
t304 = Icges(6,4) - Icges(7,5);
t303 = Icges(7,4) + Icges(6,5);
t302 = Icges(2,2) + Icges(3,3);
t301 = Icges(6,2) + Icges(7,3);
t300 = Icges(2,6) - Icges(3,6);
t299 = Icges(7,6) - Icges(6,6);
t298 = Icges(6,3) + Icges(7,2);
t297 = rSges(7,1) + pkin(5);
t296 = rSges(7,3) + qJ(6);
t273 = sin(qJ(1));
t295 = t308 * t273;
t274 = cos(qJ(1));
t294 = t308 * t274;
t271 = sin(pkin(9));
t272 = cos(pkin(9));
t187 = -t271 * t273 - t272 * t274;
t188 = t271 * t274 - t272 * t273;
t230 = cos(qJ(5));
t228 = sin(qJ(5));
t231 = cos(qJ(4));
t265 = t228 * t231;
t148 = t187 * t230 - t188 * t265;
t264 = t230 * t231;
t149 = -t187 * t228 - t188 * t264;
t229 = sin(qJ(4));
t266 = t188 * t229;
t293 = t301 * t148 - t304 * t149 - t299 * t266;
t150 = -t187 * t265 - t188 * t230;
t151 = -t187 * t264 + t188 * t228;
t267 = t187 * t229;
t292 = t301 * t150 - t304 * t151 - t299 * t267;
t291 = t299 * t148 + t303 * t149 - t298 * t266;
t290 = t299 * t150 + t303 * t151 - t298 * t267;
t289 = -t304 * t148 + t306 * t149 - t303 * t266;
t288 = -t304 * t150 + t306 * t151 - t303 * t267;
t287 = t299 * t231 + (-t301 * t228 + t304 * t230) * t229;
t286 = t298 * t231 + (-t299 * t228 - t303 * t230) * t229;
t285 = t303 * t231 + (t304 * t228 - t306 * t230) * t229;
t284 = -t302 * t274 - t295;
t283 = t302 * t273 - t294;
t282 = t307 * t273 + t294;
t281 = t307 * t274 - t295;
t270 = Icges(4,4) * t187;
t269 = Icges(5,4) * t229;
t268 = Icges(5,4) * t231;
t263 = -rSges(7,2) * t266 + t296 * t148 + t297 * t149;
t262 = -rSges(7,2) * t267 + t296 * t150 + t297 * t151;
t261 = rSges(7,2) * t231 + (-t296 * t228 - t297 * t230) * t229;
t260 = qJD(5) * t229;
t210 = pkin(1) * t273 - qJ(2) * t274;
t259 = V_base(4) * t210 + V_base(3);
t258 = V_base(5) * pkin(6) + V_base(1);
t255 = t274 * pkin(2);
t254 = t273 * pkin(2);
t179 = qJD(4) * t188 + V_base(4);
t178 = -qJD(4) * t187 + V_base(5);
t222 = V_base(6) + qJD(1);
t251 = qJD(2) * t273 + t258;
t250 = -t210 - t254;
t213 = pkin(1) * t274 + qJ(2) * t273;
t249 = -t213 - t255;
t248 = -pkin(4) * t231 - pkin(8) * t229;
t247 = V_base(4) * t254 - qJD(3) + t259;
t246 = -rSges(5,1) * t231 + rSges(5,2) * t229;
t245 = -Icges(5,1) * t231 + t269;
t244 = Icges(5,2) * t229 - t268;
t243 = -Icges(5,5) * t231 + Icges(5,6) * t229;
t164 = -pkin(3) * t188 - pkin(7) * t187;
t242 = -t164 + t250;
t241 = -qJD(2) * t274 + t222 * t213 + V_base(2);
t240 = -V_base(5) * qJ(3) + t251;
t239 = (-Icges(5,3) * t187 + t188 * t243) * t178 + (Icges(5,3) * t188 + t187 * t243) * t179 + (-Icges(5,5) * t229 - Icges(5,6) * t231) * t222;
t238 = V_base(4) * qJ(3) + t222 * t255 + t241;
t165 = -pkin(3) * t187 + pkin(7) * t188;
t237 = -V_base(4) * pkin(6) + t222 * t165 + t238;
t236 = V_base(4) * t164 + (-t165 + t249) * V_base(5) + t247;
t154 = t248 * t187;
t216 = -t229 * pkin(4) + pkin(8) * t231;
t235 = t222 * t154 - t179 * t216 + t237;
t153 = t248 * t188;
t234 = t178 * t216 + (-t153 + t242) * t222 + t240;
t233 = t179 * t153 - t178 * t154 + t236;
t137 = -Icges(5,6) * t187 + t188 * t244;
t138 = Icges(5,6) * t188 + t187 * t244;
t139 = -Icges(5,5) * t187 + t188 * t245;
t140 = Icges(5,5) * t188 + t187 * t245;
t198 = -Icges(5,2) * t231 - t269;
t203 = -Icges(5,1) * t229 - t268;
t232 = (t138 * t229 - t140 * t231) * t179 + (t137 * t229 - t139 * t231) * t178 + (t198 * t229 - t203 * t231) * t222;
t215 = rSges(2,1) * t274 - rSges(2,2) * t273;
t214 = rSges(3,1) * t274 + rSges(3,3) * t273;
t212 = rSges(2,1) * t273 + rSges(2,2) * t274;
t211 = rSges(3,1) * t273 - rSges(3,3) * t274;
t209 = -t229 * rSges(5,1) - rSges(5,2) * t231;
t208 = qJD(5) * t231 + t222;
t192 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t191 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t190 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t184 = Icges(4,4) * t188;
t177 = rSges(6,3) * t231 + (-rSges(6,1) * t230 + rSges(6,2) * t228) * t229;
t169 = V_base(5) * rSges(2,3) - t212 * t222 + t258;
t168 = t215 * t222 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t167 = t212 * V_base(4) - t215 * V_base(5) + V_base(3);
t163 = -rSges(4,1) * t187 - rSges(4,2) * t188;
t162 = -rSges(4,1) * t188 + rSges(4,2) * t187;
t161 = -Icges(4,1) * t187 - t184;
t160 = -Icges(4,1) * t188 + t270;
t159 = -Icges(4,2) * t188 - t270;
t158 = Icges(4,2) * t187 - t184;
t147 = -t187 * t260 + t179;
t146 = -t188 * t260 + t178;
t144 = V_base(5) * rSges(3,2) + (-t210 - t211) * t222 + t251;
t143 = t222 * t214 + (-rSges(3,2) - pkin(6)) * V_base(4) + t241;
t142 = t188 * rSges(5,3) + t187 * t246;
t141 = -t187 * rSges(5,3) + t188 * t246;
t134 = t211 * V_base(4) + (-t213 - t214) * V_base(5) + t259;
t130 = (-qJ(3) - rSges(4,3)) * V_base(5) + (-t162 + t250) * t222 + t251;
t129 = t222 * t163 + (rSges(4,3) - pkin(6)) * V_base(4) + t238;
t128 = rSges(6,1) * t151 - rSges(6,2) * t150 - rSges(6,3) * t267;
t126 = rSges(6,1) * t149 - rSges(6,2) * t148 - rSges(6,3) * t266;
t112 = V_base(4) * t162 + (-t163 + t249) * V_base(5) + t247;
t111 = t178 * t209 + (-t141 + t242) * t222 + t240;
t110 = t222 * t142 - t179 * t209 + t237;
t109 = t179 * t141 - t178 * t142 + t236;
t108 = -t208 * t126 + t146 * t177 + t234;
t107 = t208 * t128 - t147 * t177 + t235;
t106 = qJD(6) * t150 + t146 * t261 - t208 * t263 + t234;
t105 = qJD(6) * t148 - t147 * t261 + t208 * t262 + t235;
t104 = t147 * t126 - t146 * t128 + t233;
t103 = -qJD(6) * t229 * t228 - t146 * t262 + t147 * t263 + t233;
t1 = m(1) * (t190 ^ 2 + t191 ^ 2 + t192 ^ 2) / 0.2e1 + m(2) * (t167 ^ 2 + t168 ^ 2 + t169 ^ 2) / 0.2e1 + m(3) * (t134 ^ 2 + t143 ^ 2 + t144 ^ 2) / 0.2e1 + m(4) * (t112 ^ 2 + t129 ^ 2 + t130 ^ 2) / 0.2e1 + m(7) * (t103 ^ 2 + t105 ^ 2 + t106 ^ 2) / 0.2e1 + m(6) * (t104 ^ 2 + t107 ^ 2 + t108 ^ 2) / 0.2e1 + m(5) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + t179 * (t232 * t187 + t239 * t188) / 0.2e1 + t178 * (-t239 * t187 + t232 * t188) / 0.2e1 + ((t287 * t148 + t285 * t149 - t286 * t266) * t208 + (t292 * t148 + t288 * t149 - t290 * t266) * t147 + (t293 * t148 + t289 * t149 - t291 * t266) * t146) * t146 / 0.2e1 + ((t287 * t150 + t285 * t151 - t286 * t267) * t208 + (t292 * t150 + t288 * t151 - t290 * t267) * t147 + (t293 * t150 + t289 * t151 - t291 * t267) * t146) * t147 / 0.2e1 + ((t291 * t146 + t290 * t147 + t286 * t208) * t231 + ((-t287 * t228 - t285 * t230) * t208 + (-t292 * t228 - t288 * t230) * t147 + (-t293 * t228 - t289 * t230) * t146) * t229) * t208 / 0.2e1 + ((-t138 * t231 - t229 * t140) * t179 + (-t137 * t231 - t229 * t139) * t178 + (-t198 * t231 - t229 * t203 + Icges(3,2) + Icges(2,3) + Icges(4,3)) * t222) * t222 / 0.2e1 + ((-t158 * t188 - t160 * t187 + t284 * t273 + t282 * t274 + Icges(1,4)) * V_base(5) + (-t159 * t188 - t161 * t187 + t283 * t273 + t281 * t274 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t158 * t187 - t160 * t188 + t282 * t273 - t284 * t274 + Icges(1,2)) * V_base(5) + (t159 * t187 - t161 * t188 + t281 * t273 - t283 * t274 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t222 * (Icges(4,5) * t188 - Icges(4,6) * t187 + t305 * t273 + t300 * t274) + V_base(4) * t222 * (Icges(4,5) * t187 + Icges(4,6) * t188 - t300 * t273 + t305 * t274) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
