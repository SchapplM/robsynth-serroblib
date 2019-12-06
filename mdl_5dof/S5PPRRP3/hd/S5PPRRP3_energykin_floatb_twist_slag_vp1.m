% Calculate kinetic energy for
% S5PPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRRP3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP3_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PPRRP3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP3_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRP3_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRP3_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:10:29
% EndTime: 2019-12-05 15:10:32
% DurationCPUTime: 2.38s
% Computational Cost: add. (1149->281), mult. (2582->395), div. (0->0), fcn. (2825->8), ass. (0->135)
t296 = Icges(5,1) + Icges(6,1);
t295 = -Icges(5,4) + Icges(6,5);
t294 = Icges(6,4) + Icges(5,5);
t293 = Icges(5,2) + Icges(6,3);
t292 = Icges(6,2) + Icges(5,3);
t291 = -Icges(5,6) + Icges(6,6);
t290 = rSges(6,1) + pkin(4);
t289 = rSges(6,3) + qJ(5);
t231 = sin(pkin(7));
t233 = cos(pkin(7));
t288 = Icges(2,5) * t233 - Icges(2,6) * t231 + Icges(1,5);
t287 = Icges(2,5) * t231 + Icges(2,6) * t233 + Icges(1,6);
t232 = cos(pkin(8));
t235 = sin(qJ(3));
t264 = t233 * t235;
t236 = cos(qJ(3));
t265 = t231 * t236;
t189 = t232 * t265 - t264;
t234 = sin(qJ(4));
t230 = sin(pkin(8));
t274 = cos(qJ(4));
t253 = t230 * t274;
t167 = t189 * t234 - t231 * t253;
t268 = t230 * t234;
t168 = t189 * t274 + t231 * t268;
t263 = t233 * t236;
t266 = t231 * t235;
t188 = t232 * t266 + t263;
t286 = t293 * t167 + t295 * t168 + t291 * t188;
t191 = t232 * t263 + t266;
t169 = t191 * t234 - t233 * t253;
t170 = t191 * t274 + t233 * t268;
t190 = t232 * t264 - t265;
t285 = t293 * t169 + t295 * t170 + t291 * t190;
t284 = t291 * t167 + t294 * t168 + t292 * t188;
t283 = t291 * t169 + t294 * t170 + t292 * t190;
t282 = t295 * t167 + t296 * t168 + t294 * t188;
t281 = t295 * t169 + t296 * t170 + t294 * t190;
t193 = t232 * t274 + t236 * t268;
t194 = -t232 * t234 + t236 * t253;
t267 = t230 * t235;
t280 = t293 * t193 + t295 * t194 + t291 * t267;
t279 = t291 * t193 + t294 * t194 + t292 * t267;
t278 = t295 * t193 + t296 * t194 + t294 * t267;
t273 = Icges(2,4) * t231;
t272 = Icges(3,4) * t230;
t271 = Icges(3,4) * t232;
t270 = t230 * t231;
t269 = t230 * t233;
t262 = rSges(6,2) * t188 + t289 * t167 + t290 * t168;
t261 = rSges(6,2) * t190 + t289 * t169 + t290 * t170;
t260 = rSges(6,2) * t267 + t289 * t193 + t290 * t194;
t259 = qJD(3) * t230;
t258 = V_base(5) * qJ(1) + V_base(1);
t254 = qJD(1) + V_base(3);
t205 = t233 * t259 + V_base(4);
t204 = t231 * t259 + V_base(5);
t252 = qJD(2) * t231 + t258;
t216 = pkin(1) * t231 - qJ(2) * t233;
t251 = V_base(4) * t216 + t254;
t250 = pkin(2) * t232 + pkin(5) * t230;
t223 = -qJD(3) * t232 + V_base(6);
t249 = rSges(3,1) * t232 - rSges(3,2) * t230;
t248 = Icges(3,1) * t232 - t272;
t247 = -Icges(3,2) * t230 + t271;
t246 = Icges(3,5) * t232 - Icges(3,6) * t230;
t218 = pkin(1) * t233 + qJ(2) * t231;
t245 = -qJD(2) * t233 + V_base(6) * t218 + V_base(2);
t195 = t250 * t231;
t220 = pkin(2) * t230 - pkin(5) * t232;
t244 = V_base(5) * t220 + (-t195 - t216) * V_base(6) + t252;
t196 = t250 * t233;
t243 = V_base(4) * t195 + (-t196 - t218) * V_base(5) + t251;
t242 = (-Icges(3,3) * t233 + t231 * t246) * V_base(5) + (Icges(3,3) * t231 + t233 * t246) * V_base(4) + (Icges(3,5) * t230 + Icges(3,6) * t232) * V_base(6);
t241 = V_base(6) * t196 + (-qJ(1) - t220) * V_base(4) + t245;
t161 = pkin(3) * t189 + pkin(6) * t188;
t197 = (pkin(3) * t236 + pkin(6) * t235) * t230;
t240 = -t161 * t223 + t204 * t197 + t244;
t162 = pkin(3) * t191 + pkin(6) * t190;
t239 = t205 * t161 - t162 * t204 + t243;
t238 = t223 * t162 - t197 * t205 + t241;
t176 = -Icges(3,6) * t233 + t231 * t247;
t177 = Icges(3,6) * t231 + t233 * t247;
t178 = -Icges(3,5) * t233 + t231 * t248;
t179 = Icges(3,5) * t231 + t233 * t248;
t209 = Icges(3,2) * t232 + t272;
t212 = Icges(3,1) * t230 + t271;
t237 = (-t177 * t230 + t179 * t232) * V_base(4) + (-t176 * t230 + t178 * t232) * V_base(5) + (-t209 * t230 + t212 * t232) * V_base(6);
t228 = Icges(2,4) * t233;
t219 = rSges(2,1) * t233 - rSges(2,2) * t231;
t217 = rSges(2,1) * t231 + rSges(2,2) * t233;
t215 = rSges(3,1) * t230 + rSges(3,2) * t232;
t214 = Icges(2,1) * t233 - t273;
t213 = Icges(2,1) * t231 + t228;
t211 = -Icges(2,2) * t231 + t228;
t210 = Icges(2,2) * t233 + t273;
t203 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t202 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t201 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t192 = qJD(4) * t267 + t223;
t185 = -t232 * rSges(4,3) + (rSges(4,1) * t236 - rSges(4,2) * t235) * t230;
t184 = -Icges(4,5) * t232 + (Icges(4,1) * t236 - Icges(4,4) * t235) * t230;
t183 = -Icges(4,6) * t232 + (Icges(4,4) * t236 - Icges(4,2) * t235) * t230;
t182 = -Icges(4,3) * t232 + (Icges(4,5) * t236 - Icges(4,6) * t235) * t230;
t181 = rSges(3,3) * t231 + t233 * t249;
t180 = -rSges(3,3) * t233 + t231 * t249;
t173 = V_base(5) * rSges(2,3) - t217 * V_base(6) + t258;
t172 = t219 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t166 = qJD(4) * t190 + t205;
t165 = qJD(4) * t188 + t204;
t164 = t217 * V_base(4) - t219 * V_base(5) + t254;
t160 = rSges(5,1) * t194 - rSges(5,2) * t193 + rSges(5,3) * t267;
t152 = rSges(4,1) * t191 - rSges(4,2) * t190 + rSges(4,3) * t269;
t151 = rSges(4,1) * t189 - rSges(4,2) * t188 + rSges(4,3) * t270;
t150 = Icges(4,1) * t191 - Icges(4,4) * t190 + Icges(4,5) * t269;
t149 = Icges(4,1) * t189 - Icges(4,4) * t188 + Icges(4,5) * t270;
t148 = Icges(4,4) * t191 - Icges(4,2) * t190 + Icges(4,6) * t269;
t147 = Icges(4,4) * t189 - Icges(4,2) * t188 + Icges(4,6) * t270;
t146 = Icges(4,5) * t191 - Icges(4,6) * t190 + Icges(4,3) * t269;
t145 = Icges(4,5) * t189 - Icges(4,6) * t188 + Icges(4,3) * t270;
t140 = t215 * V_base(5) + (-t180 - t216) * V_base(6) + t252;
t139 = t181 * V_base(6) + (-qJ(1) - t215) * V_base(4) + t245;
t138 = t180 * V_base(4) + (-t181 - t218) * V_base(5) + t251;
t137 = rSges(5,1) * t170 - rSges(5,2) * t169 + rSges(5,3) * t190;
t135 = rSges(5,1) * t168 - rSges(5,2) * t167 + rSges(5,3) * t188;
t121 = -t151 * t223 + t185 * t204 + t244;
t120 = t152 * t223 - t185 * t205 + t241;
t119 = t151 * t205 - t152 * t204 + t243;
t118 = -t135 * t192 + t160 * t165 + t240;
t117 = t137 * t192 - t160 * t166 + t238;
t116 = t135 * t166 - t137 * t165 + t239;
t115 = qJD(5) * t169 + t165 * t260 - t192 * t262 + t240;
t114 = qJD(5) * t167 - t166 * t260 + t192 * t261 + t238;
t113 = qJD(5) * t193 - t165 * t261 + t166 * t262 + t239;
t1 = m(1) * (t201 ^ 2 + t202 ^ 2 + t203 ^ 2) / 0.2e1 + m(2) * (t164 ^ 2 + t172 ^ 2 + t173 ^ 2) / 0.2e1 + m(3) * (t138 ^ 2 + t139 ^ 2 + t140 ^ 2) / 0.2e1 + m(4) * (t119 ^ 2 + t120 ^ 2 + t121 ^ 2) / 0.2e1 + t205 * ((t146 * t269 - t148 * t190 + t150 * t191) * t205 + (t145 * t269 - t147 * t190 + t149 * t191) * t204 + (t182 * t269 - t183 * t190 + t184 * t191) * t223) / 0.2e1 + t204 * ((t146 * t270 - t148 * t188 + t150 * t189) * t205 + (t145 * t270 - t147 * t188 + t149 * t189) * t204 + (t182 * t270 - t183 * t188 + t184 * t189) * t223) / 0.2e1 + t223 * ((-t145 * t204 - t146 * t205 - t182 * t223) * t232 + ((-t148 * t235 + t150 * t236) * t205 + (-t147 * t235 + t149 * t236) * t204 + (-t183 * t235 + t184 * t236) * t223) * t230) / 0.2e1 + m(5) * (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + m(6) * (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + ((t280 * t167 + t278 * t168 + t279 * t188) * t192 + (t285 * t167 + t281 * t168 + t283 * t188) * t166 + (t286 * t167 + t282 * t168 + t284 * t188) * t165) * t165 / 0.2e1 + ((t280 * t169 + t278 * t170 + t279 * t190) * t192 + (t285 * t169 + t281 * t170 + t283 * t190) * t166 + (t286 * t169 + t282 * t170 + t284 * t190) * t165) * t166 / 0.2e1 + ((t280 * t193 + t278 * t194 + t279 * t267) * t192 + (t285 * t193 + t281 * t194 + t283 * t267) * t166 + (t286 * t193 + t282 * t194 + t284 * t267) * t165) * t192 / 0.2e1 + (t242 * t231 + t237 * t233 + t288 * V_base(6) + (-t210 * t231 + t213 * t233 + Icges(1,4)) * V_base(5) + (-t211 * t231 + t214 * t233 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t237 * t231 - t242 * t233 + t287 * V_base(6) + (t210 * t233 + t213 * t231 + Icges(1,2)) * V_base(5) + (t211 * t233 + t214 * t231 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t209 * t232 + t212 * t230 + Icges(1,3) + Icges(2,3)) * V_base(6) + (t176 * t232 + t178 * t230 + t287) * V_base(5) + (t177 * t232 + t179 * t230 + t288) * V_base(4)) * V_base(6) / 0.2e1;
T = t1;
