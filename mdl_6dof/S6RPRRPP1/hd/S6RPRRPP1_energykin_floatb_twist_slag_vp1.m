% Calculate kinetic energy for
% S6RPRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
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
% Datum: 2019-03-09 04:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP1_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRPP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP1_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP1_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPP1_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:27:49
% EndTime: 2019-03-09 04:27:52
% DurationCPUTime: 2.93s
% Computational Cost: add. (2169->310), mult. (2034->421), div. (0->0), fcn. (1936->10), ass. (0->150)
t317 = Icges(6,1) + Icges(7,1);
t316 = -Icges(6,4) + Icges(7,5);
t315 = Icges(7,4) + Icges(6,5);
t314 = Icges(6,2) + Icges(7,3);
t313 = -Icges(7,6) + Icges(6,6);
t312 = -Icges(6,3) - Icges(7,2) - Icges(5,3);
t311 = rSges(7,1) + pkin(5);
t310 = rSges(7,3) + qJ(6);
t241 = qJ(4) + pkin(10);
t233 = sin(t241);
t235 = cos(t241);
t242 = qJ(1) + pkin(9);
t236 = cos(t242);
t234 = sin(t242);
t248 = cos(qJ(3));
t283 = t234 * t248;
t167 = t233 * t283 + t235 * t236;
t168 = -t233 * t236 + t235 * t283;
t245 = sin(qJ(3));
t284 = t234 * t245;
t309 = t167 * t314 + t168 * t316 - t284 * t313;
t280 = t236 * t248;
t169 = t233 * t280 - t234 * t235;
t170 = t233 * t234 + t235 * t280;
t281 = t236 * t245;
t308 = t169 * t314 + t170 * t316 - t281 * t313;
t307 = t316 * t167 + t168 * t317 + t315 * t284;
t306 = t316 * t169 + t170 * t317 + t315 * t281;
t305 = t313 * t248 + (t233 * t314 + t235 * t316) * t245;
t304 = -t315 * t248 + (t316 * t233 + t235 * t317) * t245;
t247 = cos(qJ(4));
t244 = sin(qJ(4));
t279 = t244 * t248;
t183 = -t234 * t279 - t236 * t247;
t278 = t247 * t248;
t282 = t236 * t244;
t184 = t234 * t278 - t282;
t301 = Icges(5,5) * t184 + Icges(5,6) * t183 - t167 * t313 + t168 * t315 - t284 * t312;
t185 = t234 * t247 - t236 * t279;
t285 = t234 * t244;
t186 = t236 * t278 + t285;
t300 = Icges(5,5) * t186 + Icges(5,6) * t185 - t169 * t313 + t170 * t315 - t281 * t312;
t299 = t312 * t248 + (Icges(5,5) * t247 - Icges(5,6) * t244 - t233 * t313 + t235 * t315) * t245;
t246 = sin(qJ(1));
t294 = pkin(1) * t246;
t249 = cos(qJ(1));
t293 = pkin(1) * t249;
t292 = pkin(4) * t247;
t291 = -pkin(6) - qJ(2);
t289 = Icges(2,4) * t246;
t288 = Icges(3,4) * t234;
t287 = Icges(4,4) * t245;
t286 = Icges(4,4) * t248;
t277 = rSges(7,2) * t284 + t310 * t167 + t311 * t168;
t276 = rSges(7,2) * t281 + t310 * t169 + t311 * t170;
t275 = -rSges(7,2) * t248 + (t310 * t233 + t311 * t235) * t245;
t274 = qJD(4) * t245;
t273 = qJD(5) * t245;
t237 = V_base(6) + qJD(1);
t272 = t237 * t293 + V_base(2);
t271 = V_base(5) * pkin(6) + V_base(1);
t211 = qJD(3) * t234 + V_base(4);
t204 = pkin(2) * t234 - pkin(7) * t236;
t268 = -t204 - t294;
t267 = V_base(5) * qJ(2) + t271;
t266 = V_base(4) * t294 + qJD(2) + V_base(3);
t265 = pkin(3) * t248 + pkin(8) * t245;
t210 = -qJD(3) * t236 + V_base(5);
t264 = rSges(4,1) * t248 - rSges(4,2) * t245;
t263 = Icges(4,1) * t248 - t287;
t262 = -Icges(4,2) * t245 + t286;
t261 = Icges(4,5) * t248 - Icges(4,6) * t245;
t260 = qJ(5) * t245 + t248 * t292;
t259 = (-Icges(4,3) * t236 + t234 * t261) * t210 + (Icges(4,3) * t234 + t236 * t261) * t211 + (Icges(4,5) * t245 + Icges(4,6) * t248) * t237;
t205 = pkin(2) * t236 + pkin(7) * t234;
t258 = t237 * t205 + t291 * V_base(4) + t272;
t193 = t265 * t234;
t227 = pkin(3) * t245 - pkin(8) * t248;
t257 = t210 * t227 + (-t193 + t268) * t237 + t267;
t256 = V_base(4) * t204 + (-t205 - t293) * V_base(5) + t266;
t194 = t265 * t236;
t255 = t237 * t194 - t211 * t227 + t258;
t166 = -qJ(5) * t248 + t245 * t292;
t181 = t234 * t274 + t210;
t254 = t181 * t166 + t236 * t273 + t257;
t145 = pkin(4) * t285 + t236 * t260;
t223 = -qJD(4) * t248 + t237;
t253 = t223 * t145 + t234 * t273 + t255;
t252 = t211 * t193 - t210 * t194 + t256;
t144 = -pkin(4) * t282 + t234 * t260;
t182 = t236 * t274 + t211;
t251 = -qJD(5) * t248 + t182 * t144 + t252;
t158 = -Icges(4,6) * t236 + t234 * t262;
t159 = Icges(4,6) * t234 + t236 * t262;
t160 = -Icges(4,5) * t236 + t234 * t263;
t161 = Icges(4,5) * t234 + t236 * t263;
t215 = Icges(4,2) * t248 + t287;
t218 = Icges(4,1) * t245 + t286;
t250 = (-t159 * t245 + t161 * t248) * t211 + (-t158 * t245 + t160 * t248) * t210 + (-t215 * t245 + t218 * t248) * t237;
t239 = Icges(2,4) * t249;
t231 = Icges(3,4) * t236;
t226 = rSges(2,1) * t249 - t246 * rSges(2,2);
t225 = t246 * rSges(2,1) + rSges(2,2) * t249;
t224 = rSges(4,1) * t245 + rSges(4,2) * t248;
t220 = Icges(2,1) * t249 - t289;
t219 = Icges(2,1) * t246 + t239;
t217 = -Icges(2,2) * t246 + t239;
t216 = Icges(2,2) * t249 + t289;
t208 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t207 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t206 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t203 = rSges(3,1) * t236 - rSges(3,2) * t234;
t202 = rSges(3,1) * t234 + rSges(3,2) * t236;
t201 = Icges(3,1) * t236 - t288;
t200 = Icges(3,1) * t234 + t231;
t199 = -Icges(3,2) * t234 + t231;
t198 = Icges(3,2) * t236 + t288;
t191 = -rSges(5,3) * t248 + (rSges(5,1) * t247 - rSges(5,2) * t244) * t245;
t189 = -Icges(5,5) * t248 + (Icges(5,1) * t247 - Icges(5,4) * t244) * t245;
t188 = -Icges(5,6) * t248 + (Icges(5,4) * t247 - Icges(5,2) * t244) * t245;
t179 = -rSges(6,3) * t248 + (rSges(6,1) * t235 - rSges(6,2) * t233) * t245;
t165 = rSges(4,3) * t234 + t236 * t264;
t164 = -rSges(4,3) * t236 + t234 * t264;
t163 = V_base(5) * rSges(2,3) - t225 * t237 + t271;
t162 = t226 * t237 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t155 = t225 * V_base(4) - t226 * V_base(5) + V_base(3);
t153 = V_base(5) * rSges(3,3) + (-t202 - t294) * t237 + t267;
t152 = t203 * t237 + (-rSges(3,3) + t291) * V_base(4) + t272;
t150 = V_base(4) * t202 + (-t203 - t293) * V_base(5) + t266;
t147 = rSges(5,1) * t186 + rSges(5,2) * t185 + rSges(5,3) * t281;
t146 = rSges(5,1) * t184 + rSges(5,2) * t183 + rSges(5,3) * t284;
t143 = Icges(5,1) * t186 + Icges(5,4) * t185 + Icges(5,5) * t281;
t142 = Icges(5,1) * t184 + Icges(5,4) * t183 + Icges(5,5) * t284;
t141 = Icges(5,4) * t186 + Icges(5,2) * t185 + Icges(5,6) * t281;
t140 = Icges(5,4) * t184 + Icges(5,2) * t183 + Icges(5,6) * t284;
t137 = rSges(6,1) * t170 - rSges(6,2) * t169 + rSges(6,3) * t281;
t135 = rSges(6,1) * t168 - rSges(6,2) * t167 + rSges(6,3) * t284;
t119 = t210 * t224 + (-t164 + t268) * t237 + t267;
t118 = t165 * t237 - t211 * t224 + t258;
t117 = t211 * t164 - t210 * t165 + t256;
t116 = -t146 * t223 + t181 * t191 + t257;
t115 = t147 * t223 - t182 * t191 + t255;
t114 = t182 * t146 - t181 * t147 + t252;
t113 = t179 * t181 + (-t135 - t144) * t223 + t254;
t112 = t137 * t223 + (-t166 - t179) * t182 + t253;
t111 = t182 * t135 + (-t137 - t145) * t181 + t251;
t110 = qJD(6) * t169 + t275 * t181 + (-t144 - t277) * t223 + t254;
t109 = qJD(6) * t167 + t276 * t223 + (-t166 - t275) * t182 + t253;
t108 = qJD(6) * t245 * t233 + t277 * t182 + (-t145 - t276) * t181 + t251;
t1 = t211 * (t259 * t234 + t250 * t236) / 0.2e1 + t210 * (t250 * t234 - t259 * t236) / 0.2e1 + m(1) * (t206 ^ 2 + t207 ^ 2 + t208 ^ 2) / 0.2e1 + m(2) * (t155 ^ 2 + t162 ^ 2 + t163 ^ 2) / 0.2e1 + m(3) * (t150 ^ 2 + t152 ^ 2 + t153 ^ 2) / 0.2e1 + m(4) * (t117 ^ 2 + t118 ^ 2 + t119 ^ 2) / 0.2e1 + m(7) * (t108 ^ 2 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + m(6) * (t111 ^ 2 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + m(5) * (t114 ^ 2 + t115 ^ 2 + t116 ^ 2) / 0.2e1 + ((t167 * t305 + t168 * t304 + t183 * t188 + t184 * t189 + t284 * t299) * t223 + (t141 * t183 + t143 * t184 + t167 * t308 + t168 * t306 + t284 * t300) * t182 + (t183 * t140 + t184 * t142 + t309 * t167 + t307 * t168 + t301 * t284) * t181) * t181 / 0.2e1 + ((t169 * t305 + t170 * t304 + t185 * t188 + t186 * t189 + t281 * t299) * t223 + (t185 * t141 + t186 * t143 + t308 * t169 + t306 * t170 + t300 * t281) * t182 + (t140 * t185 + t142 * t186 + t169 * t309 + t307 * t170 + t301 * t281) * t181) * t182 / 0.2e1 + ((-t181 * t301 - t182 * t300 - t223 * t299) * t248 + ((-t188 * t244 + t189 * t247 + t233 * t305 + t235 * t304) * t223 + (-t141 * t244 + t143 * t247 + t233 * t308 + t235 * t306) * t182 + (-t140 * t244 + t142 * t247 + t233 * t309 + t307 * t235) * t181) * t245) * t223 / 0.2e1 + ((t159 * t248 + t161 * t245) * t211 + (t158 * t248 + t160 * t245) * t210 + (t248 * t215 + t245 * t218 + Icges(2,3) + Icges(3,3)) * t237) * t237 / 0.2e1 + ((-t198 * t234 + t200 * t236 - t246 * t216 + t219 * t249 + Icges(1,4)) * V_base(5) + (-t234 * t199 + t236 * t201 - t246 * t217 + t249 * t220 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t236 * t198 + t234 * t200 + t249 * t216 + t246 * t219 + Icges(1,2)) * V_base(5) + (t199 * t236 + t201 * t234 + t217 * t249 + t246 * t220 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t237 * (Icges(2,5) * t246 + Icges(3,5) * t234 + Icges(2,6) * t249 + Icges(3,6) * t236) + V_base(4) * t237 * (Icges(2,5) * t249 + Icges(3,5) * t236 - Icges(2,6) * t246 - Icges(3,6) * t234) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
