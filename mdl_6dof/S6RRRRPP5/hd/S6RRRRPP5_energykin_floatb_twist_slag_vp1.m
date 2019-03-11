% Calculate kinetic energy for
% S6RRRRPP5
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
% Datum: 2019-03-09 21:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPP5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP5_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRPP5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP5_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP5_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPP5_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:05:11
% EndTime: 2019-03-09 21:05:14
% DurationCPUTime: 3.07s
% Computational Cost: add. (1885->308), mult. (2687->449), div. (0->0), fcn. (2661->8), ass. (0->146)
t315 = Icges(5,1) + Icges(6,1) + Icges(7,1);
t314 = -Icges(5,4) + Icges(7,4) + Icges(6,5);
t313 = Icges(7,5) - Icges(6,4) - Icges(5,5);
t312 = Icges(5,2) + Icges(7,2) + Icges(6,3);
t311 = -Icges(6,6) + Icges(7,6) + Icges(5,6);
t310 = Icges(7,3) + Icges(5,3) + Icges(6,2);
t309 = rSges(7,1) + pkin(5);
t308 = rSges(7,3) + qJ(6);
t248 = qJ(3) + qJ(4);
t244 = sin(t248);
t245 = cos(t248);
t254 = cos(qJ(1));
t251 = sin(qJ(1));
t253 = cos(qJ(2));
t285 = t251 * t253;
t198 = t244 * t285 + t245 * t254;
t199 = -t244 * t254 + t245 * t285;
t250 = sin(qJ(2));
t287 = t250 * t251;
t307 = t311 * t198 + t313 * t199 - t310 * t287;
t284 = t253 * t254;
t200 = t244 * t284 - t251 * t245;
t201 = t244 * t251 + t245 * t284;
t286 = t250 * t254;
t306 = t311 * t200 + t313 * t201 - t310 * t286;
t305 = t312 * t198 + t314 * t199 - t311 * t287;
t304 = t312 * t200 + t314 * t201 - t311 * t286;
t303 = t314 * t198 + t315 * t199 - t313 * t287;
t302 = t314 * t200 + t315 * t201 - t313 * t286;
t301 = t310 * t253 + (t311 * t244 + t313 * t245) * t250;
t300 = t311 * t253 + (t312 * t244 + t314 * t245) * t250;
t299 = t313 * t253 + (t314 * t244 + t315 * t245) * t250;
t252 = cos(qJ(3));
t294 = pkin(3) * t252;
t292 = Icges(2,4) * t251;
t291 = Icges(3,4) * t250;
t290 = Icges(3,4) * t253;
t249 = sin(qJ(3));
t289 = t249 * t251;
t288 = t249 * t254;
t283 = rSges(7,2) * t198 + t199 * t309 - t308 * t287;
t282 = rSges(7,2) * t200 + t201 * t309 - t308 * t286;
t281 = t308 * t253 + (rSges(7,2) * t244 + t245 * t309) * t250;
t280 = qJD(3) * t250;
t279 = qJD(4) * t250;
t278 = qJD(6) * t250;
t277 = V_base(5) * pkin(6) + V_base(1);
t237 = qJD(2) * t251 + V_base(4);
t242 = V_base(6) + qJD(1);
t205 = t254 * t280 + t237;
t274 = pkin(2) * t253 + pkin(8) * t250;
t236 = -qJD(2) * t254 + V_base(5);
t273 = rSges(3,1) * t253 - rSges(3,2) * t250;
t272 = Icges(3,1) * t253 - t291;
t271 = -Icges(3,2) * t250 + t290;
t270 = Icges(3,5) * t253 - Icges(3,6) * t250;
t204 = t251 * t280 + t236;
t235 = pkin(1) * t254 + pkin(7) * t251;
t269 = -V_base(4) * pkin(6) + t242 * t235 + V_base(2);
t234 = pkin(1) * t251 - pkin(7) * t254;
t268 = V_base(4) * t234 - t235 * V_base(5) + V_base(3);
t267 = pkin(9) * t250 + t253 * t294;
t211 = t274 * t251;
t233 = t250 * pkin(2) - t253 * pkin(8);
t266 = t236 * t233 + (-t211 - t234) * t242 + t277;
t265 = (-Icges(3,3) * t254 + t251 * t270) * t236 + (Icges(3,3) * t251 + t254 * t270) * t237 + (Icges(3,5) * t250 + Icges(3,6) * t253) * t242;
t212 = t274 * t254;
t264 = t242 * t212 - t233 * t237 + t269;
t263 = t237 * t211 - t212 * t236 + t268;
t159 = -pkin(3) * t288 + t251 * t267;
t167 = -pkin(9) * t253 + t250 * t294;
t228 = -qJD(3) * t253 + t242;
t262 = -t159 * t228 + t204 * t167 + t266;
t160 = pkin(3) * t289 + t254 * t267;
t261 = t228 * t160 - t167 * t205 + t264;
t181 = t251 * t279 + t204;
t202 = (pkin(4) * t245 + qJ(5) * t244) * t250;
t260 = qJD(5) * t200 + t181 * t202 + t262;
t259 = t205 * t159 - t160 * t204 + t263;
t149 = pkin(4) * t201 + qJ(5) * t200;
t213 = (-qJD(3) - qJD(4)) * t253 + t242;
t258 = qJD(5) * t198 + t213 * t149 + t261;
t148 = pkin(4) * t199 + qJ(5) * t198;
t182 = t254 * t279 + t205;
t257 = qJD(5) * t250 * t244 + t182 * t148 + t259;
t190 = -Icges(3,6) * t254 + t251 * t271;
t191 = Icges(3,6) * t251 + t254 * t271;
t193 = -Icges(3,5) * t254 + t251 * t272;
t194 = Icges(3,5) * t251 + t254 * t272;
t222 = Icges(3,2) * t253 + t291;
t225 = Icges(3,1) * t250 + t290;
t256 = (-t191 * t250 + t194 * t253) * t237 + (-t190 * t250 + t193 * t253) * t236 + (-t222 * t250 + t225 * t253) * t242;
t246 = Icges(2,4) * t254;
t231 = rSges(2,1) * t254 - rSges(2,2) * t251;
t230 = rSges(2,1) * t251 + rSges(2,2) * t254;
t229 = rSges(3,1) * t250 + rSges(3,2) * t253;
t227 = Icges(2,1) * t254 - t292;
t226 = Icges(2,1) * t251 + t246;
t224 = -Icges(2,2) * t251 + t246;
t223 = Icges(2,2) * t254 + t292;
t218 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t217 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t216 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t209 = t252 * t284 + t289;
t208 = -t249 * t284 + t251 * t252;
t207 = t252 * t285 - t288;
t206 = -t249 * t285 - t252 * t254;
t197 = rSges(3,3) * t251 + t254 * t273;
t196 = -rSges(3,3) * t254 + t251 * t273;
t195 = -rSges(4,3) * t253 + (rSges(4,1) * t252 - rSges(4,2) * t249) * t250;
t192 = -Icges(4,5) * t253 + (Icges(4,1) * t252 - Icges(4,4) * t249) * t250;
t189 = -Icges(4,6) * t253 + (Icges(4,4) * t252 - Icges(4,2) * t249) * t250;
t186 = -Icges(4,3) * t253 + (Icges(4,5) * t252 - Icges(4,6) * t249) * t250;
t180 = -rSges(5,3) * t253 + (rSges(5,1) * t245 - rSges(5,2) * t244) * t250;
t179 = -rSges(6,2) * t253 + (rSges(6,1) * t245 + rSges(6,3) * t244) * t250;
t166 = V_base(5) * rSges(2,3) - t230 * t242 + t277;
t165 = t231 * t242 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t164 = t230 * V_base(4) - t231 * V_base(5) + V_base(3);
t158 = rSges(4,1) * t209 + rSges(4,2) * t208 + rSges(4,3) * t286;
t157 = rSges(4,1) * t207 + rSges(4,2) * t206 + rSges(4,3) * t287;
t155 = Icges(4,1) * t209 + Icges(4,4) * t208 + Icges(4,5) * t286;
t154 = Icges(4,1) * t207 + Icges(4,4) * t206 + Icges(4,5) * t287;
t153 = Icges(4,4) * t209 + Icges(4,2) * t208 + Icges(4,6) * t286;
t152 = Icges(4,4) * t207 + Icges(4,2) * t206 + Icges(4,6) * t287;
t151 = Icges(4,5) * t209 + Icges(4,6) * t208 + Icges(4,3) * t286;
t150 = Icges(4,5) * t207 + Icges(4,6) * t206 + Icges(4,3) * t287;
t147 = rSges(5,1) * t201 - rSges(5,2) * t200 + rSges(5,3) * t286;
t146 = rSges(6,1) * t201 + rSges(6,2) * t286 + rSges(6,3) * t200;
t144 = rSges(5,1) * t199 - rSges(5,2) * t198 + rSges(5,3) * t287;
t143 = rSges(6,1) * t199 + rSges(6,2) * t287 + rSges(6,3) * t198;
t120 = t229 * t236 + (-t196 - t234) * t242 + t277;
t119 = t197 * t242 - t229 * t237 + t269;
t117 = t196 * t237 - t197 * t236 + t268;
t116 = -t157 * t228 + t195 * t204 + t266;
t115 = t158 * t228 - t195 * t205 + t264;
t114 = t157 * t205 - t158 * t204 + t263;
t113 = -t144 * t213 + t180 * t181 + t262;
t112 = t147 * t213 - t180 * t182 + t261;
t111 = t144 * t182 - t147 * t181 + t259;
t110 = t179 * t181 + (-t143 - t148) * t213 + t260;
t109 = t146 * t213 + (-t179 - t202) * t182 + t258;
t108 = -t254 * t278 + t281 * t181 + (-t148 - t283) * t213 + t260;
t107 = -t251 * t278 + t282 * t213 + (-t202 - t281) * t182 + t258;
t106 = t143 * t182 + (-t146 - t149) * t181 + t257;
t105 = qJD(6) * t253 + t283 * t182 + (-t149 - t282) * t181 + t257;
t1 = t237 * (t265 * t251 + t256 * t254) / 0.2e1 + t236 * (t256 * t251 - t265 * t254) / 0.2e1 + t228 * ((-t150 * t204 - t151 * t205 - t186 * t228) * t253 + ((-t153 * t249 + t155 * t252) * t205 + (-t152 * t249 + t154 * t252) * t204 + (-t189 * t249 + t192 * t252) * t228) * t250) / 0.2e1 + t204 * ((t151 * t287 + t153 * t206 + t155 * t207) * t205 + (t150 * t287 + t206 * t152 + t207 * t154) * t204 + (t186 * t287 + t189 * t206 + t192 * t207) * t228) / 0.2e1 + t205 * ((t151 * t286 + t208 * t153 + t209 * t155) * t205 + (t150 * t286 + t152 * t208 + t154 * t209) * t204 + (t186 * t286 + t189 * t208 + t192 * t209) * t228) / 0.2e1 + m(1) * (t216 ^ 2 + t217 ^ 2 + t218 ^ 2) / 0.2e1 + m(2) * (t164 ^ 2 + t165 ^ 2 + t166 ^ 2) / 0.2e1 + m(3) * (t117 ^ 2 + t119 ^ 2 + t120 ^ 2) / 0.2e1 + m(4) * (t114 ^ 2 + t115 ^ 2 + t116 ^ 2) / 0.2e1 + m(7) * (t105 ^ 2 + t107 ^ 2 + t108 ^ 2) / 0.2e1 + m(6) * (t106 ^ 2 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + m(5) * (t111 ^ 2 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + ((t191 * t253 + t194 * t250) * t237 + (t190 * t253 + t193 * t250) * t236 + (t253 * t222 + t250 * t225 + Icges(2,3)) * t242) * t242 / 0.2e1 + ((-t223 * t251 + t226 * t254 + Icges(1,4)) * V_base(5) + (-t251 * t224 + t254 * t227 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t254 * t223 + t251 * t226 + Icges(1,2)) * V_base(5) + (t224 * t254 + t227 * t251 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t198 * t300 + t199 * t299 - t287 * t301) * t213 + (t198 * t304 + t199 * t302 - t287 * t306) * t182 + (t305 * t198 + t303 * t199 - t307 * t287) * t181) * t181 / 0.2e1 + ((t200 * t300 + t201 * t299 - t286 * t301) * t213 + (t304 * t200 + t302 * t201 - t306 * t286) * t182 + (t305 * t200 + t303 * t201 - t286 * t307) * t181) * t182 / 0.2e1 + ((t181 * t307 + t306 * t182 + t301 * t213) * t253 + ((t244 * t300 + t245 * t299) * t213 + (t244 * t304 + t245 * t302) * t182 + (t244 * t305 + t245 * t303) * t181) * t250) * t213 / 0.2e1 + V_base(4) * t242 * (Icges(2,5) * t254 - Icges(2,6) * t251) + V_base(5) * t242 * (Icges(2,5) * t251 + Icges(2,6) * t254) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
