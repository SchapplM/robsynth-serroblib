% Calculate kinetic energy for
% S6RRPRPP2
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
% Datum: 2019-03-09 09:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPP2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP2_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRPP2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP2_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPP2_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPP2_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:49:39
% EndTime: 2019-03-09 09:49:41
% DurationCPUTime: 2.58s
% Computational Cost: add. (1811->281), mult. (2227->385), div. (0->0), fcn. (2145->8), ass. (0->141)
t318 = Icges(3,3) + Icges(4,3);
t237 = qJ(2) + pkin(9);
t230 = sin(t237);
t231 = cos(t237);
t240 = sin(qJ(2));
t243 = cos(qJ(2));
t317 = Icges(3,5) * t243 + Icges(4,5) * t231 - Icges(3,6) * t240 - Icges(4,6) * t230;
t316 = Icges(5,1) + Icges(6,1) + Icges(7,1);
t315 = -Icges(5,4) + Icges(7,4) + Icges(6,5);
t314 = Icges(7,5) - Icges(6,4) - Icges(5,5);
t313 = Icges(5,2) + Icges(7,2) + Icges(6,3);
t312 = -Icges(6,6) + Icges(7,6) + Icges(5,6);
t311 = Icges(7,3) + Icges(5,3) + Icges(6,2);
t310 = rSges(7,1) + pkin(5);
t309 = rSges(7,3) + qJ(6);
t241 = sin(qJ(1));
t244 = cos(qJ(1));
t284 = Icges(4,4) * t231;
t261 = -Icges(4,2) * t230 + t284;
t172 = -Icges(4,6) * t244 + t241 * t261;
t173 = Icges(4,6) * t241 + t244 * t261;
t285 = Icges(4,4) * t230;
t263 = Icges(4,1) * t231 - t285;
t174 = -Icges(4,5) * t244 + t241 * t263;
t175 = Icges(4,5) * t241 + t244 * t263;
t286 = Icges(3,4) * t243;
t262 = -Icges(3,2) * t240 + t286;
t181 = -Icges(3,6) * t244 + t241 * t262;
t182 = Icges(3,6) * t241 + t244 * t262;
t287 = Icges(3,4) * t240;
t264 = Icges(3,1) * t243 - t287;
t183 = -Icges(3,5) * t244 + t241 * t264;
t184 = Icges(3,5) * t241 + t244 * t264;
t201 = Icges(4,2) * t231 + t285;
t202 = Icges(4,1) * t230 + t284;
t214 = Icges(3,2) * t243 + t287;
t217 = Icges(3,1) * t240 + t286;
t226 = -qJD(2) * t244 + V_base(5);
t227 = qJD(2) * t241 + V_base(4);
t232 = V_base(6) + qJD(1);
t308 = (-t201 * t230 + t202 * t231 - t214 * t240 + t217 * t243) * t232 + (-t173 * t230 + t175 * t231 - t182 * t240 + t184 * t243) * t227 + (-t172 * t230 + t174 * t231 - t181 * t240 + t183 * t243) * t226;
t307 = (Icges(3,5) * t240 + Icges(4,5) * t230 + Icges(3,6) * t243 + Icges(4,6) * t231) * t232 + (t318 * t241 + t317 * t244) * t227 + (t317 * t241 - t318 * t244) * t226;
t242 = cos(qJ(4));
t278 = t242 * t244;
t239 = sin(qJ(4));
t280 = t241 * t239;
t194 = t231 * t280 + t278;
t279 = t241 * t242;
t281 = t239 * t244;
t195 = t231 * t279 - t281;
t283 = t230 * t241;
t306 = t312 * t194 + t314 * t195 - t311 * t283;
t196 = t231 * t281 - t279;
t197 = t231 * t278 + t280;
t282 = t230 * t244;
t305 = t312 * t196 + t314 * t197 - t311 * t282;
t304 = t313 * t194 + t315 * t195 - t312 * t283;
t303 = t313 * t196 + t315 * t197 - t312 * t282;
t302 = t315 * t194 + t316 * t195 - t314 * t283;
t301 = t315 * t196 + t316 * t197 - t314 * t282;
t300 = t311 * t231 + (t312 * t239 + t314 * t242) * t230;
t299 = t312 * t231 + (t313 * t239 + t315 * t242) * t230;
t298 = t314 * t231 + (t315 * t239 + t316 * t242) * t230;
t291 = pkin(2) * t240;
t290 = pkin(2) * t243;
t288 = Icges(2,4) * t241;
t277 = rSges(7,2) * t194 + t310 * t195 - t309 * t283;
t276 = t196 * rSges(7,2) + t310 * t197 - t309 * t282;
t275 = t309 * t231 + (rSges(7,2) * t239 + t310 * t242) * t230;
t167 = -qJ(3) * t244 + t241 * t290;
t224 = t241 * pkin(1) - pkin(7) * t244;
t274 = -t167 - t224;
t273 = qJD(4) * t230;
t272 = qJD(6) * t230;
t271 = V_base(5) * pkin(6) + V_base(1);
t268 = qJD(3) * t241 + t226 * t291 + t271;
t267 = pkin(3) * t231 + pkin(8) * t230;
t266 = rSges(3,1) * t243 - rSges(3,2) * t240;
t265 = rSges(4,1) * t231 - rSges(4,2) * t230;
t225 = pkin(1) * t244 + t241 * pkin(7);
t258 = -V_base(4) * pkin(6) + t232 * t225 + V_base(2);
t257 = V_base(4) * t224 - t225 * V_base(5) + V_base(3);
t256 = t227 * t167 + t257;
t168 = qJ(3) * t241 + t244 * t290;
t253 = -qJD(3) * t244 + t232 * t168 + t258;
t190 = t267 * t241;
t204 = pkin(3) * t230 - pkin(8) * t231;
t252 = t226 * t204 + (-t190 + t274) * t232 + t268;
t191 = t267 * t244;
t251 = t227 * t190 + (-t168 - t191) * t226 + t256;
t187 = (pkin(4) * t242 + qJ(5) * t239) * t230;
t192 = t241 * t273 + t226;
t250 = qJD(5) * t196 + t192 * t187 + t252;
t144 = pkin(4) * t195 + qJ(5) * t194;
t193 = t244 * t273 + t227;
t249 = qJD(5) * t230 * t239 + t193 * t144 + t251;
t248 = t232 * t191 + (-t204 - t291) * t227 + t253;
t145 = pkin(4) * t197 + qJ(5) * t196;
t207 = -qJD(4) * t231 + t232;
t247 = qJD(5) * t194 + t207 * t145 + t248;
t235 = Icges(2,4) * t244;
t223 = rSges(2,1) * t244 - t241 * rSges(2,2);
t222 = t241 * rSges(2,1) + rSges(2,2) * t244;
t221 = rSges(3,1) * t240 + rSges(3,2) * t243;
t219 = Icges(2,1) * t244 - t288;
t218 = Icges(2,1) * t241 + t235;
t216 = -Icges(2,2) * t241 + t235;
t215 = Icges(2,2) * t244 + t288;
t210 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t209 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t208 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t203 = rSges(4,1) * t230 + rSges(4,2) * t231;
t189 = t241 * rSges(3,3) + t244 * t266;
t188 = -rSges(3,3) * t244 + t241 * t266;
t177 = t241 * rSges(4,3) + t244 * t265;
t176 = -rSges(4,3) * t244 + t241 * t265;
t166 = -rSges(5,3) * t231 + (rSges(5,1) * t242 - rSges(5,2) * t239) * t230;
t165 = -rSges(6,2) * t231 + (rSges(6,1) * t242 + rSges(6,3) * t239) * t230;
t163 = V_base(5) * rSges(2,3) - t222 * t232 + t271;
t162 = t223 * t232 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t151 = t222 * V_base(4) - t223 * V_base(5) + V_base(3);
t143 = t197 * rSges(5,1) - t196 * rSges(5,2) + rSges(5,3) * t282;
t142 = t197 * rSges(6,1) + rSges(6,2) * t282 + t196 * rSges(6,3);
t140 = rSges(5,1) * t195 - rSges(5,2) * t194 + rSges(5,3) * t283;
t139 = rSges(6,1) * t195 + rSges(6,2) * t283 + rSges(6,3) * t194;
t117 = t221 * t226 + (-t188 - t224) * t232 + t271;
t116 = t189 * t232 - t221 * t227 + t258;
t115 = t188 * t227 - t189 * t226 + t257;
t114 = t203 * t226 + (-t176 + t274) * t232 + t268;
t113 = t232 * t177 + (-t203 - t291) * t227 + t253;
t112 = t176 * t227 + (-t168 - t177) * t226 + t256;
t111 = -t140 * t207 + t166 * t192 + t252;
t110 = t207 * t143 - t193 * t166 + t248;
t109 = t140 * t193 - t143 * t192 + t251;
t108 = t165 * t192 + (-t139 - t144) * t207 + t250;
t107 = t207 * t142 + (-t165 - t187) * t193 + t247;
t106 = t139 * t193 + (-t142 - t145) * t192 + t249;
t105 = -t244 * t272 + t275 * t192 + (-t144 - t277) * t207 + t250;
t104 = -t241 * t272 + t276 * t207 + (-t187 - t275) * t193 + t247;
t103 = qJD(6) * t231 + t277 * t193 + (-t145 - t276) * t192 + t249;
t1 = m(1) * (t208 ^ 2 + t209 ^ 2 + t210 ^ 2) / 0.2e1 + m(2) * (t151 ^ 2 + t162 ^ 2 + t163 ^ 2) / 0.2e1 + m(3) * (t115 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + m(5) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + m(4) * (t112 ^ 2 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + m(7) * (t103 ^ 2 + t104 ^ 2 + t105 ^ 2) / 0.2e1 + m(6) * (t106 ^ 2 + t107 ^ 2 + t108 ^ 2) / 0.2e1 + (t308 * t241 - t307 * t244) * t226 / 0.2e1 + (t307 * t241 + t308 * t244) * t227 / 0.2e1 + ((-t241 * t215 + t218 * t244 + Icges(1,4)) * V_base(5) + (-t241 * t216 + t219 * t244 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t215 * t244 + t241 * t218 + Icges(1,2)) * V_base(5) + (t216 * t244 + t241 * t219 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t299 * t194 + t298 * t195 - t300 * t283) * t207 + (t303 * t194 + t301 * t195 - t305 * t283) * t193 + (t304 * t194 + t302 * t195 - t306 * t283) * t192) * t192 / 0.2e1 + ((t299 * t196 + t298 * t197 - t300 * t282) * t207 + (t303 * t196 + t301 * t197 - t305 * t282) * t193 + (t304 * t196 + t302 * t197 - t306 * t282) * t192) * t193 / 0.2e1 + ((t306 * t192 + t305 * t193 + t300 * t207) * t231 + ((t299 * t239 + t298 * t242) * t207 + (t303 * t239 + t301 * t242) * t193 + (t304 * t239 + t302 * t242) * t192) * t230) * t207 / 0.2e1 + ((t173 * t231 + t175 * t230 + t182 * t243 + t184 * t240) * t227 + (t172 * t231 + t174 * t230 + t181 * t243 + t183 * t240) * t226 + (t201 * t231 + t202 * t230 + t214 * t243 + t217 * t240 + Icges(2,3)) * t232) * t232 / 0.2e1 + t232 * V_base(4) * (Icges(2,5) * t244 - Icges(2,6) * t241) + t232 * V_base(5) * (Icges(2,5) * t241 + Icges(2,6) * t244) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
