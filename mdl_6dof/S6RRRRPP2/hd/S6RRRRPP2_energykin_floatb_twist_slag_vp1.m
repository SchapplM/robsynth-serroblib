% Calculate kinetic energy for
% S6RRRRPP2
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
% Datum: 2019-03-09 20:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPP2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP2_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP2_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRPP2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP2_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP2_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP2_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPP2_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:48:52
% EndTime: 2019-03-09 20:48:55
% DurationCPUTime: 2.57s
% Computational Cost: add. (1875->284), mult. (2291->407), div. (0->0), fcn. (2209->8), ass. (0->145)
t314 = Icges(5,1) + Icges(6,1) + Icges(7,1);
t313 = -Icges(5,4) + Icges(7,4) + Icges(6,5);
t312 = Icges(7,5) - Icges(6,4) - Icges(5,5);
t311 = Icges(5,2) + Icges(7,2) + Icges(6,3);
t310 = -Icges(6,6) + Icges(7,6) + Icges(5,6);
t309 = Icges(7,3) + Icges(5,3) + Icges(6,2);
t308 = rSges(7,1) + pkin(5);
t307 = rSges(7,3) + qJ(6);
t239 = qJ(2) + qJ(3);
t236 = cos(t239);
t243 = cos(qJ(4));
t245 = cos(qJ(1));
t280 = t243 * t245;
t240 = sin(qJ(4));
t242 = sin(qJ(1));
t283 = t240 * t242;
t194 = t236 * t283 + t280;
t281 = t242 * t243;
t282 = t240 * t245;
t195 = t236 * t281 - t282;
t235 = sin(t239);
t285 = t235 * t242;
t306 = t310 * t194 + t312 * t195 - t309 * t285;
t196 = t236 * t282 - t281;
t197 = t236 * t280 + t283;
t284 = t235 * t245;
t305 = t310 * t196 + t312 * t197 - t309 * t284;
t304 = t311 * t194 + t313 * t195 - t310 * t285;
t303 = t311 * t196 + t313 * t197 - t310 * t284;
t302 = t313 * t194 + t314 * t195 - t312 * t285;
t301 = t313 * t196 + t314 * t197 - t312 * t284;
t300 = t309 * t236 + (t310 * t240 + t312 * t243) * t235;
t299 = t310 * t236 + (t311 * t240 + t313 * t243) * t235;
t298 = t312 * t236 + (t313 * t240 + t314 * t243) * t235;
t241 = sin(qJ(2));
t293 = pkin(2) * t241;
t244 = cos(qJ(2));
t292 = pkin(2) * t244;
t290 = Icges(2,4) * t242;
t289 = Icges(3,4) * t241;
t288 = Icges(3,4) * t244;
t287 = Icges(4,4) * t235;
t286 = Icges(4,4) * t236;
t279 = rSges(7,2) * t194 + t308 * t195 - t307 * t285;
t278 = rSges(7,2) * t196 + t308 * t197 - t307 * t284;
t277 = t307 * t236 + (rSges(7,2) * t240 + t308 * t243) * t235;
t168 = -pkin(8) * t245 + t242 * t292;
t226 = t242 * pkin(1) - t245 * pkin(7);
t276 = -t168 - t226;
t275 = qJD(4) * t235;
t274 = qJD(6) * t235;
t273 = V_base(5) * pkin(6) + V_base(1);
t229 = qJD(2) * t242 + V_base(4);
t232 = V_base(6) + qJD(1);
t228 = -qJD(2) * t245 + V_base(5);
t270 = t228 * t293 + t273;
t206 = qJD(3) * t242 + t229;
t269 = pkin(3) * t236 + pkin(9) * t235;
t268 = rSges(3,1) * t244 - rSges(3,2) * t241;
t267 = rSges(4,1) * t236 - rSges(4,2) * t235;
t266 = Icges(3,1) * t244 - t289;
t265 = Icges(4,1) * t236 - t287;
t264 = -Icges(3,2) * t241 + t288;
t263 = -Icges(4,2) * t235 + t286;
t262 = Icges(3,5) * t244 - Icges(3,6) * t241;
t261 = Icges(4,5) * t236 - Icges(4,6) * t235;
t227 = t245 * pkin(1) + t242 * pkin(7);
t260 = -V_base(4) * pkin(6) + t232 * t227 + V_base(2);
t259 = V_base(4) * t226 - t227 * V_base(5) + V_base(3);
t205 = V_base(5) + (-qJD(2) - qJD(3)) * t245;
t258 = (-Icges(4,3) * t245 + t242 * t261) * t205 + (Icges(4,3) * t242 + t245 * t261) * t206 + (Icges(4,5) * t235 + Icges(4,6) * t236) * t232;
t257 = (-Icges(3,3) * t245 + t242 * t262) * t228 + (Icges(3,3) * t242 + t245 * t262) * t229 + (Icges(3,5) * t241 + Icges(3,6) * t244) * t232;
t192 = t269 * t242;
t204 = pkin(3) * t235 - pkin(9) * t236;
t256 = t205 * t204 + (-t192 + t276) * t232 + t270;
t169 = pkin(8) * t242 + t245 * t292;
t255 = t229 * t168 - t169 * t228 + t259;
t254 = t232 * t169 - t229 * t293 + t260;
t179 = t242 * t275 + t205;
t191 = (pkin(4) * t243 + qJ(5) * t240) * t235;
t253 = qJD(5) * t196 + t179 * t191 + t256;
t193 = t269 * t245;
t252 = t206 * t192 - t193 * t205 + t255;
t251 = t232 * t193 - t204 * t206 + t254;
t145 = pkin(4) * t195 + qJ(5) * t194;
t180 = t245 * t275 + t206;
t250 = qJD(5) * t235 * t240 + t180 * t145 + t252;
t146 = pkin(4) * t197 + qJ(5) * t196;
t212 = -qJD(4) * t236 + t232;
t249 = qJD(5) * t194 + t212 * t146 + t251;
t173 = -Icges(4,6) * t245 + t242 * t263;
t174 = Icges(4,6) * t242 + t245 * t263;
t175 = -Icges(4,5) * t245 + t242 * t265;
t176 = Icges(4,5) * t242 + t245 * t265;
t201 = Icges(4,2) * t236 + t287;
t202 = Icges(4,1) * t235 + t286;
t248 = (-t174 * t235 + t176 * t236) * t206 + (-t173 * t235 + t175 * t236) * t205 + (-t201 * t235 + t202 * t236) * t232;
t183 = -Icges(3,6) * t245 + t242 * t264;
t184 = Icges(3,6) * t242 + t245 * t264;
t185 = -Icges(3,5) * t245 + t242 * t266;
t186 = Icges(3,5) * t242 + t245 * t266;
t216 = Icges(3,2) * t244 + t289;
t219 = Icges(3,1) * t241 + t288;
t247 = (-t184 * t241 + t186 * t244) * t229 + (-t183 * t241 + t185 * t244) * t228 + (-t216 * t241 + t219 * t244) * t232;
t237 = Icges(2,4) * t245;
t224 = rSges(2,1) * t245 - rSges(2,2) * t242;
t223 = rSges(2,1) * t242 + rSges(2,2) * t245;
t222 = rSges(3,1) * t241 + rSges(3,2) * t244;
t221 = Icges(2,1) * t245 - t290;
t220 = Icges(2,1) * t242 + t237;
t218 = -Icges(2,2) * t242 + t237;
t217 = Icges(2,2) * t245 + t290;
t211 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t210 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t209 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t203 = rSges(4,1) * t235 + rSges(4,2) * t236;
t190 = rSges(3,3) * t242 + t245 * t268;
t189 = -rSges(3,3) * t245 + t242 * t268;
t178 = rSges(4,3) * t242 + t245 * t267;
t177 = -rSges(4,3) * t245 + t242 * t267;
t167 = -rSges(5,3) * t236 + (rSges(5,1) * t243 - rSges(5,2) * t240) * t235;
t166 = -rSges(6,2) * t236 + (rSges(6,1) * t243 + rSges(6,3) * t240) * t235;
t155 = V_base(5) * rSges(2,3) - t223 * t232 + t273;
t154 = t224 * t232 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t152 = t223 * V_base(4) - t224 * V_base(5) + V_base(3);
t143 = rSges(5,1) * t197 - rSges(5,2) * t196 + rSges(5,3) * t284;
t142 = rSges(6,1) * t197 + rSges(6,2) * t284 + rSges(6,3) * t196;
t140 = rSges(5,1) * t195 - rSges(5,2) * t194 + rSges(5,3) * t285;
t139 = rSges(6,1) * t195 + rSges(6,2) * t285 + rSges(6,3) * t194;
t117 = t222 * t228 + (-t189 - t226) * t232 + t273;
t116 = t190 * t232 - t222 * t229 + t260;
t115 = t189 * t229 - t190 * t228 + t259;
t114 = t203 * t205 + (-t177 + t276) * t232 + t270;
t113 = t178 * t232 - t203 * t206 + t254;
t112 = t177 * t206 - t178 * t205 + t255;
t111 = -t140 * t212 + t167 * t179 + t256;
t110 = t143 * t212 - t167 * t180 + t251;
t109 = t140 * t180 - t143 * t179 + t252;
t108 = t166 * t179 + (-t139 - t145) * t212 + t253;
t107 = t142 * t212 + (-t166 - t191) * t180 + t249;
t106 = -t245 * t274 + t277 * t179 + (-t145 - t279) * t212 + t253;
t105 = -t242 * t274 + t278 * t212 + (-t191 - t277) * t180 + t249;
t104 = t139 * t180 + (-t142 - t146) * t179 + t250;
t103 = qJD(6) * t236 + t279 * t180 + (-t146 - t278) * t179 + t250;
t1 = t229 * (t242 * t257 + t245 * t247) / 0.2e1 + t228 * (t242 * t247 - t245 * t257) / 0.2e1 + t206 * (t242 * t258 + t245 * t248) / 0.2e1 + t205 * (t242 * t248 - t258 * t245) / 0.2e1 + m(1) * (t209 ^ 2 + t210 ^ 2 + t211 ^ 2) / 0.2e1 + m(2) * (t152 ^ 2 + t154 ^ 2 + t155 ^ 2) / 0.2e1 + m(4) * (t112 ^ 2 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + m(3) * (t115 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + m(7) * (t103 ^ 2 + t105 ^ 2 + t106 ^ 2) / 0.2e1 + m(6) * (t104 ^ 2 + t107 ^ 2 + t108 ^ 2) / 0.2e1 + m(5) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + ((-t217 * t242 + t220 * t245 + Icges(1,4)) * V_base(5) + (-t218 * t242 + t221 * t245 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t217 * t245 + t220 * t242 + Icges(1,2)) * V_base(5) + (t218 * t245 + t221 * t242 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t299 * t194 + t298 * t195 - t300 * t285) * t212 + (t303 * t194 + t301 * t195 - t305 * t285) * t180 + (t304 * t194 + t302 * t195 - t306 * t285) * t179) * t179 / 0.2e1 + ((t299 * t196 + t298 * t197 - t300 * t284) * t212 + (t303 * t196 + t301 * t197 - t305 * t284) * t180 + (t304 * t196 + t302 * t197 - t306 * t284) * t179) * t180 / 0.2e1 + ((t306 * t179 + t305 * t180 + t300 * t212) * t236 + ((t299 * t240 + t298 * t243) * t212 + (t303 * t240 + t301 * t243) * t180 + (t304 * t240 + t302 * t243) * t179) * t235) * t212 / 0.2e1 + ((t184 * t244 + t186 * t241) * t229 + (t183 * t244 + t185 * t241) * t228 + (t174 * t236 + t176 * t235) * t206 + (t173 * t236 + t175 * t235) * t205 + (t201 * t236 + t202 * t235 + t216 * t244 + t219 * t241 + Icges(2,3)) * t232) * t232 / 0.2e1 + t232 * V_base(4) * (Icges(2,5) * t245 - Icges(2,6) * t242) + V_base(5) * t232 * (Icges(2,5) * t242 + Icges(2,6) * t245) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
