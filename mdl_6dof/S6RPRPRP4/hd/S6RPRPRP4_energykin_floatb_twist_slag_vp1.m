% Calculate kinetic energy for
% S6RPRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
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
% Datum: 2019-03-09 03:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRP4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP4_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPRP4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP4_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP4_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRP4_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:11:07
% EndTime: 2019-03-09 03:11:10
% DurationCPUTime: 3.16s
% Computational Cost: add. (1598->273), mult. (1715->361), div. (0->0), fcn. (1561->8), ass. (0->140)
t320 = Icges(4,4) + Icges(5,6);
t319 = Icges(4,1) + Icges(5,2);
t318 = -Icges(4,2) - Icges(5,3);
t226 = cos(qJ(3));
t317 = t320 * t226;
t223 = sin(qJ(3));
t316 = t320 * t223;
t315 = Icges(5,4) - Icges(4,5);
t314 = Icges(5,5) - Icges(4,6);
t313 = t318 * t223 + t317;
t312 = t319 * t226 - t316;
t311 = Icges(5,1) + Icges(4,3);
t310 = Icges(6,1) + Icges(7,1);
t309 = Icges(6,4) - Icges(7,5);
t308 = Icges(7,4) + Icges(6,5);
t307 = Icges(6,2) + Icges(7,3);
t306 = Icges(7,6) - Icges(6,6);
t305 = Icges(6,3) + Icges(7,2);
t221 = qJ(1) + pkin(9);
t215 = sin(t221);
t216 = cos(t221);
t304 = t313 * t215 + t314 * t216;
t303 = -t314 * t215 + t313 * t216;
t302 = t312 * t215 + t315 * t216;
t301 = -t315 * t215 + t312 * t216;
t300 = t318 * t226 - t316;
t299 = t319 * t223 + t317;
t298 = t314 * t223 - t315 * t226;
t297 = rSges(7,1) + pkin(5);
t296 = rSges(7,3) + qJ(6);
t222 = sin(qJ(5));
t225 = cos(qJ(5));
t263 = t223 * t225;
t155 = t215 * t222 - t216 * t263;
t264 = t222 * t223;
t156 = t215 * t225 + t216 * t264;
t265 = t216 * t226;
t295 = t155 * t307 - t156 * t309 + t265 * t306;
t157 = t215 * t263 + t216 * t222;
t158 = t215 * t264 - t216 * t225;
t266 = t215 * t226;
t294 = -t157 * t307 - t158 * t309 + t266 * t306;
t293 = t155 * t306 + t156 * t308 + t265 * t305;
t292 = -t157 * t306 + t158 * t308 + t266 * t305;
t291 = -t155 * t309 + t156 * t310 + t265 * t308;
t290 = t157 * t309 + t158 * t310 + t266 * t308;
t289 = (t222 * t309 + t225 * t307) * t226 + t306 * t223;
t288 = (-t222 * t308 + t225 * t306) * t226 + t305 * t223;
t287 = (-t222 * t310 - t225 * t309) * t226 + t308 * t223;
t189 = -qJD(3) * t216 + V_base(5);
t190 = qJD(3) * t215 + V_base(4);
t217 = V_base(6) + qJD(1);
t284 = (t223 * t300 + t226 * t299) * t217 + (-t223 * t303 + t226 * t301) * t190 + (-t223 * t304 + t226 * t302) * t189;
t283 = (-t315 * t223 - t314 * t226) * t217 + (t215 * t311 + t298 * t216) * t190 + (t298 * t215 - t216 * t311) * t189;
t224 = sin(qJ(1));
t276 = pkin(1) * t224;
t227 = cos(qJ(1));
t275 = pkin(1) * t227;
t274 = pkin(8) * t223;
t273 = -pkin(6) - qJ(2);
t272 = Icges(2,4) * t224;
t271 = Icges(3,4) * t215;
t262 = rSges(7,2) * t265 + t296 * t155 + t297 * t156;
t261 = rSges(7,2) * t266 - t296 * t157 + t297 * t158;
t260 = rSges(7,2) * t223 + (-t297 * t222 + t296 * t225) * t226;
t259 = qJD(4) * t223;
t258 = qJD(4) * t226;
t257 = qJD(5) * t226;
t256 = t217 * t275 + V_base(2);
t255 = V_base(5) * pkin(6) + V_base(1);
t183 = pkin(2) * t215 - pkin(7) * t216;
t252 = -t183 - t276;
t251 = V_base(5) * qJ(2) + t255;
t250 = V_base(4) * t276 + qJD(2) + V_base(3);
t246 = pkin(3) * t226 + qJ(4) * t223;
t165 = t246 * t215;
t249 = -t165 + t252;
t248 = rSges(4,1) * t226 - rSges(4,2) * t223;
t247 = -rSges(5,2) * t226 + rSges(5,3) * t223;
t206 = pkin(3) * t223 - qJ(4) * t226;
t239 = t189 * t206 + t216 * t259 + t251;
t184 = pkin(2) * t216 + pkin(7) * t215;
t236 = t217 * t184 + t273 * V_base(4) + t256;
t166 = t246 * t216;
t235 = t217 * t166 + t215 * t259 + t236;
t234 = V_base(4) * t183 + (-t184 - t275) * V_base(5) + t250;
t233 = t190 * t165 + t234;
t172 = -pkin(4) * t216 + pkin(8) * t266;
t232 = t189 * t274 + (-t172 + t249) * t217 + t239;
t171 = pkin(4) * t215 + pkin(8) * t265;
t231 = t217 * t171 + (-t206 - t274) * t190 + t235;
t230 = t190 * t172 + (-t166 - t171) * t189 + t233;
t219 = Icges(2,4) * t227;
t214 = Icges(3,4) * t216;
t210 = rSges(2,1) * t227 - t224 * rSges(2,2);
t209 = t224 * rSges(2,1) + rSges(2,2) * t227;
t208 = rSges(4,1) * t223 + rSges(4,2) * t226;
t207 = -rSges(5,2) * t223 - rSges(5,3) * t226;
t205 = qJD(5) * t223 + t217;
t202 = Icges(2,1) * t227 - t272;
t201 = Icges(2,1) * t224 + t219;
t199 = -Icges(2,2) * t224 + t219;
t198 = Icges(2,2) * t227 + t272;
t187 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t186 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t185 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t182 = rSges(3,1) * t216 - rSges(3,2) * t215;
t181 = rSges(3,1) * t215 + rSges(3,2) * t216;
t180 = Icges(3,1) * t216 - t271;
t179 = Icges(3,1) * t215 + t214;
t178 = -Icges(3,2) * t215 + t214;
t177 = Icges(3,2) * t216 + t271;
t168 = rSges(6,3) * t223 + (-rSges(6,1) * t222 - rSges(6,2) * t225) * t226;
t154 = t216 * t257 + t190;
t153 = t215 * t257 + t189;
t149 = -rSges(5,1) * t216 + t215 * t247;
t148 = rSges(5,1) * t215 + t216 * t247;
t147 = rSges(4,3) * t215 + t216 * t248;
t146 = -rSges(4,3) * t216 + t215 * t248;
t145 = V_base(5) * rSges(2,3) - t209 * t217 + t255;
t144 = t210 * t217 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t130 = t209 * V_base(4) - t210 * V_base(5) + V_base(3);
t126 = V_base(5) * rSges(3,3) + (-t181 - t276) * t217 + t251;
t125 = t182 * t217 + (-rSges(3,3) + t273) * V_base(4) + t256;
t124 = V_base(4) * t181 + (-t182 - t275) * V_base(5) + t250;
t123 = rSges(6,1) * t158 + rSges(6,2) * t157 + rSges(6,3) * t266;
t121 = rSges(6,1) * t156 - rSges(6,2) * t155 + rSges(6,3) * t265;
t107 = t189 * t208 + (-t146 + t252) * t217 + t251;
t106 = t147 * t217 - t190 * t208 + t236;
t105 = t190 * t146 - t189 * t147 + t234;
t104 = t189 * t207 + (-t149 + t249) * t217 + t239;
t103 = t148 * t217 + (-t206 - t207) * t190 + t235;
t102 = -t258 + t190 * t149 + (-t148 - t166) * t189 + t233;
t101 = -t123 * t205 + t153 * t168 + t232;
t100 = t121 * t205 - t154 * t168 + t231;
t99 = -t153 * t121 + t154 * t123 + t230 - t258;
t98 = qJD(6) * t155 + t153 * t260 - t205 * t261 + t232;
t97 = -qJD(6) * t157 - t154 * t260 + t205 * t262 + t231;
t96 = (qJD(6) * t225 - qJD(4)) * t226 + t261 * t154 - t262 * t153 + t230;
t1 = m(1) * (t185 ^ 2 + t186 ^ 2 + t187 ^ 2) / 0.2e1 + m(2) * (t130 ^ 2 + t144 ^ 2 + t145 ^ 2) / 0.2e1 + m(3) * (t124 ^ 2 + t125 ^ 2 + t126 ^ 2) / 0.2e1 + m(4) * (t105 ^ 2 + t106 ^ 2 + t107 ^ 2) / 0.2e1 + m(7) * (t96 ^ 2 + t97 ^ 2 + t98 ^ 2) / 0.2e1 + m(6) * (t100 ^ 2 + t101 ^ 2 + t99 ^ 2) / 0.2e1 + m(5) * (t102 ^ 2 + t103 ^ 2 + t104 ^ 2) / 0.2e1 + ((-t289 * t157 + t287 * t158 + t288 * t266) * t205 + (-t295 * t157 + t291 * t158 + t293 * t266) * t154 + (-t294 * t157 + t290 * t158 + t292 * t266) * t153) * t153 / 0.2e1 + ((t289 * t155 + t287 * t156 + t288 * t265) * t205 + (t295 * t155 + t291 * t156 + t293 * t265) * t154 + (t294 * t155 + t290 * t156 + t292 * t265) * t153) * t154 / 0.2e1 + (t284 * t215 - t283 * t216) * t189 / 0.2e1 + (t283 * t215 + t284 * t216) * t190 / 0.2e1 + (((-t287 * t222 + t289 * t225) * t205 + (-t291 * t222 + t295 * t225) * t154 + (-t290 * t222 + t294 * t225) * t153) * t226 + (t292 * t153 + t293 * t154 + t288 * t205) * t223) * t205 / 0.2e1 + ((-t177 * t215 + t179 * t216 - t224 * t198 + t201 * t227 + Icges(1,4)) * V_base(5) + (-t178 * t215 + t180 * t216 - t224 * t199 + t202 * t227 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t177 * t216 + t179 * t215 + t198 * t227 + t224 * t201 + Icges(1,2)) * V_base(5) + (t178 * t216 + t180 * t215 + t199 * t227 + t224 * t202 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t223 * t301 + t226 * t303) * t190 + (t223 * t302 + t226 * t304) * t189 + (t299 * t223 - t300 * t226 + Icges(2,3) + Icges(3,3)) * t217) * t217 / 0.2e1 + t217 * V_base(4) * (Icges(2,5) * t227 + Icges(3,5) * t216 - Icges(2,6) * t224 - Icges(3,6) * t215) + t217 * V_base(5) * (Icges(2,5) * t224 + Icges(3,5) * t215 + Icges(2,6) * t227 + Icges(3,6) * t216) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
