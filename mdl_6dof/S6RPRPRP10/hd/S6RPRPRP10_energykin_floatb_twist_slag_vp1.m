% Calculate kinetic energy for
% S6RPRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
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
% Datum: 2019-03-09 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRP10_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP10_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPRP10_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP10_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP10_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRP10_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:30:52
% EndTime: 2019-03-09 03:30:55
% DurationCPUTime: 3.21s
% Computational Cost: add. (931->260), mult. (1727->345), div. (0->0), fcn. (1575->6), ass. (0->135)
t332 = -Icges(4,4) - Icges(5,6);
t331 = Icges(4,1) + Icges(5,2);
t330 = Icges(4,2) + Icges(5,3);
t220 = cos(qJ(3));
t329 = t332 * t220;
t217 = sin(qJ(3));
t328 = t332 * t217;
t327 = -Icges(5,4) + Icges(4,5);
t326 = Icges(5,5) - Icges(4,6);
t325 = t330 * t220 - t328;
t324 = t331 * t217 - t329;
t323 = Icges(2,4) + Icges(3,6);
t322 = Icges(2,1) + Icges(3,2);
t321 = Icges(5,1) + Icges(4,3);
t320 = Icges(6,1) + Icges(7,1);
t319 = -Icges(3,4) + Icges(2,5);
t318 = Icges(6,4) - Icges(7,5);
t317 = Icges(7,4) + Icges(6,5);
t316 = Icges(3,5) - Icges(2,6);
t315 = Icges(2,2) + Icges(3,3);
t314 = Icges(6,2) + Icges(7,3);
t313 = Icges(7,6) - Icges(6,6);
t312 = Icges(6,3) + Icges(7,2);
t218 = sin(qJ(1));
t221 = cos(qJ(1));
t311 = -t325 * t218 + t326 * t221;
t310 = t326 * t218 + t325 * t221;
t309 = t324 * t218 + t327 * t221;
t308 = t327 * t218 - t324 * t221;
t307 = t330 * t217 + t329;
t306 = t331 * t220 + t328;
t305 = t327 * t217 - t326 * t220;
t304 = rSges(7,1) + pkin(5);
t303 = rSges(7,3) + qJ(6);
t302 = t323 * t221;
t301 = t323 * t218;
t204 = qJD(3) * t218 + V_base(5);
t205 = qJD(3) * t221 + V_base(4);
t208 = V_base(6) + qJD(1);
t300 = t204 * (-t308 * t217 + t310 * t220) + (-t309 * t217 + t311 * t220) * t205 - (t306 * t217 - t307 * t220) * t208;
t216 = sin(qJ(5));
t219 = cos(qJ(5));
t261 = t220 * t221;
t159 = t216 * t218 - t219 * t261;
t160 = t216 * t261 + t218 * t219;
t263 = t217 * t221;
t299 = t314 * t159 - t318 * t160 - t313 * t263;
t262 = t218 * t220;
t161 = t216 * t221 + t219 * t262;
t162 = -t216 * t262 + t219 * t221;
t264 = t217 * t218;
t298 = t314 * t161 - t318 * t162 + t313 * t264;
t297 = t313 * t159 + t317 * t160 - t312 * t263;
t296 = t313 * t161 + t317 * t162 + t312 * t264;
t295 = -t318 * t159 + t320 * t160 - t317 * t263;
t294 = -t318 * t161 + t320 * t162 + t317 * t264;
t293 = t313 * t220 + (-t318 * t216 - t314 * t219) * t217;
t292 = t312 * t220 + (t317 * t216 - t313 * t219) * t217;
t291 = t317 * t220 + (t320 * t216 + t318 * t219) * t217;
t290 = -t315 * t221 - t301;
t289 = t315 * t218 - t302;
t288 = t322 * t218 + t302;
t287 = t322 * t221 - t301;
t284 = (t326 * t217 + t327 * t220) * t208 + (t305 * t218 + t321 * t221) * t205 + (t321 * t218 - t305 * t221) * t204;
t274 = pkin(7) * t218;
t273 = pkin(7) * t221;
t272 = pkin(8) * t220;
t260 = -rSges(7,2) * t263 + t303 * t159 + t160 * t304;
t259 = rSges(7,2) * t264 + t303 * t161 + t162 * t304;
t258 = rSges(7,2) * t220 + (t216 * t304 - t303 * t219) * t217;
t257 = qJD(4) * t220;
t256 = qJD(5) * t217;
t193 = t218 * pkin(1) - qJ(2) * t221;
t255 = V_base(4) * t193 + V_base(3);
t254 = V_base(5) * pkin(6) + V_base(1);
t251 = -t193 - t274;
t250 = qJD(2) * t218 + t254;
t245 = pkin(3) * t217 - qJ(4) * t220;
t166 = t245 * t221;
t249 = t166 + t251;
t248 = V_base(5) * pkin(2) + t250;
t247 = rSges(4,1) * t217 + rSges(4,2) * t220;
t246 = rSges(5,2) * t217 + rSges(5,3) * t220;
t199 = pkin(1) * t221 + t218 * qJ(2);
t232 = -qJD(2) * t221 + t208 * t199 + V_base(2);
t196 = pkin(3) * t220 + qJ(4) * t217;
t229 = t204 * t196 - t218 * t257 + t248;
t228 = V_base(4) * t274 + (-t199 - t273) * V_base(5) + t255;
t227 = t208 * t273 + (-pkin(2) - pkin(6)) * V_base(4) + t232;
t226 = qJD(4) * t217 - t205 * t166 + t228;
t165 = t245 * t218;
t225 = t208 * t165 + t221 * t257 + t227;
t169 = t218 * pkin(4) - pkin(8) * t263;
t170 = pkin(4) * t221 + pkin(8) * t264;
t224 = t205 * t169 + (-t165 - t170) * t204 + t226;
t223 = t204 * t272 + (-t169 + t249) * t208 + t229;
t222 = t208 * t170 + (-t196 - t272) * t205 + t225;
t201 = rSges(2,1) * t221 - t218 * rSges(2,2);
t200 = -rSges(3,2) * t221 + t218 * rSges(3,3);
t198 = rSges(4,1) * t220 - rSges(4,2) * t217;
t197 = -rSges(5,2) * t220 + rSges(5,3) * t217;
t195 = t218 * rSges(2,1) + rSges(2,2) * t221;
t194 = -t218 * rSges(3,2) - rSges(3,3) * t221;
t192 = qJD(5) * t220 + t208;
t173 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t172 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t171 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t158 = t218 * t256 + t205;
t157 = -t221 * t256 + t204;
t154 = rSges(5,1) * t221 - t218 * t246;
t153 = t218 * rSges(5,1) + t221 * t246;
t152 = t218 * rSges(4,3) - t221 * t247;
t151 = rSges(4,3) * t221 + t218 * t247;
t150 = rSges(6,3) * t220 + (rSges(6,1) * t216 + rSges(6,2) * t219) * t217;
t127 = V_base(5) * rSges(2,3) - t195 * t208 + t254;
t126 = t201 * t208 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t125 = t195 * V_base(4) - t201 * V_base(5) + V_base(3);
t122 = rSges(6,1) * t162 - rSges(6,2) * t161 + rSges(6,3) * t264;
t120 = t160 * rSges(6,1) - t159 * rSges(6,2) - rSges(6,3) * t263;
t106 = V_base(5) * rSges(3,1) + (-t193 - t194) * t208 + t250;
t105 = t208 * t200 + (-rSges(3,1) - pkin(6)) * V_base(4) + t232;
t104 = t194 * V_base(4) + (-t199 - t200) * V_base(5) + t255;
t103 = t198 * t204 + (-t152 + t251) * t208 + t248;
t102 = t208 * t151 - t205 * t198 + t227;
t101 = -t204 * t151 + t205 * t152 + t228;
t100 = t197 * t204 + (-t153 + t249) * t208 + t229;
t99 = t208 * t154 + (-t196 - t197) * t205 + t225;
t98 = t205 * t153 + (-t154 - t165) * t204 + t226;
t97 = -t120 * t192 + t150 * t157 + t223;
t96 = t192 * t122 - t158 * t150 + t222;
t95 = t158 * t120 - t157 * t122 + t224;
t94 = qJD(6) * t161 + t157 * t258 - t192 * t260 + t223;
t93 = qJD(6) * t159 - t158 * t258 + t192 * t259 + t222;
t92 = -qJD(6) * t217 * t219 - t157 * t259 + t158 * t260 + t224;
t1 = m(1) * (t171 ^ 2 + t172 ^ 2 + t173 ^ 2) / 0.2e1 + m(2) * (t125 ^ 2 + t126 ^ 2 + t127 ^ 2) / 0.2e1 + m(5) * (t100 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(4) * (t101 ^ 2 + t102 ^ 2 + t103 ^ 2) / 0.2e1 + m(3) * (t104 ^ 2 + t105 ^ 2 + t106 ^ 2) / 0.2e1 + m(7) * (t92 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + m(6) * (t95 ^ 2 + t96 ^ 2 + t97 ^ 2) / 0.2e1 + ((t159 * t293 + t160 * t291 - t263 * t292) * t192 + (t159 * t298 + t160 * t294 - t263 * t296) * t158 + (t299 * t159 + t295 * t160 - t297 * t263) * t157) * t157 / 0.2e1 + ((t161 * t293 + t162 * t291 + t264 * t292) * t192 + (t298 * t161 + t294 * t162 + t296 * t264) * t158 + (t161 * t299 + t295 * t162 + t297 * t264) * t157) * t158 / 0.2e1 + ((t157 * t297 + t158 * t296 + t192 * t292) * t220 + ((t216 * t291 - t219 * t293) * t192 + (t216 * t294 - t219 * t298) * t158 + (t295 * t216 - t219 * t299) * t157) * t217) * t192 / 0.2e1 + (t284 * t218 + t300 * t221) * t204 / 0.2e1 + (-t300 * t218 + t284 * t221) * t205 / 0.2e1 + ((t218 * t290 + t221 * t288 + Icges(1,4)) * V_base(5) + (t289 * t218 + t287 * t221 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t288 * t218 - t290 * t221 + Icges(1,2)) * V_base(5) + (t218 * t287 - t289 * t221 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t311 * t217 + t309 * t220) * t205 + (t310 * t217 + t308 * t220) * t204 + (t307 * t217 + t306 * t220 + Icges(3,1) + Icges(2,3)) * t208) * t208 / 0.2e1 + t208 * V_base(5) * (t319 * t218 - t316 * t221) + t208 * V_base(4) * (t316 * t218 + t319 * t221) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
