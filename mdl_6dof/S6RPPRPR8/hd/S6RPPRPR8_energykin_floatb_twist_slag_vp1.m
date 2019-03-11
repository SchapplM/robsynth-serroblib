% Calculate kinetic energy for
% S6RPPRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
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
% Datum: 2019-03-09 01:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRPR8_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR8_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPPRPR8_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR8_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR8_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRPR8_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:55:14
% EndTime: 2019-03-09 01:55:17
% DurationCPUTime: 3.68s
% Computational Cost: add. (1156->293), mult. (1442->391), div. (0->0), fcn. (1224->8), ass. (0->154)
t323 = -Icges(5,4) - Icges(6,6);
t322 = Icges(5,1) + Icges(6,2);
t321 = Icges(5,2) + Icges(6,3);
t211 = pkin(9) + qJ(4);
t201 = cos(t211);
t320 = t323 * t201;
t200 = sin(t211);
t319 = t323 * t200;
t318 = -Icges(6,4) + Icges(5,5);
t317 = Icges(6,5) - Icges(5,6);
t316 = t321 * t201 - t319;
t315 = t322 * t200 - t320;
t314 = Icges(2,4) + Icges(3,6);
t313 = Icges(2,1) + Icges(3,2);
t312 = Icges(6,1) + Icges(5,3);
t311 = -Icges(3,4) + Icges(2,5);
t310 = Icges(3,5) - Icges(2,6);
t309 = Icges(2,2) + Icges(3,3);
t216 = sin(qJ(1));
t218 = cos(qJ(1));
t308 = -t316 * t216 + t317 * t218;
t307 = t317 * t216 + t316 * t218;
t306 = t315 * t216 + t318 * t218;
t305 = t318 * t216 - t315 * t218;
t304 = t321 * t200 + t320;
t303 = t322 * t201 + t319;
t302 = t318 * t200 - t317 * t201;
t301 = t314 * t218;
t300 = t314 * t216;
t195 = qJD(4) * t216 + V_base(5);
t196 = qJD(4) * t218 + V_base(4);
t202 = V_base(6) + qJD(1);
t299 = (-t305 * t200 + t307 * t201) * t195 + (-t306 * t200 + t308 * t201) * t196 - (t303 * t200 - t304 * t201) * t202;
t298 = t313 * t216 + t301;
t297 = t313 * t218 - t300;
t296 = t311 * t216 - t310 * t218;
t295 = t310 * t216 + t311 * t218;
t294 = (t317 * t200 + t318 * t201) * t202 + (t302 * t216 + t312 * t218) * t196 + (t312 * t216 - t302 * t218) * t195;
t213 = cos(pkin(9));
t212 = sin(pkin(9));
t281 = Icges(4,4) * t212;
t247 = Icges(4,2) * t213 + t281;
t139 = Icges(4,6) * t216 - t218 * t247;
t280 = Icges(4,4) * t213;
t249 = Icges(4,1) * t212 + t280;
t141 = Icges(4,5) * t216 - t218 * t249;
t290 = t139 * t213 + t141 * t212 - t309 * t218 - t300;
t138 = Icges(4,6) * t218 + t216 * t247;
t140 = Icges(4,5) * t218 + t216 * t249;
t289 = t138 * t213 + t140 * t212 + t309 * t216 - t301;
t288 = -pkin(2) - pkin(6);
t284 = pkin(3) * t212;
t283 = pkin(3) * t213;
t273 = qJ(3) * t218;
t272 = t200 * t216;
t271 = t200 * t218;
t215 = sin(qJ(6));
t270 = t215 * t218;
t269 = t216 * qJ(3);
t268 = t216 * t215;
t217 = cos(qJ(6));
t267 = t216 * t217;
t266 = t217 * t218;
t264 = qJD(5) * t216;
t263 = qJD(6) * t200;
t188 = t216 * pkin(1) - qJ(2) * t218;
t262 = V_base(4) * t188 + V_base(3);
t261 = V_base(5) * pkin(6) + V_base(1);
t258 = V_base(4) * t269 + t262;
t257 = qJD(2) * t216 + t261;
t256 = -t188 - t269;
t191 = pkin(1) * t218 + t216 * qJ(2);
t255 = -t191 - t273;
t146 = pkin(7) * t216 - t218 * t284;
t254 = -t146 + t256;
t253 = rSges(4,1) * t212 + rSges(4,2) * t213;
t252 = rSges(5,1) * t200 + rSges(5,2) * t201;
t251 = rSges(6,2) * t200 + rSges(6,3) * t201;
t250 = pkin(4) * t200 - qJ(5) * t201;
t244 = Icges(4,5) * t212 + Icges(4,6) * t213;
t172 = -Icges(4,2) * t212 + t280;
t173 = Icges(4,1) * t213 - t281;
t232 = t172 * t213 + t173 * t212;
t231 = -qJD(2) * t218 + t202 * t191 + V_base(2);
t230 = V_base(5) * pkin(2) + qJD(3) * t218 + t257;
t145 = t250 * t218;
t229 = t145 + t254;
t228 = V_base(5) * t283 + t230;
t227 = qJD(3) * t216 + t202 * t273 + t231;
t162 = pkin(4) * t201 + qJ(5) * t200;
t226 = t195 * t162 + t228;
t223 = (Icges(4,3) * t218 + t216 * t244) * V_base(4) + (Icges(4,3) * t216 - t218 * t244) * V_base(5) + (Icges(4,5) * t213 - Icges(4,6) * t212) * t202;
t147 = pkin(7) * t218 + t216 * t284;
t222 = V_base(4) * t146 + (-t147 + t255) * V_base(5) + t258;
t221 = qJD(5) * t200 - t196 * t145 + t222;
t220 = t202 * t147 + (-t283 + t288) * V_base(4) + t227;
t144 = t250 * t216;
t219 = qJD(5) * t218 * t201 + t202 * t144 + t220;
t193 = rSges(2,1) * t218 - t216 * rSges(2,2);
t192 = -rSges(3,2) * t218 + t216 * rSges(3,3);
t190 = t216 * rSges(2,1) + rSges(2,2) * t218;
t189 = -t216 * rSges(3,2) - rSges(3,3) * t218;
t174 = rSges(4,1) * t213 - rSges(4,2) * t212;
t170 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t169 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t168 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t167 = qJD(6) * t201 + t202;
t165 = pkin(5) * t218 + pkin(8) * t272;
t164 = rSges(5,1) * t201 - rSges(5,2) * t200;
t163 = -rSges(6,2) * t201 + rSges(6,3) * t200;
t161 = t216 * pkin(5) - pkin(8) * t271;
t153 = -t201 * t268 + t266;
t152 = -t201 * t267 - t270;
t151 = t201 * t270 + t267;
t150 = t201 * t266 - t268;
t149 = t216 * t263 + t196;
t148 = -t218 * t263 + t195;
t143 = t216 * rSges(4,3) - t218 * t253;
t142 = rSges(4,3) * t218 + t216 * t253;
t133 = rSges(6,1) * t218 - t216 * t251;
t132 = t216 * rSges(6,1) + t218 * t251;
t131 = t216 * rSges(5,3) - t218 * t252;
t130 = rSges(5,3) * t218 + t216 * t252;
t115 = rSges(7,3) * t201 + (rSges(7,1) * t215 + rSges(7,2) * t217) * t200;
t114 = V_base(5) * rSges(2,3) - t190 * t202 + t261;
t113 = t193 * t202 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t112 = Icges(7,5) * t201 + (Icges(7,1) * t215 + Icges(7,4) * t217) * t200;
t111 = Icges(7,6) * t201 + (Icges(7,4) * t215 + Icges(7,2) * t217) * t200;
t110 = Icges(7,3) * t201 + (Icges(7,5) * t215 + Icges(7,6) * t217) * t200;
t108 = t190 * V_base(4) - t193 * V_base(5) + V_base(3);
t107 = V_base(5) * rSges(3,1) + (-t188 - t189) * t202 + t257;
t106 = t202 * t192 + (-rSges(3,1) - pkin(6)) * V_base(4) + t231;
t105 = rSges(7,1) * t153 + rSges(7,2) * t152 + rSges(7,3) * t272;
t104 = t151 * rSges(7,1) + t150 * rSges(7,2) - rSges(7,3) * t271;
t103 = Icges(7,1) * t153 + Icges(7,4) * t152 + Icges(7,5) * t272;
t102 = Icges(7,1) * t151 + Icges(7,4) * t150 - Icges(7,5) * t271;
t101 = Icges(7,4) * t153 + Icges(7,2) * t152 + Icges(7,6) * t272;
t100 = Icges(7,4) * t151 + Icges(7,2) * t150 - Icges(7,6) * t271;
t99 = Icges(7,5) * t153 + Icges(7,6) * t152 + Icges(7,3) * t272;
t98 = Icges(7,5) * t151 + Icges(7,6) * t150 - Icges(7,3) * t271;
t97 = t189 * V_base(4) + (-t191 - t192) * V_base(5) + t262;
t96 = t174 * V_base(5) + (-t143 + t256) * t202 + t230;
t95 = t202 * t142 + (-t174 + t288) * V_base(4) + t227;
t94 = V_base(4) * t143 + (-t142 + t255) * V_base(5) + t258;
t93 = t164 * t195 + (-t131 + t254) * t202 + t228;
t92 = t202 * t130 - t196 * t164 + t220;
t91 = -t195 * t130 + t196 * t131 + t222;
t90 = -t201 * t264 + t163 * t195 + (-t132 + t229) * t202 + t226;
t89 = t202 * t133 + (-t162 - t163) * t196 + t219;
t88 = t196 * t132 + (-t133 - t144) * t195 + t221;
t87 = -t104 * t167 + t115 * t148 + (pkin(8) * t195 - t264) * t201 + (-t161 + t229) * t202 + t226;
t86 = t167 * t105 - t149 * t115 + t202 * t165 + (-pkin(8) * t201 - t162) * t196 + t219;
t85 = t149 * t104 - t148 * t105 + t196 * t161 + (-t144 - t165) * t195 + t221;
t1 = m(1) * (t168 ^ 2 + t169 ^ 2 + t170 ^ 2) / 0.2e1 + m(2) * (t108 ^ 2 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + m(3) * (t106 ^ 2 + t107 ^ 2 + t97 ^ 2) / 0.2e1 + m(6) * (t88 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + m(5) * (t91 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + m(4) * (t94 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + m(7) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + t148 * ((t150 * t101 + t151 * t103 - t271 * t99) * t149 + (t150 * t100 + t151 * t102 - t98 * t271) * t148 + (-t110 * t271 + t150 * t111 + t151 * t112) * t167) / 0.2e1 + t149 * ((t152 * t101 + t153 * t103 + t99 * t272) * t149 + (t100 * t152 + t102 * t153 + t272 * t98) * t148 + (t110 * t272 + t111 * t152 + t112 * t153) * t167) / 0.2e1 + t167 * ((t110 * t167 + t98 * t148 + t99 * t149) * t201 + ((t101 * t217 + t103 * t215) * t149 + (t100 * t217 + t102 * t215) * t148 + (t111 * t217 + t112 * t215) * t167) * t200) / 0.2e1 + (t294 * t216 + t299 * t218) * t195 / 0.2e1 + (-t299 * t216 + t294 * t218) * t196 / 0.2e1 + (t223 * t218 + (t232 * t216 + t295) * t202 + (t290 * t216 + t298 * t218 + Icges(1,4)) * V_base(5) + (t289 * t216 + t297 * t218 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t223 * t216 + (-t232 * t218 + t296) * t202 + (t298 * t216 - t290 * t218 + Icges(1,2)) * V_base(5) + (t297 * t216 - t289 * t218 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t308 * t200 + t306 * t201) * t196 + (t307 * t200 + t305 * t201) * t195 + (-t139 * t212 + t141 * t213 + t296) * V_base(5) + (-t138 * t212 + t140 * t213 + t295) * V_base(4) + (-t212 * t172 + t213 * t173 + t304 * t200 + t303 * t201 + Icges(3,1) + Icges(2,3)) * t202) * t202 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
