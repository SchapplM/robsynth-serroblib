% Calculate kinetic energy for
% S6RPRPRP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-03-09 03:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRP7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP7_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPRP7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP7_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP7_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRP7_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:21:07
% EndTime: 2019-03-09 03:21:10
% DurationCPUTime: 2.75s
% Computational Cost: add. (1391->288), mult. (1733->379), div. (0->0), fcn. (1571->8), ass. (0->148)
t323 = Icges(2,4) + Icges(3,6);
t322 = Icges(2,1) + Icges(3,2);
t321 = Icges(6,1) + Icges(7,1);
t320 = -Icges(3,4) + Icges(2,5);
t319 = Icges(6,4) + Icges(7,4);
t318 = Icges(3,5) - Icges(2,6);
t317 = Icges(7,5) + Icges(6,5);
t316 = Icges(2,2) + Icges(3,3);
t315 = Icges(6,2) + Icges(7,2);
t314 = Icges(7,6) + Icges(6,6);
t313 = Icges(4,3) + Icges(5,3);
t312 = Icges(7,3) + Icges(6,3);
t212 = qJ(3) + pkin(9);
t201 = sin(t212);
t202 = cos(t212);
t216 = sin(qJ(3));
t219 = cos(qJ(3));
t311 = Icges(4,5) * t216 + Icges(5,5) * t201 + Icges(4,6) * t219 + Icges(5,6) * t202;
t220 = cos(qJ(1));
t310 = t323 * t220;
t217 = sin(qJ(1));
t309 = t323 * t217;
t269 = Icges(5,4) * t201;
t239 = Icges(5,2) * t202 + t269;
t136 = Icges(5,6) * t220 + t217 * t239;
t137 = Icges(5,6) * t217 - t220 * t239;
t268 = Icges(5,4) * t202;
t241 = Icges(5,1) * t201 + t268;
t138 = Icges(5,5) * t220 + t217 * t241;
t139 = Icges(5,5) * t217 - t220 * t241;
t271 = Icges(4,4) * t216;
t240 = Icges(4,2) * t219 + t271;
t146 = Icges(4,6) * t220 + t217 * t240;
t147 = Icges(4,6) * t217 - t220 * t240;
t270 = Icges(4,4) * t219;
t242 = Icges(4,1) * t216 + t270;
t148 = Icges(4,5) * t220 + t217 * t242;
t149 = Icges(4,5) * t217 - t220 * t242;
t164 = -Icges(5,2) * t201 + t268;
t165 = Icges(5,1) * t202 - t269;
t181 = -Icges(4,2) * t216 + t270;
t186 = Icges(4,1) * t219 - t271;
t197 = qJD(3) * t217 + V_base(5);
t198 = qJD(3) * t220 + V_base(4);
t203 = V_base(6) + qJD(1);
t308 = (t137 * t202 + t139 * t201 + t147 * t219 + t149 * t216) * t197 + (t136 * t202 + t138 * t201 + t146 * t219 + t148 * t216) * t198 + (t164 * t202 + t165 * t201 + t181 * t219 + t186 * t216) * t203;
t218 = cos(qJ(5));
t260 = t218 * t220;
t215 = sin(qJ(5));
t262 = t217 * t215;
t158 = -t201 * t262 + t260;
t261 = t217 * t218;
t263 = t215 * t220;
t159 = t201 * t261 + t263;
t265 = t202 * t217;
t307 = t314 * t158 + t317 * t159 - t312 * t265;
t160 = t201 * t263 + t261;
t161 = -t201 * t260 + t262;
t264 = t202 * t220;
t306 = t314 * t160 + t317 * t161 + t312 * t264;
t305 = t315 * t158 + t319 * t159 - t314 * t265;
t304 = t315 * t160 + t319 * t161 + t314 * t264;
t303 = t319 * t158 + t321 * t159 - t317 * t265;
t302 = t319 * t160 + t321 * t161 + t317 * t264;
t301 = (-t314 * t215 + t317 * t218) * t202 + t312 * t201;
t300 = (-t315 * t215 + t319 * t218) * t202 + t314 * t201;
t299 = (-t319 * t215 + t321 * t218) * t202 + t317 * t201;
t280 = pkin(3) * t216;
t155 = qJ(4) * t220 + t217 * t280;
t298 = qJD(4) * t217 + t203 * t155;
t297 = -t316 * t220 - t309;
t296 = t316 * t217 - t310;
t295 = t322 * t217 + t310;
t294 = t322 * t220 - t309;
t291 = (Icges(4,5) * t219 + Icges(5,5) * t202 - Icges(4,6) * t216 - Icges(5,6) * t201) * t203 + (t311 * t217 + t313 * t220) * t198 + (t313 * t217 - t311 * t220) * t197;
t276 = pkin(5) * t218;
t287 = -qJ(6) * t202 + t201 * t276;
t279 = pkin(3) * t219;
t278 = pkin(7) * t220;
t277 = t217 * pkin(7);
t273 = rSges(7,1) * t159 + rSges(7,2) * t158 - rSges(7,3) * t265 + pkin(5) * t263 + t287 * t217;
t259 = t161 * rSges(7,1) + t160 * rSges(7,2) + rSges(7,3) * t264 + pkin(5) * t262 - t287 * t220;
t258 = (rSges(7,1) * t218 - rSges(7,2) * t215 + t276) * t202 + (qJ(6) + rSges(7,3)) * t201;
t257 = qJD(2) * t220;
t256 = qJD(5) * t202;
t255 = qJD(6) * t202;
t193 = pkin(1) * t220 + t217 * qJ(2);
t254 = t203 * t193 + V_base(2);
t189 = t217 * pkin(1) - qJ(2) * t220;
t253 = V_base(4) * t189 + V_base(3);
t252 = V_base(5) * pkin(6) + V_base(1);
t249 = -t189 - t277;
t248 = qJD(2) * t217 + t252;
t154 = qJ(4) * t217 - t220 * t280;
t247 = -t154 + t249;
t246 = V_base(5) * pkin(2) + t248;
t245 = pkin(4) * t201 - pkin(8) * t202;
t244 = rSges(4,1) * t216 + rSges(4,2) * t219;
t243 = rSges(5,1) * t201 + rSges(5,2) * t202;
t230 = qJD(4) * t220 + t197 * t279 + t246;
t227 = t203 * t278 + (-pkin(2) - pkin(6)) * V_base(4) + t254;
t226 = V_base(4) * t277 + (-t193 - t278) * V_base(5) + t253;
t225 = t227 - t257;
t224 = t198 * t154 + t226;
t153 = t245 * t220;
t167 = pkin(4) * t202 + pkin(8) * t201;
t223 = t197 * t167 + (t153 + t247) * t203 + t230;
t152 = t245 * t217;
t222 = -t198 * t153 + (-t152 - t155) * t197 + t224;
t221 = t203 * t152 + (-t167 - t279) * t198 + t227 + t298;
t195 = rSges(2,1) * t220 - t217 * rSges(2,2);
t194 = -rSges(3,2) * t220 + t217 * rSges(3,3);
t192 = rSges(4,1) * t219 - rSges(4,2) * t216;
t191 = t217 * rSges(2,1) + rSges(2,2) * t220;
t190 = -t217 * rSges(3,2) - rSges(3,3) * t220;
t173 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t172 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t171 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t170 = qJD(5) * t201 + t203;
t166 = rSges(5,1) * t202 - rSges(5,2) * t201;
t157 = -t217 * t256 + t198;
t156 = t220 * t256 + t197;
t151 = t217 * rSges(4,3) - t220 * t244;
t150 = rSges(4,3) * t220 + t217 * t244;
t141 = t217 * rSges(5,3) - t220 * t243;
t140 = rSges(5,3) * t220 + t217 * t243;
t132 = rSges(6,3) * t201 + (rSges(6,1) * t218 - rSges(6,2) * t215) * t202;
t130 = V_base(5) * rSges(2,3) - t191 * t203 + t252;
t129 = t195 * t203 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t120 = t191 * V_base(4) - t195 * V_base(5) + V_base(3);
t118 = V_base(5) * rSges(3,1) + (-t189 - t190) * t203 + t248;
t117 = -t257 + t203 * t194 + (-rSges(3,1) - pkin(6)) * V_base(4) + t254;
t116 = t161 * rSges(6,1) + t160 * rSges(6,2) + rSges(6,3) * t264;
t114 = rSges(6,1) * t159 + rSges(6,2) * t158 - rSges(6,3) * t265;
t98 = t190 * V_base(4) + (-t193 - t194) * V_base(5) + t253;
t97 = t192 * t197 + (-t151 + t249) * t203 + t246;
t96 = t203 * t150 - t198 * t192 + t225;
t95 = -t197 * t150 + t198 * t151 + t226;
t94 = t166 * t197 + (-t141 + t247) * t203 + t230;
t93 = t203 * t140 + (-t166 - t279) * t198 + t225 + t298;
t92 = t198 * t141 + (-t140 - t155) * t197 + t224;
t91 = -t116 * t170 + t132 * t156 + t223;
t90 = t170 * t114 - t157 * t132 + t221 - t257;
t89 = -t156 * t114 + t157 * t116 + t222;
t88 = t156 * t258 - t170 * t259 - t217 * t255 + t223;
t87 = (-qJD(2) + t255) * t220 + t273 * t170 - t258 * t157 + t221;
t86 = qJD(6) * t201 - t156 * t273 + t157 * t259 + t222;
t1 = m(7) * (t86 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(6) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + m(5) * (t92 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + m(1) * (t171 ^ 2 + t172 ^ 2 + t173 ^ 2) / 0.2e1 + m(2) * (t120 ^ 2 + t129 ^ 2 + t130 ^ 2) / 0.2e1 + m(3) * (t117 ^ 2 + t118 ^ 2 + t98 ^ 2) / 0.2e1 + m(4) * (t95 ^ 2 + t96 ^ 2 + t97 ^ 2) / 0.2e1 + ((t300 * t160 + t299 * t161 + t301 * t264) * t170 + (t305 * t160 + t303 * t161 + t307 * t264) * t157 + (t304 * t160 + t302 * t161 + t306 * t264) * t156) * t156 / 0.2e1 + ((t300 * t158 + t299 * t159 - t301 * t265) * t170 + (t305 * t158 + t303 * t159 - t307 * t265) * t157 + (t304 * t158 + t302 * t159 - t306 * t265) * t156) * t157 / 0.2e1 + (((-t300 * t215 + t299 * t218) * t170 + (-t305 * t215 + t303 * t218) * t157 + (-t304 * t215 + t302 * t218) * t156) * t202 + (t306 * t156 + t307 * t157 + t301 * t170) * t201) * t170 / 0.2e1 + (t291 * t217 - t308 * t220) * t197 / 0.2e1 + (t308 * t217 + t291 * t220) * t198 / 0.2e1 + ((t297 * t217 + t295 * t220 + Icges(1,4)) * V_base(5) + (t296 * t217 + t294 * t220 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t295 * t217 - t297 * t220 + Icges(1,2)) * V_base(5) + (t294 * t217 - t296 * t220 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t201 * t136 + t202 * t138 - t146 * t216 + t148 * t219) * t198 + (-t137 * t201 + t139 * t202 - t147 * t216 + t149 * t219) * t197 + (-t164 * t201 + t165 * t202 - t181 * t216 + t186 * t219 + Icges(3,1) + Icges(2,3)) * t203) * t203 / 0.2e1 + t203 * V_base(5) * (t320 * t217 - t318 * t220) + t203 * V_base(4) * (t318 * t217 + t320 * t220) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
