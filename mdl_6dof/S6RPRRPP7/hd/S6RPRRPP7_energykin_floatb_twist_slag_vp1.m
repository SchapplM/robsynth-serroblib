% Calculate kinetic energy for
% S6RPRRPP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
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
% Datum: 2019-03-09 04:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPP7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP7_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRPP7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP7_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP7_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPP7_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:50:02
% EndTime: 2019-03-09 04:50:04
% DurationCPUTime: 2.94s
% Computational Cost: add. (1051->261), mult. (2013->347), div. (0->0), fcn. (1931->6), ass. (0->131)
t302 = Icges(2,4) + Icges(3,6);
t301 = Icges(2,1) + Icges(3,2);
t300 = -Icges(3,4) + Icges(2,5);
t299 = Icges(3,5) - Icges(2,6);
t298 = Icges(2,2) + Icges(3,3);
t297 = Icges(5,1) + Icges(6,1) + Icges(7,1);
t296 = Icges(5,4) - Icges(7,4) - Icges(6,5);
t295 = -Icges(7,5) + Icges(6,4) + Icges(5,5);
t294 = Icges(5,2) + Icges(7,2) + Icges(6,3);
t293 = Icges(6,6) - Icges(7,6) - Icges(5,6);
t292 = -Icges(7,3) - Icges(5,3) - Icges(6,2);
t291 = rSges(7,1) + pkin(5);
t290 = rSges(7,3) + qJ(6);
t222 = cos(qJ(1));
t289 = t302 * t222;
t219 = sin(qJ(1));
t288 = t302 * t219;
t218 = sin(qJ(3));
t220 = cos(qJ(4));
t254 = t222 * t220;
t217 = sin(qJ(4));
t258 = t219 * t217;
t167 = t218 * t258 - t254;
t257 = t219 * t220;
t259 = t217 * t222;
t168 = t218 * t257 + t259;
t132 = pkin(4) * t168 + qJ(5) * t167;
t169 = t218 * t259 + t257;
t210 = V_base(6) + qJD(1);
t195 = qJD(4) * t218 + t210;
t287 = -qJD(5) * t169 + t195 * t132;
t286 = -t298 * t222 - t288;
t285 = t298 * t219 - t289;
t284 = t301 * t219 + t289;
t283 = t301 * t222 - t288;
t221 = cos(qJ(3));
t256 = t219 * t221;
t280 = -t293 * t167 - t295 * t168 - t292 * t256;
t170 = -t218 * t254 + t258;
t255 = t221 * t222;
t279 = t293 * t169 - t295 * t170 + t292 * t255;
t278 = t294 * t167 - t296 * t168 - t293 * t256;
t277 = -t294 * t169 - t296 * t170 + t293 * t255;
t276 = -t296 * t167 + t297 * t168 - t295 * t256;
t275 = t296 * t169 + t297 * t170 + t295 * t255;
t274 = (-t293 * t217 - t295 * t220) * t221 + t292 * t218;
t273 = (t294 * t217 - t296 * t220) * t221 + t293 * t218;
t272 = (-t296 * t217 + t297 * t220) * t221 + t295 * t218;
t263 = Icges(4,4) * t218;
t236 = Icges(4,2) * t221 + t263;
t150 = Icges(4,6) * t222 + t219 * t236;
t151 = Icges(4,6) * t219 - t222 * t236;
t262 = Icges(4,4) * t221;
t237 = Icges(4,1) * t218 + t262;
t155 = Icges(4,5) * t222 + t219 * t237;
t156 = Icges(4,5) * t219 - t222 * t237;
t187 = -Icges(4,2) * t218 + t262;
t192 = Icges(4,1) * t221 - t263;
t205 = qJD(3) * t219 + V_base(5);
t206 = qJD(3) * t222 + V_base(4);
t271 = (t150 * t221 + t155 * t218) * t206 + (t151 * t221 + t156 * t218) * t205 + (t187 * t221 + t192 * t218) * t210;
t266 = pkin(7) * t219;
t265 = pkin(7) * t222;
t253 = rSges(7,2) * t167 + t168 * t291 + t290 * t256;
t252 = -t169 * rSges(7,2) + t170 * t291 - t290 * t255;
t251 = (rSges(7,2) * t217 + t220 * t291) * t221 - t290 * t218;
t250 = qJD(2) * t222;
t249 = qJD(4) * t221;
t248 = qJD(6) * t221;
t200 = pkin(1) * t222 + t219 * qJ(2);
t247 = t210 * t200 + V_base(2);
t196 = t219 * pkin(1) - qJ(2) * t222;
t246 = V_base(4) * t196 + V_base(3);
t245 = V_base(5) * pkin(6) + V_base(1);
t242 = -t196 - t266;
t241 = qJD(2) * t219 + t245;
t240 = V_base(5) * pkin(2) + t241;
t239 = pkin(3) * t218 - pkin(8) * t221;
t238 = rSges(4,1) * t218 + rSges(4,2) * t221;
t235 = Icges(4,5) * t218 + Icges(4,6) * t221;
t231 = (Icges(4,3) * t222 + t219 * t235) * t206 + (Icges(4,3) * t219 - t222 * t235) * t205 + (Icges(4,5) * t221 - Icges(4,6) * t218) * t210;
t230 = t210 * t265 + (-pkin(2) - pkin(6)) * V_base(4) + t247;
t229 = V_base(4) * t266 + (-t200 - t265) * V_base(5) + t246;
t173 = t239 * t219;
t203 = pkin(3) * t221 + pkin(8) * t218;
t228 = t210 * t173 - t206 * t203 + t230;
t174 = t239 * t222;
t227 = t205 * t203 + (t174 + t242) * t210 + t240;
t226 = -t205 * t173 - t206 * t174 + t229;
t165 = t222 * t249 + t205;
t172 = (pkin(4) * t220 + qJ(5) * t217) * t221;
t225 = qJD(5) * t167 + t165 * t172 + t227;
t224 = t228 - t250;
t133 = pkin(4) * t170 - qJ(5) * t169;
t166 = -t219 * t249 + t206;
t223 = qJD(5) * t221 * t217 + t166 * t133 + t226;
t202 = rSges(2,1) * t222 - t219 * rSges(2,2);
t201 = -rSges(3,2) * t222 + t219 * rSges(3,3);
t199 = rSges(4,1) * t221 - rSges(4,2) * t218;
t198 = t219 * rSges(2,1) + rSges(2,2) * t222;
t197 = -t219 * rSges(3,2) - rSges(3,3) * t222;
t179 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t178 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t177 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t161 = t219 * rSges(4,3) - t222 * t238;
t160 = rSges(5,3) * t218 + (rSges(5,1) * t220 - rSges(5,2) * t217) * t221;
t159 = rSges(6,2) * t218 + (rSges(6,1) * t220 + rSges(6,3) * t217) * t221;
t157 = rSges(4,3) * t222 + t219 * t238;
t137 = V_base(5) * rSges(2,3) - t198 * t210 + t245;
t136 = t202 * t210 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t135 = t198 * V_base(4) - t202 * V_base(5) + V_base(3);
t131 = t170 * rSges(5,1) + t169 * rSges(5,2) + rSges(5,3) * t255;
t130 = t170 * rSges(6,1) + rSges(6,2) * t255 - t169 * rSges(6,3);
t128 = rSges(5,1) * t168 - rSges(5,2) * t167 - rSges(5,3) * t256;
t127 = rSges(6,1) * t168 - rSges(6,2) * t256 + rSges(6,3) * t167;
t106 = V_base(5) * rSges(3,1) + (-t196 - t197) * t210 + t241;
t105 = -t250 + t210 * t201 + (-rSges(3,1) - pkin(6)) * V_base(4) + t247;
t103 = t197 * V_base(4) + (-t200 - t201) * V_base(5) + t246;
t102 = t199 * t205 + (-t161 + t242) * t210 + t240;
t101 = t210 * t157 - t206 * t199 + t230 - t250;
t100 = -t205 * t157 + t206 * t161 + t229;
t99 = -t131 * t195 + t160 * t165 + t227;
t98 = t195 * t128 - t166 * t160 + t224;
t97 = -t165 * t128 + t166 * t131 + t226;
t96 = t159 * t165 + (-t130 - t133) * t195 + t225;
t95 = t195 * t127 + (-t159 - t172) * t166 + t224 + t287;
t94 = t166 * t130 + (-t127 - t132) * t165 + t223;
t93 = t219 * t248 + t251 * t165 + (-t133 - t252) * t195 + t225;
t92 = (-qJD(2) - t248) * t222 + t253 * t195 + (-t172 - t251) * t166 + t228 + t287;
t91 = -qJD(6) * t218 + t252 * t166 + (-t132 - t253) * t165 + t223;
t1 = m(2) * (t135 ^ 2 + t136 ^ 2 + t137 ^ 2) / 0.2e1 + m(3) * (t103 ^ 2 + t105 ^ 2 + t106 ^ 2) / 0.2e1 + m(5) * (t97 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(4) * (t100 ^ 2 + t101 ^ 2 + t102 ^ 2) / 0.2e1 + m(7) * (t91 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + m(6) * (t94 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + m(1) * (t177 ^ 2 + t178 ^ 2 + t179 ^ 2) / 0.2e1 + t206 * (t271 * t219 + t231 * t222) / 0.2e1 + t205 * (t231 * t219 - t271 * t222) / 0.2e1 + ((-t273 * t169 + t272 * t170 - t274 * t255) * t195 + (-t278 * t169 + t276 * t170 - t280 * t255) * t166 + (-t277 * t169 + t275 * t170 - t279 * t255) * t165) * t165 / 0.2e1 + ((t273 * t167 + t272 * t168 + t274 * t256) * t195 + (t278 * t167 + t276 * t168 + t280 * t256) * t166 + (t277 * t167 + t275 * t168 + t279 * t256) * t165) * t166 / 0.2e1 + (((t273 * t217 + t272 * t220) * t195 + (t278 * t217 + t276 * t220) * t166 + (t277 * t217 + t275 * t220) * t165) * t221 + (-t279 * t165 - t280 * t166 - t274 * t195) * t218) * t195 / 0.2e1 + ((-t150 * t218 + t155 * t221) * t206 + (-t151 * t218 + t156 * t221) * t205 + (-t187 * t218 + t192 * t221 + Icges(3,1) + Icges(2,3)) * t210) * t210 / 0.2e1 + ((t286 * t219 + t284 * t222 + Icges(1,4)) * V_base(5) + (t285 * t219 + t283 * t222 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t284 * t219 - t286 * t222 + Icges(1,2)) * V_base(5) + (t283 * t219 - t285 * t222 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t210 * (t300 * t219 - t299 * t222) + V_base(4) * t210 * (t299 * t219 + t300 * t222) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
