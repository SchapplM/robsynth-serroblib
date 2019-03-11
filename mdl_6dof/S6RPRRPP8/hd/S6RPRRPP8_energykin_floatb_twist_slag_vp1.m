% Calculate kinetic energy for
% S6RPRRPP8
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
% Datum: 2019-03-09 04:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPP8_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP8_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRPP8_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP8_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP8_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPP8_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:53:49
% EndTime: 2019-03-09 04:53:52
% DurationCPUTime: 2.96s
% Computational Cost: add. (1053->258), mult. (2018->348), div. (0->0), fcn. (1938->6), ass. (0->129)
t303 = Icges(2,4) + Icges(3,6);
t302 = Icges(2,1) + Icges(3,2);
t301 = -Icges(3,4) + Icges(2,5);
t300 = Icges(3,5) - Icges(2,6);
t299 = Icges(2,2) + Icges(3,3);
t298 = Icges(5,1) + Icges(6,2) + Icges(7,3);
t297 = Icges(5,4) + Icges(6,6) - Icges(7,6);
t296 = Icges(5,5) + Icges(7,5) - Icges(6,4);
t295 = Icges(5,2) + Icges(7,2) + Icges(6,3);
t294 = Icges(5,6) - Icges(6,5) - Icges(7,4);
t293 = Icges(5,3) + Icges(7,1) + Icges(6,1);
t292 = rSges(7,1) + pkin(5);
t291 = rSges(7,3) + qJ(6);
t225 = cos(qJ(1));
t290 = t303 * t225;
t222 = sin(qJ(1));
t289 = t303 * t222;
t288 = -t299 * t225 - t289;
t287 = t299 * t222 - t290;
t286 = t302 * t222 + t290;
t285 = t302 * t225 - t289;
t221 = sin(qJ(3));
t223 = cos(qJ(4));
t255 = t225 * t223;
t220 = sin(qJ(4));
t260 = t222 * t220;
t169 = t221 * t260 - t255;
t259 = t222 * t223;
t261 = t220 * t225;
t170 = t221 * t259 + t261;
t224 = cos(qJ(3));
t258 = t222 * t224;
t282 = -t297 * t169 + t298 * t170 - t296 * t258;
t171 = t221 * t261 + t259;
t172 = t221 * t255 - t260;
t256 = t224 * t225;
t281 = t297 * t171 - t298 * t172 + t296 * t256;
t280 = t295 * t169 - t297 * t170 + t294 * t258;
t279 = -t295 * t171 + t297 * t172 - t294 * t256;
t278 = -t294 * t169 + t296 * t170 - t293 * t258;
t277 = t294 * t171 - t296 * t172 + t293 * t256;
t276 = (-t294 * t220 + t296 * t223) * t224 + t293 * t221;
t275 = (-t295 * t220 + t297 * t223) * t224 + t294 * t221;
t274 = (-t297 * t220 + t298 * t223) * t224 + t296 * t221;
t265 = Icges(4,4) * t221;
t240 = Icges(4,2) * t224 + t265;
t148 = Icges(4,6) * t225 + t222 * t240;
t149 = Icges(4,6) * t222 - t225 * t240;
t264 = Icges(4,4) * t224;
t241 = Icges(4,1) * t221 + t264;
t151 = Icges(4,5) * t225 + t222 * t241;
t152 = Icges(4,5) * t222 - t225 * t241;
t189 = -Icges(4,2) * t221 + t264;
t194 = Icges(4,1) * t224 - t265;
t207 = qJD(3) * t222 + V_base(5);
t208 = qJD(3) * t225 + V_base(4);
t213 = V_base(6) + qJD(1);
t273 = (t148 * t224 + t151 * t221) * t208 + (t149 * t224 + t152 * t221) * t207 + (t189 * t224 + t194 * t221) * t213;
t268 = pkin(7) * t222;
t267 = pkin(7) * t225;
t257 = t223 * t224;
t254 = rSges(7,2) * t169 + t291 * t170 - t258 * t292;
t253 = -t171 * rSges(7,2) - t291 * t172 + t256 * t292;
t252 = (rSges(7,2) * t220 + rSges(7,3) * t223) * t224 + qJ(6) * t257 + t292 * t221;
t251 = qJD(4) * t224;
t198 = t222 * pkin(1) - qJ(2) * t225;
t250 = V_base(4) * t198 + V_base(3);
t249 = V_base(5) * pkin(6) + V_base(1);
t246 = -t198 - t268;
t245 = qJD(2) * t222 + t249;
t244 = V_base(5) * pkin(2) + t245;
t243 = pkin(3) * t221 - pkin(8) * t224;
t242 = rSges(4,1) * t221 + rSges(4,2) * t224;
t239 = Icges(4,5) * t221 + Icges(4,6) * t224;
t202 = pkin(1) * t225 + t222 * qJ(2);
t235 = -qJD(2) * t225 + t213 * t202 + V_base(2);
t234 = (Icges(4,3) * t225 + t222 * t239) * t208 + (Icges(4,3) * t222 - t225 * t239) * t207 + (Icges(4,5) * t224 - Icges(4,6) * t221) * t213;
t233 = V_base(4) * t268 + (-t202 - t267) * V_base(5) + t250;
t232 = t213 * t267 + (-pkin(2) - pkin(6)) * V_base(4) + t235;
t176 = t243 * t225;
t205 = pkin(3) * t224 + pkin(8) * t221;
t231 = t207 * t205 + (t176 + t246) * t213 + t244;
t175 = t243 * t222;
t230 = -t207 * t175 - t208 * t176 + t233;
t167 = t225 * t251 + t207;
t174 = (pkin(4) * t223 + qJ(5) * t220) * t224;
t229 = qJD(5) * t169 + t167 * t174 + t231;
t228 = t213 * t175 - t208 * t205 + t232;
t135 = -pkin(4) * t172 - qJ(5) * t171;
t168 = -t222 * t251 + t208;
t227 = qJD(5) * t224 * t220 + t168 * t135 + t230;
t134 = pkin(4) * t170 + qJ(5) * t169;
t197 = qJD(4) * t221 + t213;
t226 = -qJD(5) * t171 + t197 * t134 + t228;
t204 = rSges(2,1) * t225 - t222 * rSges(2,2);
t203 = -rSges(3,2) * t225 + t222 * rSges(3,3);
t201 = rSges(4,1) * t224 - rSges(4,2) * t221;
t200 = t222 * rSges(2,1) + rSges(2,2) * t225;
t199 = -t222 * rSges(3,2) - rSges(3,3) * t225;
t181 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t180 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t179 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t163 = rSges(6,1) * t221 + (-rSges(6,2) * t223 + rSges(6,3) * t220) * t224;
t161 = t222 * rSges(4,3) - t225 * t242;
t160 = rSges(5,3) * t221 + (rSges(5,1) * t223 - rSges(5,2) * t220) * t224;
t159 = rSges(4,3) * t225 + t222 * t242;
t139 = V_base(5) * rSges(2,3) - t200 * t213 + t249;
t138 = t204 * t213 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t137 = t200 * V_base(4) - t204 * V_base(5) + V_base(3);
t133 = -t172 * rSges(5,1) + t171 * rSges(5,2) + rSges(5,3) * t256;
t132 = rSges(5,1) * t170 - rSges(5,2) * t169 - rSges(5,3) * t258;
t131 = rSges(6,1) * t256 + t172 * rSges(6,2) - t171 * rSges(6,3);
t129 = -rSges(6,1) * t258 - rSges(6,2) * t170 + rSges(6,3) * t169;
t108 = V_base(5) * rSges(3,1) + (-t198 - t199) * t213 + t245;
t107 = t213 * t203 + (-rSges(3,1) - pkin(6)) * V_base(4) + t235;
t105 = t199 * V_base(4) + (-t202 - t203) * V_base(5) + t250;
t104 = t201 * t207 + (-t161 + t246) * t213 + t244;
t103 = t213 * t159 - t208 * t201 + t232;
t102 = -t207 * t159 + t208 * t161 + t233;
t101 = -t133 * t197 + t160 * t167 + t231;
t100 = t197 * t132 - t168 * t160 + t228;
t99 = -t167 * t132 + t168 * t133 + t230;
t98 = t163 * t167 + (-t131 - t135) * t197 + t229;
t97 = t197 * t129 + (-t163 - t174) * t168 + t226;
t96 = t168 * t131 + (-t129 - t134) * t167 + t227;
t95 = qJD(6) * t170 + t252 * t167 + (-t135 - t253) * t197 + t229;
t94 = -qJD(6) * t172 + t254 * t197 + (-t174 - t252) * t168 + t226;
t93 = qJD(6) * t257 + t253 * t168 + (-t134 - t254) * t167 + t227;
t1 = m(7) * (t93 ^ 2 + t94 ^ 2 + t95 ^ 2) / 0.2e1 + m(6) * (t96 ^ 2 + t97 ^ 2 + t98 ^ 2) / 0.2e1 + m(5) * (t100 ^ 2 + t101 ^ 2 + t99 ^ 2) / 0.2e1 + m(4) * (t102 ^ 2 + t103 ^ 2 + t104 ^ 2) / 0.2e1 + m(3) * (t105 ^ 2 + t107 ^ 2 + t108 ^ 2) / 0.2e1 + t208 * (t222 * t273 + t234 * t225) / 0.2e1 + t207 * (t234 * t222 - t225 * t273) / 0.2e1 + m(2) * (t137 ^ 2 + t138 ^ 2 + t139 ^ 2) / 0.2e1 + m(1) * (t179 ^ 2 + t180 ^ 2 + t181 ^ 2) / 0.2e1 + ((t171 * t275 - t172 * t274 + t256 * t276) * t197 + (-t171 * t280 - t172 * t282 + t256 * t278) * t168 + (-t171 * t279 - t172 * t281 + t256 * t277) * t167) * t167 / 0.2e1 + ((-t169 * t275 + t170 * t274 - t258 * t276) * t197 + (t169 * t280 + t170 * t282 - t258 * t278) * t168 + (t169 * t279 + t170 * t281 - t258 * t277) * t167) * t168 / 0.2e1 + (((-t220 * t275 + t274 * t223) * t197 + (t220 * t280 + t282 * t223) * t168 + (t220 * t279 + t281 * t223) * t167) * t224 + (t167 * t277 + t168 * t278 + t197 * t276) * t221) * t197 / 0.2e1 + ((-t148 * t221 + t151 * t224) * t208 + (-t149 * t221 + t152 * t224) * t207 + (-t189 * t221 + t194 * t224 + Icges(3,1) + Icges(2,3)) * t213) * t213 / 0.2e1 + ((t222 * t288 + t286 * t225 + Icges(1,4)) * V_base(5) + (t222 * t287 + t225 * t285 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t286 * t222 - t225 * t288 + Icges(1,2)) * V_base(5) + (t222 * t285 - t225 * t287 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t213 * (t301 * t222 - t300 * t225) + V_base(4) * t213 * (t300 * t222 + t301 * t225) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
