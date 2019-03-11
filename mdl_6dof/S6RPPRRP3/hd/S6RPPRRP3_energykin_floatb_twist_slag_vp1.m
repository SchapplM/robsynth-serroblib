% Calculate kinetic energy for
% S6RPPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
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
% Datum: 2019-03-09 02:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRP3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP3_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPPRRP3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP3_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP3_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRP3_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:02:39
% EndTime: 2019-03-09 02:02:41
% DurationCPUTime: 2.50s
% Computational Cost: add. (1442->266), mult. (1532->352), div. (0->0), fcn. (1378->8), ass. (0->136)
t298 = Icges(3,4) + Icges(4,6);
t297 = Icges(3,1) + Icges(4,2);
t296 = Icges(6,1) + Icges(7,1);
t295 = -Icges(4,4) + Icges(3,5);
t294 = Icges(6,4) - Icges(7,5);
t293 = Icges(7,4) + Icges(6,5);
t292 = Icges(4,5) - Icges(3,6);
t291 = Icges(3,2) + Icges(4,3);
t290 = Icges(6,2) + Icges(7,3);
t289 = Icges(7,6) - Icges(6,6);
t288 = Icges(6,3) + Icges(7,2);
t287 = rSges(7,1) + pkin(5);
t286 = rSges(7,3) + qJ(6);
t211 = qJ(1) + pkin(9);
t205 = cos(t211);
t285 = t298 * t205;
t204 = sin(t211);
t284 = t298 * t204;
t215 = cos(qJ(5));
t212 = sin(qJ(5));
t213 = sin(qJ(4));
t250 = t212 * t213;
t141 = t204 * t250 - t205 * t215;
t249 = t213 * t215;
t142 = t204 * t249 + t205 * t212;
t216 = cos(qJ(4));
t252 = t204 * t216;
t283 = t141 * t290 - t142 * t294 - t252 * t289;
t143 = t204 * t215 + t205 * t250;
t144 = t204 * t212 - t205 * t249;
t251 = t205 * t216;
t282 = -t143 * t290 - t144 * t294 + t251 * t289;
t281 = t141 * t289 + t142 * t293 - t252 * t288;
t280 = -t143 * t289 + t144 * t293 + t251 * t288;
t279 = -t141 * t294 + t142 * t296 - t252 * t293;
t278 = t143 * t294 + t144 * t296 + t251 * t293;
t277 = (t212 * t290 - t215 * t294) * t216 + t289 * t213;
t276 = (t212 * t289 + t215 * t293) * t216 + t288 * t213;
t275 = (-t212 * t294 + t215 * t296) * t216 + t293 * t213;
t274 = -t205 * t291 - t284;
t273 = t204 * t291 - t285;
t272 = t204 * t297 + t285;
t271 = t205 * t297 - t284;
t256 = Icges(5,4) * t213;
t231 = Icges(5,2) * t216 + t256;
t129 = Icges(5,6) * t205 + t204 * t231;
t130 = Icges(5,6) * t204 - t205 * t231;
t255 = Icges(5,4) * t216;
t232 = Icges(5,1) * t213 + t255;
t131 = Icges(5,5) * t205 + t204 * t232;
t132 = Icges(5,5) * t204 - t205 * t232;
t180 = qJD(4) * t204 + V_base(5);
t181 = qJD(4) * t205 + V_base(4);
t185 = -Icges(5,2) * t213 + t255;
t188 = Icges(5,1) * t216 - t256;
t206 = V_base(6) + qJD(1);
t268 = (t129 * t216 + t131 * t213) * t181 + (t130 * t216 + t132 * t213) * t180 + (t185 * t216 + t188 * t213) * t206;
t214 = sin(qJ(1));
t263 = pkin(1) * t214;
t217 = cos(qJ(1));
t262 = pkin(1) * t217;
t261 = pkin(7) * t204;
t260 = pkin(7) * t205;
t259 = -pkin(6) - qJ(2);
t258 = Icges(2,4) * t214;
t248 = -rSges(7,2) * t252 + t286 * t141 + t287 * t142;
t247 = rSges(7,2) * t251 - t286 * t143 + t287 * t144;
t246 = rSges(7,2) * t213 + (t286 * t212 + t287 * t215) * t216;
t245 = qJD(5) * t216;
t244 = t206 * t262 + V_base(2);
t243 = V_base(5) * pkin(6) + V_base(1);
t170 = pkin(2) * t204 - qJ(3) * t205;
t240 = -t170 - t263;
t173 = pkin(2) * t205 + qJ(3) * t204;
t239 = -t173 - t262;
t238 = V_base(5) * qJ(2) + t243;
t237 = V_base(4) * t263 + qJD(2) + V_base(3);
t236 = qJD(3) * t204 + t238;
t235 = pkin(4) * t213 - pkin(8) * t216;
t234 = V_base(4) * t170 + t237;
t233 = rSges(5,1) * t213 + rSges(5,2) * t216;
t230 = Icges(5,5) * t213 + Icges(5,6) * t216;
t226 = V_base(5) * pkin(3) + t236;
t225 = t240 - t261;
t224 = -qJD(3) * t205 + t206 * t173 + t244;
t223 = (Icges(5,3) * t205 + t204 * t230) * t181 + (Icges(5,3) * t204 - t205 * t230) * t180 + (Icges(5,5) * t216 - Icges(5,6) * t213) * t206;
t222 = t206 * t260 + (-pkin(3) + t259) * V_base(4) + t224;
t221 = V_base(4) * t261 + (t239 - t260) * V_base(5) + t234;
t155 = t235 * t205;
t195 = pkin(4) * t216 + pkin(8) * t213;
t220 = t180 * t195 + (t155 + t225) * t206 + t226;
t154 = t235 * t204;
t219 = t206 * t154 - t181 * t195 + t222;
t218 = -t180 * t154 - t181 * t155 + t221;
t208 = Icges(2,4) * t217;
t194 = rSges(2,1) * t217 - t214 * rSges(2,2);
t193 = rSges(5,1) * t216 - rSges(5,2) * t213;
t192 = t214 * rSges(2,1) + rSges(2,2) * t217;
t191 = qJD(5) * t213 + t206;
t190 = Icges(2,1) * t217 - t258;
t189 = Icges(2,1) * t214 + t208;
t187 = -Icges(2,2) * t214 + t208;
t186 = Icges(2,2) * t217 + t258;
t178 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t177 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t176 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t175 = rSges(3,1) * t205 - rSges(3,2) * t204;
t174 = -rSges(4,2) * t205 + rSges(4,3) * t204;
t172 = rSges(3,1) * t204 + rSges(3,2) * t205;
t171 = -rSges(4,2) * t204 - rSges(4,3) * t205;
t153 = rSges(6,3) * t213 + (rSges(6,1) * t215 - rSges(6,2) * t212) * t216;
t140 = -t204 * t245 + t181;
t139 = t205 * t245 + t180;
t136 = rSges(5,3) * t204 - t205 * t233;
t135 = rSges(5,3) * t205 + t204 * t233;
t134 = V_base(5) * rSges(2,3) - t192 * t206 + t243;
t133 = t194 * t206 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t126 = t192 * V_base(4) - t194 * V_base(5) + V_base(3);
t122 = V_base(5) * rSges(3,3) + (-t172 - t263) * t206 + t238;
t121 = t175 * t206 + (-rSges(3,3) + t259) * V_base(4) + t244;
t120 = V_base(4) * t172 + (-t175 - t262) * V_base(5) + t237;
t119 = rSges(6,1) * t144 + rSges(6,2) * t143 + rSges(6,3) * t251;
t117 = rSges(6,1) * t142 - rSges(6,2) * t141 - rSges(6,3) * t252;
t103 = V_base(5) * rSges(4,1) + (-t171 + t240) * t206 + t236;
t102 = t174 * t206 + (-rSges(4,1) + t259) * V_base(4) + t224;
t101 = V_base(4) * t171 + (-t174 + t239) * V_base(5) + t234;
t100 = t180 * t193 + (-t136 + t225) * t206 + t226;
t99 = t135 * t206 - t181 * t193 + t222;
t98 = -t180 * t135 + t181 * t136 + t221;
t97 = -t119 * t191 + t139 * t153 + t220;
t96 = t117 * t191 - t140 * t153 + t219;
t95 = -t139 * t117 + t140 * t119 + t218;
t94 = qJD(6) * t141 + t139 * t246 - t191 * t247 + t220;
t93 = -qJD(6) * t143 - t140 * t246 + t191 * t248 + t219;
t92 = qJD(6) * t216 * t212 - t139 * t248 + t140 * t247 + t218;
t1 = m(1) * (t176 ^ 2 + t177 ^ 2 + t178 ^ 2) / 0.2e1 + m(2) * (t126 ^ 2 + t133 ^ 2 + t134 ^ 2) / 0.2e1 + m(3) * (t120 ^ 2 + t121 ^ 2 + t122 ^ 2) / 0.2e1 + m(4) * (t101 ^ 2 + t102 ^ 2 + t103 ^ 2) / 0.2e1 + m(7) * (t92 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + m(6) * (t95 ^ 2 + t96 ^ 2 + t97 ^ 2) / 0.2e1 + m(5) * (t100 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + t181 * (t204 * t268 + t223 * t205) / 0.2e1 + t180 * (t223 * t204 - t205 * t268) / 0.2e1 + ((-t143 * t277 + t144 * t275 + t251 * t276) * t191 + (-t143 * t283 + t279 * t144 + t281 * t251) * t140 + (-t143 * t282 + t144 * t278 + t251 * t280) * t139) * t139 / 0.2e1 + ((t141 * t277 + t142 * t275 - t252 * t276) * t191 + (t141 * t283 + t279 * t142 - t281 * t252) * t140 + (t141 * t282 + t142 * t278 - t252 * t280) * t139) * t140 / 0.2e1 + (((t212 * t277 + t215 * t275) * t191 + (t212 * t283 + t279 * t215) * t140 + (t212 * t282 + t215 * t278) * t139) * t216 + (t139 * t280 + t140 * t281 + t191 * t276) * t213) * t191 / 0.2e1 + ((-t129 * t213 + t131 * t216) * t181 + (-t130 * t213 + t132 * t216) * t180 + (-t185 * t213 + t188 * t216 + Icges(4,1) + Icges(2,3) + Icges(3,3)) * t206) * t206 / 0.2e1 + ((-t214 * t186 + t189 * t217 + t204 * t274 + t205 * t272 + Icges(1,4)) * V_base(5) + (-t214 * t187 + t190 * t217 + t204 * t273 + t205 * t271 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t186 * t217 + t214 * t189 + t204 * t272 - t205 * t274 + Icges(1,2)) * V_base(5) + (t187 * t217 + t214 * t190 + t204 * t271 - t205 * t273 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t206 * (Icges(2,5) * t214 + Icges(2,6) * t217 + t204 * t295 - t205 * t292) + V_base(4) * t206 * (Icges(2,5) * t217 - Icges(2,6) * t214 + t204 * t292 + t205 * t295) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
