% Calculate kinetic energy for
% S6RPRPPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
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
% Datum: 2019-03-09 03:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPPR8_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR8_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPPR8_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR8_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR8_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPPR8_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:58:41
% EndTime: 2019-03-09 02:58:44
% DurationCPUTime: 3.15s
% Computational Cost: add. (823->260), mult. (1467->350), div. (0->0), fcn. (1249->6), ass. (0->135)
t309 = -Icges(4,4) - Icges(6,4) + Icges(5,5);
t308 = Icges(4,1) + Icges(5,1) + Icges(6,2);
t307 = Icges(6,1) + Icges(4,2) + Icges(5,3);
t206 = cos(qJ(3));
t306 = t309 * t206;
t203 = sin(qJ(3));
t305 = t309 * t203;
t304 = Icges(5,4) + Icges(4,5) + Icges(6,6);
t303 = Icges(6,5) + Icges(4,6) - Icges(5,6);
t302 = -t307 * t206 + t305;
t301 = t308 * t203 - t306;
t300 = Icges(2,4) + Icges(3,6);
t299 = Icges(2,1) + Icges(3,2);
t298 = -Icges(3,4) + Icges(2,5);
t297 = Icges(3,5) - Icges(2,6);
t296 = Icges(2,2) + Icges(3,3);
t295 = Icges(5,2) + Icges(4,3) + Icges(6,3);
t204 = sin(qJ(1));
t207 = cos(qJ(1));
t294 = t302 * t204 - t303 * t207;
t293 = -t303 * t204 - t302 * t207;
t292 = t301 * t204 + t304 * t207;
t291 = t304 * t204 - t301 * t207;
t290 = t307 * t203 + t306;
t289 = t308 * t206 + t305;
t288 = t304 * t203 + t303 * t206;
t287 = t300 * t207;
t286 = t300 * t204;
t190 = qJD(3) * t204 + V_base(5);
t191 = qJD(3) * t207 + V_base(4);
t194 = V_base(6) + qJD(1);
t285 = t190 * (t291 * t203 - t293 * t206) + t191 * (t292 * t203 - t294 * t206) + t194 * (t289 * t203 - t290 * t206);
t284 = -t296 * t207 - t286;
t283 = t296 * t204 - t287;
t282 = t299 * t204 + t287;
t281 = t299 * t207 - t286;
t278 = (-t303 * t203 + t304 * t206) * t194 + (t288 * t204 + t295 * t207) * t191 + (t295 * t204 - t288 * t207) * t190;
t271 = pkin(4) * t206;
t270 = pkin(7) * t204;
t269 = pkin(7) * t207;
t259 = t203 * t204;
t258 = t203 * t207;
t257 = t204 * t206;
t256 = t206 * t207;
t239 = pkin(3) * t203 - qJ(4) * t206;
t145 = t239 * t204;
t151 = pkin(4) * t259 - qJ(5) * t207;
t255 = -t145 - t151;
t254 = qJD(4) * t206;
t253 = qJD(6) * t203;
t179 = t204 * pkin(1) - qJ(2) * t207;
t252 = V_base(4) * t179 + V_base(3);
t251 = V_base(5) * pkin(6) + V_base(1);
t182 = pkin(3) * t206 + qJ(4) * t203;
t248 = -t182 - t271;
t247 = -t179 - t270;
t246 = qJD(2) * t204 + t251;
t146 = t239 * t207;
t245 = t146 + t247;
t244 = V_base(5) * pkin(2) + t246;
t243 = pkin(5) * t206 - pkin(8) * t203;
t242 = rSges(4,1) * t203 + rSges(4,2) * t206;
t241 = rSges(5,1) * t203 - rSges(5,3) * t206;
t240 = rSges(6,1) * t206 + rSges(6,2) * t203;
t185 = pkin(1) * t207 + t204 * qJ(2);
t220 = -qJD(2) * t207 + t194 * t185 + V_base(2);
t152 = -pkin(4) * t258 - t204 * qJ(5);
t219 = -t152 + t245;
t215 = t190 * t182 - t204 * t254 + t244;
t214 = V_base(4) * t270 + (-t185 - t269) * V_base(5) + t252;
t213 = t194 * t269 + (-pkin(2) - pkin(6)) * V_base(4) + t220;
t212 = qJD(4) * t203 - t191 * t146 + t214;
t211 = t194 * t145 + t207 * t254 + t213;
t210 = -qJD(5) * t207 + t190 * t271 + t215;
t209 = t191 * t152 + t212;
t208 = -qJD(5) * t204 + t194 * t151 + t211;
t205 = cos(qJ(6));
t202 = sin(qJ(6));
t188 = pkin(5) * t203 + pkin(8) * t206;
t187 = rSges(2,1) * t207 - t204 * rSges(2,2);
t186 = -rSges(3,2) * t207 + t204 * rSges(3,3);
t184 = rSges(4,1) * t206 - rSges(4,2) * t203;
t183 = rSges(5,1) * t206 + rSges(5,3) * t203;
t181 = t204 * rSges(2,1) + rSges(2,2) * t207;
t180 = -t204 * rSges(3,2) - rSges(3,3) * t207;
t178 = rSges(6,1) * t203 - rSges(6,2) * t206;
t177 = qJD(6) * t206 + t194;
t155 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t154 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t153 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t148 = t243 * t207;
t147 = t243 * t204;
t143 = -t204 * t202 + t205 * t256;
t142 = -t202 * t256 - t204 * t205;
t141 = -t202 * t207 - t205 * t257;
t140 = t202 * t257 - t205 * t207;
t139 = t204 * t253 + t191;
t138 = -t207 * t253 + t190;
t135 = -t204 * rSges(6,3) + t207 * t240;
t134 = t204 * rSges(4,3) - t207 * t242;
t133 = t204 * rSges(5,2) - t207 * t241;
t132 = -rSges(6,3) * t207 - t204 * t240;
t131 = rSges(4,3) * t207 + t204 * t242;
t130 = rSges(5,2) * t207 + t204 * t241;
t129 = rSges(7,3) * t206 + (rSges(7,1) * t205 - rSges(7,2) * t202) * t203;
t122 = Icges(7,5) * t206 + (Icges(7,1) * t205 - Icges(7,4) * t202) * t203;
t115 = Icges(7,6) * t206 + (Icges(7,4) * t205 - Icges(7,2) * t202) * t203;
t108 = Icges(7,3) * t206 + (Icges(7,5) * t205 - Icges(7,6) * t202) * t203;
t104 = V_base(5) * rSges(2,3) - t181 * t194 + t251;
t103 = t187 * t194 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t102 = t181 * V_base(4) - t187 * V_base(5) + V_base(3);
t101 = t143 * rSges(7,1) + t142 * rSges(7,2) - rSges(7,3) * t258;
t100 = rSges(7,1) * t141 + rSges(7,2) * t140 + rSges(7,3) * t259;
t99 = Icges(7,1) * t143 + Icges(7,4) * t142 - Icges(7,5) * t258;
t98 = Icges(7,1) * t141 + Icges(7,4) * t140 + Icges(7,5) * t259;
t97 = Icges(7,4) * t143 + Icges(7,2) * t142 - Icges(7,6) * t258;
t96 = Icges(7,4) * t141 + Icges(7,2) * t140 + Icges(7,6) * t259;
t95 = Icges(7,5) * t143 + Icges(7,6) * t142 - Icges(7,3) * t258;
t94 = Icges(7,5) * t141 + Icges(7,6) * t140 + Icges(7,3) * t259;
t93 = V_base(5) * rSges(3,1) + (-t179 - t180) * t194 + t246;
t92 = t194 * t186 + (-rSges(3,1) - pkin(6)) * V_base(4) + t220;
t91 = t180 * V_base(4) + (-t185 - t186) * V_base(5) + t252;
t90 = t184 * t190 + (-t134 + t247) * t194 + t244;
t89 = t194 * t131 - t191 * t184 + t213;
t88 = -t190 * t131 + t191 * t134 + t214;
t87 = t183 * t190 + (-t133 + t245) * t194 + t215;
t86 = t194 * t130 + (-t182 - t183) * t191 + t211;
t85 = t191 * t133 + (-t130 - t145) * t190 + t212;
t84 = t190 * t178 + (-t135 + t219) * t194 + t210;
t83 = t194 * t132 + (-t178 + t248) * t191 + t208;
t82 = t191 * t135 + (-t132 + t255) * t190 + t209;
t81 = -t177 * t101 + t138 * t129 + t190 * t188 + (-t148 + t219) * t194 + t210;
t80 = t177 * t100 - t139 * t129 - t194 * t147 + (-t188 + t248) * t191 + t208;
t79 = -t138 * t100 + t139 * t101 + t191 * t148 + (t147 + t255) * t190 + t209;
t1 = m(1) * (t153 ^ 2 + t154 ^ 2 + t155 ^ 2) / 0.2e1 + m(2) * (t102 ^ 2 + t103 ^ 2 + t104 ^ 2) / 0.2e1 + m(4) * (t88 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + m(3) * (t91 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + m(7) * (t79 ^ 2 + t80 ^ 2 + t81 ^ 2) / 0.2e1 + m(6) * (t82 ^ 2 + t83 ^ 2 + t84 ^ 2) / 0.2e1 + m(5) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + t177 * ((t108 * t177 + t95 * t138 + t94 * t139) * t206 + ((-t202 * t96 + t205 * t98) * t139 + (-t202 * t97 + t205 * t99) * t138 + (-t115 * t202 + t122 * t205) * t177) * t203) / 0.2e1 + t138 * ((t142 * t96 + t143 * t98 - t94 * t258) * t139 + (t142 * t97 + t143 * t99 - t95 * t258) * t138 + (-t108 * t258 + t142 * t115 + t143 * t122) * t177) / 0.2e1 + t139 * ((t140 * t96 + t141 * t98 + t94 * t259) * t139 + (t140 * t97 + t141 * t99 + t95 * t259) * t138 + (t108 * t259 + t115 * t140 + t122 * t141) * t177) / 0.2e1 + (t278 * t204 - t285 * t207) * t190 / 0.2e1 + (t285 * t204 + t278 * t207) * t191 / 0.2e1 + ((t204 * t284 + t282 * t207 + Icges(1,4)) * V_base(5) + (t283 * t204 + t281 * t207 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t282 * t204 - t284 * t207 + Icges(1,2)) * V_base(5) + (t204 * t281 - t207 * t283 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t294 * t203 + t292 * t206) * t191 + (t293 * t203 + t291 * t206) * t190 + (t290 * t203 + t289 * t206 + Icges(3,1) + Icges(2,3)) * t194) * t194 / 0.2e1 + t194 * V_base(5) * (t298 * t204 - t297 * t207) + t194 * V_base(4) * (t297 * t204 + t298 * t207) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
