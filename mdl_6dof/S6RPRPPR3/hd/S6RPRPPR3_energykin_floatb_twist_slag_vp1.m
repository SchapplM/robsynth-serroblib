% Calculate kinetic energy for
% S6RPRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
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
% Datum: 2019-03-09 02:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPPR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR3_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPPR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR3_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR3_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPPR3_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:44:02
% EndTime: 2019-03-09 02:44:05
% DurationCPUTime: 2.96s
% Computational Cost: add. (1404->271), mult. (1455->366), div. (0->0), fcn. (1235->8), ass. (0->139)
t298 = Icges(4,4) + Icges(6,4) - Icges(5,5);
t297 = Icges(4,1) + Icges(5,1) + Icges(6,2);
t296 = Icges(6,1) + Icges(4,2) + Icges(5,3);
t213 = sin(qJ(3));
t295 = t298 * t213;
t216 = cos(qJ(3));
t294 = t298 * t216;
t293 = Icges(5,4) + Icges(4,5) + Icges(6,6);
t292 = Icges(6,5) + Icges(4,6) - Icges(5,6);
t291 = t296 * t213 - t294;
t290 = t297 * t216 - t295;
t289 = Icges(5,2) + Icges(4,3) + Icges(6,3);
t211 = qJ(1) + pkin(9);
t205 = sin(t211);
t206 = cos(t211);
t288 = t205 * t291 + t206 * t292;
t287 = -t205 * t292 + t206 * t291;
t286 = t290 * t205 - t206 * t293;
t285 = t205 * t293 + t290 * t206;
t284 = -t296 * t216 - t295;
t283 = t297 * t213 + t294;
t282 = -t292 * t213 + t216 * t293;
t173 = -qJD(3) * t206 + V_base(5);
t174 = qJD(3) * t205 + V_base(4);
t207 = V_base(6) + qJD(1);
t279 = (t213 * t284 + t216 * t283) * t207 + (t213 * t287 + t216 * t285) * t174 + (t213 * t288 + t216 * t286) * t173;
t278 = (t213 * t293 + t292 * t216) * t207 + (t205 * t289 + t282 * t206) * t174 + (t282 * t205 - t206 * t289) * t173;
t214 = sin(qJ(1));
t274 = pkin(1) * t214;
t217 = cos(qJ(1));
t273 = pkin(1) * t217;
t272 = pkin(4) * t213;
t271 = -pkin(6) - qJ(2);
t270 = Icges(2,4) * t214;
t269 = Icges(3,4) * t205;
t262 = t205 * t216;
t261 = t206 * t216;
t212 = sin(qJ(6));
t260 = t212 * t213;
t215 = cos(qJ(6));
t259 = t213 * t215;
t242 = pkin(3) * t216 + qJ(4) * t213;
t151 = t242 * t206;
t158 = pkin(4) * t261 - qJ(5) * t205;
t258 = -t151 - t158;
t257 = qJD(4) * t213;
t256 = qJD(6) * t216;
t255 = t207 * t273 + V_base(2);
t254 = V_base(5) * pkin(6) + V_base(1);
t168 = pkin(2) * t205 - pkin(7) * t206;
t251 = -t168 - t274;
t193 = pkin(3) * t213 - qJ(4) * t216;
t250 = -t193 - t272;
t249 = V_base(5) * qJ(2) + t254;
t248 = V_base(4) * t274 + qJD(2) + V_base(3);
t150 = t242 * t205;
t247 = -t150 + t251;
t246 = pkin(5) * t213 + pkin(8) * t216;
t245 = rSges(4,1) * t216 - rSges(4,2) * t213;
t244 = rSges(5,1) * t216 + rSges(5,3) * t213;
t243 = rSges(6,1) * t213 - rSges(6,2) * t216;
t157 = pkin(4) * t262 + qJ(5) * t206;
t232 = -t157 + t247;
t231 = t173 * t193 + t206 * t257 + t249;
t169 = pkin(2) * t206 + pkin(7) * t205;
t227 = t207 * t169 + t271 * V_base(4) + t255;
t226 = -qJD(5) * t205 + t173 * t272 + t231;
t225 = t207 * t151 + t205 * t257 + t227;
t224 = V_base(4) * t168 + (-t169 - t273) * V_base(5) + t248;
t223 = qJD(5) * t206 + t207 * t158 + t225;
t222 = -qJD(4) * t216 + t174 * t150 + t224;
t221 = t174 * t157 + t222;
t209 = Icges(2,4) * t217;
t204 = Icges(3,4) * t206;
t199 = -pkin(5) * t216 + pkin(8) * t213;
t198 = rSges(2,1) * t217 - t214 * rSges(2,2);
t197 = -rSges(6,1) * t216 - rSges(6,2) * t213;
t196 = t214 * rSges(2,1) + rSges(2,2) * t217;
t195 = rSges(4,1) * t213 + rSges(4,2) * t216;
t194 = rSges(5,1) * t213 - rSges(5,3) * t216;
t192 = qJD(6) * t213 + t207;
t189 = Icges(2,1) * t217 - t270;
t188 = Icges(2,1) * t214 + t209;
t184 = -Icges(2,2) * t214 + t209;
t183 = Icges(2,2) * t217 + t270;
t172 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t171 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t170 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t167 = rSges(3,1) * t206 - rSges(3,2) * t205;
t166 = rSges(3,1) * t205 + rSges(3,2) * t206;
t165 = Icges(3,1) * t206 - t269;
t164 = Icges(3,1) * t205 + t204;
t163 = -Icges(3,2) * t205 + t204;
t162 = Icges(3,2) * t206 + t269;
t155 = t246 * t206;
t154 = t246 * t205;
t152 = rSges(7,3) * t213 + (-rSges(7,1) * t215 + rSges(7,2) * t212) * t216;
t149 = Icges(7,5) * t213 + (-Icges(7,1) * t215 + Icges(7,4) * t212) * t216;
t148 = Icges(7,6) * t213 + (-Icges(7,4) * t215 + Icges(7,2) * t212) * t216;
t147 = Icges(7,3) * t213 + (-Icges(7,5) * t215 + Icges(7,6) * t212) * t216;
t146 = -t205 * t212 + t206 * t259;
t145 = -t205 * t215 - t206 * t260;
t144 = t205 * t259 + t206 * t212;
t143 = -t205 * t260 + t206 * t215;
t142 = t206 * t256 + t174;
t141 = t205 * t256 + t173;
t137 = rSges(4,3) * t205 + t206 * t245;
t136 = rSges(5,2) * t205 + t206 * t244;
t135 = -rSges(6,3) * t205 + t206 * t243;
t134 = -rSges(4,3) * t206 + t205 * t245;
t133 = -rSges(5,2) * t206 + t205 * t244;
t132 = rSges(6,3) * t206 + t205 * t243;
t131 = V_base(5) * rSges(2,3) - t196 * t207 + t254;
t130 = t198 * t207 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t110 = t196 * V_base(4) - t198 * V_base(5) + V_base(3);
t108 = V_base(5) * rSges(3,3) + (-t166 - t274) * t207 + t249;
t107 = t167 * t207 + (-rSges(3,3) + t271) * V_base(4) + t255;
t106 = V_base(4) * t166 + (-t167 - t273) * V_base(5) + t248;
t105 = rSges(7,1) * t146 + rSges(7,2) * t145 + rSges(7,3) * t261;
t104 = rSges(7,1) * t144 + rSges(7,2) * t143 + rSges(7,3) * t262;
t103 = Icges(7,1) * t146 + Icges(7,4) * t145 + Icges(7,5) * t261;
t102 = Icges(7,1) * t144 + Icges(7,4) * t143 + Icges(7,5) * t262;
t101 = Icges(7,4) * t146 + Icges(7,2) * t145 + Icges(7,6) * t261;
t100 = Icges(7,4) * t144 + Icges(7,2) * t143 + Icges(7,6) * t262;
t99 = Icges(7,5) * t146 + Icges(7,6) * t145 + Icges(7,3) * t261;
t98 = Icges(7,5) * t144 + Icges(7,6) * t143 + Icges(7,3) * t262;
t97 = t173 * t195 + (-t134 + t251) * t207 + t249;
t96 = t137 * t207 - t174 * t195 + t227;
t95 = t174 * t134 - t173 * t137 + t224;
t94 = t173 * t194 + (-t133 + t247) * t207 + t231;
t93 = t136 * t207 + (-t193 - t194) * t174 + t225;
t92 = t173 * t197 + (-t132 + t232) * t207 + t226;
t91 = t135 * t207 + (-t197 + t250) * t174 + t223;
t90 = t174 * t133 + (-t136 - t151) * t173 + t222;
t89 = t174 * t132 + (-t135 + t258) * t173 + t221;
t88 = -t104 * t192 + t141 * t152 + t173 * t199 + (-t154 + t232) * t207 + t226;
t87 = t105 * t192 - t142 * t152 + t155 * t207 + (-t199 + t250) * t174 + t223;
t86 = t142 * t104 - t141 * t105 + t174 * t154 + (-t155 + t258) * t173 + t221;
t1 = t142 * ((t145 * t101 + t146 * t103 + t99 * t261) * t142 + (t100 * t145 + t102 * t146 + t261 * t98) * t141 + (t145 * t148 + t146 * t149 + t147 * t261) * t192) / 0.2e1 + t141 * ((t101 * t143 + t103 * t144 + t262 * t99) * t142 + (t143 * t100 + t144 * t102 + t98 * t262) * t141 + (t143 * t148 + t144 * t149 + t147 * t262) * t192) / 0.2e1 + t192 * ((t141 * t98 + t142 * t99 + t147 * t192) * t213 + ((t101 * t212 - t103 * t215) * t142 + (t100 * t212 - t102 * t215) * t141 + (t148 * t212 - t149 * t215) * t192) * t216) / 0.2e1 + m(1) * (t170 ^ 2 + t171 ^ 2 + t172 ^ 2) / 0.2e1 + m(2) * (t110 ^ 2 + t130 ^ 2 + t131 ^ 2) / 0.2e1 + m(3) * (t106 ^ 2 + t107 ^ 2 + t108 ^ 2) / 0.2e1 + m(5) * (t90 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + m(4) * (t95 ^ 2 + t96 ^ 2 + t97 ^ 2) / 0.2e1 + m(7) * (t86 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(6) * (t89 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + (t279 * t205 - t278 * t206) * t173 / 0.2e1 + (t278 * t205 + t279 * t206) * t174 / 0.2e1 + ((-t162 * t205 + t164 * t206 - t214 * t183 + t188 * t217 + Icges(1,4)) * V_base(5) + (-t163 * t205 + t165 * t206 - t214 * t184 + t189 * t217 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t162 * t206 + t164 * t205 + t183 * t217 + t214 * t188 + Icges(1,2)) * V_base(5) + (t163 * t206 + t165 * t205 + t184 * t217 + t214 * t189 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t213 * t285 - t216 * t287) * t174 + (t213 * t286 - t216 * t288) * t173 + (t213 * t283 - t216 * t284 + Icges(2,3) + Icges(3,3)) * t207) * t207 / 0.2e1 + t207 * V_base(5) * (Icges(2,5) * t214 + Icges(3,5) * t205 + Icges(2,6) * t217 + Icges(3,6) * t206) + t207 * V_base(4) * (Icges(2,5) * t217 + Icges(3,5) * t206 - Icges(2,6) * t214 - Icges(3,6) * t205) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
