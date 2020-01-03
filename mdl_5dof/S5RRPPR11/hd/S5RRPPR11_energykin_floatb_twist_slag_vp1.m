% Calculate kinetic energy for
% S5RRPPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPPR11_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR11_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPPR11_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR11_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR11_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR11_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:46:20
% EndTime: 2019-12-31 19:46:23
% DurationCPUTime: 3.36s
% Computational Cost: add. (966->287), mult. (1565->403), div. (0->0), fcn. (1459->8), ass. (0->139)
t279 = Icges(3,4) + Icges(4,6);
t278 = Icges(3,1) + Icges(4,2);
t277 = Icges(3,2) + Icges(4,3);
t208 = cos(qJ(2));
t276 = t279 * t208;
t206 = sin(qJ(2));
t275 = t279 * t206;
t274 = Icges(4,4) - Icges(3,5);
t273 = Icges(4,5) - Icges(3,6);
t272 = t277 * t206 - t276;
t271 = t278 * t208 - t275;
t270 = Icges(4,1) + Icges(3,3);
t207 = sin(qJ(1));
t209 = cos(qJ(1));
t269 = t272 * t207 - t273 * t209;
t268 = t273 * t207 + t272 * t209;
t267 = t271 * t207 + t274 * t209;
t266 = -t274 * t207 + t271 * t209;
t265 = -t277 * t208 - t275;
t264 = t278 * t206 + t276;
t263 = t273 * t206 - t274 * t208;
t188 = -qJD(2) * t209 + V_base(5);
t189 = qJD(2) * t207 + V_base(4);
t197 = V_base(6) + qJD(1);
t262 = (t265 * t206 + t264 * t208) * t197 + (t268 * t206 + t266 * t208) * t189 + (t269 * t206 + t267 * t208) * t188;
t261 = (-t274 * t206 - t273 * t208) * t197 + (t270 * t207 + t263 * t209) * t189 + (t263 * t207 - t270 * t209) * t188;
t203 = sin(pkin(8));
t257 = pkin(4) * t203;
t204 = cos(pkin(8));
t256 = pkin(4) * t204;
t255 = Icges(2,4) * t207;
t250 = qJ(4) * t206;
t249 = t206 * t209;
t202 = pkin(8) + qJ(5);
t195 = sin(t202);
t248 = t207 * t195;
t196 = cos(t202);
t247 = t207 * t196;
t246 = t207 * t203;
t245 = t207 * t204;
t244 = t207 * t208;
t243 = t208 * t209;
t228 = pkin(2) * t208 + qJ(3) * t206;
t158 = t228 * t207;
t186 = t207 * pkin(1) - pkin(6) * t209;
t241 = -t158 - t186;
t159 = t228 * t209;
t163 = t207 * pkin(3) + qJ(4) * t243;
t240 = -t159 - t163;
t239 = qJD(3) * t206;
t238 = qJD(4) * t208;
t237 = qJD(5) * t208;
t236 = V_base(5) * pkin(5) + V_base(1);
t164 = -pkin(3) * t209 + qJ(4) * t244;
t233 = -t164 + t241;
t181 = pkin(2) * t206 - qJ(3) * t208;
t232 = -t181 - t250;
t231 = t188 * t181 + t209 * t239 + t236;
t230 = rSges(3,1) * t208 - rSges(3,2) * t206;
t229 = -rSges(4,2) * t208 + rSges(4,3) * t206;
t187 = pkin(1) * t209 + t207 * pkin(6);
t221 = -V_base(4) * pkin(5) + t197 * t187 + V_base(2);
t220 = V_base(4) * t186 - t187 * V_base(5) + V_base(3);
t219 = t188 * t250 + t209 * t238 + t231;
t218 = pkin(7) * t208 + t206 * t257;
t215 = t197 * t159 + t207 * t239 + t221;
t214 = -qJD(3) * t208 + t189 * t158 + t220;
t213 = t197 * t163 + t207 * t238 + t215;
t212 = qJD(4) * t206 + t189 * t164 + t214;
t200 = Icges(2,4) * t209;
t185 = rSges(2,1) * t209 - t207 * rSges(2,2);
t184 = t207 * rSges(2,1) + rSges(2,2) * t209;
t183 = rSges(3,1) * t206 + rSges(3,2) * t208;
t182 = -rSges(4,2) * t206 - rSges(4,3) * t208;
t180 = qJD(5) * t206 + t197;
t179 = Icges(2,1) * t209 - t255;
t178 = Icges(2,1) * t207 + t200;
t176 = -Icges(2,2) * t207 + t200;
t175 = Icges(2,2) * t209 + t255;
t167 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t166 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t165 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t157 = t209 * t237 + t189;
t156 = t207 * t237 + t188;
t155 = -t204 * t209 + t206 * t246;
t154 = t203 * t209 + t206 * t245;
t153 = t203 * t249 + t245;
t152 = t204 * t249 - t246;
t149 = pkin(7) * t206 - t208 * t257;
t148 = -rSges(4,1) * t209 + t207 * t229;
t147 = t207 * rSges(4,1) + t209 * t229;
t146 = t207 * rSges(3,3) + t209 * t230;
t145 = -rSges(3,3) * t209 + t207 * t230;
t132 = -t196 * t209 + t206 * t248;
t131 = t195 * t209 + t206 * t247;
t130 = t195 * t249 + t247;
t129 = t196 * t249 - t248;
t128 = rSges(5,3) * t206 + (-rSges(5,1) * t203 - rSges(5,2) * t204) * t208;
t126 = Icges(5,5) * t206 + (-Icges(5,1) * t203 - Icges(5,4) * t204) * t208;
t125 = Icges(5,6) * t206 + (-Icges(5,4) * t203 - Icges(5,2) * t204) * t208;
t124 = Icges(5,3) * t206 + (-Icges(5,5) * t203 - Icges(5,6) * t204) * t208;
t121 = rSges(6,3) * t206 + (-rSges(6,1) * t195 - rSges(6,2) * t196) * t208;
t120 = Icges(6,5) * t206 + (-Icges(6,1) * t195 - Icges(6,4) * t196) * t208;
t119 = Icges(6,6) * t206 + (-Icges(6,4) * t195 - Icges(6,2) * t196) * t208;
t118 = Icges(6,3) * t206 + (-Icges(6,5) * t195 - Icges(6,6) * t196) * t208;
t117 = V_base(5) * rSges(2,3) - t184 * t197 + t236;
t116 = t185 * t197 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t115 = t184 * V_base(4) - t185 * V_base(5) + V_base(3);
t114 = t207 * t218 - t209 * t256;
t113 = t207 * t256 + t209 * t218;
t112 = rSges(5,1) * t155 + rSges(5,2) * t154 + rSges(5,3) * t244;
t111 = t153 * rSges(5,1) + t152 * rSges(5,2) + rSges(5,3) * t243;
t110 = Icges(5,1) * t155 + Icges(5,4) * t154 + Icges(5,5) * t244;
t109 = Icges(5,1) * t153 + Icges(5,4) * t152 + Icges(5,5) * t243;
t108 = Icges(5,4) * t155 + Icges(5,2) * t154 + Icges(5,6) * t244;
t107 = Icges(5,4) * t153 + Icges(5,2) * t152 + Icges(5,6) * t243;
t106 = Icges(5,5) * t155 + Icges(5,6) * t154 + Icges(5,3) * t244;
t105 = Icges(5,5) * t153 + Icges(5,6) * t152 + Icges(5,3) * t243;
t104 = rSges(6,1) * t132 + rSges(6,2) * t131 + rSges(6,3) * t244;
t103 = t130 * rSges(6,1) + t129 * rSges(6,2) + rSges(6,3) * t243;
t102 = Icges(6,1) * t132 + Icges(6,4) * t131 + Icges(6,5) * t244;
t101 = Icges(6,1) * t130 + Icges(6,4) * t129 + Icges(6,5) * t243;
t100 = Icges(6,4) * t132 + Icges(6,2) * t131 + Icges(6,6) * t244;
t99 = Icges(6,4) * t130 + Icges(6,2) * t129 + Icges(6,6) * t243;
t98 = Icges(6,5) * t132 + Icges(6,6) * t131 + Icges(6,3) * t244;
t97 = Icges(6,5) * t130 + Icges(6,6) * t129 + Icges(6,3) * t243;
t96 = t183 * t188 + (-t145 - t186) * t197 + t236;
t95 = t146 * t197 - t183 * t189 + t221;
t94 = t145 * t189 - t146 * t188 + t220;
t93 = t182 * t188 + (-t148 + t241) * t197 + t231;
t92 = t147 * t197 + (-t181 - t182) * t189 + t215;
t91 = t148 * t189 + (-t147 - t159) * t188 + t214;
t90 = t128 * t188 + (-t112 + t233) * t197 + t219;
t89 = t111 * t197 + (-t128 + t232) * t189 + t213;
t88 = t112 * t189 + (-t111 + t240) * t188 + t212;
t87 = -t104 * t180 + t121 * t156 + t149 * t188 + (-t114 + t233) * t197 + t219;
t86 = t103 * t180 + t113 * t197 - t121 * t157 + (-t149 + t232) * t189 + t213;
t85 = -t103 * t156 + t104 * t157 + t114 * t189 + (-t113 + t240) * t188 + t212;
t1 = m(1) * (t165 ^ 2 + t166 ^ 2 + t167 ^ 2) / 0.2e1 + m(2) * (t115 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + m(3) * (t94 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + m(4) * (t91 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + m(5) * (t88 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + m(6) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + t157 * ((t130 * t101 + t129 * t99 + t97 * t243) * t157 + (t129 * t100 + t130 * t102 + t243 * t98) * t156 + (t118 * t243 + t129 * t119 + t130 * t120) * t180) / 0.2e1 + t156 * ((t101 * t132 + t131 * t99 + t244 * t97) * t157 + (t100 * t131 + t102 * t132 + t98 * t244) * t156 + (t118 * t244 + t119 * t131 + t120 * t132) * t180) / 0.2e1 + t180 * ((t118 * t180 + t98 * t156 + t97 * t157) * t206 + ((-t101 * t195 - t196 * t99) * t157 + (-t100 * t196 - t102 * t195) * t156 + (-t119 * t196 - t120 * t195) * t180) * t208) / 0.2e1 + ((-t207 * t175 + t178 * t209 + Icges(1,4)) * V_base(5) + (-t207 * t176 + t179 * t209 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t175 * t209 + t207 * t178 + Icges(1,2)) * V_base(5) + (t176 * t209 + t207 * t179 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t105 * t244 + t107 * t154 + t109 * t155) * t189 + (t106 * t244 + t108 * t154 + t155 * t110) * t188 + (t124 * t244 + t125 * t154 + t126 * t155) * t197 - t261 * t209 + t262 * t207) * t188 / 0.2e1 + ((t105 * t243 + t152 * t107 + t153 * t109) * t189 + (t106 * t243 + t152 * t108 + t153 * t110) * t188 + (t124 * t243 + t152 * t125 + t153 * t126) * t197 + t262 * t209 + t261 * t207) * t189 / 0.2e1 + (((-t107 * t204 - t109 * t203 - t268) * t208 + (t105 + t266) * t206) * t189 + ((-t108 * t204 - t110 * t203 - t269) * t208 + (t106 + t267) * t206) * t188 + (Icges(2,3) + (-t125 * t204 - t126 * t203 - t265) * t208 + (t124 + t264) * t206) * t197) * t197 / 0.2e1 + t197 * V_base(4) * (Icges(2,5) * t209 - Icges(2,6) * t207) + V_base(5) * t197 * (Icges(2,5) * t207 + Icges(2,6) * t209) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
