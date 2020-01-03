% Calculate kinetic energy for
% S5RRPRR15
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR15_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR15_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPRR15_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR15_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR15_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR15_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:41:05
% EndTime: 2019-12-31 20:41:09
% DurationCPUTime: 3.77s
% Computational Cost: add. (1002->284), mult. (1625->421), div. (0->0), fcn. (1519->8), ass. (0->138)
t276 = Icges(3,4) + Icges(4,6);
t275 = Icges(3,1) + Icges(4,2);
t274 = -Icges(3,2) - Icges(4,3);
t208 = cos(qJ(2));
t273 = t276 * t208;
t205 = sin(qJ(2));
t272 = t276 * t205;
t271 = Icges(4,4) - Icges(3,5);
t270 = Icges(4,5) - Icges(3,6);
t269 = t274 * t205 + t273;
t268 = t275 * t208 - t272;
t267 = Icges(4,1) + Icges(3,3);
t206 = sin(qJ(1));
t209 = cos(qJ(1));
t266 = t269 * t206 + t270 * t209;
t265 = -t270 * t206 + t269 * t209;
t264 = t268 * t206 + t271 * t209;
t263 = -t271 * t206 + t268 * t209;
t262 = t274 * t208 - t272;
t261 = t275 * t205 + t273;
t260 = t270 * t205 - t271 * t208;
t189 = -qJD(2) * t209 + V_base(5);
t190 = qJD(2) * t206 + V_base(4);
t196 = V_base(6) + qJD(1);
t259 = (t262 * t205 + t261 * t208) * t196 + (-t265 * t205 + t263 * t208) * t190 + (-t266 * t205 + t264 * t208) * t189;
t258 = (-t271 * t205 - t270 * t208) * t196 + (t267 * t206 + t260 * t209) * t190 + (t260 * t206 - t267 * t209) * t189;
t204 = sin(qJ(4));
t254 = pkin(4) * t204;
t253 = t205 * pkin(7);
t207 = cos(qJ(4));
t252 = pkin(4) * t207;
t250 = Icges(2,4) * t206;
t245 = t205 * t206;
t244 = t205 * t209;
t243 = t206 * t207;
t242 = t206 * t208;
t241 = t207 * t209;
t240 = t208 * t209;
t229 = pkin(2) * t208 + qJ(3) * t205;
t158 = t229 * t206;
t187 = pkin(1) * t206 - pkin(6) * t209;
t239 = -t158 - t187;
t238 = qJD(3) * t205;
t237 = qJD(4) * t208;
t236 = qJD(5) * t208;
t235 = V_base(5) * pkin(5) + V_base(1);
t153 = t209 * t237 + t190;
t181 = qJD(4) * t205 + t196;
t182 = pkin(2) * t205 - qJ(3) * t208;
t232 = t189 * t182 + t209 * t238 + t235;
t231 = rSges(3,1) * t208 - rSges(3,2) * t205;
t230 = -rSges(4,2) * t208 + rSges(4,3) * t205;
t152 = t206 * t237 + t189;
t188 = pkin(1) * t209 + pkin(6) * t206;
t222 = -V_base(4) * pkin(5) + t196 * t188 + V_base(2);
t221 = V_base(4) * t187 - t188 * V_base(5) + V_base(3);
t220 = pkin(8) * t208 + t205 * t254;
t159 = t229 * t209;
t217 = t196 * t159 + t206 * t238 + t222;
t216 = -qJD(3) * t208 + t190 * t158 + t221;
t165 = -t209 * pkin(3) + pkin(7) * t242;
t215 = t189 * t253 + (-t165 + t239) * t196 + t232;
t164 = t206 * pkin(3) + pkin(7) * t240;
t214 = t196 * t164 + (-t182 - t253) * t190 + t217;
t213 = t190 * t165 + (-t159 - t164) * t189 + t216;
t203 = qJ(4) + qJ(5);
t201 = Icges(2,4) * t209;
t200 = cos(t203);
t199 = sin(t203);
t186 = rSges(2,1) * t209 - rSges(2,2) * t206;
t185 = rSges(2,1) * t206 + rSges(2,2) * t209;
t184 = rSges(3,1) * t205 + rSges(3,2) * t208;
t183 = -rSges(4,2) * t205 - rSges(4,3) * t208;
t180 = Icges(2,1) * t209 - t250;
t179 = Icges(2,1) * t206 + t201;
t177 = -Icges(2,2) * t206 + t201;
t176 = Icges(2,2) * t209 + t250;
t168 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t167 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t166 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t161 = qJD(5) * t205 + t181;
t157 = t204 * t245 - t241;
t156 = t204 * t209 + t205 * t243;
t155 = t204 * t244 + t243;
t154 = -t204 * t206 + t205 * t241;
t150 = pkin(8) * t205 - t208 * t254;
t148 = t199 * t245 - t200 * t209;
t147 = t199 * t209 + t200 * t245;
t146 = t199 * t244 + t200 * t206;
t145 = -t199 * t206 + t200 * t244;
t144 = -rSges(4,1) * t209 + t206 * t230;
t143 = rSges(4,1) * t206 + t209 * t230;
t142 = rSges(3,3) * t206 + t209 * t231;
t141 = rSges(5,3) * t205 + (-rSges(5,1) * t204 - rSges(5,2) * t207) * t208;
t140 = -rSges(3,3) * t209 + t206 * t231;
t131 = Icges(5,5) * t205 + (-Icges(5,1) * t204 - Icges(5,4) * t207) * t208;
t128 = Icges(5,6) * t205 + (-Icges(5,4) * t204 - Icges(5,2) * t207) * t208;
t125 = Icges(5,3) * t205 + (-Icges(5,5) * t204 - Icges(5,6) * t207) * t208;
t122 = t209 * t236 + t153;
t121 = t206 * t236 + t152;
t120 = rSges(6,3) * t205 + (-rSges(6,1) * t199 - rSges(6,2) * t200) * t208;
t118 = Icges(6,5) * t205 + (-Icges(6,1) * t199 - Icges(6,4) * t200) * t208;
t117 = Icges(6,6) * t205 + (-Icges(6,4) * t199 - Icges(6,2) * t200) * t208;
t116 = Icges(6,3) * t205 + (-Icges(6,5) * t199 - Icges(6,6) * t200) * t208;
t115 = V_base(5) * rSges(2,3) - t185 * t196 + t235;
t114 = t186 * t196 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t113 = t185 * V_base(4) - t186 * V_base(5) + V_base(3);
t112 = t206 * t220 - t209 * t252;
t111 = t206 * t252 + t209 * t220;
t110 = rSges(5,1) * t157 + rSges(5,2) * t156 + rSges(5,3) * t242;
t109 = rSges(5,1) * t155 + rSges(5,2) * t154 + rSges(5,3) * t240;
t108 = Icges(5,1) * t157 + Icges(5,4) * t156 + Icges(5,5) * t242;
t107 = Icges(5,1) * t155 + Icges(5,4) * t154 + Icges(5,5) * t240;
t106 = Icges(5,4) * t157 + Icges(5,2) * t156 + Icges(5,6) * t242;
t105 = Icges(5,4) * t155 + Icges(5,2) * t154 + Icges(5,6) * t240;
t104 = Icges(5,5) * t157 + Icges(5,6) * t156 + Icges(5,3) * t242;
t103 = Icges(5,5) * t155 + Icges(5,6) * t154 + Icges(5,3) * t240;
t102 = rSges(6,1) * t148 + rSges(6,2) * t147 + rSges(6,3) * t242;
t101 = rSges(6,1) * t146 + rSges(6,2) * t145 + rSges(6,3) * t240;
t100 = Icges(6,1) * t148 + Icges(6,4) * t147 + Icges(6,5) * t242;
t99 = Icges(6,1) * t146 + Icges(6,4) * t145 + Icges(6,5) * t240;
t98 = Icges(6,4) * t148 + Icges(6,2) * t147 + Icges(6,6) * t242;
t97 = Icges(6,4) * t146 + Icges(6,2) * t145 + Icges(6,6) * t240;
t96 = Icges(6,5) * t148 + Icges(6,6) * t147 + Icges(6,3) * t242;
t95 = Icges(6,5) * t146 + Icges(6,6) * t145 + Icges(6,3) * t240;
t94 = t184 * t189 + (-t140 - t187) * t196 + t235;
t93 = t142 * t196 - t184 * t190 + t222;
t92 = t140 * t190 - t142 * t189 + t221;
t91 = t183 * t189 + (-t144 + t239) * t196 + t232;
t90 = t143 * t196 + (-t182 - t183) * t190 + t217;
t89 = t144 * t190 + (-t143 - t159) * t189 + t216;
t88 = -t110 * t181 + t141 * t152 + t215;
t87 = t109 * t181 - t141 * t153 + t214;
t86 = -t109 * t152 + t110 * t153 + t213;
t85 = -t102 * t161 - t112 * t181 + t120 * t121 + t150 * t152 + t215;
t84 = t101 * t161 + t111 * t181 - t120 * t122 - t150 * t153 + t214;
t83 = -t101 * t121 + t102 * t122 - t111 * t152 + t112 * t153 + t213;
t1 = m(1) * (t166 ^ 2 + t167 ^ 2 + t168 ^ 2) / 0.2e1 + m(2) * (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + m(3) * (t92 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + m(4) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + m(5) * (t86 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + t153 * ((t103 * t240 + t105 * t154 + t107 * t155) * t153 + (t104 * t240 + t106 * t154 + t108 * t155) * t152 + (t125 * t240 + t128 * t154 + t131 * t155) * t181) / 0.2e1 + t152 * ((t103 * t242 + t105 * t156 + t107 * t157) * t153 + (t104 * t242 + t106 * t156 + t108 * t157) * t152 + (t125 * t242 + t128 * t156 + t131 * t157) * t181) / 0.2e1 + t181 * ((t103 * t153 + t104 * t152 + t125 * t181) * t205 + ((-t105 * t207 - t107 * t204) * t153 + (-t106 * t207 - t108 * t204) * t152 + (-t128 * t207 - t131 * t204) * t181) * t208) / 0.2e1 + m(6) * (t83 ^ 2 + t84 ^ 2 + t85 ^ 2) / 0.2e1 + t122 * ((t145 * t97 + t146 * t99 + t240 * t95) * t122 + (t100 * t146 + t145 * t98 + t240 * t96) * t121 + (t116 * t240 + t117 * t145 + t118 * t146) * t161) / 0.2e1 + t121 * ((t147 * t97 + t148 * t99 + t242 * t95) * t122 + (t100 * t148 + t147 * t98 + t242 * t96) * t121 + (t116 * t242 + t117 * t147 + t118 * t148) * t161) / 0.2e1 + t161 * ((t116 * t161 + t96 * t121 + t95 * t122) * t205 + ((-t199 * t99 - t200 * t97) * t122 + (-t100 * t199 - t200 * t98) * t121 + (-t117 * t200 - t118 * t199) * t161) * t208) / 0.2e1 + (t259 * t206 - t258 * t209) * t189 / 0.2e1 + (t258 * t206 + t259 * t209) * t190 / 0.2e1 + ((-t176 * t206 + t179 * t209 + Icges(1,4)) * V_base(5) + (-t177 * t206 + t180 * t209 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t176 * t209 + t179 * t206 + Icges(1,2)) * V_base(5) + (t177 * t209 + t180 * t206 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t263 * t205 + t265 * t208) * t190 + (t264 * t205 + t266 * t208) * t189 + (t261 * t205 - t262 * t208 + Icges(2,3)) * t196) * t196 / 0.2e1 + t196 * V_base(4) * (Icges(2,5) * t209 - Icges(2,6) * t206) + V_base(5) * t196 * (Icges(2,5) * t206 + Icges(2,6) * t209) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
