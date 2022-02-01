% Calculate kinetic energy for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% m [6x1]
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
% Datum: 2022-01-20 11:18
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPRR6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR6_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR6_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR6_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:16:45
% EndTime: 2022-01-20 11:16:48
% DurationCPUTime: 2.50s
% Computational Cost: add. (1410->294), mult. (1396->427), div. (0->0), fcn. (1288->10), ass. (0->141)
t243 = -pkin(5) - pkin(6);
t199 = sin(qJ(1));
t241 = pkin(1) * t199;
t201 = cos(qJ(1));
t240 = pkin(1) * t201;
t200 = cos(qJ(4));
t239 = pkin(4) * t200;
t237 = Icges(2,4) * t199;
t195 = qJ(1) + qJ(2);
t188 = sin(t195);
t236 = Icges(3,4) * t188;
t196 = sin(pkin(9));
t235 = Icges(4,4) * t196;
t197 = cos(pkin(9));
t234 = Icges(4,4) * t197;
t233 = t188 * t196;
t232 = t188 * t197;
t198 = sin(qJ(4));
t231 = t188 * t198;
t190 = cos(t195);
t230 = t190 * t196;
t229 = t190 * t197;
t228 = t190 * t198;
t227 = t197 * t198;
t226 = t197 * t200;
t225 = qJD(4) * t196;
t224 = qJD(5) * t196;
t186 = V_base(6) + qJD(1);
t223 = t186 * t240 + V_base(2);
t222 = V_base(4) * t241 + V_base(3);
t221 = V_base(5) * pkin(5) + V_base(1);
t160 = t190 * t225 + V_base(4);
t159 = t188 * t225 + V_base(5);
t156 = pkin(2) * t190 + qJ(3) * t188;
t218 = -t156 - t240;
t154 = pkin(2) * t188 - qJ(3) * t190;
t217 = V_base(4) * t154 + t222;
t184 = qJD(2) + t186;
t216 = pkin(3) * t197 + pkin(7) * t196;
t215 = rSges(4,1) * t197 - rSges(4,2) * t196;
t214 = Icges(4,1) * t197 - t235;
t213 = -Icges(4,2) * t196 + t234;
t212 = Icges(4,5) * t197 - Icges(4,6) * t196;
t211 = -qJD(3) * t190 + t184 * t156 + t223;
t210 = V_base(5) * pkin(6) - t186 * t241 + t221;
t209 = pkin(8) * t196 + t197 * t239;
t208 = qJD(3) * t188 + t210;
t207 = (-Icges(4,3) * t190 + t188 * t212) * V_base(5) + (Icges(4,3) * t188 + t190 * t212) * V_base(4) + (Icges(4,5) * t196 + Icges(4,6) * t197) * t184;
t144 = t216 * t188;
t145 = t216 * t190;
t206 = V_base(4) * t144 + (-t145 + t218) * V_base(5) + t217;
t169 = t196 * pkin(3) - t197 * pkin(7);
t205 = t184 * t145 + (-t169 + t243) * V_base(4) + t211;
t204 = V_base(5) * t169 + (-t144 - t154) * t184 + t208;
t116 = -Icges(4,6) * t190 + t188 * t213;
t117 = Icges(4,6) * t188 + t190 * t213;
t118 = -Icges(4,5) * t190 + t188 * t214;
t119 = Icges(4,5) * t188 + t190 * t214;
t166 = Icges(4,2) * t197 + t235;
t167 = Icges(4,1) * t196 + t234;
t203 = (-t117 * t196 + t119 * t197) * V_base(4) + (-t116 * t196 + t118 * t197) * V_base(5) + (-t166 * t196 + t167 * t197) * t184;
t194 = qJ(4) + qJ(5);
t191 = Icges(2,4) * t201;
t189 = cos(t194);
t187 = sin(t194);
t183 = Icges(3,4) * t190;
t179 = rSges(2,1) * t201 - rSges(2,2) * t199;
t178 = rSges(2,1) * t199 + rSges(2,2) * t201;
t175 = Icges(2,1) * t201 - t237;
t174 = Icges(2,1) * t199 + t191;
t173 = -Icges(2,2) * t199 + t191;
t172 = Icges(2,2) * t201 + t237;
t168 = rSges(4,1) * t196 + rSges(4,2) * t197;
t164 = -qJD(4) * t197 + t184;
t163 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t162 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t161 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t157 = rSges(3,1) * t190 - rSges(3,2) * t188;
t155 = rSges(3,1) * t188 + rSges(3,2) * t190;
t153 = Icges(3,1) * t190 - t236;
t152 = Icges(3,1) * t188 + t183;
t151 = -Icges(3,2) * t188 + t183;
t150 = Icges(3,2) * t190 + t236;
t149 = Icges(3,5) * t190 - Icges(3,6) * t188;
t148 = Icges(3,5) * t188 + Icges(3,6) * t190;
t147 = (-qJD(4) - qJD(5)) * t197 + t184;
t143 = t190 * t226 + t231;
t142 = t188 * t200 - t190 * t227;
t141 = t188 * t226 - t228;
t140 = -t188 * t227 - t190 * t200;
t138 = -rSges(5,3) * t197 + (rSges(5,1) * t200 - rSges(5,2) * t198) * t196;
t137 = -Icges(5,5) * t197 + (Icges(5,1) * t200 - Icges(5,4) * t198) * t196;
t136 = -Icges(5,6) * t197 + (Icges(5,4) * t200 - Icges(5,2) * t198) * t196;
t135 = -Icges(5,3) * t197 + (Icges(5,5) * t200 - Icges(5,6) * t198) * t196;
t133 = t190 * t224 + t160;
t132 = t188 * t224 + t159;
t131 = t187 * t188 + t189 * t229;
t130 = -t187 * t229 + t188 * t189;
t129 = -t187 * t190 + t189 * t232;
t128 = -t187 * t232 - t189 * t190;
t127 = -rSges(6,3) * t197 + (rSges(6,1) * t189 - rSges(6,2) * t187) * t196;
t126 = -Icges(6,5) * t197 + (Icges(6,1) * t189 - Icges(6,4) * t187) * t196;
t125 = -Icges(6,6) * t197 + (Icges(6,4) * t189 - Icges(6,2) * t187) * t196;
t124 = -Icges(6,3) * t197 + (Icges(6,5) * t189 - Icges(6,6) * t187) * t196;
t123 = -pkin(8) * t197 + t196 * t239;
t122 = rSges(4,3) * t188 + t190 * t215;
t121 = -rSges(4,3) * t190 + t188 * t215;
t113 = V_base(5) * rSges(2,3) - t178 * t186 + t221;
t112 = t179 * t186 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t111 = t178 * V_base(4) - t179 * V_base(5) + V_base(3);
t110 = V_base(5) * rSges(3,3) - t155 * t184 + t210;
t109 = t157 * t184 + (-rSges(3,3) + t243) * V_base(4) + t223;
t108 = t155 * V_base(4) + (-t157 - t240) * V_base(5) + t222;
t107 = rSges(5,1) * t143 + rSges(5,2) * t142 + rSges(5,3) * t230;
t106 = rSges(5,1) * t141 + rSges(5,2) * t140 + rSges(5,3) * t233;
t105 = pkin(4) * t231 + t190 * t209;
t104 = -pkin(4) * t228 + t188 * t209;
t103 = Icges(5,1) * t143 + Icges(5,4) * t142 + Icges(5,5) * t230;
t102 = Icges(5,1) * t141 + Icges(5,4) * t140 + Icges(5,5) * t233;
t101 = Icges(5,4) * t143 + Icges(5,2) * t142 + Icges(5,6) * t230;
t100 = Icges(5,4) * t141 + Icges(5,2) * t140 + Icges(5,6) * t233;
t99 = Icges(5,5) * t143 + Icges(5,6) * t142 + Icges(5,3) * t230;
t98 = Icges(5,5) * t141 + Icges(5,6) * t140 + Icges(5,3) * t233;
t97 = rSges(6,1) * t131 + rSges(6,2) * t130 + rSges(6,3) * t230;
t96 = rSges(6,1) * t129 + rSges(6,2) * t128 + rSges(6,3) * t233;
t95 = Icges(6,1) * t131 + Icges(6,4) * t130 + Icges(6,5) * t230;
t94 = Icges(6,1) * t129 + Icges(6,4) * t128 + Icges(6,5) * t233;
t93 = Icges(6,4) * t131 + Icges(6,2) * t130 + Icges(6,6) * t230;
t92 = Icges(6,4) * t129 + Icges(6,2) * t128 + Icges(6,6) * t233;
t91 = Icges(6,5) * t131 + Icges(6,6) * t130 + Icges(6,3) * t230;
t90 = Icges(6,5) * t129 + Icges(6,6) * t128 + Icges(6,3) * t233;
t89 = t168 * V_base(5) + (-t121 - t154) * t184 + t208;
t88 = t122 * t184 + (-t168 + t243) * V_base(4) + t211;
t87 = t121 * V_base(4) + (-t122 + t218) * V_base(5) + t217;
t86 = -t106 * t164 + t138 * t159 + t204;
t85 = t107 * t164 - t138 * t160 + t205;
t84 = t106 * t160 - t107 * t159 + t206;
t83 = -t104 * t164 + t123 * t159 + t127 * t132 - t147 * t96 + t204;
t82 = t105 * t164 - t123 * t160 - t127 * t133 + t147 * t97 + t205;
t81 = t104 * t160 - t105 * t159 - t132 * t97 + t133 * t96 + t206;
t1 = m(1) * (t161 ^ 2 + t162 ^ 2 + t163 ^ 2) / 0.2e1 + m(2) * (t111 ^ 2 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + m(3) * (t108 ^ 2 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + m(4) * (t87 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + m(5) * (t84 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + t160 * ((t142 * t101 + t143 * t103 + t99 * t230) * t160 + (t100 * t142 + t102 * t143 + t230 * t98) * t159 + (t135 * t230 + t136 * t142 + t137 * t143) * t164) / 0.2e1 + t159 * ((t101 * t140 + t103 * t141 + t233 * t99) * t160 + (t140 * t100 + t141 * t102 + t98 * t233) * t159 + (t135 * t233 + t136 * t140 + t137 * t141) * t164) / 0.2e1 + t164 * ((-t135 * t164 - t98 * t159 - t99 * t160) * t197 + ((-t101 * t198 + t103 * t200) * t160 + (-t100 * t198 + t102 * t200) * t159 + (-t136 * t198 + t137 * t200) * t164) * t196) / 0.2e1 + m(6) * (t81 ^ 2 + t82 ^ 2 + t83 ^ 2) / 0.2e1 + t133 * ((t130 * t93 + t131 * t95 + t91 * t230) * t133 + (t130 * t92 + t131 * t94 + t230 * t90) * t132 + (t124 * t230 + t125 * t130 + t126 * t131) * t147) / 0.2e1 + t132 * ((t128 * t93 + t129 * t95 + t233 * t91) * t133 + (t128 * t92 + t129 * t94 + t90 * t233) * t132 + (t124 * t233 + t125 * t128 + t126 * t129) * t147) / 0.2e1 + t147 * ((-t124 * t147 - t90 * t132 - t91 * t133) * t197 + ((-t187 * t93 + t189 * t95) * t133 + (-t187 * t92 + t189 * t94) * t132 + (-t125 * t187 + t126 * t189) * t147) * t196) / 0.2e1 + ((t116 * t197 + t118 * t196 + t148) * V_base(5) + (t117 * t197 + t119 * t196 + t149) * V_base(4) + (t197 * t166 + t196 * t167 + Icges(3,3)) * t184) * t184 / 0.2e1 + (t149 * t184 + t207 * t188 + t203 * t190 + (-t150 * t188 + t152 * t190 - t172 * t199 + t174 * t201 + Icges(1,4)) * V_base(5) + (-t188 * t151 + t190 * t153 - t199 * t173 + t201 * t175 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t148 * t184 + t203 * t188 - t207 * t190 + (t190 * t150 + t188 * t152 + t201 * t172 + t199 * t174 + Icges(1,2)) * V_base(5) + (t151 * t190 + t153 * t188 + t173 * t201 + t175 * t199 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t199 + Icges(2,6) * t201) * V_base(5) + (Icges(2,5) * t201 - Icges(2,6) * t199) * V_base(4) + Icges(2,3) * t186 / 0.2e1) * t186;
T = t1;
