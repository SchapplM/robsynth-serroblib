% Calculate kinetic energy for
% S5PRRPR8
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRPR8_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR8_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRRPR8_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR8_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR8_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR8_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:42:15
% EndTime: 2019-12-31 17:42:17
% DurationCPUTime: 2.07s
% Computational Cost: add. (1383->280), mult. (1408->402), div. (0->0), fcn. (1246->10), ass. (0->144)
t263 = Icges(4,3) + Icges(5,3);
t194 = qJ(2) + qJ(3);
t188 = pkin(9) + t194;
t181 = sin(t188);
t182 = cos(t188);
t189 = sin(t194);
t190 = cos(t194);
t262 = Icges(4,5) * t190 + Icges(5,5) * t182 - Icges(4,6) * t189 - Icges(5,6) * t181;
t195 = sin(pkin(8));
t196 = cos(pkin(8));
t245 = Icges(5,4) * t182;
t218 = -Icges(5,2) * t181 + t245;
t111 = -Icges(5,6) * t196 + t218 * t195;
t112 = Icges(5,6) * t195 + t218 * t196;
t246 = Icges(5,4) * t181;
t221 = Icges(5,1) * t182 - t246;
t113 = -Icges(5,5) * t196 + t221 * t195;
t114 = Icges(5,5) * t195 + t221 * t196;
t247 = Icges(4,4) * t190;
t219 = -Icges(4,2) * t189 + t247;
t125 = -Icges(4,6) * t196 + t219 * t195;
t126 = Icges(4,6) * t195 + t219 * t196;
t248 = Icges(4,4) * t189;
t222 = Icges(4,1) * t190 - t248;
t127 = -Icges(4,5) * t196 + t222 * t195;
t128 = Icges(4,5) * t195 + t222 * t196;
t147 = Icges(5,2) * t182 + t246;
t148 = Icges(5,1) * t181 + t245;
t152 = Icges(4,2) * t190 + t248;
t153 = Icges(4,1) * t189 + t247;
t154 = V_base(5) + (-qJD(2) - qJD(3)) * t196;
t180 = qJD(2) * t195 + V_base(4);
t155 = qJD(3) * t195 + t180;
t259 = (-t147 * t181 + t148 * t182 - t152 * t189 + t153 * t190) * V_base(6) + (-t112 * t181 + t114 * t182 - t126 * t189 + t128 * t190) * t155 + (-t111 * t181 + t113 * t182 - t125 * t189 + t127 * t190) * t154;
t258 = (Icges(4,5) * t189 + Icges(5,5) * t181 + Icges(4,6) * t190 + Icges(5,6) * t182) * V_base(6) + (t263 * t195 + t262 * t196) * t155 + (t262 * t195 - t263 * t196) * t154;
t198 = sin(qJ(2));
t255 = pkin(2) * t198;
t254 = pkin(3) * t189;
t200 = cos(qJ(2));
t253 = t200 * pkin(2);
t251 = Icges(2,4) * t195;
t250 = Icges(3,4) * t198;
t249 = Icges(3,4) * t200;
t244 = t181 * t195;
t243 = t181 * t196;
t197 = sin(qJ(5));
t242 = t195 * t197;
t199 = cos(qJ(5));
t241 = t195 * t199;
t240 = t196 * t197;
t239 = t196 * t199;
t117 = -pkin(6) * t196 + t253 * t195;
t173 = t195 * pkin(1) - t196 * pkin(5);
t238 = -t117 - t173;
t237 = pkin(3) * t190;
t235 = qJD(5) * t181;
t234 = V_base(5) * qJ(1) + V_base(1);
t230 = qJD(1) + V_base(3);
t100 = -qJ(4) * t196 + t237 * t195;
t229 = -t100 + t238;
t179 = -qJD(2) * t196 + V_base(5);
t228 = t179 * t255 + t234;
t227 = pkin(4) * t182 + pkin(7) * t181;
t226 = rSges(3,1) * t200 - rSges(3,2) * t198;
t225 = rSges(4,1) * t190 - rSges(4,2) * t189;
t224 = rSges(5,1) * t182 - rSges(5,2) * t181;
t223 = Icges(3,1) * t200 - t250;
t220 = -Icges(3,2) * t198 + t249;
t217 = Icges(3,5) * t200 - Icges(3,6) * t198;
t214 = qJD(4) * t195 + t154 * t254 + t228;
t174 = t196 * pkin(1) + t195 * pkin(5);
t213 = -V_base(4) * qJ(1) + V_base(6) * t174 + V_base(2);
t212 = V_base(4) * t173 - t174 * V_base(5) + t230;
t209 = (-Icges(3,3) * t196 + t217 * t195) * t179 + (Icges(3,3) * t195 + t217 * t196) * t180 + (Icges(3,5) * t198 + Icges(3,6) * t200) * V_base(6);
t118 = pkin(6) * t195 + t253 * t196;
t208 = V_base(6) * t118 - t180 * t255 + t213;
t207 = t180 * t117 - t118 * t179 + t212;
t206 = t155 * t100 + t207;
t101 = qJ(4) * t195 + t237 * t196;
t205 = -qJD(4) * t196 + V_base(6) * t101 + t208;
t135 = -Icges(3,6) * t196 + t220 * t195;
t136 = Icges(3,6) * t195 + t220 * t196;
t137 = -Icges(3,5) * t196 + t223 * t195;
t138 = Icges(3,5) * t195 + t223 * t196;
t176 = Icges(3,2) * t200 + t250;
t177 = Icges(3,1) * t198 + t249;
t202 = (-t136 * t198 + t138 * t200) * t180 + (-t135 * t198 + t137 * t200) * t179 + (-t176 * t198 + t177 * t200) * V_base(6);
t187 = Icges(2,4) * t196;
t178 = rSges(3,1) * t198 + rSges(3,2) * t200;
t172 = rSges(2,1) * t196 - rSges(2,2) * t195;
t171 = rSges(2,1) * t195 + rSges(2,2) * t196;
t170 = Icges(2,1) * t196 - t251;
t169 = Icges(2,1) * t195 + t187;
t168 = -Icges(2,2) * t195 + t187;
t167 = Icges(2,2) * t196 + t251;
t164 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t163 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t162 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t161 = -qJD(5) * t182 + V_base(6);
t156 = rSges(4,1) * t189 + rSges(4,2) * t190;
t150 = pkin(4) * t181 - pkin(7) * t182;
t149 = rSges(5,1) * t181 + rSges(5,2) * t182;
t145 = t182 * t239 + t242;
t144 = -t182 * t240 + t241;
t143 = t182 * t241 - t240;
t142 = -t182 * t242 - t239;
t140 = rSges(3,3) * t195 + t226 * t196;
t139 = -rSges(3,3) * t196 + t226 * t195;
t132 = t227 * t196;
t131 = t227 * t195;
t130 = rSges(4,3) * t195 + t225 * t196;
t129 = -rSges(4,3) * t196 + t225 * t195;
t122 = t196 * t235 + t155;
t121 = t195 * t235 + t154;
t120 = V_base(5) * rSges(2,3) - t171 * V_base(6) + t234;
t119 = t172 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t116 = rSges(5,3) * t195 + t224 * t196;
t115 = -rSges(5,3) * t196 + t224 * t195;
t107 = -rSges(6,3) * t182 + (rSges(6,1) * t199 - rSges(6,2) * t197) * t181;
t106 = -Icges(6,5) * t182 + (Icges(6,1) * t199 - Icges(6,4) * t197) * t181;
t105 = -Icges(6,6) * t182 + (Icges(6,4) * t199 - Icges(6,2) * t197) * t181;
t104 = -Icges(6,3) * t182 + (Icges(6,5) * t199 - Icges(6,6) * t197) * t181;
t103 = t171 * V_base(4) - t172 * V_base(5) + t230;
t97 = rSges(6,1) * t145 + rSges(6,2) * t144 + rSges(6,3) * t243;
t96 = rSges(6,1) * t143 + rSges(6,2) * t142 + rSges(6,3) * t244;
t95 = Icges(6,1) * t145 + Icges(6,4) * t144 + Icges(6,5) * t243;
t94 = Icges(6,1) * t143 + Icges(6,4) * t142 + Icges(6,5) * t244;
t93 = Icges(6,4) * t145 + Icges(6,2) * t144 + Icges(6,6) * t243;
t92 = Icges(6,4) * t143 + Icges(6,2) * t142 + Icges(6,6) * t244;
t91 = Icges(6,5) * t145 + Icges(6,6) * t144 + Icges(6,3) * t243;
t90 = Icges(6,5) * t143 + Icges(6,6) * t142 + Icges(6,3) * t244;
t89 = t178 * t179 + (-t139 - t173) * V_base(6) + t234;
t88 = t140 * V_base(6) - t178 * t180 + t213;
t87 = t139 * t180 - t140 * t179 + t212;
t86 = t154 * t156 + (-t129 + t238) * V_base(6) + t228;
t85 = t130 * V_base(6) - t155 * t156 + t208;
t84 = t129 * t155 - t130 * t154 + t207;
t83 = t149 * t154 + (-t115 + t229) * V_base(6) + t214;
t82 = t116 * V_base(6) + (-t149 - t254) * t155 + t205;
t81 = t115 * t155 + (-t101 - t116) * t154 + t206;
t80 = t107 * t121 + t150 * t154 - t161 * t96 + (-t131 + t229) * V_base(6) + t214;
t79 = -t107 * t122 + t132 * V_base(6) + t161 * t97 + (-t150 - t254) * t155 + t205;
t78 = -t121 * t97 + t122 * t96 + t131 * t155 + (-t101 - t132) * t154 + t206;
t1 = m(1) * (t162 ^ 2 + t163 ^ 2 + t164 ^ 2) / 0.2e1 + m(2) * (t103 ^ 2 + t119 ^ 2 + t120 ^ 2) / 0.2e1 + m(3) * (t87 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + t180 * (t209 * t195 + t202 * t196) / 0.2e1 + t179 * (t202 * t195 - t209 * t196) / 0.2e1 + m(4) * (t84 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + m(5) * (t81 ^ 2 + t82 ^ 2 + t83 ^ 2) / 0.2e1 + m(6) * (t78 ^ 2 + t79 ^ 2 + t80 ^ 2) / 0.2e1 + t122 * ((t144 * t93 + t145 * t95 + t91 * t243) * t122 + (t144 * t92 + t145 * t94 + t90 * t243) * t121 + (t104 * t243 + t105 * t144 + t106 * t145) * t161) / 0.2e1 + t121 * ((t142 * t93 + t143 * t95 + t91 * t244) * t122 + (t142 * t92 + t143 * t94 + t90 * t244) * t121 + (t104 * t244 + t105 * t142 + t106 * t143) * t161) / 0.2e1 + t161 * ((-t104 * t161 - t90 * t121 - t91 * t122) * t182 + ((-t197 * t93 + t199 * t95) * t122 + (-t197 * t92 + t199 * t94) * t121 + (-t105 * t197 + t106 * t199) * t161) * t181) / 0.2e1 + (t259 * t195 - t258 * t196) * t154 / 0.2e1 + (t258 * t195 + t259 * t196) * t155 / 0.2e1 + ((-t167 * t195 + t169 * t196 + Icges(1,4)) * V_base(5) + (-t168 * t195 + t170 * t196 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t167 * t196 + t169 * t195 + Icges(1,2)) * V_base(5) + (t168 * t196 + t170 * t195 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t136 * t200 + t138 * t198) * t180 + (t135 * t200 + t137 * t198) * t179 + (t112 * t182 + t114 * t181 + t126 * t190 + t128 * t189) * t155 + (t111 * t182 + t113 * t181 + t125 * t190 + t127 * t189) * t154 + (t147 * t182 + t148 * t181 + t152 * t190 + t153 * t189 + t176 * t200 + t177 * t198 + Icges(1,3) + Icges(2,3)) * V_base(6)) * V_base(6) / 0.2e1 + V_base(6) * V_base(4) * (Icges(2,5) * t196 - Icges(2,6) * t195 + Icges(1,5)) + V_base(6) * V_base(5) * (Icges(2,5) * t195 + Icges(2,6) * t196 + Icges(1,6));
T = t1;
