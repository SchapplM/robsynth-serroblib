% Calculate kinetic energy for
% S5PRPRR2
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR2_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRPRR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR2_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR2_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR2_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:44:43
% EndTime: 2019-12-05 15:44:45
% DurationCPUTime: 2.14s
% Computational Cost: add. (1365->280), mult. (1390->402), div. (0->0), fcn. (1228->10), ass. (0->144)
t263 = Icges(3,3) + Icges(4,3);
t194 = qJ(2) + pkin(9);
t187 = sin(t194);
t188 = cos(t194);
t199 = sin(qJ(2));
t201 = cos(qJ(2));
t262 = Icges(3,5) * t201 + Icges(4,5) * t188 - Icges(3,6) * t199 - Icges(4,6) * t187;
t195 = sin(pkin(8));
t196 = cos(pkin(8));
t247 = Icges(4,4) * t188;
t219 = -Icges(4,2) * t187 + t247;
t121 = -Icges(4,6) * t196 + t219 * t195;
t122 = Icges(4,6) * t195 + t219 * t196;
t248 = Icges(4,4) * t187;
t222 = Icges(4,1) * t188 - t248;
t123 = -Icges(4,5) * t196 + t222 * t195;
t124 = Icges(4,5) * t195 + t222 * t196;
t249 = Icges(3,4) * t201;
t220 = -Icges(3,2) * t199 + t249;
t135 = -Icges(3,6) * t196 + t220 * t195;
t136 = Icges(3,6) * t195 + t220 * t196;
t250 = Icges(3,4) * t199;
t223 = Icges(3,1) * t201 - t250;
t137 = -Icges(3,5) * t196 + t223 * t195;
t138 = Icges(3,5) * t195 + t223 * t196;
t152 = Icges(4,2) * t188 + t248;
t153 = Icges(4,1) * t187 + t247;
t176 = Icges(3,2) * t201 + t250;
t177 = Icges(3,1) * t199 + t249;
t179 = -qJD(2) * t196 + V_base(5);
t180 = qJD(2) * t195 + V_base(4);
t259 = (-t152 * t187 + t153 * t188 - t176 * t199 + t177 * t201) * V_base(6) + (-t122 * t187 + t124 * t188 - t136 * t199 + t138 * t201) * t180 + (-t121 * t187 + t123 * t188 - t135 * t199 + t137 * t201) * t179;
t258 = (Icges(3,5) * t199 + Icges(4,5) * t187 + Icges(3,6) * t201 + Icges(4,6) * t188) * V_base(6) + (t263 * t195 + t262 * t196) * t180 + (t262 * t195 - t263 * t196) * t179;
t255 = pkin(2) * t199;
t254 = pkin(3) * t187;
t253 = t201 * pkin(2);
t251 = Icges(2,4) * t195;
t190 = qJ(4) + t194;
t181 = sin(t190);
t246 = Icges(5,4) * t181;
t182 = cos(t190);
t245 = Icges(5,4) * t182;
t244 = t181 * t195;
t243 = t181 * t196;
t198 = sin(qJ(5));
t242 = t195 * t198;
t200 = cos(qJ(5));
t241 = t195 * t200;
t240 = t196 * t198;
t239 = t196 * t200;
t117 = -qJ(3) * t196 + t253 * t195;
t173 = pkin(1) * t195 - pkin(5) * t196;
t238 = -t117 - t173;
t237 = pkin(3) * t188;
t235 = qJD(5) * t181;
t234 = V_base(5) * qJ(1) + V_base(1);
t230 = qJD(1) + V_base(3);
t100 = -pkin(6) * t196 + t237 * t195;
t229 = -t100 + t238;
t156 = qJD(4) * t195 + t180;
t228 = qJD(3) * t195 + t179 * t255 + t234;
t227 = pkin(4) * t182 + pkin(7) * t181;
t226 = rSges(3,1) * t201 - rSges(3,2) * t199;
t225 = rSges(4,1) * t188 - rSges(4,2) * t187;
t224 = rSges(5,1) * t182 - rSges(5,2) * t181;
t221 = Icges(5,1) * t182 - t246;
t218 = -Icges(5,2) * t181 + t245;
t215 = Icges(5,5) * t182 - Icges(5,6) * t181;
t214 = t179 * t254 + t228;
t174 = pkin(1) * t196 + pkin(5) * t195;
t213 = -V_base(4) * qJ(1) + V_base(6) * t174 + V_base(2);
t155 = V_base(5) + (-qJD(2) - qJD(4)) * t196;
t212 = V_base(4) * t173 - t174 * V_base(5) + t230;
t211 = t180 * t117 + t212;
t210 = (-Icges(5,3) * t196 + t215 * t195) * t155 + (Icges(5,3) * t195 + t215 * t196) * t156 + (Icges(5,5) * t181 + Icges(5,6) * t182) * V_base(6);
t118 = qJ(3) * t195 + t253 * t196;
t207 = -qJD(3) * t196 + V_base(6) * t118 + t213;
t101 = pkin(6) * t195 + t237 * t196;
t206 = t180 * t100 + (-t101 - t118) * t179 + t211;
t205 = V_base(6) * t101 + (-t254 - t255) * t180 + t207;
t111 = -Icges(5,6) * t196 + t218 * t195;
t112 = Icges(5,6) * t195 + t218 * t196;
t113 = -Icges(5,5) * t196 + t221 * t195;
t114 = Icges(5,5) * t195 + t221 * t196;
t146 = Icges(5,2) * t182 + t246;
t147 = Icges(5,1) * t181 + t245;
t204 = (-t112 * t181 + t114 * t182) * t156 + (-t111 * t181 + t113 * t182) * t155 + (-t146 * t181 + t147 * t182) * V_base(6);
t189 = Icges(2,4) * t196;
t178 = t199 * rSges(3,1) + rSges(3,2) * t201;
t172 = rSges(2,1) * t196 - rSges(2,2) * t195;
t171 = rSges(2,1) * t195 + rSges(2,2) * t196;
t170 = Icges(2,1) * t196 - t251;
t169 = Icges(2,1) * t195 + t189;
t168 = -Icges(2,2) * t195 + t189;
t167 = Icges(2,2) * t196 + t251;
t164 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t163 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t162 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t161 = -qJD(5) * t182 + V_base(6);
t154 = rSges(4,1) * t187 + rSges(4,2) * t188;
t149 = pkin(4) * t181 - pkin(7) * t182;
t148 = rSges(5,1) * t181 + rSges(5,2) * t182;
t144 = t182 * t239 + t242;
t143 = -t182 * t240 + t241;
t142 = t182 * t241 - t240;
t141 = -t182 * t242 - t239;
t140 = t195 * rSges(3,3) + t226 * t196;
t139 = -t196 * rSges(3,3) + t226 * t195;
t132 = t227 * t196;
t131 = t227 * t195;
t130 = rSges(4,3) * t195 + t225 * t196;
t129 = -rSges(4,3) * t196 + t225 * t195;
t128 = t196 * t235 + t156;
t127 = t195 * t235 + t155;
t126 = V_base(5) * rSges(2,3) - t171 * V_base(6) + t234;
t125 = t172 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t116 = rSges(5,3) * t195 + t224 * t196;
t115 = -rSges(5,3) * t196 + t224 * t195;
t107 = -rSges(6,3) * t182 + (rSges(6,1) * t200 - rSges(6,2) * t198) * t181;
t106 = -Icges(6,5) * t182 + (Icges(6,1) * t200 - Icges(6,4) * t198) * t181;
t105 = -Icges(6,6) * t182 + (Icges(6,4) * t200 - Icges(6,2) * t198) * t181;
t104 = -Icges(6,3) * t182 + (Icges(6,5) * t200 - Icges(6,6) * t198) * t181;
t103 = t171 * V_base(4) - t172 * V_base(5) + t230;
t97 = rSges(6,1) * t144 + rSges(6,2) * t143 + rSges(6,3) * t243;
t96 = rSges(6,1) * t142 + rSges(6,2) * t141 + rSges(6,3) * t244;
t95 = Icges(6,1) * t144 + Icges(6,4) * t143 + Icges(6,5) * t243;
t94 = Icges(6,1) * t142 + Icges(6,4) * t141 + Icges(6,5) * t244;
t93 = Icges(6,4) * t144 + Icges(6,2) * t143 + Icges(6,6) * t243;
t92 = Icges(6,4) * t142 + Icges(6,2) * t141 + Icges(6,6) * t244;
t91 = Icges(6,5) * t144 + Icges(6,6) * t143 + Icges(6,3) * t243;
t90 = Icges(6,5) * t142 + Icges(6,6) * t141 + Icges(6,3) * t244;
t89 = t178 * t179 + (-t139 - t173) * V_base(6) + t234;
t88 = t140 * V_base(6) - t178 * t180 + t213;
t87 = t139 * t180 - t140 * t179 + t212;
t86 = t154 * t179 + (-t129 + t238) * V_base(6) + t228;
t85 = t130 * V_base(6) + (-t154 - t255) * t180 + t207;
t84 = t129 * t180 + (-t118 - t130) * t179 + t211;
t83 = t148 * t155 + (-t115 + t229) * V_base(6) + t214;
t82 = t116 * V_base(6) - t148 * t156 + t205;
t81 = t115 * t156 - t116 * t155 + t206;
t80 = t107 * t127 + t149 * t155 - t161 * t96 + (-t131 + t229) * V_base(6) + t214;
t79 = -t107 * t128 + t132 * V_base(6) - t149 * t156 + t161 * t97 + t205;
t78 = -t127 * t97 + t128 * t96 + t131 * t156 - t132 * t155 + t206;
t1 = m(1) * (t162 ^ 2 + t163 ^ 2 + t164 ^ 2) / 0.2e1 + m(2) * (t103 ^ 2 + t125 ^ 2 + t126 ^ 2) / 0.2e1 + m(3) * (t87 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + m(4) * (t84 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + m(5) * (t81 ^ 2 + t82 ^ 2 + t83 ^ 2) / 0.2e1 + t156 * (t210 * t195 + t204 * t196) / 0.2e1 + t155 * (t204 * t195 - t210 * t196) / 0.2e1 + m(6) * (t78 ^ 2 + t79 ^ 2 + t80 ^ 2) / 0.2e1 + t128 * ((t143 * t93 + t144 * t95 + t91 * t243) * t128 + (t143 * t92 + t144 * t94 + t90 * t243) * t127 + (t104 * t243 + t105 * t143 + t106 * t144) * t161) / 0.2e1 + t127 * ((t141 * t93 + t142 * t95 + t91 * t244) * t128 + (t141 * t92 + t142 * t94 + t90 * t244) * t127 + (t104 * t244 + t105 * t141 + t106 * t142) * t161) / 0.2e1 + t161 * ((-t104 * t161 - t90 * t127 - t91 * t128) * t182 + ((-t198 * t93 + t200 * t95) * t128 + (-t198 * t92 + t200 * t94) * t127 + (-t105 * t198 + t106 * t200) * t161) * t181) / 0.2e1 + (t259 * t195 - t258 * t196) * t179 / 0.2e1 + (t258 * t195 + t259 * t196) * t180 / 0.2e1 + ((-t167 * t195 + t169 * t196 + Icges(1,4)) * V_base(5) + (-t168 * t195 + t170 * t196 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t167 * t196 + t169 * t195 + Icges(1,2)) * V_base(5) + (t168 * t196 + t170 * t195 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t112 * t182 + t114 * t181) * t156 + (t111 * t182 + t113 * t181) * t155 + (t122 * t188 + t124 * t187 + t136 * t201 + t199 * t138) * t180 + (t121 * t188 + t123 * t187 + t135 * t201 + t199 * t137) * t179 + (t146 * t182 + t147 * t181 + t152 * t188 + t153 * t187 + t176 * t201 + t199 * t177 + Icges(1,3) + Icges(2,3)) * V_base(6)) * V_base(6) / 0.2e1 + V_base(6) * V_base(4) * (Icges(2,5) * t196 - Icges(2,6) * t195 + Icges(1,5)) + V_base(6) * V_base(5) * (Icges(2,5) * t195 + Icges(2,6) * t196 + Icges(1,6));
T = t1;
