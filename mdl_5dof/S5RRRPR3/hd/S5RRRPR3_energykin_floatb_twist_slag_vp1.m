% Calculate kinetic energy for
% S5RRRPR3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2020-01-03 12:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRPR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR3_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR3_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR3_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:08:47
% EndTime: 2020-01-03 12:08:49
% DurationCPUTime: 1.80s
% Computational Cost: add. (1294->241), mult. (946->330), div. (0->0), fcn. (726->10), ass. (0->126)
t237 = Icges(4,3) + Icges(5,3);
t167 = qJ(3) + pkin(9);
t156 = sin(t167);
t157 = cos(t167);
t170 = sin(qJ(3));
t172 = cos(qJ(3));
t236 = Icges(4,5) * t172 + Icges(5,5) * t157 - Icges(4,6) * t170 - Icges(5,6) * t156;
t168 = qJ(1) + qJ(2);
t160 = sin(t168);
t161 = cos(t168);
t216 = Icges(4,4) * t172;
t193 = -Icges(4,2) * t170 + t216;
t102 = -Icges(4,6) * t161 + t193 * t160;
t103 = -Icges(4,6) * t160 - t193 * t161;
t217 = Icges(4,4) * t170;
t196 = Icges(4,1) * t172 - t217;
t104 = -Icges(4,5) * t161 + t196 * t160;
t105 = -Icges(4,5) * t160 - t196 * t161;
t215 = Icges(5,4) * t156;
t117 = Icges(5,2) * t157 + t215;
t214 = Icges(5,4) * t157;
t118 = Icges(5,1) * t156 + t214;
t135 = -qJD(3) * t160 + V_base(6);
t136 = -qJD(3) * t161 + V_base(5);
t140 = Icges(4,2) * t172 + t217;
t143 = Icges(4,1) * t170 + t216;
t159 = V_base(4) + qJD(1);
t154 = qJD(2) + t159;
t192 = -Icges(5,2) * t156 + t214;
t91 = -Icges(5,6) * t161 + t192 * t160;
t92 = -Icges(5,6) * t160 - t192 * t161;
t195 = Icges(5,1) * t157 - t215;
t93 = -Icges(5,5) * t161 + t195 * t160;
t94 = -Icges(5,5) * t160 - t195 * t161;
t235 = t135 * (t103 * t170 - t105 * t172 + t156 * t92 - t157 * t94) + t136 * (t102 * t170 - t104 * t172 + t156 * t91 - t157 * t93) + t154 * (t117 * t156 - t118 * t157 + t140 * t170 - t143 * t172);
t234 = (-Icges(4,5) * t170 - Icges(5,5) * t156 - Icges(4,6) * t172 - Icges(5,6) * t157) * t154 + (-t236 * t160 + t237 * t161) * t136 + (t237 * t160 + t236 * t161) * t135;
t158 = qJ(5) + t167;
t153 = cos(t158);
t152 = sin(t158);
t213 = Icges(6,4) * t152;
t110 = Icges(6,2) * t153 + t213;
t212 = Icges(6,4) * t153;
t111 = Icges(6,1) * t152 + t212;
t209 = -qJD(3) - qJD(5);
t112 = t209 * t160 + V_base(6);
t113 = t209 * t161 + V_base(5);
t191 = -Icges(6,2) * t152 + t212;
t81 = -Icges(6,6) * t161 + t191 * t160;
t82 = -Icges(6,6) * t160 - t191 * t161;
t194 = Icges(6,1) * t153 - t213;
t83 = -Icges(6,5) * t161 + t194 * t160;
t84 = -Icges(6,5) * t160 - t194 * t161;
t230 = (t110 * t152 - t111 * t153) * t154 + (t152 * t81 - t153 * t83) * t113 + (t152 * t82 - t153 * t84) * t112;
t229 = -pkin(5) - pkin(6);
t225 = pkin(1) * t159;
t224 = pkin(3) * t170;
t223 = pkin(4) * t156;
t222 = t172 * pkin(3);
t130 = -pkin(2) * t161 - pkin(7) * t160;
t86 = -qJ(4) * t160 - t222 * t161;
t220 = -t130 - t86;
t173 = cos(qJ(1));
t219 = Icges(2,4) * t173;
t218 = Icges(3,4) * t161;
t211 = pkin(4) * t157;
t171 = sin(qJ(1));
t208 = t171 * t225 + V_base(3);
t207 = V_base(6) * pkin(5) + V_base(2);
t204 = V_base(6) * pkin(6) + t173 * t225 + t207;
t203 = rSges(4,1) * t172 - rSges(4,2) * t170;
t202 = rSges(5,1) * t157 - rSges(5,2) * t156;
t201 = rSges(6,1) * t153 - rSges(6,2) * t152;
t188 = Icges(6,5) * t153 - Icges(6,6) * t152;
t182 = -(Icges(6,5) * t152 + Icges(6,6) * t153) * t154 - (-Icges(6,3) * t160 - t188 * t161) * t112 - (-Icges(6,3) * t161 + t188 * t160) * t113;
t129 = pkin(2) * t160 - pkin(7) * t161;
t179 = t154 * t129 + t229 * V_base(5) + t208;
t178 = -qJD(4) * t161 + t135 * t224 + t204;
t177 = V_base(1) + (-V_base(6) * t171 - t173 * V_base(5)) * pkin(1);
t85 = -qJ(4) * t161 + t222 * t160;
t176 = -qJD(4) * t160 + t154 * t85 + t179;
t175 = -V_base(6) * t129 + V_base(5) * t130 + t177;
t174 = t136 * t86 + t175;
t162 = Icges(2,4) * t171;
t151 = Icges(3,4) * t160;
t148 = -rSges(2,1) * t173 + t171 * rSges(2,2);
t147 = t171 * rSges(2,1) + rSges(2,2) * t173;
t146 = rSges(4,1) * t170 + rSges(4,2) * t172;
t145 = -Icges(2,1) * t173 + t162;
t144 = Icges(2,1) * t171 + t219;
t142 = Icges(2,2) * t171 - t219;
t141 = Icges(2,2) * t173 + t162;
t134 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t133 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t132 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t128 = -rSges(3,1) * t161 + rSges(3,2) * t160;
t127 = rSges(3,1) * t160 + rSges(3,2) * t161;
t126 = -Icges(3,1) * t161 + t151;
t125 = Icges(3,1) * t160 + t218;
t124 = Icges(3,2) * t160 - t218;
t123 = Icges(3,2) * t161 + t151;
t120 = rSges(5,1) * t156 + rSges(5,2) * t157;
t114 = rSges(6,1) * t152 + rSges(6,2) * t153;
t107 = -rSges(4,3) * t160 - t203 * t161;
t106 = -rSges(4,3) * t161 + t203 * t160;
t99 = V_base(6) * rSges(2,3) - t148 * t159 + t207;
t98 = t147 * t159 + V_base(3) + (-rSges(2,3) - pkin(5)) * V_base(5);
t97 = -t147 * V_base(6) + t148 * V_base(5) + V_base(1);
t96 = -rSges(5,3) * t160 - t202 * t161;
t95 = -rSges(5,3) * t161 + t202 * t160;
t88 = -rSges(6,3) * t160 - t201 * t161;
t87 = -rSges(6,3) * t161 + t201 * t160;
t76 = V_base(6) * rSges(3,3) - t128 * t154 + t204;
t75 = t127 * t154 + (-rSges(3,3) + t229) * V_base(5) + t208;
t74 = -pkin(8) * t160 - t211 * t161;
t73 = -pkin(8) * t161 + t211 * t160;
t72 = -V_base(6) * t127 + V_base(5) * t128 + t177;
t71 = t135 * t146 + (-t107 - t130) * t154 + t204;
t70 = t106 * t154 - t136 * t146 + t179;
t69 = -t135 * t106 + t136 * t107 + t175;
t68 = t120 * t135 + (-t96 + t220) * t154 + t178;
t67 = t154 * t95 + (-t120 - t224) * t136 + t176;
t66 = t136 * t96 + (-t85 - t95) * t135 + t174;
t65 = t135 * t223 + t112 * t114 + (-t74 - t88 + t220) * t154 + t178;
t64 = -t113 * t114 + (t73 + t87) * t154 + (-t223 - t224) * t136 + t176;
t63 = -t112 * t87 + t113 * t88 + t136 * t74 + (-t73 - t85) * t135 + t174;
t1 = m(1) * (t132 ^ 2 + t133 ^ 2 + t134 ^ 2) / 0.2e1 + m(2) * (t97 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(3) * (t72 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + m(4) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + m(5) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + m(6) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + t113 * (-t230 * t160 + t182 * t161) / 0.2e1 + t112 * (t182 * t160 + t230 * t161) / 0.2e1 + (t234 * t160 + t235 * t161) * t135 / 0.2e1 + (-t235 * t160 + t234 * t161) * t136 / 0.2e1 + ((t124 * t161 + t126 * t160 + t142 * t173 + t171 * t145 + Icges(1,6)) * V_base(6) + (t161 * t123 + t160 * t125 + t173 * t141 + t171 * t144 + Icges(1,2)) * V_base(5)) * V_base(5) / 0.2e1 + ((t160 * t124 - t161 * t126 + t171 * t142 - t173 * t145 + Icges(1,3)) * V_base(6) + (t123 * t160 - t125 * t161 + t171 * t141 - t144 * t173 + Icges(1,6)) * V_base(5)) * V_base(6) / 0.2e1 + ((t152 * t83 + t153 * t81) * t113 + (t152 * t84 + t153 * t82) * t112 + (t102 * t172 + t104 * t170 + t156 * t93 + t157 * t91) * t136 + (t103 * t172 + t105 * t170 + t156 * t94 + t157 * t92) * t135 + (t153 * t110 + t152 * t111 + t157 * t117 + t156 * t118 + t172 * t140 + t170 * t143 + Icges(3,3)) * t154) * t154 / 0.2e1 + t154 * V_base(6) * (-Icges(3,5) * t161 + Icges(3,6) * t160) + t154 * V_base(5) * (Icges(3,5) * t160 + Icges(3,6) * t161) + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + Icges(1,1) * V_base(4) / 0.2e1) * V_base(4) + ((Icges(2,5) * t171 + Icges(2,6) * t173) * V_base(5) + (-Icges(2,5) * t173 + Icges(2,6) * t171) * V_base(6) + Icges(2,3) * t159 / 0.2e1) * t159;
T = t1;
