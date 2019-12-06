% Calculate kinetic energy for
% S5RPRRR2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-05 18:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR2_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR2_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRRR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR2_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR2_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR2_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR2_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:11:13
% EndTime: 2019-12-05 18:11:15
% DurationCPUTime: 1.78s
% Computational Cost: add. (1295->265), mult. (1174->381), div. (0->0), fcn. (956->10), ass. (0->144)
t187 = sin(pkin(9));
t248 = pkin(2) * t187;
t186 = pkin(9) + qJ(3);
t174 = sin(t186);
t247 = pkin(3) * t174;
t176 = qJ(4) + t186;
t170 = sin(t176);
t246 = pkin(4) * t170;
t188 = cos(pkin(9));
t245 = t188 * pkin(2);
t190 = sin(qJ(1));
t244 = Icges(2,4) * t190;
t243 = Icges(3,4) * t187;
t242 = Icges(3,4) * t188;
t241 = Icges(4,4) * t174;
t175 = cos(t186);
t240 = Icges(4,4) * t175;
t239 = Icges(5,4) * t170;
t171 = cos(t176);
t238 = Icges(5,4) * t171;
t173 = qJ(5) + t176;
t166 = sin(t173);
t237 = Icges(6,4) * t166;
t167 = cos(t173);
t236 = Icges(6,4) * t167;
t191 = cos(qJ(1));
t108 = -pkin(6) * t191 + t190 * t245;
t160 = pkin(1) * t190 - qJ(2) * t191;
t234 = -t108 - t160;
t233 = pkin(4) * t171;
t232 = pkin(3) * t175;
t229 = -qJD(3) - qJD(4);
t228 = V_base(4) * t160 + V_base(3);
t227 = V_base(5) * pkin(5) + V_base(1);
t85 = -pkin(7) * t191 + t190 * t232;
t224 = -t85 + t234;
t165 = qJD(3) * t190 + V_base(4);
t223 = qJD(2) * t190 + t227;
t144 = qJD(4) * t190 + t165;
t222 = V_base(5) * t248 + t223;
t221 = rSges(3,1) * t188 - rSges(3,2) * t187;
t220 = rSges(4,1) * t175 - rSges(4,2) * t174;
t219 = rSges(5,1) * t171 - rSges(5,2) * t170;
t218 = rSges(6,1) * t167 - rSges(6,2) * t166;
t217 = Icges(3,1) * t188 - t243;
t216 = Icges(4,1) * t175 - t241;
t215 = Icges(5,1) * t171 - t239;
t214 = Icges(6,1) * t167 - t237;
t213 = -Icges(3,2) * t187 + t242;
t212 = -Icges(4,2) * t174 + t240;
t211 = -Icges(5,2) * t170 + t238;
t210 = -Icges(6,2) * t166 + t236;
t209 = Icges(3,5) * t188 - Icges(3,6) * t187;
t208 = Icges(4,5) * t175 - Icges(4,6) * t174;
t207 = Icges(5,5) * t171 - Icges(5,6) * t170;
t206 = Icges(6,5) * t167 - Icges(6,6) * t166;
t162 = pkin(1) * t191 + t190 * qJ(2);
t177 = V_base(6) + qJD(1);
t205 = -qJD(2) * t191 + t162 * t177 + V_base(2);
t164 = -qJD(3) * t191 + V_base(5);
t204 = t164 * t247 + t222;
t131 = V_base(5) + (-qJD(5) + t229) * t191;
t132 = qJD(5) * t190 + t144;
t203 = (Icges(6,5) * t166 + Icges(6,6) * t167) * t177 + (-Icges(6,3) * t191 + t190 * t206) * t131 + (Icges(6,3) * t190 + t191 * t206) * t132;
t143 = t191 * t229 + V_base(5);
t202 = (Icges(5,5) * t170 + Icges(5,6) * t171) * t177 + (-Icges(5,3) * t191 + t190 * t207) * t143 + (Icges(5,3) * t190 + t191 * t207) * t144;
t201 = (-Icges(4,3) * t191 + t190 * t208) * t164 + (Icges(4,3) * t190 + t191 * t208) * t165 + (Icges(4,5) * t174 + Icges(4,6) * t175) * t177;
t109 = pkin(6) * t190 + t191 * t245;
t200 = V_base(4) * t108 + (-t109 - t162) * V_base(5) + t228;
t199 = (-Icges(3,3) * t191 + t190 * t209) * V_base(5) + (Icges(3,3) * t190 + t191 * t209) * V_base(4) + (Icges(3,5) * t187 + Icges(3,6) * t188) * t177;
t86 = pkin(7) * t190 + t191 * t232;
t198 = -t164 * t86 + t165 * t85 + t200;
t197 = t177 * t109 + (-pkin(5) - t248) * V_base(4) + t205;
t196 = -t165 * t247 + t177 * t86 + t197;
t128 = Icges(6,2) * t167 + t237;
t129 = Icges(6,1) * t166 + t236;
t90 = -Icges(6,6) * t191 + t190 * t210;
t91 = Icges(6,6) * t190 + t191 * t210;
t92 = -Icges(6,5) * t191 + t190 * t214;
t93 = Icges(6,5) * t190 + t191 * t214;
t195 = (-t166 * t91 + t167 * t93) * t132 + (-t166 * t90 + t167 * t92) * t131 + (-t128 * t166 + t129 * t167) * t177;
t100 = -Icges(5,6) * t191 + t190 * t211;
t101 = Icges(5,6) * t190 + t191 * t211;
t102 = -Icges(5,5) * t191 + t190 * t215;
t103 = Icges(5,5) * t190 + t191 * t215;
t134 = Icges(5,2) * t171 + t239;
t135 = Icges(5,1) * t170 + t238;
t194 = (-t101 * t170 + t103 * t171) * t144 + (-t100 * t170 + t102 * t171) * t143 + (-t134 * t170 + t135 * t171) * t177;
t112 = -Icges(4,6) * t191 + t190 * t212;
t113 = Icges(4,6) * t190 + t191 * t212;
t114 = -Icges(4,5) * t191 + t190 * t216;
t115 = Icges(4,5) * t190 + t191 * t216;
t140 = Icges(4,2) * t175 + t241;
t141 = Icges(4,1) * t174 + t240;
t193 = (-t113 * t174 + t115 * t175) * t165 + (-t112 * t174 + t114 * t175) * t164 + (-t140 * t174 + t141 * t175) * t177;
t121 = -Icges(3,6) * t191 + t190 * t213;
t122 = Icges(3,6) * t190 + t191 * t213;
t123 = -Icges(3,5) * t191 + t190 * t217;
t124 = Icges(3,5) * t190 + t191 * t217;
t151 = Icges(3,2) * t188 + t243;
t152 = Icges(3,1) * t187 + t242;
t192 = (-t122 * t187 + t124 * t188) * V_base(4) + (-t121 * t187 + t123 * t188) * V_base(5) + (-t151 * t187 + t152 * t188) * t177;
t183 = Icges(2,4) * t191;
t163 = rSges(2,1) * t191 - rSges(2,2) * t190;
t161 = rSges(2,1) * t190 + rSges(2,2) * t191;
t159 = Icges(2,1) * t191 - t244;
t158 = Icges(2,1) * t190 + t183;
t157 = -Icges(2,2) * t190 + t183;
t156 = Icges(2,2) * t191 + t244;
t155 = Icges(2,5) * t191 - Icges(2,6) * t190;
t154 = Icges(2,5) * t190 + Icges(2,6) * t191;
t153 = rSges(3,1) * t187 + rSges(3,2) * t188;
t149 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t148 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t147 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t142 = rSges(4,1) * t174 + rSges(4,2) * t175;
t136 = rSges(5,1) * t170 + rSges(5,2) * t171;
t130 = rSges(6,1) * t166 + rSges(6,2) * t167;
t126 = t190 * rSges(3,3) + t191 * t221;
t125 = -rSges(3,3) * t191 + t190 * t221;
t117 = t190 * rSges(4,3) + t191 * t220;
t116 = -rSges(4,3) * t191 + t190 * t220;
t107 = t190 * rSges(5,3) + t191 * t219;
t106 = -rSges(5,3) * t191 + t190 * t219;
t105 = V_base(5) * rSges(2,3) - t161 * t177 + t227;
t104 = t163 * t177 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t97 = t161 * V_base(4) - t163 * V_base(5) + V_base(3);
t95 = t190 * rSges(6,3) + t191 * t218;
t94 = -rSges(6,3) * t191 + t190 * t218;
t82 = pkin(8) * t190 + t191 * t233;
t81 = -pkin(8) * t191 + t190 * t233;
t80 = t153 * V_base(5) + (-t125 - t160) * t177 + t223;
t79 = t177 * t126 + (-pkin(5) - t153) * V_base(4) + t205;
t78 = t125 * V_base(4) + (-t126 - t162) * V_base(5) + t228;
t77 = t142 * t164 + (-t116 + t234) * t177 + t222;
t76 = t117 * t177 - t142 * t165 + t197;
t75 = t116 * t165 - t117 * t164 + t200;
t74 = t136 * t143 + (-t106 + t224) * t177 + t204;
t73 = t107 * t177 - t136 * t144 + t196;
t72 = t106 * t144 - t107 * t143 + t198;
t71 = t143 * t246 + t130 * t131 + (-t81 - t94 + t224) * t177 + t204;
t70 = -t144 * t246 - t132 * t130 + (t82 + t95) * t177 + t196;
t69 = -t131 * t95 + t132 * t94 - t143 * t82 + t144 * t81 + t198;
t1 = m(1) * (t147 ^ 2 + t148 ^ 2 + t149 ^ 2) / 0.2e1 + m(2) * (t104 ^ 2 + t105 ^ 2 + t97 ^ 2) / 0.2e1 + m(3) * (t78 ^ 2 + t79 ^ 2 + t80 ^ 2) / 0.2e1 + m(4) * (t75 ^ 2 + t76 ^ 2 + t77 ^ 2) / 0.2e1 + t165 * (t201 * t190 + t193 * t191) / 0.2e1 + t164 * (t193 * t190 - t201 * t191) / 0.2e1 + m(5) * (t72 ^ 2 + t73 ^ 2 + t74 ^ 2) / 0.2e1 + t144 * (t202 * t190 + t194 * t191) / 0.2e1 + t143 * (t194 * t190 - t202 * t191) / 0.2e1 + m(6) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + t132 * (t203 * t190 + t195 * t191) / 0.2e1 + t131 * (t195 * t190 - t203 * t191) / 0.2e1 + (t155 * t177 + t199 * t190 + t192 * t191 + (-t156 * t190 + t158 * t191 + Icges(1,4)) * V_base(5) + (-t190 * t157 + t159 * t191 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t154 * t177 + t192 * t190 - t199 * t191 + (t156 * t191 + t190 * t158 + Icges(1,2)) * V_base(5) + (t157 * t191 + t159 * t190 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t113 * t175 + t115 * t174) * t165 + (t112 * t175 + t114 * t174) * t164 + (t101 * t171 + t103 * t170) * t144 + (t100 * t171 + t102 * t170) * t143 + (t166 * t93 + t167 * t91) * t132 + (t166 * t92 + t167 * t90) * t131 + (t121 * t188 + t123 * t187 + t154) * V_base(5) + (t122 * t188 + t124 * t187 + t155) * V_base(4) + (t128 * t167 + t129 * t166 + t134 * t171 + t135 * t170 + t140 * t175 + t141 * t174 + t151 * t188 + t152 * t187 + Icges(2,3)) * t177) * t177 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
