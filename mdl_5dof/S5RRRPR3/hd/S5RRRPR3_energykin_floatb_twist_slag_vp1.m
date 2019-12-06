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
% Datum: 2019-12-05 18:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:42:15
% EndTime: 2019-12-05 18:42:17
% DurationCPUTime: 1.74s
% Computational Cost: add. (1294->240), mult. (946->330), div. (0->0), fcn. (726->10), ass. (0->126)
t241 = Icges(4,3) + Icges(5,3);
t169 = qJ(3) + pkin(9);
t159 = sin(t169);
t160 = cos(t169);
t172 = sin(qJ(3));
t174 = cos(qJ(3));
t240 = Icges(4,5) * t174 + Icges(5,5) * t160 - Icges(4,6) * t172 - Icges(5,6) * t159;
t170 = qJ(1) + qJ(2);
t163 = sin(t170);
t164 = cos(t170);
t217 = Icges(4,4) * t174;
t195 = -Icges(4,2) * t172 + t217;
t102 = Icges(4,6) * t164 - t195 * t163;
t103 = Icges(4,6) * t163 + t195 * t164;
t218 = Icges(4,4) * t172;
t198 = Icges(4,1) * t174 - t218;
t104 = Icges(4,5) * t164 - t198 * t163;
t105 = Icges(4,5) * t163 + t198 * t164;
t216 = Icges(5,4) * t159;
t117 = Icges(5,2) * t160 + t216;
t215 = Icges(5,4) * t160;
t118 = Icges(5,1) * t159 + t215;
t135 = qJD(3) * t163 + V_base(6);
t136 = qJD(3) * t164 + V_base(5);
t140 = Icges(4,2) * t174 + t218;
t143 = Icges(4,1) * t172 + t217;
t162 = V_base(4) + qJD(1);
t157 = qJD(2) + t162;
t194 = -Icges(5,2) * t159 + t215;
t91 = Icges(5,6) * t164 - t194 * t163;
t92 = Icges(5,6) * t163 + t194 * t164;
t197 = Icges(5,1) * t160 - t216;
t93 = Icges(5,5) * t164 - t197 * t163;
t94 = Icges(5,5) * t163 + t197 * t164;
t239 = (t103 * t172 - t105 * t174 + t159 * t92 - t160 * t94) * t135 + (t102 * t172 - t104 * t174 + t159 * t91 - t160 * t93) * t136 + (t117 * t159 - t118 * t160 + t140 * t172 - t143 * t174) * t157;
t238 = (Icges(4,5) * t172 + Icges(5,5) * t159 + Icges(4,6) * t174 + Icges(5,6) * t160) * t157 + (-t240 * t163 + t241 * t164) * t136 + (t241 * t163 + t240 * t164) * t135;
t161 = qJ(5) + t169;
t156 = cos(t161);
t155 = sin(t161);
t214 = Icges(6,4) * t155;
t110 = Icges(6,2) * t156 + t214;
t213 = Icges(6,4) * t156;
t111 = Icges(6,1) * t155 + t213;
t112 = qJD(5) * t163 + t135;
t113 = qJD(5) * t164 + t136;
t193 = -Icges(6,2) * t155 + t213;
t81 = Icges(6,6) * t164 - t193 * t163;
t82 = Icges(6,6) * t163 + t193 * t164;
t196 = Icges(6,1) * t156 - t214;
t83 = Icges(6,5) * t164 - t196 * t163;
t84 = Icges(6,5) * t163 + t196 * t164;
t234 = (t110 * t155 - t111 * t156) * t157 + (t155 * t81 - t156 * t83) * t113 + (t155 * t82 - t156 * t84) * t112;
t233 = -pkin(5) - pkin(6);
t173 = sin(qJ(1));
t229 = pkin(1) * t173;
t175 = cos(qJ(1));
t228 = pkin(1) * t175;
t227 = pkin(3) * t172;
t226 = pkin(4) * t159;
t225 = t174 * pkin(3);
t130 = pkin(2) * t164 + pkin(7) * t163;
t86 = qJ(4) * t163 + t225 * t164;
t223 = -t130 - t86;
t222 = Icges(2,4) * t173;
t221 = Icges(2,4) * t175;
t220 = Icges(3,4) * t163;
t219 = Icges(3,4) * t164;
t212 = pkin(4) * t160;
t210 = V_base(6) * pkin(5) + V_base(2);
t207 = V_base(5) * t228 + V_base(6) * t229 + V_base(1);
t206 = rSges(4,1) * t174 - rSges(4,2) * t172;
t205 = rSges(5,1) * t160 - rSges(5,2) * t159;
t204 = rSges(6,1) * t156 - rSges(6,2) * t155;
t203 = -t162 * t229 + V_base(3);
t190 = Icges(6,5) * t156 - Icges(6,6) * t155;
t184 = V_base(6) * pkin(6) - t162 * t228 + t210;
t183 = (Icges(6,5) * t155 + Icges(6,6) * t156) * t157 + (Icges(6,3) * t163 + t190 * t164) * t112 + (Icges(6,3) * t164 - t190 * t163) * t113;
t129 = -pkin(2) * t163 + pkin(7) * t164;
t180 = -t129 * V_base(6) + V_base(5) * t130 + t207;
t179 = qJD(4) * t164 + t135 * t227 + t184;
t178 = t136 * t86 + t180;
t177 = t157 * t129 + t233 * V_base(5) + t203;
t85 = qJ(4) * t164 - t225 * t163;
t176 = qJD(4) * t163 + t157 * t85 + t177;
t148 = rSges(2,1) * t175 - t173 * rSges(2,2);
t147 = -t173 * rSges(2,1) - rSges(2,2) * t175;
t146 = rSges(4,1) * t172 + rSges(4,2) * t174;
t145 = Icges(2,1) * t175 - t222;
t144 = -Icges(2,1) * t173 - t221;
t142 = -Icges(2,2) * t173 + t221;
t141 = -Icges(2,2) * t175 - t222;
t134 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t133 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t132 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t128 = rSges(3,1) * t164 - rSges(3,2) * t163;
t127 = -rSges(3,1) * t163 - rSges(3,2) * t164;
t126 = Icges(3,1) * t164 - t220;
t125 = -Icges(3,1) * t163 - t219;
t124 = -Icges(3,2) * t163 + t219;
t123 = -Icges(3,2) * t164 - t220;
t120 = rSges(5,1) * t159 + rSges(5,2) * t160;
t114 = rSges(6,1) * t155 + rSges(6,2) * t156;
t107 = rSges(4,3) * t163 + t206 * t164;
t106 = rSges(4,3) * t164 - t206 * t163;
t99 = V_base(6) * rSges(2,3) - t148 * t162 + t210;
t98 = t147 * t162 + V_base(3) + (-rSges(2,3) - pkin(5)) * V_base(5);
t97 = -t147 * V_base(6) + t148 * V_base(5) + V_base(1);
t96 = rSges(5,3) * t163 + t205 * t164;
t95 = rSges(5,3) * t164 - t205 * t163;
t88 = rSges(6,3) * t163 + t204 * t164;
t87 = rSges(6,3) * t164 - t204 * t163;
t76 = V_base(6) * rSges(3,3) - t157 * t128 + t184;
t75 = t127 * t157 + (-rSges(3,3) + t233) * V_base(5) + t203;
t74 = pkin(8) * t163 + t212 * t164;
t73 = pkin(8) * t164 - t212 * t163;
t72 = -t127 * V_base(6) + t128 * V_base(5) + t207;
t71 = t135 * t146 + (-t107 - t130) * t157 + t184;
t70 = t106 * t157 - t136 * t146 + t177;
t69 = -t106 * t135 + t107 * t136 + t180;
t68 = t135 * t120 + (-t96 + t223) * t157 + t179;
t67 = t157 * t95 + (-t120 - t227) * t136 + t176;
t66 = t136 * t96 + (-t85 - t95) * t135 + t178;
t65 = t135 * t226 + t112 * t114 + (-t74 - t88 + t223) * t157 + t179;
t64 = -t113 * t114 + (t73 + t87) * t157 + (-t226 - t227) * t136 + t176;
t63 = -t112 * t87 + t113 * t88 + t136 * t74 + (-t73 - t85) * t135 + t178;
t1 = m(1) * (t132 ^ 2 + t133 ^ 2 + t134 ^ 2) / 0.2e1 + m(2) * (t97 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(3) * (t72 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + m(4) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + m(5) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + m(6) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + t113 * (t234 * t163 + t183 * t164) / 0.2e1 + t112 * (t183 * t163 - t234 * t164) / 0.2e1 + (t238 * t163 - t239 * t164) * t135 / 0.2e1 + (t239 * t163 + t238 * t164) * t136 / 0.2e1 + ((-t124 * t164 - t126 * t163 - t142 * t175 - t173 * t145 + Icges(1,6)) * V_base(6) + (-t164 * t123 - t163 * t125 - t175 * t141 - t173 * t144 + Icges(1,2)) * V_base(5)) * V_base(5) / 0.2e1 + ((-t124 * t163 + t126 * t164 - t173 * t142 + t145 * t175 + Icges(1,3)) * V_base(6) + (-t123 * t163 + t125 * t164 - t173 * t141 + t144 * t175 + Icges(1,6)) * V_base(5)) * V_base(6) / 0.2e1 + ((t155 * t83 + t156 * t81) * t113 + (t155 * t84 + t156 * t82) * t112 + (t102 * t174 + t104 * t172 + t159 * t93 + t160 * t91) * t136 + (t103 * t174 + t105 * t172 + t159 * t94 + t160 * t92) * t135 + (t156 * t110 + t155 * t111 + t160 * t117 + t159 * t118 + t174 * t140 + t172 * t143 + Icges(3,3)) * t157) * t157 / 0.2e1 + t157 * V_base(6) * (Icges(3,5) * t164 - Icges(3,6) * t163) + t157 * V_base(5) * (-Icges(3,5) * t163 - Icges(3,6) * t164) + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + Icges(1,1) * V_base(4) / 0.2e1) * V_base(4) + ((-Icges(2,5) * t173 - Icges(2,6) * t175) * V_base(5) + (Icges(2,5) * t175 - Icges(2,6) * t173) * V_base(6) + Icges(2,3) * t162 / 0.2e1) * t162;
T = t1;
