% Calculate kinetic energy for
% S5RPRPR4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-05 17:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRPR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR4_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR4_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR4_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:53:05
% EndTime: 2019-12-05 17:53:07
% DurationCPUTime: 1.82s
% Computational Cost: add. (1262->244), mult. (946->323), div. (0->0), fcn. (726->10), ass. (0->125)
t242 = Icges(4,3) + Icges(5,3);
t168 = qJ(3) + pkin(9);
t158 = sin(t168);
t160 = cos(t168);
t171 = sin(qJ(3));
t173 = cos(qJ(3));
t241 = Icges(4,5) * t173 + Icges(5,5) * t160 - Icges(4,6) * t171 - Icges(5,6) * t158;
t169 = qJ(1) + pkin(8);
t159 = sin(t169);
t161 = cos(t169);
t217 = Icges(4,4) * t173;
t193 = -Icges(4,2) * t171 + t217;
t100 = Icges(4,6) * t161 - t159 * t193;
t101 = Icges(4,6) * t159 + t161 * t193;
t218 = Icges(4,4) * t171;
t196 = Icges(4,1) * t173 - t218;
t102 = Icges(4,5) * t161 - t159 * t196;
t103 = Icges(4,5) * t159 + t161 * t196;
t216 = Icges(5,4) * t158;
t120 = Icges(5,2) * t160 + t216;
t215 = Icges(5,4) * t160;
t123 = Icges(5,1) * t158 + t215;
t135 = qJD(3) * t159 + V_base(6);
t136 = qJD(3) * t161 + V_base(5);
t140 = Icges(4,2) * t173 + t218;
t143 = Icges(4,1) * t171 + t217;
t163 = V_base(4) + qJD(1);
t192 = -Icges(5,2) * t158 + t215;
t91 = Icges(5,6) * t161 - t159 * t192;
t92 = Icges(5,6) * t159 + t161 * t192;
t195 = Icges(5,1) * t160 - t216;
t93 = Icges(5,5) * t161 - t159 * t195;
t94 = Icges(5,5) * t159 + t161 * t195;
t240 = t135 * (t101 * t171 - t103 * t173 + t158 * t92 - t160 * t94) + t136 * (t100 * t171 - t102 * t173 + t158 * t91 - t160 * t93) + t163 * (t120 * t158 - t123 * t160 + t140 * t171 - t143 * t173);
t237 = (Icges(4,5) * t171 + Icges(5,5) * t158 + Icges(4,6) * t173 + Icges(5,6) * t160) * t163 + (-t241 * t159 + t242 * t161) * t136 + (t242 * t159 + t241 * t161) * t135;
t109 = qJD(5) * t159 + t135;
t110 = qJD(5) * t161 + t136;
t162 = qJ(5) + t168;
t156 = cos(t162);
t155 = sin(t162);
t214 = Icges(6,4) * t155;
t113 = Icges(6,2) * t156 + t214;
t213 = Icges(6,4) * t156;
t114 = Icges(6,1) * t155 + t213;
t191 = -Icges(6,2) * t155 + t213;
t83 = Icges(6,6) * t161 - t159 * t191;
t84 = Icges(6,6) * t159 + t161 * t191;
t194 = Icges(6,1) * t156 - t214;
t85 = Icges(6,5) * t161 - t159 * t194;
t86 = Icges(6,5) * t159 + t161 * t194;
t233 = (t113 * t155 - t114 * t156) * t163 + (t155 * t83 - t156 * t85) * t110 + (t155 * t84 - t156 * t86) * t109;
t172 = sin(qJ(1));
t229 = pkin(1) * t172;
t174 = cos(qJ(1));
t228 = pkin(1) * t174;
t227 = pkin(3) * t171;
t226 = pkin(4) * t158;
t225 = t173 * pkin(3);
t224 = -pkin(5) - qJ(2);
t222 = Icges(2,4) * t172;
t221 = Icges(2,4) * t174;
t220 = Icges(3,4) * t159;
t219 = Icges(3,4) * t161;
t212 = pkin(4) * t160;
t210 = V_base(6) * pkin(5) + V_base(2);
t130 = pkin(2) * t161 + pkin(6) * t159;
t207 = -t130 - t228;
t206 = V_base(6) * qJ(2) + t210;
t80 = qJ(4) * t159 + t161 * t225;
t205 = t207 - t80;
t204 = V_base(5) * t228 + V_base(6) * t229 + qJD(2) + V_base(1);
t203 = rSges(4,1) * t173 - rSges(4,2) * t171;
t202 = rSges(5,1) * t160 - rSges(5,2) * t158;
t201 = rSges(6,1) * t156 - rSges(6,2) * t155;
t188 = Icges(6,5) * t156 - Icges(6,6) * t155;
t182 = qJD(4) * t161 + t135 * t227 + t206;
t181 = (Icges(6,3) * t159 + t161 * t188) * t109 + (Icges(6,3) * t161 - t159 * t188) * t110 + (Icges(6,5) * t155 + Icges(6,6) * t156) * t163;
t129 = -pkin(2) * t159 + pkin(6) * t161;
t178 = t163 * t129 + t224 * V_base(5) + V_base(3);
t79 = qJ(4) * t161 - t159 * t225;
t177 = qJD(4) * t159 + t163 * t79 + t178;
t176 = -t129 * V_base(6) + V_base(5) * t130 + t204;
t175 = t136 * t80 + t176;
t148 = rSges(2,1) * t174 - t172 * rSges(2,2);
t147 = -t172 * rSges(2,1) - rSges(2,2) * t174;
t146 = rSges(4,1) * t171 + rSges(4,2) * t173;
t145 = Icges(2,1) * t174 - t222;
t144 = -Icges(2,1) * t172 - t221;
t142 = -Icges(2,2) * t172 + t221;
t141 = -Icges(2,2) * t174 - t222;
t134 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t133 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t132 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t128 = rSges(3,1) * t161 - rSges(3,2) * t159;
t127 = -rSges(3,1) * t159 - rSges(3,2) * t161;
t126 = rSges(5,1) * t158 + rSges(5,2) * t160;
t125 = Icges(3,1) * t161 - t220;
t124 = -Icges(3,1) * t159 - t219;
t122 = -Icges(3,2) * t159 + t219;
t121 = -Icges(3,2) * t161 - t220;
t115 = rSges(6,1) * t155 + rSges(6,2) * t156;
t107 = rSges(4,3) * t159 + t161 * t203;
t106 = rSges(4,3) * t161 - t159 * t203;
t105 = V_base(6) * rSges(2,3) - t148 * t163 + t210;
t104 = t147 * t163 + V_base(3) + (-rSges(2,3) - pkin(5)) * V_base(5);
t97 = -t147 * V_base(6) + t148 * V_base(5) + V_base(1);
t96 = rSges(5,3) * t159 + t161 * t202;
t95 = rSges(5,3) * t161 - t159 * t202;
t88 = rSges(6,3) * t159 + t161 * t201;
t87 = rSges(6,3) * t161 - t159 * t201;
t77 = V_base(6) * rSges(3,3) + (-t128 - t228) * t163 + t206;
t76 = V_base(3) + (t127 - t229) * t163 + (-rSges(3,3) + t224) * V_base(5);
t74 = pkin(7) * t159 + t161 * t212;
t73 = pkin(7) * t161 - t159 * t212;
t72 = -t127 * V_base(6) + t128 * V_base(5) + t204;
t71 = t135 * t146 + (-t107 + t207) * t163 + t206;
t70 = -t136 * t146 + (t106 - t229) * t163 + t178;
t69 = -t106 * t135 + t107 * t136 + t176;
t68 = t135 * t126 + (t205 - t96) * t163 + t182;
t67 = (t95 - t229) * t163 + (-t126 - t227) * t136 + t177;
t66 = t136 * t96 + (-t79 - t95) * t135 + t175;
t65 = t135 * t226 + t109 * t115 + (t205 - t74 - t88) * t163 + t182;
t64 = -t110 * t115 + (-t226 - t227) * t136 + (t73 + t87 - t229) * t163 + t177;
t63 = -t109 * t87 + t110 * t88 + t136 * t74 + (-t73 - t79) * t135 + t175;
t1 = m(1) * (t132 ^ 2 + t133 ^ 2 + t134 ^ 2) / 0.2e1 + m(2) * (t104 ^ 2 + t105 ^ 2 + t97 ^ 2) / 0.2e1 + m(3) * (t72 ^ 2 + t76 ^ 2 + t77 ^ 2) / 0.2e1 + m(4) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + m(5) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + m(6) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + t110 * (t233 * t159 + t181 * t161) / 0.2e1 + t109 * (t181 * t159 - t233 * t161) / 0.2e1 + (t237 * t159 - t240 * t161) * t135 / 0.2e1 + (t240 * t159 + t237 * t161) * t136 / 0.2e1 + ((-t122 * t161 - t125 * t159 - t142 * t174 - t172 * t145 + Icges(1,6)) * V_base(6) + (-t121 * t161 - t124 * t159 - t141 * t174 - t172 * t144 + Icges(1,2)) * V_base(5)) * V_base(5) / 0.2e1 + ((-t122 * t159 + t125 * t161 - t172 * t142 + t145 * t174 + Icges(1,3)) * V_base(6) + (-t121 * t159 + t124 * t161 - t172 * t141 + t144 * t174 + Icges(1,6)) * V_base(5)) * V_base(6) / 0.2e1 + ((t155 * t85 + t156 * t83) * t110 + (t155 * t86 + t156 * t84) * t109 + (t100 * t173 + t102 * t171 + t158 * t93 + t160 * t91) * t136 + (t101 * t173 + t103 * t171 + t158 * t94 + t160 * t92) * t135 + (t113 * t156 + t114 * t155 + t120 * t160 + t123 * t158 + t140 * t173 + t143 * t171 + Icges(2,3) + Icges(3,3)) * t163) * t163 / 0.2e1 + t163 * V_base(6) * (Icges(2,5) * t174 + Icges(3,5) * t161 - Icges(2,6) * t172 - Icges(3,6) * t159) + t163 * V_base(5) * (-Icges(2,5) * t172 - Icges(3,5) * t159 - Icges(2,6) * t174 - Icges(3,6) * t161) + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + Icges(1,1) * V_base(4) / 0.2e1) * V_base(4);
T = t1;
