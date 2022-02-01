% Calculate kinetic energy for
% S5RPRPR2
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
% Datum: 2022-01-23 09:19
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRPR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR2_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR2_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR2_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:18:45
% EndTime: 2022-01-23 09:18:47
% DurationCPUTime: 1.44s
% Computational Cost: add. (1103->232), mult. (736->300), div. (0->0), fcn. (514->10), ass. (0->120)
t160 = sin(qJ(1));
t203 = pkin(1) * t160;
t161 = cos(qJ(1));
t202 = pkin(1) * t161;
t156 = qJ(1) + pkin(8);
t146 = sin(t156);
t201 = pkin(2) * t146;
t148 = cos(t156);
t200 = pkin(2) * t148;
t157 = sin(pkin(9));
t199 = pkin(4) * t157;
t158 = cos(pkin(9));
t198 = pkin(4) * t158;
t197 = -pkin(5) - qJ(2);
t196 = Icges(2,4) * t160;
t195 = Icges(3,4) * t146;
t149 = qJ(3) + t156;
t141 = sin(t149);
t194 = Icges(4,4) * t141;
t193 = Icges(5,4) * t157;
t192 = Icges(5,4) * t158;
t155 = pkin(9) + qJ(5);
t145 = sin(t155);
t191 = Icges(6,4) * t145;
t147 = cos(t155);
t190 = Icges(6,4) * t147;
t188 = -pkin(6) + t197;
t150 = V_base(6) + qJD(1);
t187 = t150 * t202 + V_base(2);
t186 = V_base(5) * pkin(5) + V_base(1);
t183 = t150 * t200 + t187;
t182 = V_base(5) * qJ(2) + t186;
t181 = V_base(4) * t203 + qJD(2) + V_base(3);
t180 = -t200 - t202;
t179 = V_base(4) * t201 + t181;
t178 = rSges(5,1) * t158 - rSges(5,2) * t157;
t177 = rSges(6,1) * t147 - rSges(6,2) * t145;
t176 = Icges(5,1) * t158 - t193;
t175 = Icges(6,1) * t147 - t191;
t174 = -Icges(5,2) * t157 + t192;
t173 = -Icges(6,2) * t145 + t190;
t172 = Icges(5,5) * t158 - Icges(5,6) * t157;
t171 = Icges(6,5) * t147 - Icges(6,6) * t145;
t142 = cos(t149);
t101 = pkin(3) * t141 - qJ(4) * t142;
t170 = V_base(4) * t101 + t179;
t103 = pkin(3) * t142 + qJ(4) * t141;
t169 = -t103 + t180;
t143 = qJD(3) + t150;
t168 = -qJD(4) * t142 + t143 * t103 + t183;
t117 = -qJD(5) * t142 + V_base(5);
t118 = qJD(5) * t141 + V_base(4);
t167 = (Icges(6,5) * t145 + Icges(6,6) * t147) * t143 + t117 * (-Icges(6,3) * t142 + t171 * t141) + t118 * (Icges(6,3) * t141 + t171 * t142);
t166 = (Icges(5,5) * t157 + Icges(5,6) * t158) * t143 + (-Icges(5,3) * t142 + t172 * t141) * V_base(5) + (Icges(5,3) * t141 + t172 * t142) * V_base(4);
t165 = V_base(5) * pkin(6) + (-t201 - t203) * t150 + t182;
t164 = qJD(4) * t141 + t165;
t108 = Icges(6,2) * t147 + t191;
t111 = Icges(6,1) * t145 + t190;
t76 = -Icges(6,6) * t142 + t141 * t173;
t77 = Icges(6,6) * t141 + t142 * t173;
t78 = -Icges(6,5) * t142 + t141 * t175;
t79 = Icges(6,5) * t141 + t142 * t175;
t163 = (-t145 * t77 + t147 * t79) * t118 + (-t145 * t76 + t147 * t78) * t117 + (-t108 * t145 + t111 * t147) * t143;
t124 = Icges(5,2) * t158 + t193;
t125 = Icges(5,1) * t157 + t192;
t84 = -Icges(5,6) * t142 + t141 * t174;
t85 = Icges(5,6) * t141 + t142 * t174;
t86 = -Icges(5,5) * t142 + t141 * t176;
t87 = Icges(5,5) * t141 + t142 * t176;
t162 = (-t157 * t85 + t158 * t87) * V_base(4) + (-t157 * t84 + t158 * t86) * V_base(5) + (-t124 * t157 + t125 * t158) * t143;
t152 = Icges(2,4) * t161;
t140 = Icges(3,4) * t148;
t138 = Icges(4,4) * t142;
t134 = rSges(2,1) * t161 - t160 * rSges(2,2);
t133 = t160 * rSges(2,1) + rSges(2,2) * t161;
t132 = Icges(2,1) * t161 - t196;
t131 = Icges(2,1) * t160 + t152;
t130 = -Icges(2,2) * t160 + t152;
t129 = Icges(2,2) * t161 + t196;
t126 = rSges(5,1) * t157 + rSges(5,2) * t158;
t121 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t120 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t119 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t116 = rSges(3,1) * t148 - rSges(3,2) * t146;
t115 = rSges(3,1) * t146 + rSges(3,2) * t148;
t114 = rSges(6,1) * t145 + rSges(6,2) * t147;
t113 = Icges(3,1) * t148 - t195;
t112 = Icges(3,1) * t146 + t140;
t110 = -Icges(3,2) * t146 + t140;
t109 = Icges(3,2) * t148 + t195;
t104 = rSges(4,1) * t142 - rSges(4,2) * t141;
t102 = rSges(4,1) * t141 + rSges(4,2) * t142;
t100 = Icges(4,1) * t142 - t194;
t99 = Icges(4,1) * t141 + t138;
t98 = -Icges(4,2) * t141 + t138;
t97 = Icges(4,2) * t142 + t194;
t96 = Icges(4,5) * t142 - Icges(4,6) * t141;
t95 = Icges(4,5) * t141 + Icges(4,6) * t142;
t92 = V_base(5) * rSges(2,3) - t133 * t150 + t186;
t91 = t134 * t150 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t90 = t133 * V_base(4) - t134 * V_base(5) + V_base(3);
t89 = rSges(5,3) * t141 + t142 * t178;
t88 = -rSges(5,3) * t142 + t141 * t178;
t81 = rSges(6,3) * t141 + t142 * t177;
t80 = -rSges(6,3) * t142 + t141 * t177;
t73 = pkin(7) * t141 + t142 * t198;
t72 = -pkin(7) * t142 + t141 * t198;
t71 = V_base(5) * rSges(3,3) + (-t115 - t203) * t150 + t182;
t70 = t116 * t150 + (-rSges(3,3) + t197) * V_base(4) + t187;
t69 = V_base(4) * t115 + (-t116 - t202) * V_base(5) + t181;
t68 = V_base(5) * rSges(4,3) - t102 * t143 + t165;
t67 = t104 * t143 + (-rSges(4,3) + t188) * V_base(4) + t183;
t66 = V_base(4) * t102 + (-t104 + t180) * V_base(5) + t179;
t65 = t126 * V_base(5) + (-t101 - t88) * t143 + t164;
t64 = t143 * t89 + (-t126 + t188) * V_base(4) + t168;
t63 = V_base(4) * t88 + (t169 - t89) * V_base(5) + t170;
t62 = V_base(5) * t199 + t114 * t117 + (-t101 - t72 - t80) * t143 + t164;
t61 = -t114 * t118 + (t73 + t81) * t143 + (t188 - t199) * V_base(4) + t168;
t60 = -t117 * t81 + t118 * t80 + V_base(4) * t72 + (t169 - t73) * V_base(5) + t170;
t1 = m(1) * (t119 ^ 2 + t120 ^ 2 + t121 ^ 2) / 0.2e1 + m(2) * (t90 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(3) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + m(4) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + m(5) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + m(6) * (t60 ^ 2 + t61 ^ 2 + t62 ^ 2) / 0.2e1 + t118 * (t167 * t141 + t163 * t142) / 0.2e1 + t117 * (t163 * t141 - t167 * t142) / 0.2e1 + ((t145 * t79 + t147 * t77) * t118 + (t145 * t78 + t147 * t76) * t117 + (t157 * t86 + t158 * t84 + t95) * V_base(5) + (t157 * t87 + t158 * t85 + t96) * V_base(4) + (t108 * t147 + t111 * t145 + t124 * t158 + t125 * t157 + Icges(4,3)) * t143) * t143 / 0.2e1 + (t166 * t141 + t162 * t142 + t96 * t143 + (-t109 * t146 + t112 * t148 - t160 * t129 + t131 * t161 - t141 * t97 + t142 * t99 + Icges(1,4)) * V_base(5) + (t100 * t142 - t146 * t110 + t113 * t148 - t160 * t130 + t132 * t161 - t141 * t98 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t162 * t141 - t166 * t142 + t95 * t143 + (t109 * t148 + t112 * t146 + t129 * t161 + t160 * t131 + t141 * t99 + t142 * t97 + Icges(1,2)) * V_base(5) + (t100 * t141 + t110 * t148 + t113 * t146 + t130 * t161 + t160 * t132 + t142 * t98 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t161 + Icges(3,5) * t148 - Icges(2,6) * t160 - Icges(3,6) * t146) * V_base(4) + (Icges(2,5) * t160 + Icges(3,5) * t146 + Icges(2,6) * t161 + Icges(3,6) * t148) * V_base(5) + (Icges(2,3) / 0.2e1 + Icges(3,3) / 0.2e1) * t150) * t150 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
