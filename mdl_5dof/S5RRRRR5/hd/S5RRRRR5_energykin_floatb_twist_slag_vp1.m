% Calculate kinetic energy for
% S5RRRRR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-05 18:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRRR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR5_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR5_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR5_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:58:06
% EndTime: 2019-12-05 18:58:08
% DurationCPUTime: 1.39s
% Computational Cost: add. (1207->225), mult. (760->315), div. (0->0), fcn. (538->10), ass. (0->118)
t153 = qJ(4) + qJ(5);
t147 = cos(t153);
t145 = sin(t153);
t189 = Icges(6,4) * t145;
t107 = Icges(6,2) * t147 + t189;
t188 = Icges(6,4) * t147;
t110 = Icges(6,1) * t145 + t188;
t144 = V_base(4) + qJD(1);
t140 = qJD(2) + t144;
t137 = qJD(3) + t140;
t154 = qJ(1) + qJ(2);
t149 = qJ(3) + t154;
t141 = sin(t149);
t142 = cos(t149);
t172 = -Icges(6,2) * t145 + t188;
t73 = Icges(6,6) * t142 - t141 * t172;
t74 = Icges(6,6) * t141 + t142 * t172;
t174 = Icges(6,1) * t147 - t189;
t75 = Icges(6,5) * t142 - t141 * t174;
t76 = Icges(6,5) * t141 + t142 * t174;
t119 = qJD(4) * t141 + V_base(6);
t91 = qJD(5) * t141 + t119;
t120 = qJD(4) * t142 + V_base(5);
t92 = qJD(5) * t142 + t120;
t208 = (t107 * t145 - t110 * t147) * t137 + (t145 * t73 - t147 * t75) * t92 + (t145 * t74 - t147 * t76) * t91;
t157 = cos(qJ(4));
t155 = sin(qJ(4));
t191 = Icges(5,4) * t155;
t124 = Icges(5,2) * t157 + t191;
t190 = Icges(5,4) * t157;
t127 = Icges(5,1) * t155 + t190;
t173 = -Icges(5,2) * t155 + t190;
t81 = Icges(5,6) * t142 - t141 * t173;
t82 = Icges(5,6) * t141 + t142 * t173;
t175 = Icges(5,1) * t157 - t191;
t83 = Icges(5,5) * t142 - t141 * t175;
t84 = Icges(5,5) * t141 + t142 * t175;
t207 = (t124 * t155 - t127 * t157) * t137 + (t155 * t81 - t157 * t83) * t120 + (t155 * t82 - t157 * t84) * t119;
t206 = -pkin(5) - pkin(6);
t156 = sin(qJ(1));
t204 = pkin(1) * t156;
t158 = cos(qJ(1));
t203 = pkin(1) * t158;
t146 = sin(t154);
t202 = pkin(2) * t146;
t148 = cos(t154);
t201 = pkin(2) * t148;
t200 = pkin(4) * t155;
t199 = pkin(4) * t157;
t197 = Icges(2,4) * t156;
t196 = Icges(2,4) * t158;
t195 = Icges(3,4) * t146;
t194 = Icges(3,4) * t148;
t193 = Icges(4,4) * t141;
t192 = Icges(4,4) * t142;
t187 = -pkin(7) + t206;
t186 = V_base(6) * pkin(5) + V_base(2);
t183 = V_base(5) * t203 + V_base(6) * t204 + V_base(1);
t182 = rSges(5,1) * t157 - rSges(5,2) * t155;
t181 = rSges(6,1) * t147 - rSges(6,2) * t145;
t180 = -t144 * t204 + V_base(3);
t171 = Icges(5,5) * t157 - Icges(5,6) * t155;
t170 = Icges(6,5) * t147 - Icges(6,6) * t145;
t167 = V_base(5) * t201 + V_base(6) * t202 + t183;
t166 = V_base(6) * pkin(6) - t144 * t203 + t186;
t165 = (Icges(6,5) * t145 + Icges(6,6) * t147) * t137 + (Icges(6,3) * t142 - t141 * t170) * t92 + (Icges(6,3) * t141 + t142 * t170) * t91;
t164 = t119 * (Icges(5,3) * t141 + t142 * t171) + t120 * (Icges(5,3) * t142 - t141 * t171) + (Icges(5,5) * t155 + Icges(5,6) * t157) * t137;
t163 = -t140 * t202 + t180;
t102 = -t141 * pkin(3) + t142 * pkin(8);
t103 = t142 * pkin(3) + t141 * pkin(8);
t162 = -t102 * V_base(6) + V_base(5) * t103 + t167;
t161 = V_base(6) * pkin(7) - t140 * t201 + t166;
t160 = t137 * t102 + t187 * V_base(5) + t163;
t132 = rSges(2,1) * t158 - rSges(2,2) * t156;
t131 = -rSges(2,1) * t156 - rSges(2,2) * t158;
t130 = rSges(5,1) * t155 + rSges(5,2) * t157;
t129 = Icges(2,1) * t158 - t197;
t128 = -Icges(2,1) * t156 - t196;
t126 = -Icges(2,2) * t156 + t196;
t125 = -Icges(2,2) * t158 - t197;
t118 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t117 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t116 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t115 = rSges(3,1) * t148 - rSges(3,2) * t146;
t114 = -rSges(3,1) * t146 - rSges(3,2) * t148;
t113 = rSges(6,1) * t145 + rSges(6,2) * t147;
t112 = Icges(3,1) * t148 - t195;
t111 = -Icges(3,1) * t146 - t194;
t109 = -Icges(3,2) * t146 + t194;
t108 = -Icges(3,2) * t148 - t195;
t101 = rSges(4,1) * t142 - rSges(4,2) * t141;
t100 = -rSges(4,1) * t141 - rSges(4,2) * t142;
t99 = Icges(4,1) * t142 - t193;
t98 = -Icges(4,1) * t141 - t192;
t97 = -Icges(4,2) * t141 + t192;
t96 = -Icges(4,2) * t142 - t193;
t89 = V_base(6) * rSges(2,3) - t132 * t144 + t186;
t88 = t131 * t144 + V_base(3) + (-rSges(2,3) - pkin(5)) * V_base(5);
t87 = -t131 * V_base(6) + t132 * V_base(5) + V_base(1);
t86 = rSges(5,3) * t141 + t142 * t182;
t85 = rSges(5,3) * t142 - t141 * t182;
t78 = rSges(6,3) * t141 + t142 * t181;
t77 = rSges(6,3) * t142 - t141 * t181;
t70 = pkin(9) * t141 + t142 * t199;
t69 = pkin(9) * t142 - t141 * t199;
t68 = V_base(6) * rSges(3,3) - t115 * t140 + t166;
t67 = t114 * t140 + (-rSges(3,3) + t206) * V_base(5) + t180;
t66 = -t114 * V_base(6) + t115 * V_base(5) + t183;
t65 = V_base(6) * rSges(4,3) - t101 * t137 + t161;
t64 = t100 * t137 + (-rSges(4,3) + t187) * V_base(5) + t163;
t63 = -t100 * V_base(6) + t101 * V_base(5) + t167;
t62 = t119 * t130 + (-t103 - t86) * t137 + t161;
t61 = -t120 * t130 + t137 * t85 + t160;
t60 = -t119 * t85 + t120 * t86 + t162;
t59 = t119 * t200 + t113 * t91 + (-t103 - t70 - t78) * t137 + t161;
t58 = -t120 * t200 - t113 * t92 + (t69 + t77) * t137 + t160;
t57 = -t119 * t69 + t120 * t70 - t77 * t91 + t78 * t92 + t162;
t1 = m(1) * (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + m(2) * (t87 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + m(3) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + m(4) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + m(5) * (t60 ^ 2 + t61 ^ 2 + t62 ^ 2) / 0.2e1 + t120 * (t207 * t141 + t164 * t142) / 0.2e1 + t119 * (t164 * t141 - t207 * t142) / 0.2e1 + m(6) * (t57 ^ 2 + t58 ^ 2 + t59 ^ 2) / 0.2e1 + t92 * (t208 * t141 + t165 * t142) / 0.2e1 + t91 * (t165 * t141 - t208 * t142) / 0.2e1 + ((t155 * t83 + t157 * t81) * t120 + (t155 * t84 + t157 * t82) * t119 + (t145 * t75 + t147 * t73) * t92 + (t145 * t76 + t147 * t74) * t91 + (t107 * t147 + t110 * t145 + t124 * t157 + t127 * t155 + Icges(4,3)) * t137) * t137 / 0.2e1 + ((-t109 * t148 - t112 * t146 - t126 * t158 - t129 * t156 - t141 * t99 - t142 * t97 + Icges(1,6)) * V_base(6) + (-t108 * t148 - t111 * t146 - t125 * t158 - t128 * t156 - t141 * t98 - t142 * t96 + Icges(1,2)) * V_base(5)) * V_base(5) / 0.2e1 + ((-t109 * t146 + t112 * t148 - t126 * t156 + t129 * t158 - t141 * t97 + t142 * t99 + Icges(1,3)) * V_base(6) + (-t108 * t146 + t111 * t148 - t125 * t156 + t128 * t158 - t141 * t96 + t142 * t98 + Icges(1,6)) * V_base(5)) * V_base(6) / 0.2e1 + t137 * V_base(6) * (Icges(4,5) * t142 - Icges(4,6) * t141) + t137 * V_base(5) * (-Icges(4,5) * t141 - Icges(4,6) * t142) + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + Icges(1,1) * V_base(4) / 0.2e1) * V_base(4) + ((-Icges(2,5) * t156 - Icges(2,6) * t158) * V_base(5) + (Icges(2,5) * t158 - Icges(2,6) * t156) * V_base(6) + Icges(2,3) * t144 / 0.2e1) * t144 + ((-Icges(3,5) * t146 - Icges(3,6) * t148) * V_base(5) + (Icges(3,5) * t148 - Icges(3,6) * t146) * V_base(6) + Icges(3,3) * t140 / 0.2e1) * t140;
T = t1;
