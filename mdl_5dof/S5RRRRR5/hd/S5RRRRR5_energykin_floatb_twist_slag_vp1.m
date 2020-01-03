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
% Datum: 2020-01-03 12:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:13:01
% EndTime: 2020-01-03 12:13:02
% DurationCPUTime: 1.31s
% Computational Cost: add. (1207->230), mult. (760->314), div. (0->0), fcn. (538->10), ass. (0->119)
t154 = qJ(4) + qJ(5);
t147 = cos(t154);
t145 = sin(t154);
t191 = Icges(6,4) * t145;
t107 = Icges(6,2) * t147 + t191;
t190 = Icges(6,4) * t147;
t110 = Icges(6,1) * t145 + t190;
t144 = V_base(4) + qJD(1);
t140 = qJD(2) + t144;
t138 = qJD(3) + t140;
t155 = qJ(1) + qJ(2);
t150 = qJ(3) + t155;
t141 = sin(t150);
t142 = cos(t150);
t170 = -Icges(6,2) * t145 + t190;
t73 = -Icges(6,6) * t142 + t141 * t170;
t74 = -Icges(6,6) * t141 - t142 * t170;
t172 = Icges(6,1) * t147 - t191;
t75 = -Icges(6,5) * t142 + t141 * t172;
t76 = -Icges(6,5) * t141 - t142 * t172;
t188 = -qJD(4) - qJD(5);
t91 = t141 * t188 + V_base(6);
t92 = t142 * t188 + V_base(5);
t207 = (t107 * t145 - t110 * t147) * t138 + (t145 * t73 - t147 * t75) * t92 + (t145 * t74 - t147 * t76) * t91;
t121 = -qJD(4) * t141 + V_base(6);
t122 = -qJD(4) * t142 + V_base(5);
t158 = cos(qJ(4));
t156 = sin(qJ(4));
t193 = Icges(5,4) * t156;
t126 = Icges(5,2) * t158 + t193;
t192 = Icges(5,4) * t158;
t129 = Icges(5,1) * t156 + t192;
t171 = -Icges(5,2) * t156 + t192;
t81 = -Icges(5,6) * t142 + t141 * t171;
t82 = -Icges(5,6) * t141 - t142 * t171;
t173 = Icges(5,1) * t158 - t193;
t83 = -Icges(5,5) * t142 + t141 * t173;
t84 = -Icges(5,5) * t141 - t142 * t173;
t206 = (t126 * t156 - t129 * t158) * t138 + (t156 * t81 - t158 * t83) * t122 + (t156 * t82 - t158 * t84) * t121;
t205 = -pkin(5) - pkin(6);
t157 = sin(qJ(1));
t203 = pkin(1) * t157;
t159 = cos(qJ(1));
t202 = pkin(1) * t159;
t146 = sin(t155);
t201 = pkin(2) * t146;
t148 = cos(t155);
t200 = pkin(2) * t148;
t199 = pkin(4) * t156;
t198 = pkin(4) * t158;
t196 = Icges(2,4) * t159;
t195 = Icges(3,4) * t148;
t194 = Icges(4,4) * t142;
t189 = -pkin(7) + t205;
t187 = t144 * t203 + V_base(3);
t186 = V_base(6) * pkin(5) + V_base(2);
t183 = t140 * t201 + t187;
t182 = V_base(6) * pkin(6) + t144 * t202 + t186;
t181 = -t201 - t203;
t180 = -t200 - t202;
t179 = rSges(5,1) * t158 - rSges(5,2) * t156;
t178 = rSges(6,1) * t147 - rSges(6,2) * t145;
t169 = Icges(5,5) * t158 - Icges(5,6) * t156;
t168 = Icges(6,5) * t147 - Icges(6,6) * t145;
t165 = V_base(6) * pkin(7) + t140 * t200 + t182;
t164 = -(Icges(6,5) * t145 + Icges(6,6) * t147) * t138 - (-Icges(6,3) * t142 + t141 * t168) * t92 - (-Icges(6,3) * t141 - t142 * t168) * t91;
t163 = -t121 * (-Icges(5,3) * t141 - t142 * t169) - t122 * (-Icges(5,3) * t142 + t141 * t169) - (Icges(5,5) * t156 + Icges(5,6) * t158) * t138;
t102 = t141 * pkin(3) - t142 * pkin(8);
t162 = t138 * t102 + t189 * V_base(5) + t183;
t103 = -t142 * pkin(3) - t141 * pkin(8);
t161 = V_base(1) + (-t102 + t181) * V_base(6) + (t103 + t180) * V_base(5);
t149 = Icges(2,4) * t157;
t139 = Icges(3,4) * t146;
t137 = Icges(4,4) * t141;
t134 = -rSges(2,1) * t159 + rSges(2,2) * t157;
t133 = rSges(2,1) * t157 + rSges(2,2) * t159;
t132 = rSges(5,1) * t156 + rSges(5,2) * t158;
t131 = -Icges(2,1) * t159 + t149;
t130 = Icges(2,1) * t157 + t196;
t128 = Icges(2,2) * t157 - t196;
t127 = Icges(2,2) * t159 + t149;
t120 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t119 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t118 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t115 = -rSges(3,1) * t148 + rSges(3,2) * t146;
t114 = rSges(3,1) * t146 + rSges(3,2) * t148;
t113 = rSges(6,1) * t145 + rSges(6,2) * t147;
t112 = -Icges(3,1) * t148 + t139;
t111 = Icges(3,1) * t146 + t195;
t109 = Icges(3,2) * t146 - t195;
t108 = Icges(3,2) * t148 + t139;
t101 = -rSges(4,1) * t142 + rSges(4,2) * t141;
t100 = rSges(4,1) * t141 + rSges(4,2) * t142;
t99 = -Icges(4,1) * t142 + t137;
t98 = Icges(4,1) * t141 + t194;
t97 = Icges(4,2) * t141 - t194;
t96 = Icges(4,2) * t142 + t137;
t89 = V_base(6) * rSges(2,3) - t134 * t144 + t186;
t88 = t133 * t144 + V_base(3) + (-rSges(2,3) - pkin(5)) * V_base(5);
t87 = -t133 * V_base(6) + t134 * V_base(5) + V_base(1);
t86 = -rSges(5,3) * t141 - t142 * t179;
t85 = -rSges(5,3) * t142 + t141 * t179;
t78 = -rSges(6,3) * t141 - t142 * t178;
t77 = -rSges(6,3) * t142 + t141 * t178;
t70 = -pkin(9) * t141 - t142 * t198;
t69 = -pkin(9) * t142 + t141 * t198;
t68 = V_base(6) * rSges(3,3) - t115 * t140 + t182;
t67 = t114 * t140 + (-rSges(3,3) + t205) * V_base(5) + t187;
t66 = -t114 * V_base(6) + t115 * V_base(5) + V_base(1) + (-V_base(6) * t157 - t159 * V_base(5)) * pkin(1);
t65 = V_base(6) * rSges(4,3) - t101 * t138 + t165;
t64 = t100 * t138 + (-rSges(4,3) + t189) * V_base(5) + t183;
t63 = V_base(1) + (-t100 + t181) * V_base(6) + (t101 + t180) * V_base(5);
t62 = t121 * t132 + (-t103 - t86) * t138 + t165;
t61 = -t122 * t132 + t138 * t85 + t162;
t60 = -t121 * t85 + t122 * t86 + t161;
t59 = t121 * t199 + t113 * t91 + (-t103 - t70 - t78) * t138 + t165;
t58 = -t122 * t199 - t113 * t92 + (t69 + t77) * t138 + t162;
t57 = -t121 * t69 + t122 * t70 - t77 * t91 + t78 * t92 + t161;
t1 = m(1) * (t118 ^ 2 + t119 ^ 2 + t120 ^ 2) / 0.2e1 + m(2) * (t87 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + m(3) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + m(4) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + m(5) * (t60 ^ 2 + t61 ^ 2 + t62 ^ 2) / 0.2e1 + t122 * (-t206 * t141 + t163 * t142) / 0.2e1 + t121 * (t163 * t141 + t206 * t142) / 0.2e1 + m(6) * (t57 ^ 2 + t58 ^ 2 + t59 ^ 2) / 0.2e1 + t92 * (-t207 * t141 + t164 * t142) / 0.2e1 + t91 * (t164 * t141 + t207 * t142) / 0.2e1 + ((t156 * t83 + t158 * t81) * t122 + (t156 * t84 + t158 * t82) * t121 + (t145 * t75 + t147 * t73) * t92 + (t145 * t76 + t147 * t74) * t91 + (t107 * t147 + t110 * t145 + t126 * t158 + t129 * t156 + Icges(4,3)) * t138) * t138 / 0.2e1 + ((t109 * t148 + t112 * t146 + t128 * t159 + t131 * t157 + t141 * t99 + t142 * t97 + Icges(1,6)) * V_base(6) + (t108 * t148 + t111 * t146 + t127 * t159 + t130 * t157 + t141 * t98 + t142 * t96 + Icges(1,2)) * V_base(5)) * V_base(5) / 0.2e1 + ((t109 * t146 - t112 * t148 + t128 * t157 - t131 * t159 + t141 * t97 - t142 * t99 + Icges(1,3)) * V_base(6) + (t108 * t146 - t111 * t148 + t127 * t157 - t130 * t159 + t141 * t96 - t142 * t98 + Icges(1,6)) * V_base(5)) * V_base(6) / 0.2e1 + t138 * V_base(6) * (-Icges(4,5) * t142 + Icges(4,6) * t141) + V_base(5) * t138 * (Icges(4,5) * t141 + Icges(4,6) * t142) + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + Icges(1,1) * V_base(4) / 0.2e1) * V_base(4) + ((Icges(2,5) * t157 + Icges(2,6) * t159) * V_base(5) + (-Icges(2,5) * t159 + Icges(2,6) * t157) * V_base(6) + Icges(2,3) * t144 / 0.2e1) * t144 + ((Icges(3,5) * t146 + Icges(3,6) * t148) * V_base(5) + (-Icges(3,5) * t148 + Icges(3,6) * t146) * V_base(6) + Icges(3,3) * t140 / 0.2e1) * t140;
T = t1;
