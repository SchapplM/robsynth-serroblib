% Calculate kinetic energy for
% S5RRPRR4
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
% Datum: 2019-12-05 18:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPRR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR4_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR4_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR4_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:31:47
% EndTime: 2019-12-05 18:31:48
% DurationCPUTime: 1.30s
% Computational Cost: add. (1186->228), mult. (760->308), div. (0->0), fcn. (538->10), ass. (0->117)
t152 = qJ(4) + qJ(5);
t147 = cos(t152);
t145 = sin(t152);
t188 = Icges(6,4) * t145;
t107 = Icges(6,2) * t147 + t188;
t187 = Icges(6,4) * t147;
t110 = Icges(6,1) * t145 + t187;
t144 = V_base(4) + qJD(1);
t141 = qJD(2) + t144;
t153 = qJ(1) + qJ(2);
t143 = pkin(9) + t153;
t139 = sin(t143);
t140 = cos(t143);
t170 = -Icges(6,2) * t145 + t187;
t73 = Icges(6,6) * t140 - t139 * t170;
t74 = Icges(6,6) * t139 + t140 * t170;
t172 = Icges(6,1) * t147 - t188;
t75 = Icges(6,5) * t140 - t139 * t172;
t76 = Icges(6,5) * t139 + t140 * t172;
t116 = qJD(4) * t139 + V_base(6);
t91 = qJD(5) * t139 + t116;
t117 = qJD(4) * t140 + V_base(5);
t92 = qJD(5) * t140 + t117;
t207 = (t107 * t145 - t110 * t147) * t141 + (t145 * t73 - t147 * t75) * t92 + (t145 * t74 - t147 * t76) * t91;
t156 = cos(qJ(4));
t154 = sin(qJ(4));
t190 = Icges(5,4) * t154;
t124 = Icges(5,2) * t156 + t190;
t189 = Icges(5,4) * t156;
t127 = Icges(5,1) * t154 + t189;
t171 = -Icges(5,2) * t154 + t189;
t81 = Icges(5,6) * t140 - t139 * t171;
t82 = Icges(5,6) * t139 + t140 * t171;
t173 = Icges(5,1) * t156 - t190;
t83 = Icges(5,5) * t140 - t139 * t173;
t84 = Icges(5,5) * t139 + t140 * t173;
t206 = (t124 * t154 - t127 * t156) * t141 + (t154 * t81 - t156 * t83) * t117 + (t154 * t82 - t156 * t84) * t116;
t205 = -pkin(5) - pkin(6);
t155 = sin(qJ(1));
t203 = pkin(1) * t155;
t157 = cos(qJ(1));
t202 = pkin(1) * t157;
t146 = sin(t153);
t201 = pkin(2) * t146;
t148 = cos(t153);
t200 = pkin(2) * t148;
t199 = pkin(4) * t154;
t198 = pkin(4) * t156;
t196 = Icges(2,4) * t155;
t195 = Icges(2,4) * t157;
t194 = Icges(3,4) * t146;
t193 = Icges(3,4) * t148;
t192 = Icges(4,4) * t139;
t191 = Icges(4,4) * t140;
t186 = -qJ(3) + t205;
t185 = V_base(6) * pkin(5) + V_base(2);
t103 = t140 * pkin(3) + t139 * pkin(7);
t182 = -t103 - t200;
t181 = V_base(5) * t202 + V_base(6) * t203 + V_base(1);
t180 = rSges(5,1) * t156 - rSges(5,2) * t154;
t179 = rSges(6,1) * t147 - rSges(6,2) * t145;
t178 = -t144 * t203 + V_base(3);
t169 = Icges(5,5) * t156 - Icges(5,6) * t154;
t168 = Icges(6,5) * t147 - Icges(6,6) * t145;
t165 = V_base(5) * t200 + V_base(6) * t201 + qJD(3) + t181;
t164 = V_base(6) * pkin(6) - t144 * t202 + t185;
t163 = (Icges(6,5) * t145 + Icges(6,6) * t147) * t141 + (Icges(6,3) * t140 - t139 * t168) * t92 + (Icges(6,3) * t139 + t140 * t168) * t91;
t162 = t116 * (Icges(5,3) * t139 + t140 * t169) + t117 * (Icges(5,3) * t140 - t139 * t169) + (Icges(5,5) * t154 + Icges(5,6) * t156) * t141;
t161 = V_base(6) * qJ(3) + t164;
t102 = -t139 * pkin(3) + t140 * pkin(7);
t160 = -t102 * V_base(6) + V_base(5) * t103 + t165;
t159 = t141 * t102 + t186 * V_base(5) + t178;
t132 = rSges(2,1) * t157 - rSges(2,2) * t155;
t131 = -rSges(2,1) * t155 - rSges(2,2) * t157;
t130 = rSges(5,1) * t154 + rSges(5,2) * t156;
t129 = Icges(2,1) * t157 - t196;
t128 = -Icges(2,1) * t155 - t195;
t126 = -Icges(2,2) * t155 + t195;
t125 = -Icges(2,2) * t157 - t196;
t120 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t119 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t118 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t115 = rSges(3,1) * t148 - rSges(3,2) * t146;
t114 = -rSges(3,1) * t146 - rSges(3,2) * t148;
t113 = rSges(6,1) * t145 + rSges(6,2) * t147;
t112 = Icges(3,1) * t148 - t194;
t111 = -Icges(3,1) * t146 - t193;
t109 = -Icges(3,2) * t146 + t193;
t108 = -Icges(3,2) * t148 - t194;
t101 = rSges(4,1) * t140 - rSges(4,2) * t139;
t100 = -rSges(4,1) * t139 - rSges(4,2) * t140;
t99 = Icges(4,1) * t140 - t192;
t98 = -Icges(4,1) * t139 - t191;
t97 = -Icges(4,2) * t139 + t191;
t96 = -Icges(4,2) * t140 - t192;
t89 = V_base(6) * rSges(2,3) - t132 * t144 + t185;
t88 = t131 * t144 + V_base(3) + (-rSges(2,3) - pkin(5)) * V_base(5);
t87 = -t131 * V_base(6) + t132 * V_base(5) + V_base(1);
t86 = rSges(5,3) * t139 + t140 * t180;
t85 = rSges(5,3) * t140 - t139 * t180;
t78 = rSges(6,3) * t139 + t140 * t179;
t77 = rSges(6,3) * t140 - t139 * t179;
t70 = pkin(8) * t139 + t140 * t198;
t69 = pkin(8) * t140 - t139 * t198;
t68 = V_base(6) * rSges(3,3) - t115 * t141 + t164;
t67 = t114 * t141 + (-rSges(3,3) + t205) * V_base(5) + t178;
t66 = -t114 * V_base(6) + t115 * V_base(5) + t181;
t65 = V_base(6) * rSges(4,3) + (-t101 - t200) * t141 + t161;
t64 = (t100 - t201) * t141 + (-rSges(4,3) + t186) * V_base(5) + t178;
t63 = -t100 * V_base(6) + t101 * V_base(5) + t165;
t62 = t116 * t130 + (t182 - t86) * t141 + t161;
t61 = -t117 * t130 + (t85 - t201) * t141 + t159;
t60 = -t116 * t85 + t117 * t86 + t160;
t59 = t116 * t199 + t113 * t91 + (t182 - t70 - t78) * t141 + t161;
t58 = -t117 * t199 - t113 * t92 + (t69 + t77 - t201) * t141 + t159;
t57 = -t116 * t69 + t117 * t70 - t77 * t91 + t78 * t92 + t160;
t1 = m(1) * (t118 ^ 2 + t119 ^ 2 + t120 ^ 2) / 0.2e1 + m(2) * (t87 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + m(3) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + m(4) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + m(5) * (t60 ^ 2 + t61 ^ 2 + t62 ^ 2) / 0.2e1 + t117 * (t206 * t139 + t162 * t140) / 0.2e1 + t116 * (t162 * t139 - t206 * t140) / 0.2e1 + m(6) * (t57 ^ 2 + t58 ^ 2 + t59 ^ 2) / 0.2e1 + t92 * (t207 * t139 + t163 * t140) / 0.2e1 + t91 * (t163 * t139 - t207 * t140) / 0.2e1 + ((t154 * t83 + t156 * t81) * t117 + (t154 * t84 + t156 * t82) * t116 + (t145 * t75 + t147 * t73) * t92 + (t145 * t76 + t147 * t74) * t91 + (t107 * t147 + t110 * t145 + t124 * t156 + t127 * t154 + Icges(3,3) + Icges(4,3)) * t141) * t141 / 0.2e1 + ((-t109 * t148 - t112 * t146 - t126 * t157 - t129 * t155 - t139 * t99 - t140 * t97 + Icges(1,6)) * V_base(6) + (-t108 * t148 - t111 * t146 - t125 * t157 - t128 * t155 - t139 * t98 - t140 * t96 + Icges(1,2)) * V_base(5)) * V_base(5) / 0.2e1 + ((-t109 * t146 + t112 * t148 - t126 * t155 + t129 * t157 - t139 * t97 + t140 * t99 + Icges(1,3)) * V_base(6) + (-t108 * t146 + t111 * t148 - t125 * t155 + t128 * t157 - t139 * t96 + t140 * t98 + Icges(1,6)) * V_base(5)) * V_base(6) / 0.2e1 + V_base(5) * t141 * (-Icges(3,5) * t146 - Icges(4,5) * t139 - Icges(3,6) * t148 - Icges(4,6) * t140) + V_base(6) * t141 * (Icges(3,5) * t148 + Icges(4,5) * t140 - Icges(3,6) * t146 - Icges(4,6) * t139) + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + Icges(1,1) * V_base(4) / 0.2e1) * V_base(4) + ((-Icges(2,5) * t155 - Icges(2,6) * t157) * V_base(5) + (Icges(2,5) * t157 - Icges(2,6) * t155) * V_base(6) + Icges(2,3) * t144 / 0.2e1) * t144;
T = t1;
