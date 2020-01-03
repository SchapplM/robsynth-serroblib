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
% Datum: 2020-01-03 12:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:01:46
% EndTime: 2020-01-03 12:01:47
% DurationCPUTime: 1.27s
% Computational Cost: add. (1186->230), mult. (760->309), div. (0->0), fcn. (538->10), ass. (0->119)
t153 = qJ(4) + qJ(5);
t147 = cos(t153);
t145 = sin(t153);
t191 = Icges(6,4) * t145;
t107 = Icges(6,2) * t147 + t191;
t190 = Icges(6,4) * t147;
t110 = Icges(6,1) * t145 + t190;
t144 = V_base(4) + qJD(1);
t141 = qJD(2) + t144;
t154 = qJ(1) + qJ(2);
t143 = pkin(9) + t154;
t139 = sin(t143);
t140 = cos(t143);
t169 = -Icges(6,2) * t145 + t190;
t73 = -Icges(6,6) * t140 + t139 * t169;
t74 = -Icges(6,6) * t139 - t140 * t169;
t171 = Icges(6,1) * t147 - t191;
t75 = -Icges(6,5) * t140 + t139 * t171;
t76 = -Icges(6,5) * t139 - t140 * t171;
t189 = -qJD(4) - qJD(5);
t91 = t139 * t189 + V_base(6);
t92 = t140 * t189 + V_base(5);
t207 = (t107 * t145 - t110 * t147) * t141 + (t145 * t73 - t147 * t75) * t92 + (t145 * t74 - t147 * t76) * t91;
t118 = -qJD(4) * t139 + V_base(6);
t119 = -qJD(4) * t140 + V_base(5);
t157 = cos(qJ(4));
t155 = sin(qJ(4));
t193 = Icges(5,4) * t155;
t126 = Icges(5,2) * t157 + t193;
t192 = Icges(5,4) * t157;
t129 = Icges(5,1) * t155 + t192;
t170 = -Icges(5,2) * t155 + t192;
t81 = -Icges(5,6) * t140 + t139 * t170;
t82 = -Icges(5,6) * t139 - t140 * t170;
t172 = Icges(5,1) * t157 - t193;
t83 = -Icges(5,5) * t140 + t139 * t172;
t84 = -Icges(5,5) * t139 - t140 * t172;
t206 = (t126 * t155 - t129 * t157) * t141 + (t155 * t81 - t157 * t83) * t119 + (t155 * t82 - t157 * t84) * t118;
t205 = -pkin(5) - pkin(6);
t156 = sin(qJ(1));
t203 = pkin(1) * t156;
t158 = cos(qJ(1));
t202 = pkin(1) * t158;
t146 = sin(t154);
t201 = pkin(2) * t146;
t148 = cos(t154);
t200 = pkin(2) * t148;
t199 = pkin(4) * t155;
t198 = pkin(4) * t157;
t196 = Icges(2,4) * t158;
t195 = Icges(3,4) * t148;
t194 = Icges(4,4) * t140;
t188 = -qJ(3) + t205;
t187 = t144 * t203 + V_base(3);
t186 = V_base(6) * pkin(5) + V_base(2);
t183 = qJD(3) + V_base(1);
t182 = t141 * t201 + t187;
t181 = V_base(6) * pkin(6) + t144 * t202 + t186;
t180 = -t201 - t203;
t179 = -t200 - t202;
t178 = rSges(5,1) * t157 - rSges(5,2) * t155;
t177 = rSges(6,1) * t147 - rSges(6,2) * t145;
t168 = Icges(5,5) * t157 - Icges(5,6) * t155;
t167 = Icges(6,5) * t147 - Icges(6,6) * t145;
t164 = V_base(6) * qJ(3) + t141 * t200 + t181;
t163 = -(Icges(6,5) * t145 + Icges(6,6) * t147) * t141 - (-Icges(6,3) * t140 + t139 * t167) * t92 - (-Icges(6,3) * t139 - t140 * t167) * t91;
t162 = -t118 * (-Icges(5,3) * t139 - t140 * t168) - t119 * (-Icges(5,3) * t140 + t139 * t168) - (Icges(5,5) * t155 + Icges(5,6) * t157) * t141;
t102 = t139 * pkin(3) - t140 * pkin(7);
t161 = t141 * t102 + t188 * V_base(5) + t182;
t103 = -t140 * pkin(3) - t139 * pkin(7);
t160 = (-t102 + t180) * V_base(6) + t183 + (t103 + t179) * V_base(5);
t150 = Icges(2,4) * t156;
t138 = Icges(3,4) * t146;
t137 = Icges(4,4) * t139;
t134 = -rSges(2,1) * t158 + rSges(2,2) * t156;
t133 = rSges(2,1) * t156 + rSges(2,2) * t158;
t132 = rSges(5,1) * t155 + rSges(5,2) * t157;
t131 = -Icges(2,1) * t158 + t150;
t130 = Icges(2,1) * t156 + t196;
t128 = Icges(2,2) * t156 - t196;
t127 = Icges(2,2) * t158 + t150;
t122 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t121 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t120 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t115 = -rSges(3,1) * t148 + rSges(3,2) * t146;
t114 = rSges(3,1) * t146 + rSges(3,2) * t148;
t113 = rSges(6,1) * t145 + rSges(6,2) * t147;
t112 = -Icges(3,1) * t148 + t138;
t111 = Icges(3,1) * t146 + t195;
t109 = Icges(3,2) * t146 - t195;
t108 = Icges(3,2) * t148 + t138;
t101 = -rSges(4,1) * t140 + rSges(4,2) * t139;
t100 = rSges(4,1) * t139 + rSges(4,2) * t140;
t99 = -Icges(4,1) * t140 + t137;
t98 = Icges(4,1) * t139 + t194;
t97 = Icges(4,2) * t139 - t194;
t96 = Icges(4,2) * t140 + t137;
t89 = V_base(6) * rSges(2,3) - t134 * t144 + t186;
t88 = t133 * t144 + V_base(3) + (-rSges(2,3) - pkin(5)) * V_base(5);
t87 = -t133 * V_base(6) + t134 * V_base(5) + V_base(1);
t86 = -rSges(5,3) * t139 - t140 * t178;
t85 = -rSges(5,3) * t140 + t139 * t178;
t78 = -rSges(6,3) * t139 - t140 * t177;
t77 = -rSges(6,3) * t140 + t139 * t177;
t70 = -pkin(8) * t139 - t140 * t198;
t69 = -pkin(8) * t140 + t139 * t198;
t68 = V_base(6) * rSges(3,3) - t115 * t141 + t181;
t67 = t114 * t141 + (-rSges(3,3) + t205) * V_base(5) + t187;
t66 = -t114 * V_base(6) + t115 * V_base(5) + V_base(1) + (-V_base(6) * t156 - t158 * V_base(5)) * pkin(1);
t65 = V_base(6) * rSges(4,3) - t101 * t141 + t164;
t64 = t100 * t141 + (-rSges(4,3) + t188) * V_base(5) + t182;
t63 = (-t100 + t180) * V_base(6) + (t101 + t179) * V_base(5) + t183;
t62 = t118 * t132 + (-t103 - t86) * t141 + t164;
t61 = -t119 * t132 + t141 * t85 + t161;
t60 = -t118 * t85 + t119 * t86 + t160;
t59 = t118 * t199 + t113 * t91 + (-t103 - t70 - t78) * t141 + t164;
t58 = -t119 * t199 - t113 * t92 + (t69 + t77) * t141 + t161;
t57 = -t118 * t69 + t119 * t70 - t77 * t91 + t78 * t92 + t160;
t1 = m(1) * (t120 ^ 2 + t121 ^ 2 + t122 ^ 2) / 0.2e1 + m(2) * (t87 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + m(3) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + m(4) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + m(5) * (t60 ^ 2 + t61 ^ 2 + t62 ^ 2) / 0.2e1 + t119 * (-t206 * t139 + t162 * t140) / 0.2e1 + t118 * (t162 * t139 + t206 * t140) / 0.2e1 + m(6) * (t57 ^ 2 + t58 ^ 2 + t59 ^ 2) / 0.2e1 + t92 * (-t207 * t139 + t163 * t140) / 0.2e1 + t91 * (t163 * t139 + t207 * t140) / 0.2e1 + ((t155 * t83 + t157 * t81) * t119 + (t155 * t84 + t157 * t82) * t118 + (t145 * t75 + t147 * t73) * t92 + (t145 * t76 + t147 * t74) * t91 + (t107 * t147 + t110 * t145 + t126 * t157 + t129 * t155 + Icges(3,3) + Icges(4,3)) * t141) * t141 / 0.2e1 + ((t109 * t148 + t112 * t146 + t128 * t158 + t131 * t156 + t139 * t99 + t140 * t97 + Icges(1,6)) * V_base(6) + (t108 * t148 + t111 * t146 + t127 * t158 + t130 * t156 + t139 * t98 + t140 * t96 + Icges(1,2)) * V_base(5)) * V_base(5) / 0.2e1 + ((t109 * t146 - t112 * t148 + t128 * t156 - t131 * t158 + t139 * t97 - t140 * t99 + Icges(1,3)) * V_base(6) + (t108 * t146 - t111 * t148 + t127 * t156 - t130 * t158 + t139 * t96 - t140 * t98 + Icges(1,6)) * V_base(5)) * V_base(6) / 0.2e1 + V_base(5) * t141 * (Icges(3,5) * t146 + Icges(4,5) * t139 + Icges(3,6) * t148 + Icges(4,6) * t140) + V_base(6) * t141 * (-Icges(3,5) * t148 - Icges(4,5) * t140 + Icges(3,6) * t146 + Icges(4,6) * t139) + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + Icges(1,1) * V_base(4) / 0.2e1) * V_base(4) + ((Icges(2,5) * t156 + Icges(2,6) * t158) * V_base(5) + (-Icges(2,5) * t158 + Icges(2,6) * t156) * V_base(6) + Icges(2,3) * t144 / 0.2e1) * t144;
T = t1;
