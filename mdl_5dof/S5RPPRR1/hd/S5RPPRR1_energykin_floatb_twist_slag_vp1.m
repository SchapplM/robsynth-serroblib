% Calculate kinetic energy for
% S5RPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR1_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPRR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR1_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR1_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR1_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:37:56
% EndTime: 2019-12-05 17:37:58
% DurationCPUTime: 1.53s
% Computational Cost: add. (569->204), mult. (774->268), div. (0->0), fcn. (556->6), ass. (0->106)
t203 = Icges(2,4) + Icges(3,6) - Icges(4,6);
t202 = Icges(2,1) + Icges(3,2) + Icges(4,3);
t201 = -Icges(3,4) + Icges(2,5) + Icges(4,5);
t200 = Icges(4,4) + Icges(3,5) - Icges(2,6);
t199 = Icges(2,2) + Icges(4,2) + Icges(3,3);
t143 = sin(qJ(1));
t198 = t203 * t143;
t145 = cos(qJ(1));
t197 = t203 * t145;
t196 = t199 * t145 + t198;
t195 = -t199 * t143 + t197;
t194 = t202 * t143 + t197;
t193 = t202 * t145 - t198;
t190 = -pkin(2) - pkin(5);
t142 = sin(qJ(4));
t188 = pkin(4) * t142;
t144 = cos(qJ(4));
t187 = pkin(4) * t144;
t186 = t143 * pkin(6);
t185 = t145 * pkin(6);
t182 = Icges(5,4) * t142;
t181 = Icges(5,4) * t144;
t141 = qJ(4) + qJ(5);
t134 = sin(t141);
t180 = Icges(6,4) * t134;
t135 = cos(t141);
t179 = Icges(6,4) * t135;
t175 = qJ(3) * t143;
t174 = qJ(3) * t145;
t114 = pkin(1) * t143 - qJ(2) * t145;
t173 = V_base(4) * t114 + V_base(3);
t172 = V_base(5) * pkin(5) + V_base(1);
t124 = qJD(4) * t145 + V_base(5);
t169 = V_base(4) * t175 + t173;
t168 = qJD(2) * t143 + t172;
t167 = -t114 - t175;
t119 = pkin(1) * t145 + qJ(2) * t143;
t166 = -t119 - t174;
t165 = rSges(5,1) * t142 + rSges(5,2) * t144;
t164 = rSges(6,1) * t134 + rSges(6,2) * t135;
t163 = Icges(5,1) * t142 + t181;
t162 = Icges(6,1) * t134 + t179;
t161 = Icges(5,2) * t144 + t182;
t160 = Icges(6,2) * t135 + t180;
t159 = Icges(5,5) * t142 + Icges(5,6) * t144;
t158 = Icges(6,5) * t134 + Icges(6,6) * t135;
t129 = V_base(6) + qJD(1);
t157 = -qJD(2) * t145 + t129 * t119 + V_base(2);
t156 = V_base(5) * pkin(2) + qJD(3) * t145 + t168;
t155 = t167 - t185;
t154 = V_base(5) * pkin(3) + t156;
t87 = qJD(5) * t145 + t124;
t88 = V_base(4) + (-qJD(4) - qJD(5)) * t143;
t153 = (Icges(6,5) * t135 - Icges(6,6) * t134) * t129 + (Icges(6,3) * t145 + t158 * t143) * t87 + (-Icges(6,3) * t143 + t158 * t145) * t88;
t152 = qJD(3) * t143 + t129 * t174 + t157;
t125 = -qJD(4) * t143 + V_base(4);
t151 = (Icges(5,5) * t144 - Icges(5,6) * t142) * t129 + (Icges(5,3) * t145 + t159 * t143) * t124 + (-Icges(5,3) * t143 + t159 * t145) * t125;
t150 = V_base(4) * t185 + t169 + (t166 + t186) * V_base(5);
t149 = (-pkin(3) + t190) * V_base(4) + t152;
t66 = Icges(6,6) * t145 + t160 * t143;
t67 = -Icges(6,6) * t143 + t160 * t145;
t68 = Icges(6,5) * t145 + t162 * t143;
t69 = -Icges(6,5) * t143 + t162 * t145;
t84 = -Icges(6,2) * t134 + t179;
t85 = Icges(6,1) * t135 - t180;
t148 = (t134 * t69 + t135 * t67) * t88 + (t134 * t68 + t135 * t66) * t87 + (t134 * t85 + t135 * t84) * t129;
t104 = -Icges(5,2) * t142 + t181;
t111 = Icges(5,1) * t144 - t182;
t74 = Icges(5,6) * t145 + t161 * t143;
t75 = -Icges(5,6) * t143 + t161 * t145;
t76 = Icges(5,5) * t145 + t163 * t143;
t77 = -Icges(5,5) * t143 + t163 * t145;
t147 = (t142 * t77 + t144 * t75) * t125 + (t142 * t76 + t144 * t74) * t124 + (t104 * t144 + t111 * t142) * t129;
t122 = rSges(2,1) * t145 - rSges(2,2) * t143;
t121 = -rSges(3,2) * t145 + rSges(3,3) * t143;
t120 = -rSges(4,2) * t145 + rSges(4,3) * t143;
t118 = rSges(5,1) * t144 - rSges(5,2) * t142;
t117 = rSges(2,1) * t143 + rSges(2,2) * t145;
t116 = -rSges(3,2) * t143 - rSges(3,3) * t145;
t115 = rSges(4,2) * t143 + rSges(4,3) * t145;
t92 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t91 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t90 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t86 = rSges(6,1) * t135 - rSges(6,2) * t134;
t81 = pkin(7) * t145 + t143 * t188;
t80 = -pkin(7) * t143 + t145 * t188;
t79 = -rSges(5,3) * t143 + t165 * t145;
t78 = rSges(5,3) * t145 + t165 * t143;
t71 = -rSges(6,3) * t143 + t164 * t145;
t70 = rSges(6,3) * t145 + t164 * t143;
t63 = V_base(5) * rSges(2,3) - t117 * t129 + t172;
t62 = t122 * t129 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t61 = t117 * V_base(4) - t122 * V_base(5) + V_base(3);
t60 = V_base(5) * rSges(3,1) + (-t114 - t116) * t129 + t168;
t59 = t121 * t129 + (-rSges(3,1) - pkin(5)) * V_base(4) + t157;
t58 = t116 * V_base(4) + (-t119 - t121) * V_base(5) + t173;
t57 = V_base(5) * rSges(4,1) + (-t120 + t167) * t129 + t156;
t56 = t115 * t129 + (-rSges(4,1) + t190) * V_base(4) + t152;
t55 = t120 * V_base(4) + (-t115 + t166) * V_base(5) + t169;
t54 = t118 * t124 + (t155 - t78) * t129 + t154;
t53 = -t118 * t125 + (t79 - t186) * t129 + t149;
t52 = -t124 * t79 + t125 * t78 + t150;
t51 = t124 * t187 + t86 * t87 + (t155 - t70 - t81) * t129 + t154;
t50 = -t125 * t187 - t86 * t88 + (t71 + t80 - t186) * t129 + t149;
t49 = -t124 * t80 + t125 * t81 + t70 * t88 - t71 * t87 + t150;
t1 = m(1) * (t90 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(2) * (t61 ^ 2 + t62 ^ 2 + t63 ^ 2) / 0.2e1 + m(3) * (t58 ^ 2 + t59 ^ 2 + t60 ^ 2) / 0.2e1 + m(4) * (t55 ^ 2 + t56 ^ 2 + t57 ^ 2) / 0.2e1 + m(5) * (t52 ^ 2 + t53 ^ 2 + t54 ^ 2) / 0.2e1 + t125 * (-t151 * t143 + t147 * t145) / 0.2e1 + t124 * (t147 * t143 + t151 * t145) / 0.2e1 + m(6) * (t49 ^ 2 + t50 ^ 2 + t51 ^ 2) / 0.2e1 + t88 * (-t153 * t143 + t148 * t145) / 0.2e1 + t87 * (t148 * t143 + t153 * t145) / 0.2e1 + ((-t196 * t143 + t194 * t145 + Icges(1,4)) * V_base(5) + (-t195 * t143 + t193 * t145 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t194 * t143 + t196 * t145 + Icges(1,2)) * V_base(5) + (t193 * t143 + t195 * t145 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t142 * t75 + t144 * t77) * t125 + (-t142 * t74 + t144 * t76) * t124 + (-t134 * t67 + t135 * t69) * t88 + (-t134 * t66 + t135 * t68) * t87 + (-t104 * t142 + t111 * t144 - t134 * t84 + t135 * t85 + Icges(3,1) + Icges(4,1) + Icges(2,3)) * t129) * t129 / 0.2e1 + t129 * V_base(5) * (t201 * t143 - t200 * t145) + t129 * V_base(4) * (t200 * t143 + t201 * t145) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
