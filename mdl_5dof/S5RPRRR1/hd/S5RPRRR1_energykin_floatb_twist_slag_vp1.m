% Calculate kinetic energy for
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
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
% Datum: 2019-12-05 18:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(1,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR1_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRRR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_energykin_floatb_twist_slag_vp1: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR1_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR1_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR1_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:08:33
% EndTime: 2019-12-05 18:08:36
% DurationCPUTime: 2.52s
% Computational Cost: add. (768->247), mult. (1582->375), div. (0->0), fcn. (1646->8), ass. (0->113)
t187 = Icges(2,4) - Icges(3,5);
t186 = Icges(2,1) + Icges(3,1);
t185 = Icges(3,4) + Icges(2,5);
t184 = Icges(2,2) + Icges(3,3);
t183 = Icges(2,6) - Icges(3,6);
t146 = sin(qJ(1));
t182 = t187 * t146;
t150 = cos(qJ(1));
t181 = t187 * t150;
t180 = -t184 * t150 - t182;
t179 = t184 * t146 - t181;
t176 = t186 * t146 + t181;
t175 = t186 * t150 - t182;
t145 = sin(qJ(3));
t172 = Icges(4,4) * t145;
t149 = cos(qJ(3));
t171 = Icges(4,4) * t149;
t138 = V_base(6) + qJD(1);
t169 = qJ(2) * t138;
t144 = sin(qJ(4));
t168 = t144 * t145;
t167 = t145 * t146;
t148 = cos(qJ(4));
t166 = t145 * t148;
t165 = t145 * t150;
t164 = t146 * t149;
t163 = t149 * t150;
t162 = qJD(4) * t145;
t134 = qJD(3) * t146 + V_base(4);
t159 = qJD(2) * t146 + t150 * t169 + V_base(1);
t100 = t150 * t162 + t134;
t133 = -qJD(3) * t150 + V_base(5);
t158 = rSges(4,1) * t149 - rSges(4,2) * t145;
t157 = Icges(4,1) * t149 - t172;
t156 = -Icges(4,2) * t145 + t171;
t155 = Icges(4,5) * t149 - Icges(4,6) * t145;
t154 = -qJD(2) * t150 + t146 * t169 + V_base(2);
t99 = t146 * t162 + t133;
t125 = -qJD(4) * t149 + t138;
t153 = (Icges(4,5) * t145 + Icges(4,6) * t149) * t138 + (-Icges(4,3) * t150 + t146 * t155) * t133 + (Icges(4,3) * t146 + t150 * t155) * t134;
t152 = V_base(3) + (-t146 * V_base(5) - t150 * V_base(4)) * qJ(2);
t115 = Icges(4,2) * t149 + t172;
t120 = Icges(4,1) * t145 + t171;
t90 = -Icges(4,6) * t150 + t146 * t156;
t91 = Icges(4,6) * t146 + t150 * t156;
t93 = -Icges(4,5) * t150 + t146 * t157;
t94 = Icges(4,5) * t146 + t150 * t157;
t151 = (-t145 * t91 + t149 * t94) * t134 + (-t145 * t90 + t149 * t93) * t133 + (-t115 * t145 + t120 * t149) * t138;
t147 = cos(qJ(5));
t143 = sin(qJ(5));
t130 = rSges(2,1) * t150 - t146 * rSges(2,2);
t129 = rSges(3,1) * t150 + t146 * rSges(3,3);
t128 = t146 * rSges(2,1) + rSges(2,2) * t150;
t127 = t146 * rSges(3,1) - rSges(3,3) * t150;
t126 = rSges(4,1) * t145 + rSges(4,2) * t149;
t109 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t108 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t107 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t106 = t146 * t144 + t148 * t163;
t105 = t144 * t163 - t146 * t148;
t104 = -t144 * t150 + t148 * t164;
t103 = t144 * t164 + t148 * t150;
t102 = -t143 * t149 + t147 * t166;
t101 = -t143 * t166 - t147 * t149;
t98 = qJD(5) * t168 + t125;
t97 = t146 * rSges(4,3) + t150 * t158;
t96 = -rSges(4,3) * t150 + t146 * t158;
t95 = -rSges(5,3) * t149 + (rSges(5,1) * t148 - rSges(5,2) * t144) * t145;
t92 = -Icges(5,5) * t149 + (Icges(5,1) * t148 - Icges(5,4) * t144) * t145;
t89 = -Icges(5,6) * t149 + (Icges(5,4) * t148 - Icges(5,2) * t144) * t145;
t86 = -Icges(5,3) * t149 + (Icges(5,5) * t148 - Icges(5,6) * t144) * t145;
t85 = V_base(5) * rSges(2,3) - t128 * t138 + V_base(1);
t84 = -V_base(4) * rSges(2,3) + t130 * t138 + V_base(2);
t83 = t106 * t147 + t143 * t165;
t82 = -t106 * t143 + t147 * t165;
t81 = t104 * t147 + t143 * t167;
t80 = -t104 * t143 + t147 * t167;
t79 = t128 * V_base(4) - t130 * V_base(5) + V_base(3);
t78 = qJD(5) * t105 + t100;
t77 = qJD(5) * t103 + t99;
t76 = V_base(5) * rSges(3,2) - t127 * t138 + t159;
t75 = -V_base(4) * rSges(3,2) + t138 * t129 + t154;
t74 = t106 * rSges(5,1) - t105 * rSges(5,2) + rSges(5,3) * t165;
t73 = rSges(5,1) * t104 - rSges(5,2) * t103 + rSges(5,3) * t167;
t72 = rSges(6,1) * t102 + rSges(6,2) * t101 + rSges(6,3) * t168;
t71 = Icges(5,1) * t106 - Icges(5,4) * t105 + Icges(5,5) * t165;
t70 = Icges(5,1) * t104 - Icges(5,4) * t103 + Icges(5,5) * t167;
t69 = Icges(6,1) * t102 + Icges(6,4) * t101 + Icges(6,5) * t168;
t68 = Icges(5,4) * t106 - Icges(5,2) * t105 + Icges(5,6) * t165;
t67 = Icges(5,4) * t104 - Icges(5,2) * t103 + Icges(5,6) * t167;
t66 = Icges(6,4) * t102 + Icges(6,2) * t101 + Icges(6,6) * t168;
t65 = Icges(5,5) * t106 - Icges(5,6) * t105 + Icges(5,3) * t165;
t64 = Icges(5,5) * t104 - Icges(5,6) * t103 + Icges(5,3) * t167;
t63 = Icges(6,5) * t102 + Icges(6,6) * t101 + Icges(6,3) * t168;
t62 = V_base(4) * t127 - V_base(5) * t129 + t152;
t61 = t126 * t133 - t138 * t96 + t159;
t60 = -t134 * t126 + t138 * t97 + t154;
t59 = -t133 * t97 + t134 * t96 + t152;
t58 = rSges(6,1) * t83 + rSges(6,2) * t82 + rSges(6,3) * t105;
t57 = rSges(6,1) * t81 + rSges(6,2) * t80 + rSges(6,3) * t103;
t56 = Icges(6,1) * t83 + Icges(6,4) * t82 + Icges(6,5) * t105;
t55 = Icges(6,1) * t81 + Icges(6,4) * t80 + Icges(6,5) * t103;
t54 = Icges(6,4) * t83 + Icges(6,2) * t82 + Icges(6,6) * t105;
t53 = Icges(6,4) * t81 + Icges(6,2) * t80 + Icges(6,6) * t103;
t52 = Icges(6,5) * t83 + Icges(6,6) * t82 + Icges(6,3) * t105;
t51 = Icges(6,5) * t81 + Icges(6,6) * t80 + Icges(6,3) * t103;
t50 = -t125 * t73 + t95 * t99 + t159;
t49 = -t100 * t95 + t125 * t74 + t154;
t48 = t100 * t73 - t99 * t74 + t152;
t47 = -t57 * t98 + t72 * t77 + t159;
t46 = t98 * t58 - t78 * t72 + t154;
t45 = t78 * t57 - t77 * t58 + t152;
t1 = m(1) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + m(2) * (t79 ^ 2 + t84 ^ 2 + t85 ^ 2) / 0.2e1 + m(3) * (t62 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + m(4) * (t59 ^ 2 + t60 ^ 2 + t61 ^ 2) / 0.2e1 + t134 * (t153 * t146 + t151 * t150) / 0.2e1 + t133 * (t151 * t146 - t153 * t150) / 0.2e1 + m(5) * (t48 ^ 2 + t49 ^ 2 + t50 ^ 2) / 0.2e1 + t100 * ((-t105 * t68 + t106 * t71 + t65 * t165) * t100 + (-t105 * t67 + t106 * t70 + t165 * t64) * t99 + (-t105 * t89 + t106 * t92 + t165 * t86) * t125) / 0.2e1 + t99 * ((-t103 * t68 + t104 * t71 + t167 * t65) * t100 + (-t103 * t67 + t104 * t70 + t64 * t167) * t99 + (-t103 * t89 + t104 * t92 + t167 * t86) * t125) / 0.2e1 + t125 * ((-t65 * t100 - t86 * t125 - t64 * t99) * t149 + ((-t144 * t68 + t148 * t71) * t100 + (-t144 * t67 + t148 * t70) * t99 + (-t144 * t89 + t148 * t92) * t125) * t145) / 0.2e1 + m(6) * (t45 ^ 2 + t46 ^ 2 + t47 ^ 2) / 0.2e1 + t78 * ((t105 * t52 + t54 * t82 + t56 * t83) * t78 + (t105 * t51 + t53 * t82 + t55 * t83) * t77 + (t105 * t63 + t66 * t82 + t69 * t83) * t98) / 0.2e1 + t77 * ((t103 * t52 + t54 * t80 + t56 * t81) * t78 + (t103 * t51 + t53 * t80 + t55 * t81) * t77 + (t103 * t63 + t66 * t80 + t69 * t81) * t98) / 0.2e1 + t98 * ((t101 * t54 + t102 * t56 + t168 * t52) * t78 + (t101 * t53 + t102 * t55 + t168 * t51) * t77 + (t101 * t66 + t102 * t69 + t63 * t168) * t98) / 0.2e1 + ((t145 * t94 + t149 * t91) * t134 + (t145 * t93 + t149 * t90) * t133 + (t115 * t149 + t120 * t145 + Icges(3,2) + Icges(2,3)) * t138) * t138 / 0.2e1 + ((t146 * t180 + t176 * t150 + Icges(1,4)) * V_base(5) + (t179 * t146 + t175 * t150 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t176 * t146 - t150 * t180 + Icges(1,2)) * V_base(5) + (t146 * t175 - t150 * t179 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t138 * (t185 * t146 + t183 * t150) + V_base(4) * t138 * (-t183 * t146 + t185 * t150) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
