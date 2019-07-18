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
% Datum: 2019-07-18 13:26
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
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
% StartTime: 2019-07-18 13:24:50
% EndTime: 2019-07-18 13:24:52
% DurationCPUTime: 2.48s
% Computational Cost: add. (768->247), mult. (1582->375), div. (0->0), fcn. (1646->8), ass. (0->113)
t188 = Icges(2,4) - Icges(3,5);
t187 = Icges(2,1) + Icges(3,1);
t186 = Icges(3,4) + Icges(2,5);
t185 = Icges(2,2) + Icges(3,3);
t184 = Icges(2,6) - Icges(3,6);
t147 = sin(qJ(1));
t183 = t188 * t147;
t151 = cos(qJ(1));
t182 = t188 * t151;
t181 = -t151 * t185 - t183;
t180 = t147 * t185 - t182;
t177 = t147 * t187 + t182;
t176 = t151 * t187 - t183;
t146 = sin(qJ(3));
t173 = Icges(4,4) * t146;
t150 = cos(qJ(3));
t172 = Icges(4,4) * t150;
t139 = V_base(6) + qJD(1);
t170 = qJ(2) * t139;
t145 = sin(qJ(4));
t169 = t145 * t146;
t168 = t146 * t147;
t149 = cos(qJ(4));
t167 = t146 * t149;
t166 = t146 * t151;
t165 = t147 * t150;
t164 = t150 * t151;
t163 = qJD(4) * t146;
t135 = qJD(3) * t147 + V_base(4);
t160 = qJD(2) * t147 + t151 * t170 + V_base(1);
t101 = t151 * t163 + t135;
t134 = -qJD(3) * t151 + V_base(5);
t159 = rSges(4,1) * t150 - rSges(4,2) * t146;
t158 = Icges(4,1) * t150 - t173;
t157 = -Icges(4,2) * t146 + t172;
t156 = Icges(4,5) * t150 - Icges(4,6) * t146;
t155 = -qJD(2) * t151 + t147 * t170 + V_base(2);
t100 = t147 * t163 + t134;
t126 = -qJD(4) * t150 + t139;
t154 = (Icges(4,5) * t146 + Icges(4,6) * t150) * t139 + (-Icges(4,3) * t151 + t147 * t156) * t134 + (Icges(4,3) * t147 + t151 * t156) * t135;
t153 = V_base(3) + (-t147 * V_base(5) - t151 * V_base(4)) * qJ(2);
t116 = Icges(4,2) * t150 + t173;
t121 = Icges(4,1) * t146 + t172;
t91 = -Icges(4,6) * t151 + t147 * t157;
t92 = Icges(4,6) * t147 + t151 * t157;
t94 = -Icges(4,5) * t151 + t147 * t158;
t95 = Icges(4,5) * t147 + t151 * t158;
t152 = (-t146 * t92 + t150 * t95) * t135 + (-t146 * t91 + t150 * t94) * t134 + (-t116 * t146 + t121 * t150) * t139;
t148 = cos(qJ(5));
t144 = sin(qJ(5));
t131 = rSges(2,1) * t151 - t147 * rSges(2,2);
t130 = rSges(3,1) * t151 + t147 * rSges(3,3);
t129 = t147 * rSges(2,1) + rSges(2,2) * t151;
t128 = t147 * rSges(3,1) - rSges(3,3) * t151;
t127 = rSges(4,1) * t146 + rSges(4,2) * t150;
t110 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t109 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t108 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t107 = t147 * t145 + t149 * t164;
t106 = t145 * t164 - t147 * t149;
t105 = -t145 * t151 + t149 * t165;
t104 = t145 * t165 + t149 * t151;
t103 = -t144 * t150 + t148 * t167;
t102 = -t144 * t167 - t148 * t150;
t99 = qJD(5) * t169 + t126;
t98 = t147 * rSges(4,3) + t151 * t159;
t97 = -rSges(4,3) * t151 + t147 * t159;
t96 = -rSges(5,3) * t150 + (rSges(5,1) * t149 - rSges(5,2) * t145) * t146;
t93 = -Icges(5,5) * t150 + (Icges(5,1) * t149 - Icges(5,4) * t145) * t146;
t90 = -Icges(5,6) * t150 + (Icges(5,4) * t149 - Icges(5,2) * t145) * t146;
t87 = -Icges(5,3) * t150 + (Icges(5,5) * t149 - Icges(5,6) * t145) * t146;
t86 = V_base(5) * rSges(2,3) - t129 * t139 + V_base(1);
t85 = -V_base(4) * rSges(2,3) + t131 * t139 + V_base(2);
t84 = t107 * t148 + t144 * t166;
t83 = -t107 * t144 + t148 * t166;
t82 = t105 * t148 + t144 * t168;
t81 = -t105 * t144 + t148 * t168;
t80 = t129 * V_base(4) - t131 * V_base(5) + V_base(3);
t79 = qJD(5) * t106 + t101;
t78 = qJD(5) * t104 + t100;
t77 = V_base(5) * rSges(3,2) - t128 * t139 + t160;
t76 = -V_base(4) * rSges(3,2) + t139 * t130 + t155;
t75 = t107 * rSges(5,1) - t106 * rSges(5,2) + rSges(5,3) * t166;
t74 = rSges(5,1) * t105 - rSges(5,2) * t104 + rSges(5,3) * t168;
t73 = rSges(6,1) * t103 + rSges(6,2) * t102 + rSges(6,3) * t169;
t72 = Icges(5,1) * t107 - Icges(5,4) * t106 + Icges(5,5) * t166;
t71 = Icges(5,1) * t105 - Icges(5,4) * t104 + Icges(5,5) * t168;
t70 = Icges(6,1) * t103 + Icges(6,4) * t102 + Icges(6,5) * t169;
t69 = Icges(5,4) * t107 - Icges(5,2) * t106 + Icges(5,6) * t166;
t68 = Icges(5,4) * t105 - Icges(5,2) * t104 + Icges(5,6) * t168;
t67 = Icges(6,4) * t103 + Icges(6,2) * t102 + Icges(6,6) * t169;
t66 = Icges(5,5) * t107 - Icges(5,6) * t106 + Icges(5,3) * t166;
t65 = Icges(5,5) * t105 - Icges(5,6) * t104 + Icges(5,3) * t168;
t64 = Icges(6,5) * t103 + Icges(6,6) * t102 + Icges(6,3) * t169;
t63 = V_base(4) * t128 - V_base(5) * t130 + t153;
t62 = t127 * t134 - t139 * t97 + t160;
t61 = -t135 * t127 + t139 * t98 + t155;
t60 = -t134 * t98 + t135 * t97 + t153;
t59 = rSges(6,1) * t84 + rSges(6,2) * t83 + rSges(6,3) * t106;
t58 = rSges(6,1) * t82 + rSges(6,2) * t81 + rSges(6,3) * t104;
t57 = Icges(6,1) * t84 + Icges(6,4) * t83 + Icges(6,5) * t106;
t56 = Icges(6,1) * t82 + Icges(6,4) * t81 + Icges(6,5) * t104;
t55 = Icges(6,4) * t84 + Icges(6,2) * t83 + Icges(6,6) * t106;
t54 = Icges(6,4) * t82 + Icges(6,2) * t81 + Icges(6,6) * t104;
t53 = Icges(6,5) * t84 + Icges(6,6) * t83 + Icges(6,3) * t106;
t52 = Icges(6,5) * t82 + Icges(6,6) * t81 + Icges(6,3) * t104;
t51 = t100 * t96 - t126 * t74 + t160;
t50 = -t101 * t96 + t126 * t75 + t155;
t49 = -t100 * t75 + t101 * t74 + t153;
t48 = -t58 * t99 + t73 * t78 + t160;
t47 = t99 * t59 - t79 * t73 + t155;
t46 = t79 * t58 - t78 * t59 + t153;
t1 = m(1) * (t108 ^ 2 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + m(2) * (t80 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + m(3) * (t63 ^ 2 + t76 ^ 2 + t77 ^ 2) / 0.2e1 + m(4) * (t60 ^ 2 + t61 ^ 2 + t62 ^ 2) / 0.2e1 + t135 * (t154 * t147 + t152 * t151) / 0.2e1 + t134 * (t152 * t147 - t154 * t151) / 0.2e1 + m(5) * (t49 ^ 2 + t50 ^ 2 + t51 ^ 2) / 0.2e1 + t101 * ((-t106 * t69 + t107 * t72 + t66 * t166) * t101 + (-t106 * t68 + t107 * t71 + t166 * t65) * t100 + (-t106 * t90 + t107 * t93 + t166 * t87) * t126) / 0.2e1 + t100 * ((-t104 * t69 + t105 * t72 + t168 * t66) * t101 + (-t104 * t68 + t105 * t71 + t65 * t168) * t100 + (-t104 * t90 + t105 * t93 + t168 * t87) * t126) / 0.2e1 + t126 * ((-t65 * t100 - t66 * t101 - t87 * t126) * t150 + ((-t145 * t69 + t149 * t72) * t101 + (-t145 * t68 + t149 * t71) * t100 + (-t145 * t90 + t149 * t93) * t126) * t146) / 0.2e1 + m(6) * (t46 ^ 2 + t47 ^ 2 + t48 ^ 2) / 0.2e1 + t79 * ((t106 * t53 + t55 * t83 + t57 * t84) * t79 + (t106 * t52 + t54 * t83 + t56 * t84) * t78 + (t106 * t64 + t67 * t83 + t70 * t84) * t99) / 0.2e1 + t78 * ((t104 * t53 + t55 * t81 + t57 * t82) * t79 + (t104 * t52 + t54 * t81 + t56 * t82) * t78 + (t104 * t64 + t67 * t81 + t70 * t82) * t99) / 0.2e1 + t99 * ((t102 * t55 + t103 * t57 + t169 * t53) * t79 + (t102 * t54 + t103 * t56 + t169 * t52) * t78 + (t102 * t67 + t103 * t70 + t169 * t64) * t99) / 0.2e1 + ((t146 * t95 + t150 * t92) * t135 + (t146 * t94 + t150 * t91) * t134 + (t116 * t150 + t121 * t146 + Icges(3,2) + Icges(2,3)) * t139) * t139 / 0.2e1 + ((t147 * t181 + t177 * t151 + Icges(1,4)) * V_base(5) + (t180 * t147 + t176 * t151 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t177 * t147 - t181 * t151 + Icges(1,2)) * V_base(5) + (t147 * t176 - t151 * t180 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t139 * (t147 * t186 + t151 * t184) + V_base(4) * t139 * (-t184 * t147 + t186 * t151) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
