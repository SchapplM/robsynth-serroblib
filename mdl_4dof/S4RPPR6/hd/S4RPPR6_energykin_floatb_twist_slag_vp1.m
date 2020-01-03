% Calculate kinetic energy for
% S4RPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPPR6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR6_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RPPR6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR6_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR6_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR6_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:35
% EndTime: 2019-12-31 16:40:37
% DurationCPUTime: 1.78s
% Computational Cost: add. (508->193), mult. (990->257), div. (0->0), fcn. (900->6), ass. (0->97)
t197 = Icges(3,4) - Icges(4,5);
t196 = Icges(3,1) + Icges(4,1);
t195 = Icges(3,2) + Icges(4,3);
t138 = sin(pkin(6));
t194 = t197 * t138;
t139 = cos(pkin(6));
t193 = t197 * t139;
t192 = Icges(4,4) + Icges(3,5);
t191 = Icges(3,6) - Icges(4,6);
t190 = t195 * t138 - t193;
t189 = t196 * t139 - t194;
t141 = sin(qJ(1));
t143 = cos(qJ(1));
t188 = t190 * t141 + t191 * t143;
t187 = -t191 * t141 + t190 * t143;
t186 = t189 * t141 - t192 * t143;
t185 = t192 * t141 + t189 * t143;
t184 = -t195 * t139 - t194;
t183 = t196 * t138 + t193;
t182 = t191 * t138 - t192 * t139;
t181 = Icges(2,2) + Icges(4,2) + Icges(3,3);
t134 = V_base(6) + qJD(1);
t136 = Icges(2,4) * t143;
t173 = Icges(2,4) * t141;
t180 = (-t192 * t138 - t191 * t139) * t134 + (t182 * t141 + t181 * t143 + t173) * V_base(5) + (-t181 * t141 + t182 * t143 + t136) * V_base(4);
t179 = (t184 * t138 + t183 * t139) * t134 + (Icges(2,1) * t141 + t188 * t138 + t186 * t139 + t136) * V_base(5) + (Icges(2,1) * t143 + t187 * t138 + t185 * t139 - t173) * V_base(4);
t177 = pkin(3) * t138;
t176 = pkin(3) * t139;
t117 = pkin(2) * t138 - qJ(3) * t139;
t175 = -pkin(4) - t117;
t126 = t141 * pkin(1) - qJ(2) * t143;
t159 = pkin(2) * t139 + qJ(3) * t138;
t99 = t159 * t141;
t174 = -t126 - t99;
t100 = t159 * t143;
t128 = pkin(1) * t143 + t141 * qJ(2);
t168 = -t100 - t128;
t167 = qJD(3) * t138;
t166 = V_base(4) * t126 + V_base(3);
t165 = V_base(5) * pkin(4) + V_base(1);
t162 = qJD(2) * t141 + t165;
t161 = rSges(3,1) * t139 - rSges(3,2) * t138;
t160 = rSges(4,1) * t139 + rSges(4,3) * t138;
t140 = sin(qJ(4));
t142 = cos(qJ(4));
t104 = t138 * t142 - t139 * t140;
t152 = t138 * t140 + t139 * t142;
t151 = -qJD(2) * t143 + t134 * t128 + V_base(2);
t150 = V_base(5) * t117 + t143 * t167 + t162;
t149 = -qJD(3) * t139 + V_base(4) * t99 + t166;
t148 = t134 * t100 + t141 * t167 + t151;
t131 = -qJD(4) * t141 + V_base(4);
t130 = qJD(4) * t143 + V_base(5);
t129 = rSges(2,1) * t143 - t141 * rSges(2,2);
t127 = t141 * rSges(2,1) + rSges(2,2) * t143;
t121 = Icges(2,5) * t143 - Icges(2,6) * t141;
t120 = Icges(2,5) * t141 + Icges(2,6) * t143;
t119 = rSges(3,1) * t138 + rSges(3,2) * t139;
t118 = rSges(4,1) * t138 - rSges(4,3) * t139;
t110 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t109 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t108 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t107 = -t141 * pkin(5) + t143 * t176;
t106 = pkin(5) * t143 + t141 * t176;
t98 = t152 * t143;
t97 = t104 * t143;
t96 = t152 * t141;
t95 = t104 * t141;
t93 = t141 * rSges(3,3) + t143 * t161;
t92 = t141 * rSges(4,2) + t143 * t160;
t91 = -rSges(3,3) * t143 + t141 * t161;
t90 = -rSges(4,2) * t143 + t141 * t160;
t76 = V_base(5) * rSges(2,3) - t127 * t134 + t165;
t75 = t129 * t134 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t74 = t127 * V_base(4) - t129 * V_base(5) + V_base(3);
t73 = rSges(5,1) * t104 - rSges(5,2) * t152;
t72 = Icges(5,1) * t104 - Icges(5,4) * t152;
t71 = Icges(5,4) * t104 - Icges(5,2) * t152;
t70 = Icges(5,5) * t104 - Icges(5,6) * t152;
t69 = rSges(5,1) * t98 + rSges(5,2) * t97 - rSges(5,3) * t141;
t68 = t96 * rSges(5,1) + t95 * rSges(5,2) + rSges(5,3) * t143;
t67 = Icges(5,1) * t98 + Icges(5,4) * t97 - Icges(5,5) * t141;
t66 = Icges(5,1) * t96 + Icges(5,4) * t95 + Icges(5,5) * t143;
t65 = Icges(5,4) * t98 + Icges(5,2) * t97 - Icges(5,6) * t141;
t64 = Icges(5,4) * t96 + Icges(5,2) * t95 + Icges(5,6) * t143;
t63 = Icges(5,5) * t98 + Icges(5,6) * t97 - Icges(5,3) * t141;
t62 = Icges(5,5) * t96 + Icges(5,6) * t95 + Icges(5,3) * t143;
t61 = t119 * V_base(5) + (-t126 - t91) * t134 + t162;
t60 = t134 * t93 + (-pkin(4) - t119) * V_base(4) + t151;
t59 = t91 * V_base(4) + (-t128 - t93) * V_base(5) + t166;
t58 = t118 * V_base(5) + (-t90 + t174) * t134 + t150;
t57 = t134 * t92 + (-t118 + t175) * V_base(4) + t148;
t56 = t90 * V_base(4) + (-t92 + t168) * V_base(5) + t149;
t55 = V_base(5) * t177 + t130 * t73 + (-t106 - t68 + t174) * t134 + t150;
t54 = -t131 * t73 + (t107 + t69) * t134 + (t175 - t177) * V_base(4) + t148;
t53 = t106 * V_base(4) - t130 * t69 + t131 * t68 + (-t107 + t168) * V_base(5) + t149;
t1 = m(1) * (t108 ^ 2 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + m(2) * (t74 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + m(3) * (t59 ^ 2 + t60 ^ 2 + t61 ^ 2) / 0.2e1 + m(4) * (t56 ^ 2 + t57 ^ 2 + t58 ^ 2) / 0.2e1 + m(5) * (t53 ^ 2 + t54 ^ 2 + t55 ^ 2) / 0.2e1 + t131 * ((-t141 * t63 + t65 * t97 + t67 * t98) * t131 + (-t141 * t62 + t64 * t97 + t66 * t98) * t130 + (-t141 * t70 + t71 * t97 + t72 * t98) * t134) / 0.2e1 + t130 * ((t143 * t63 + t95 * t65 + t96 * t67) * t131 + (t143 * t62 + t95 * t64 + t96 * t66) * t130 + (t143 * t70 + t95 * t71 + t96 * t72) * t134) / 0.2e1 + ((t104 * t67 - t152 * t65) * t131 + (t104 * t66 - t152 * t64) * t130 + (t186 * t138 - t188 * t139 + t120) * V_base(5) + (t185 * t138 - t187 * t139 + t121) * V_base(4) + (t104 * t72 + t183 * t138 - t184 * t139 - t152 * t71 + Icges(2,3)) * t134) * t134 / 0.2e1 + (Icges(1,1) * V_base(4) + t121 * t134 - t180 * t141 + t179 * t143) * V_base(4) / 0.2e1 + (Icges(1,2) * V_base(5) + t120 * t134 + t179 * t141 + t180 * t143) * V_base(5) / 0.2e1 + V_base(4) * V_base(5) * Icges(1,4) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
