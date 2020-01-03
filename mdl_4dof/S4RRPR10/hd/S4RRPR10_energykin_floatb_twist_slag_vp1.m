% Calculate kinetic energy for
% S4RRPR10
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
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRPR10_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR10_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RRPR10_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR10_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR10_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR10_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:11:06
% EndTime: 2019-12-31 17:11:08
% DurationCPUTime: 2.05s
% Computational Cost: add. (553->206), mult. (1036->296), div. (0->0), fcn. (924->6), ass. (0->103)
t205 = Icges(3,4) + Icges(4,6);
t204 = Icges(3,1) + Icges(4,2);
t203 = -Icges(3,2) - Icges(4,3);
t146 = cos(qJ(2));
t202 = t205 * t146;
t143 = sin(qJ(2));
t201 = t205 * t143;
t200 = Icges(4,4) - Icges(3,5);
t199 = Icges(4,5) - Icges(3,6);
t198 = t203 * t143 + t202;
t197 = t204 * t146 - t201;
t144 = sin(qJ(1));
t147 = cos(qJ(1));
t196 = t198 * t144 + t199 * t147;
t195 = -t199 * t144 + t198 * t147;
t194 = t197 * t144 + t200 * t147;
t193 = -t200 * t144 + t197 * t147;
t192 = Icges(4,1) + Icges(3,3);
t191 = t203 * t146 - t201;
t190 = t204 * t143 + t202;
t189 = t199 * t143 - t200 * t146;
t134 = -qJD(2) * t147 + V_base(5);
t135 = qJD(2) * t144 + V_base(4);
t138 = V_base(6) + qJD(1);
t188 = (t191 * t143 + t190 * t146) * t138 + (-t195 * t143 + t193 * t146) * t135 + (-t196 * t143 + t194 * t146) * t134;
t187 = (-t200 * t143 - t199 * t146) * t138 + (t192 * t144 + t189 * t147) * t135 + (t189 * t144 - t192 * t147) * t134;
t183 = pkin(6) * t143;
t182 = Icges(2,4) * t144;
t142 = sin(qJ(4));
t177 = t142 * t147;
t176 = t144 * t142;
t145 = cos(qJ(4));
t175 = t144 * t145;
t174 = t144 * t146;
t173 = t145 * t147;
t172 = t146 * t147;
t162 = pkin(2) * t146 + qJ(3) * t143;
t105 = t162 * t144;
t132 = t144 * pkin(1) - pkin(5) * t147;
t171 = -t105 - t132;
t170 = qJD(3) * t143;
t169 = qJD(4) * t146;
t168 = V_base(5) * pkin(4) + V_base(1);
t127 = pkin(2) * t143 - qJ(3) * t146;
t165 = t134 * t127 + t147 * t170 + t168;
t164 = rSges(3,1) * t146 - rSges(3,2) * t143;
t163 = -rSges(4,2) * t146 + rSges(4,3) * t143;
t133 = pkin(1) * t147 + t144 * pkin(5);
t155 = -V_base(4) * pkin(4) + t138 * t133 + V_base(2);
t154 = V_base(4) * t132 - t133 * V_base(5) + V_base(3);
t106 = t162 * t147;
t151 = t138 * t106 + t144 * t170 + t155;
t150 = -qJD(3) * t146 + t135 * t105 + t154;
t140 = Icges(2,4) * t147;
t131 = rSges(2,1) * t147 - t144 * rSges(2,2);
t130 = t144 * rSges(2,1) + rSges(2,2) * t147;
t129 = rSges(3,1) * t143 + rSges(3,2) * t146;
t128 = -rSges(4,2) * t143 - rSges(4,3) * t146;
t126 = qJD(4) * t143 + t138;
t125 = Icges(2,1) * t147 - t182;
t124 = Icges(2,1) * t144 + t140;
t122 = -Icges(2,2) * t144 + t140;
t121 = Icges(2,2) * t147 + t182;
t113 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t112 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t111 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t110 = -pkin(3) * t147 + pkin(6) * t174;
t109 = t144 * pkin(3) + pkin(6) * t172;
t104 = t143 * t176 - t173;
t103 = t143 * t175 + t177;
t102 = t143 * t177 + t175;
t101 = t143 * t173 - t176;
t100 = t147 * t169 + t135;
t99 = t144 * t169 + t134;
t97 = -rSges(4,1) * t147 + t144 * t163;
t96 = t144 * rSges(4,1) + t147 * t163;
t95 = t144 * rSges(3,3) + t147 * t164;
t94 = rSges(5,3) * t143 + (-rSges(5,1) * t142 - rSges(5,2) * t145) * t146;
t93 = -rSges(3,3) * t147 + t144 * t164;
t84 = Icges(5,5) * t143 + (-Icges(5,1) * t142 - Icges(5,4) * t145) * t146;
t81 = Icges(5,6) * t143 + (-Icges(5,4) * t142 - Icges(5,2) * t145) * t146;
t78 = Icges(5,3) * t143 + (-Icges(5,5) * t142 - Icges(5,6) * t145) * t146;
t75 = V_base(5) * rSges(2,3) - t130 * t138 + t168;
t74 = t131 * t138 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t73 = t130 * V_base(4) - t131 * V_base(5) + V_base(3);
t72 = rSges(5,1) * t104 + rSges(5,2) * t103 + rSges(5,3) * t174;
t71 = t102 * rSges(5,1) + t101 * rSges(5,2) + rSges(5,3) * t172;
t70 = Icges(5,1) * t104 + Icges(5,4) * t103 + Icges(5,5) * t174;
t69 = Icges(5,1) * t102 + Icges(5,4) * t101 + Icges(5,5) * t172;
t68 = Icges(5,4) * t104 + Icges(5,2) * t103 + Icges(5,6) * t174;
t67 = Icges(5,4) * t102 + Icges(5,2) * t101 + Icges(5,6) * t172;
t66 = Icges(5,5) * t104 + Icges(5,6) * t103 + Icges(5,3) * t174;
t65 = Icges(5,5) * t102 + Icges(5,6) * t101 + Icges(5,3) * t172;
t64 = t129 * t134 + (-t132 - t93) * t138 + t168;
t63 = -t129 * t135 + t138 * t95 + t155;
t62 = -t134 * t95 + t135 * t93 + t154;
t61 = t128 * t134 + (-t97 + t171) * t138 + t165;
t60 = t138 * t96 + (-t127 - t128) * t135 + t151;
t59 = t135 * t97 + (-t106 - t96) * t134 + t150;
t58 = t134 * t183 - t126 * t72 + t94 * t99 + (-t110 + t171) * t138 + t165;
t57 = -t100 * t94 + t109 * t138 + t126 * t71 + (-t127 - t183) * t135 + t151;
t56 = t100 * t72 + t110 * t135 - t71 * t99 + (-t106 - t109) * t134 + t150;
t1 = m(1) * (t111 ^ 2 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + m(2) * (t73 ^ 2 + t74 ^ 2 + t75 ^ 2) / 0.2e1 + m(3) * (t62 ^ 2 + t63 ^ 2 + t64 ^ 2) / 0.2e1 + m(4) * (t59 ^ 2 + t60 ^ 2 + t61 ^ 2) / 0.2e1 + m(5) * (t56 ^ 2 + t57 ^ 2 + t58 ^ 2) / 0.2e1 + t100 * ((t101 * t67 + t102 * t69 + t65 * t172) * t100 + (t101 * t68 + t102 * t70 + t172 * t66) * t99 + (t101 * t81 + t102 * t84 + t172 * t78) * t126) / 0.2e1 + t99 * ((t103 * t67 + t104 * t69 + t174 * t65) * t100 + (t103 * t68 + t104 * t70 + t66 * t174) * t99 + (t103 * t81 + t104 * t84 + t174 * t78) * t126) / 0.2e1 + t126 * ((t65 * t100 + t78 * t126 + t66 * t99) * t143 + ((-t142 * t69 - t145 * t67) * t100 + (-t142 * t70 - t145 * t68) * t99 + (-t142 * t84 - t145 * t81) * t126) * t146) / 0.2e1 + (t188 * t144 - t187 * t147) * t134 / 0.2e1 + (t187 * t144 + t188 * t147) * t135 / 0.2e1 + ((-t144 * t121 + t124 * t147 + Icges(1,4)) * V_base(5) + (-t144 * t122 + t147 * t125 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t147 * t121 + t144 * t124 + Icges(1,2)) * V_base(5) + (t122 * t147 + t144 * t125 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t193 * t143 + t195 * t146) * t135 + (t194 * t143 + t196 * t146) * t134 + (t190 * t143 - t191 * t146 + Icges(2,3)) * t138) * t138 / 0.2e1 + t138 * V_base(4) * (Icges(2,5) * t147 - Icges(2,6) * t144) + V_base(5) * t138 * (Icges(2,5) * t144 + Icges(2,6) * t147) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
