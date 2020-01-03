% Calculate kinetic energy for
% S4PRPR7
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
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRPR7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR7_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PRPR7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR7_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR7_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR7_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:25:34
% EndTime: 2019-12-31 16:25:36
% DurationCPUTime: 1.96s
% Computational Cost: add. (521->206), mult. (1036->293), div. (0->0), fcn. (924->6), ass. (0->101)
t208 = Icges(3,4) + Icges(4,6);
t207 = Icges(3,1) + Icges(4,2);
t206 = -Icges(3,2) - Icges(4,3);
t149 = cos(qJ(2));
t205 = t208 * t149;
t147 = sin(qJ(2));
t204 = t208 * t147;
t203 = Icges(4,4) - Icges(3,5);
t202 = Icges(4,5) - Icges(3,6);
t201 = t206 * t147 + t205;
t200 = t207 * t149 - t204;
t144 = sin(pkin(6));
t145 = cos(pkin(6));
t199 = t201 * t144 + t202 * t145;
t198 = -t202 * t144 + t201 * t145;
t197 = t200 * t144 + t203 * t145;
t196 = -t203 * t144 + t200 * t145;
t195 = Icges(4,1) + Icges(3,3);
t194 = t206 * t149 - t204;
t193 = t207 * t147 + t205;
t192 = t202 * t147 - t203 * t149;
t136 = -qJD(2) * t145 + V_base(5);
t137 = qJD(2) * t144 + V_base(4);
t189 = (t194 * t147 + t193 * t149) * V_base(6) + (-t198 * t147 + t196 * t149) * t137 + (-t199 * t147 + t197 * t149) * t136;
t188 = (-t203 * t147 - t202 * t149) * V_base(6) + (t195 * t144 + t192 * t145) * t137 + (t192 * t144 - t195 * t145) * t136;
t185 = pkin(5) * t147;
t184 = Icges(2,4) * t144;
t179 = t144 * t149;
t178 = t145 * t149;
t146 = sin(qJ(4));
t177 = t146 * t147;
t148 = cos(qJ(4));
t176 = t147 * t148;
t164 = pkin(2) * t149 + qJ(3) * t147;
t108 = t164 * t144;
t125 = pkin(1) * t144 - pkin(4) * t145;
t175 = -t108 - t125;
t174 = qJD(3) * t147;
t173 = qJD(4) * t149;
t172 = V_base(5) * qJ(1) + V_base(1);
t168 = qJD(1) + V_base(3);
t133 = t147 * pkin(2) - qJ(3) * t149;
t167 = t136 * t133 + t145 * t174 + t172;
t166 = rSges(3,1) * t149 - rSges(3,2) * t147;
t165 = -rSges(4,2) * t149 + rSges(4,3) * t147;
t126 = pkin(1) * t145 + pkin(4) * t144;
t157 = -V_base(4) * qJ(1) + V_base(6) * t126 + V_base(2);
t156 = V_base(4) * t125 - V_base(5) * t126 + t168;
t109 = t164 * t145;
t155 = V_base(6) * t109 + t144 * t174 + t157;
t152 = -qJD(3) * t149 + t137 * t108 + t156;
t142 = Icges(2,4) * t145;
t138 = qJD(4) * t147 + V_base(6);
t135 = t147 * rSges(3,1) + rSges(3,2) * t149;
t134 = -t147 * rSges(4,2) - rSges(4,3) * t149;
t124 = rSges(2,1) * t145 - rSges(2,2) * t144;
t123 = rSges(2,1) * t144 + rSges(2,2) * t145;
t122 = Icges(2,1) * t145 - t184;
t121 = Icges(2,1) * t144 + t142;
t120 = -Icges(2,2) * t144 + t142;
t119 = Icges(2,2) * t145 + t184;
t116 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t115 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t114 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t113 = -t145 * pkin(3) + pkin(5) * t179;
t112 = t144 * pkin(3) + pkin(5) * t178;
t107 = t144 * t177 - t145 * t148;
t106 = t144 * t176 + t145 * t146;
t105 = t144 * t148 + t145 * t177;
t104 = -t144 * t146 + t145 * t176;
t103 = t145 * t173 + t137;
t102 = t144 * t173 + t136;
t99 = t147 * rSges(5,3) + (-rSges(5,1) * t146 - rSges(5,2) * t148) * t149;
t98 = Icges(5,5) * t147 + (-Icges(5,1) * t146 - Icges(5,4) * t148) * t149;
t97 = Icges(5,6) * t147 + (-Icges(5,4) * t146 - Icges(5,2) * t148) * t149;
t96 = Icges(5,3) * t147 + (-Icges(5,5) * t146 - Icges(5,6) * t148) * t149;
t95 = -t145 * rSges(4,1) + t144 * t165;
t94 = t144 * rSges(4,1) + t145 * t165;
t93 = t144 * rSges(3,3) + t145 * t166;
t92 = -t145 * rSges(3,3) + t144 * t166;
t78 = V_base(5) * rSges(2,3) - t123 * V_base(6) + t172;
t77 = t124 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t76 = t123 * V_base(4) - t124 * V_base(5) + t168;
t75 = t107 * rSges(5,1) + t106 * rSges(5,2) + rSges(5,3) * t179;
t74 = t105 * rSges(5,1) + t104 * rSges(5,2) + rSges(5,3) * t178;
t73 = Icges(5,1) * t107 + Icges(5,4) * t106 + Icges(5,5) * t179;
t72 = Icges(5,1) * t105 + Icges(5,4) * t104 + Icges(5,5) * t178;
t71 = Icges(5,4) * t107 + Icges(5,2) * t106 + Icges(5,6) * t179;
t70 = Icges(5,4) * t105 + Icges(5,2) * t104 + Icges(5,6) * t178;
t69 = Icges(5,5) * t107 + Icges(5,6) * t106 + Icges(5,3) * t179;
t68 = Icges(5,5) * t105 + Icges(5,6) * t104 + Icges(5,3) * t178;
t67 = t135 * t136 + (-t125 - t92) * V_base(6) + t172;
t66 = -t135 * t137 + t93 * V_base(6) + t157;
t65 = -t136 * t93 + t137 * t92 + t156;
t64 = t134 * t136 + (-t95 + t175) * V_base(6) + t167;
t63 = t94 * V_base(6) + (-t133 - t134) * t137 + t155;
t62 = t137 * t95 + (-t109 - t94) * t136 + t152;
t61 = t136 * t185 + t102 * t99 - t138 * t75 + (-t113 + t175) * V_base(6) + t167;
t60 = -t103 * t99 + t112 * V_base(6) + t138 * t74 + (-t133 - t185) * t137 + t155;
t59 = -t102 * t74 + t103 * t75 + t137 * t113 + (-t109 - t112) * t136 + t152;
t1 = m(1) * (t114 ^ 2 + t115 ^ 2 + t116 ^ 2) / 0.2e1 + m(2) * (t76 ^ 2 + t77 ^ 2 + t78 ^ 2) / 0.2e1 + m(3) * (t65 ^ 2 + t66 ^ 2 + t67 ^ 2) / 0.2e1 + m(4) * (t62 ^ 2 + t63 ^ 2 + t64 ^ 2) / 0.2e1 + m(5) * (t59 ^ 2 + t60 ^ 2 + t61 ^ 2) / 0.2e1 + t103 * ((t104 * t70 + t105 * t72 + t68 * t178) * t103 + (t104 * t71 + t105 * t73 + t178 * t69) * t102 + (t104 * t97 + t105 * t98 + t178 * t96) * t138) / 0.2e1 + t102 * ((t106 * t70 + t107 * t72 + t179 * t68) * t103 + (t106 * t71 + t107 * t73 + t69 * t179) * t102 + (t106 * t97 + t107 * t98 + t179 * t96) * t138) / 0.2e1 + t138 * ((t69 * t102 + t68 * t103 + t96 * t138) * t147 + ((-t146 * t72 - t148 * t70) * t103 + (-t146 * t73 - t148 * t71) * t102 + (-t146 * t98 - t148 * t97) * t138) * t149) / 0.2e1 + (t189 * t144 - t188 * t145) * t136 / 0.2e1 + (t188 * t144 + t189 * t145) * t137 / 0.2e1 + ((-t119 * t144 + t121 * t145 + Icges(1,4)) * V_base(5) + (-t120 * t144 + t122 * t145 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t119 * t145 + t121 * t144 + Icges(1,2)) * V_base(5) + (t120 * t145 + t122 * t144 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t196 * t147 + t198 * t149) * t137 + (t197 * t147 + t199 * t149) * t136 + (t193 * t147 - t194 * t149 + Icges(1,3) + Icges(2,3)) * V_base(6)) * V_base(6) / 0.2e1 + V_base(6) * V_base(4) * (Icges(2,5) * t145 - Icges(2,6) * t144 + Icges(1,5)) + V_base(6) * V_base(5) * (Icges(2,5) * t144 + Icges(2,6) * t145 + Icges(1,6));
T = t1;
