% Calculate kinetic energy for
% S4PRRP5
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
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRRP5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP5_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PRRP5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP5_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP5_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP5_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:28:48
% EndTime: 2019-12-31 16:28:50
% DurationCPUTime: 1.55s
% Computational Cost: add. (647->204), mult. (1305->294), div. (0->0), fcn. (1249->6), ass. (0->100)
t212 = Icges(4,1) + Icges(5,1);
t211 = Icges(4,4) + Icges(5,4);
t210 = -Icges(5,5) - Icges(4,5);
t209 = Icges(4,2) + Icges(5,2);
t208 = -Icges(5,6) - Icges(4,6);
t207 = -Icges(5,3) - Icges(4,3);
t152 = sin(pkin(6));
t153 = cos(pkin(6));
t157 = cos(qJ(3));
t155 = sin(qJ(3));
t158 = cos(qJ(2));
t180 = t155 * t158;
t119 = -t152 * t180 - t153 * t157;
t179 = t157 * t158;
t182 = t153 * t155;
t120 = t152 * t179 - t182;
t156 = sin(qJ(2));
t183 = t152 * t156;
t206 = -t208 * t119 - t210 * t120 - t207 * t183;
t121 = t152 * t157 - t153 * t180;
t184 = t152 * t155;
t122 = t153 * t179 + t184;
t181 = t153 * t156;
t205 = -t208 * t121 - t210 * t122 - t207 * t181;
t204 = t209 * t119 + t211 * t120 - t208 * t183;
t203 = t209 * t121 + t211 * t122 - t208 * t181;
t202 = t211 * t119 + t212 * t120 - t210 * t183;
t201 = t211 * t121 + t212 * t122 - t210 * t181;
t198 = t207 * t158 + (t208 * t155 - t210 * t157) * t156;
t197 = t208 * t158 + (-t209 * t155 + t211 * t157) * t156;
t196 = t210 * t158 + (-t211 * t155 + t212 * t157) * t156;
t192 = pkin(3) * t157;
t164 = qJ(4) * t156 + t158 * t192;
t190 = rSges(5,1) * t120 + rSges(5,2) * t119 + rSges(5,3) * t183 - pkin(3) * t182 + t164 * t152;
t189 = rSges(5,1) * t122 + rSges(5,2) * t121 + rSges(5,3) * t181 + pkin(3) * t184 + t164 * t153;
t188 = (-qJ(4) - rSges(5,3)) * t158 + (rSges(5,1) * t157 - rSges(5,2) * t155 + t192) * t156;
t187 = Icges(2,4) * t152;
t186 = Icges(3,4) * t156;
t185 = Icges(3,4) * t158;
t178 = qJD(3) * t156;
t177 = qJD(4) * t156;
t176 = V_base(5) * qJ(1) + V_base(1);
t172 = qJD(1) + V_base(3);
t146 = qJD(2) * t152 + V_base(4);
t171 = pkin(2) * t158 + pkin(5) * t156;
t145 = -qJD(2) * t153 + V_base(5);
t170 = rSges(3,1) * t158 - rSges(3,2) * t156;
t169 = Icges(3,1) * t158 - t186;
t168 = -Icges(3,2) * t156 + t185;
t167 = Icges(3,5) * t158 - Icges(3,6) * t156;
t139 = pkin(1) * t153 + pkin(4) * t152;
t166 = -V_base(4) * qJ(1) + V_base(6) * t139 + V_base(2);
t138 = pkin(1) * t152 - pkin(4) * t153;
t165 = V_base(4) * t138 - V_base(5) * t139 + t172;
t163 = (Icges(3,3) * t152 + t153 * t167) * t146 + (Icges(3,5) * t156 + Icges(3,6) * t158) * V_base(6) + (-Icges(3,3) * t153 + t152 * t167) * t145;
t123 = t171 * t152;
t144 = t156 * pkin(2) - pkin(5) * t158;
t162 = t145 * t144 + (-t123 - t138) * V_base(6) + t176;
t124 = t171 * t153;
t161 = V_base(6) * t124 - t144 * t146 + t166;
t160 = t146 * t123 - t145 * t124 + t165;
t101 = -Icges(3,6) * t153 + t152 * t168;
t102 = Icges(3,6) * t152 + t153 * t168;
t103 = -Icges(3,5) * t153 + t152 * t169;
t104 = Icges(3,5) * t152 + t153 * t169;
t141 = Icges(3,2) * t158 + t186;
t142 = Icges(3,1) * t156 + t185;
t159 = (-t102 * t156 + t104 * t158) * t146 + (-t101 * t156 + t103 * t158) * t145 + (-t141 * t156 + t142 * t158) * V_base(6);
t150 = Icges(2,4) * t153;
t147 = -qJD(3) * t158 + V_base(6);
t143 = t156 * rSges(3,1) + rSges(3,2) * t158;
t137 = rSges(2,1) * t153 - rSges(2,2) * t152;
t136 = rSges(2,1) * t152 + rSges(2,2) * t153;
t135 = Icges(2,1) * t153 - t187;
t134 = Icges(2,1) * t152 + t150;
t133 = -Icges(2,2) * t152 + t150;
t132 = Icges(2,2) * t153 + t187;
t129 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t128 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t127 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t118 = t153 * t178 + t146;
t117 = t152 * t178 + t145;
t114 = -rSges(4,3) * t158 + (rSges(4,1) * t157 - rSges(4,2) * t155) * t156;
t106 = t152 * rSges(3,3) + t153 * t170;
t105 = -t153 * rSges(3,3) + t152 * t170;
t96 = V_base(5) * rSges(2,3) - t136 * V_base(6) + t176;
t95 = t137 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t94 = t136 * V_base(4) - t137 * V_base(5) + t172;
t93 = rSges(4,1) * t122 + rSges(4,2) * t121 + rSges(4,3) * t181;
t91 = rSges(4,1) * t120 + rSges(4,2) * t119 + rSges(4,3) * t183;
t75 = t143 * t145 + (-t105 - t138) * V_base(6) + t176;
t74 = t106 * V_base(6) - t143 * t146 + t166;
t73 = t105 * t146 - t106 * t145 + t165;
t72 = t114 * t117 - t147 * t91 + t162;
t71 = -t114 * t118 + t147 * t93 + t161;
t70 = -t117 * t93 + t118 * t91 + t160;
t69 = t117 * t188 - t147 * t190 + t153 * t177 + t162;
t68 = -t118 * t188 + t147 * t189 + t152 * t177 + t161;
t67 = -qJD(4) * t158 - t117 * t189 + t118 * t190 + t160;
t1 = m(1) * (t127 ^ 2 + t128 ^ 2 + t129 ^ 2) / 0.2e1 + m(2) * (t94 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + m(3) * (t73 ^ 2 + t74 ^ 2 + t75 ^ 2) / 0.2e1 + t146 * (t163 * t152 + t159 * t153) / 0.2e1 + t145 * (t152 * t159 - t153 * t163) / 0.2e1 + m(4) * (t70 ^ 2 + t71 ^ 2 + t72 ^ 2) / 0.2e1 + m(5) * (t67 ^ 2 + t68 ^ 2 + t69 ^ 2) / 0.2e1 + ((t197 * t119 + t196 * t120 + t198 * t183) * t147 + (t203 * t119 + t201 * t120 + t205 * t183) * t118 + (t204 * t119 + t202 * t120 + t206 * t183) * t117) * t117 / 0.2e1 + ((t197 * t121 + t196 * t122 + t198 * t181) * t147 + (t203 * t121 + t201 * t122 + t205 * t181) * t118 + (t204 * t121 + t202 * t122 + t206 * t181) * t117) * t118 / 0.2e1 + ((-t206 * t117 - t205 * t118 - t198 * t147) * t158 + ((-t197 * t155 + t196 * t157) * t147 + (-t203 * t155 + t201 * t157) * t118 + (-t204 * t155 + t202 * t157) * t117) * t156) * t147 / 0.2e1 + ((-t152 * t132 + t134 * t153 + Icges(1,4)) * V_base(5) + (-t152 * t133 + t153 * t135 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t153 * t132 + t152 * t134 + Icges(1,2)) * V_base(5) + (t133 * t153 + t135 * t152 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t102 * t158 + t156 * t104) * t146 + (t101 * t158 + t156 * t103) * t145 + (t158 * t141 + t156 * t142 + Icges(1,3) + Icges(2,3)) * V_base(6)) * V_base(6) / 0.2e1 + V_base(6) * V_base(4) * (Icges(2,5) * t153 - Icges(2,6) * t152 + Icges(1,5)) + V_base(6) * V_base(5) * (Icges(2,5) * t152 + Icges(2,6) * t153 + Icges(1,6));
T = t1;
