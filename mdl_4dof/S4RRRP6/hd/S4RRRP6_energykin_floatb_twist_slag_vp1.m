% Calculate kinetic energy for
% S4RRRP6
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
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRRP6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP6_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RRRP6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP6_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP6_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP6_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:18:06
% EndTime: 2019-12-31 17:18:08
% DurationCPUTime: 1.58s
% Computational Cost: add. (679->204), mult. (1305->299), div. (0->0), fcn. (1249->6), ass. (0->100)
t205 = Icges(4,1) + Icges(5,1);
t204 = Icges(4,4) + Icges(5,4);
t203 = -Icges(5,5) - Icges(4,5);
t202 = Icges(4,2) + Icges(5,2);
t201 = Icges(4,6) + Icges(5,6);
t200 = -Icges(5,3) - Icges(4,3);
t149 = sin(qJ(3));
t152 = cos(qJ(3));
t154 = cos(qJ(1));
t151 = sin(qJ(1));
t153 = cos(qJ(2));
t174 = t151 * t153;
t114 = -t149 * t174 - t152 * t154;
t178 = t149 * t154;
t115 = t152 * t174 - t178;
t150 = sin(qJ(2));
t177 = t150 * t151;
t199 = t201 * t114 - t203 * t115 - t200 * t177;
t173 = t153 * t154;
t116 = -t149 * t173 + t151 * t152;
t175 = t151 * t149;
t117 = t152 * t173 + t175;
t176 = t150 * t154;
t198 = t201 * t116 - t203 * t117 - t200 * t176;
t197 = t202 * t114 + t204 * t115 + t201 * t177;
t196 = t202 * t116 + t204 * t117 + t201 * t176;
t195 = t204 * t114 + t205 * t115 - t203 * t177;
t194 = t204 * t116 + t205 * t117 - t203 * t176;
t193 = t200 * t153 + (-t201 * t149 - t203 * t152) * t150;
t192 = -t201 * t153 + (-t202 * t149 + t204 * t152) * t150;
t191 = t203 * t153 + (-t204 * t149 + t205 * t152) * t150;
t186 = pkin(3) * t152;
t159 = qJ(4) * t150 + t153 * t186;
t184 = rSges(5,1) * t115 + rSges(5,2) * t114 + rSges(5,3) * t177 - pkin(3) * t178 + t151 * t159;
t183 = t117 * rSges(5,1) + t116 * rSges(5,2) + rSges(5,3) * t176 + pkin(3) * t175 + t154 * t159;
t182 = (-qJ(4) - rSges(5,3)) * t153 + (rSges(5,1) * t152 - rSges(5,2) * t149 + t186) * t150;
t181 = Icges(2,4) * t151;
t180 = Icges(3,4) * t150;
t179 = Icges(3,4) * t153;
t172 = qJD(3) * t150;
t171 = qJD(4) * t150;
t170 = V_base(5) * pkin(4) + V_base(1);
t142 = qJD(2) * t151 + V_base(4);
t144 = V_base(6) + qJD(1);
t167 = pkin(2) * t153 + pkin(6) * t150;
t141 = -qJD(2) * t154 + V_base(5);
t166 = rSges(3,1) * t153 - rSges(3,2) * t150;
t165 = Icges(3,1) * t153 - t180;
t164 = -Icges(3,2) * t150 + t179;
t163 = Icges(3,5) * t153 - Icges(3,6) * t150;
t140 = pkin(1) * t154 + t151 * pkin(5);
t162 = -V_base(4) * pkin(4) + t144 * t140 + V_base(2);
t139 = t151 * pkin(1) - pkin(5) * t154;
t161 = V_base(4) * t139 - t140 * V_base(5) + V_base(3);
t160 = (Icges(3,5) * t150 + Icges(3,6) * t153) * t144 + t141 * (-Icges(3,3) * t154 + t151 * t163) + t142 * (Icges(3,3) * t151 + t154 * t163);
t119 = t167 * t151;
t138 = pkin(2) * t150 - pkin(6) * t153;
t158 = t141 * t138 + (-t119 - t139) * t144 + t170;
t120 = t167 * t154;
t157 = t144 * t120 - t138 * t142 + t162;
t156 = t142 * t119 - t120 * t141 + t161;
t101 = -Icges(3,6) * t154 + t151 * t164;
t102 = Icges(3,6) * t151 + t154 * t164;
t105 = -Icges(3,5) * t154 + t151 * t165;
t106 = Icges(3,5) * t151 + t154 * t165;
t128 = Icges(3,2) * t153 + t180;
t131 = Icges(3,1) * t150 + t179;
t155 = (-t102 * t150 + t106 * t153) * t142 + (-t101 * t150 + t105 * t153) * t141 + (-t128 * t150 + t131 * t153) * t144;
t146 = Icges(2,4) * t154;
t137 = rSges(2,1) * t154 - t151 * rSges(2,2);
t136 = t151 * rSges(2,1) + rSges(2,2) * t154;
t135 = rSges(3,1) * t150 + rSges(3,2) * t153;
t134 = -qJD(3) * t153 + t144;
t133 = Icges(2,1) * t154 - t181;
t132 = Icges(2,1) * t151 + t146;
t130 = -Icges(2,2) * t151 + t146;
t129 = Icges(2,2) * t154 + t181;
t124 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t123 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t122 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t113 = t154 * t172 + t142;
t112 = t151 * t172 + t141;
t110 = t151 * rSges(3,3) + t154 * t166;
t109 = -rSges(3,3) * t154 + t151 * t166;
t108 = -rSges(4,3) * t153 + (rSges(4,1) * t152 - rSges(4,2) * t149) * t150;
t91 = V_base(5) * rSges(2,3) - t136 * t144 + t170;
t90 = t137 * t144 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t89 = t136 * V_base(4) - t137 * V_base(5) + V_base(3);
t88 = t117 * rSges(4,1) + t116 * rSges(4,2) + rSges(4,3) * t176;
t86 = rSges(4,1) * t115 + rSges(4,2) * t114 + rSges(4,3) * t177;
t70 = t135 * t141 + (-t109 - t139) * t144 + t170;
t69 = t110 * t144 - t135 * t142 + t162;
t68 = t109 * t142 - t110 * t141 + t161;
t67 = t108 * t112 - t134 * t86 + t158;
t66 = -t108 * t113 + t134 * t88 + t157;
t65 = -t112 * t88 + t113 * t86 + t156;
t64 = t112 * t182 - t134 * t184 + t154 * t171 + t158;
t63 = -t113 * t182 + t134 * t183 + t151 * t171 + t157;
t62 = -qJD(4) * t153 - t112 * t183 + t113 * t184 + t156;
t1 = m(1) * (t122 ^ 2 + t123 ^ 2 + t124 ^ 2) / 0.2e1 + m(2) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + m(3) * (t68 ^ 2 + t69 ^ 2 + t70 ^ 2) / 0.2e1 + t142 * (t160 * t151 + t155 * t154) / 0.2e1 + t141 * (t155 * t151 - t160 * t154) / 0.2e1 + m(4) * (t65 ^ 2 + t66 ^ 2 + t67 ^ 2) / 0.2e1 + m(5) * (t62 ^ 2 + t63 ^ 2 + t64 ^ 2) / 0.2e1 + ((t192 * t114 + t191 * t115 + t193 * t177) * t134 + (t196 * t114 + t194 * t115 + t198 * t177) * t113 + (t197 * t114 + t195 * t115 + t199 * t177) * t112) * t112 / 0.2e1 + ((t192 * t116 + t191 * t117 + t193 * t176) * t134 + (t196 * t116 + t194 * t117 + t198 * t176) * t113 + (t197 * t116 + t195 * t117 + t199 * t176) * t112) * t113 / 0.2e1 + ((-t199 * t112 - t198 * t113 - t193 * t134) * t153 + ((-t192 * t149 + t191 * t152) * t134 + (-t196 * t149 + t194 * t152) * t113 + (-t197 * t149 + t195 * t152) * t112) * t150) * t134 / 0.2e1 + ((t102 * t153 + t106 * t150) * t142 + (t101 * t153 + t105 * t150) * t141 + (t153 * t128 + t150 * t131 + Icges(2,3)) * t144) * t144 / 0.2e1 + ((-t151 * t129 + t132 * t154 + Icges(1,4)) * V_base(5) + (-t151 * t130 + t154 * t133 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t154 * t129 + t151 * t132 + Icges(1,2)) * V_base(5) + (t130 * t154 + t151 * t133 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t144 * (Icges(2,5) * t154 - Icges(2,6) * t151) + V_base(5) * t144 * (Icges(2,5) * t151 + Icges(2,6) * t154) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
