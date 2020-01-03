% Calculate kinetic energy for
% S4PRRP6
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
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRRP6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP6_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PRRP6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP6_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP6_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP6_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:10
% EndTime: 2019-12-31 16:30:11
% DurationCPUTime: 1.53s
% Computational Cost: add. (629->199), mult. (1296->287), div. (0->0), fcn. (1250->6), ass. (0->97)
t203 = Icges(4,1) + Icges(5,1);
t202 = -Icges(4,4) + Icges(5,5);
t201 = Icges(5,4) + Icges(4,5);
t200 = Icges(4,2) + Icges(5,3);
t199 = -Icges(5,6) + Icges(4,6);
t198 = -Icges(4,3) - Icges(5,2);
t197 = rSges(5,1) + pkin(3);
t196 = rSges(5,3) + qJ(4);
t148 = sin(pkin(6));
t149 = cos(pkin(6));
t152 = cos(qJ(3));
t150 = sin(qJ(3));
t153 = cos(qJ(2));
t174 = t150 * t153;
t114 = t148 * t174 + t149 * t152;
t173 = t152 * t153;
t115 = t148 * t173 - t149 * t150;
t151 = sin(qJ(2));
t176 = t148 * t151;
t195 = t200 * t114 + t202 * t115 - t199 * t176;
t116 = -t148 * t152 + t149 * t174;
t117 = t148 * t150 + t149 * t173;
t175 = t149 * t151;
t194 = t200 * t116 + t202 * t117 - t199 * t175;
t193 = -t199 * t114 + t201 * t115 - t198 * t176;
t192 = -t199 * t116 + t201 * t117 - t198 * t175;
t191 = t202 * t114 + t203 * t115 + t201 * t176;
t190 = t202 * t116 + t203 * t117 + t201 * t175;
t187 = t199 * t153 + (t200 * t150 + t202 * t152) * t151;
t186 = t198 * t153 + (-t199 * t150 + t201 * t152) * t151;
t185 = -t201 * t153 + (t202 * t150 + t203 * t152) * t151;
t181 = rSges(5,2) * t176 + t196 * t114 + t115 * t197;
t180 = rSges(5,2) * t175 + t196 * t116 + t117 * t197;
t179 = Icges(2,4) * t148;
t178 = Icges(3,4) * t151;
t177 = Icges(3,4) * t153;
t172 = -rSges(5,2) * t153 + (t196 * t150 + t152 * t197) * t151;
t171 = qJD(3) * t151;
t170 = V_base(5) * qJ(1) + V_base(1);
t166 = qJD(1) + V_base(3);
t142 = qJD(2) * t148 + V_base(4);
t165 = pkin(2) * t153 + pkin(5) * t151;
t141 = -qJD(2) * t149 + V_base(5);
t164 = rSges(3,1) * t153 - rSges(3,2) * t151;
t163 = Icges(3,1) * t153 - t178;
t162 = -Icges(3,2) * t151 + t177;
t161 = Icges(3,5) * t153 - Icges(3,6) * t151;
t135 = pkin(1) * t149 + pkin(4) * t148;
t160 = -V_base(4) * qJ(1) + V_base(6) * t135 + V_base(2);
t134 = pkin(1) * t148 - pkin(4) * t149;
t159 = V_base(4) * t134 - t135 * V_base(5) + t166;
t158 = (Icges(3,5) * t151 + Icges(3,6) * t153) * V_base(6) + (-Icges(3,3) * t149 + t148 * t161) * t141 + (Icges(3,3) * t148 + t149 * t161) * t142;
t118 = t165 * t148;
t140 = t151 * pkin(2) - pkin(5) * t153;
t157 = t141 * t140 + (-t118 - t134) * V_base(6) + t170;
t119 = t165 * t149;
t156 = V_base(6) * t119 - t140 * t142 + t160;
t155 = t142 * t118 - t119 * t141 + t159;
t137 = Icges(3,2) * t153 + t178;
t138 = Icges(3,1) * t151 + t177;
t96 = -Icges(3,6) * t149 + t148 * t162;
t97 = Icges(3,6) * t148 + t149 * t162;
t98 = -Icges(3,5) * t149 + t148 * t163;
t99 = Icges(3,5) * t148 + t149 * t163;
t154 = (-t151 * t97 + t153 * t99) * t142 + (-t151 * t96 + t153 * t98) * t141 + (-t137 * t151 + t138 * t153) * V_base(6);
t146 = Icges(2,4) * t149;
t143 = -qJD(3) * t153 + V_base(6);
t139 = t151 * rSges(3,1) + rSges(3,2) * t153;
t133 = rSges(2,1) * t149 - rSges(2,2) * t148;
t132 = rSges(2,1) * t148 + rSges(2,2) * t149;
t131 = Icges(2,1) * t149 - t179;
t130 = Icges(2,1) * t148 + t146;
t129 = -Icges(2,2) * t148 + t146;
t128 = Icges(2,2) * t149 + t179;
t125 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t124 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t123 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t113 = t149 * t171 + t142;
t112 = t148 * t171 + t141;
t109 = -rSges(4,3) * t153 + (rSges(4,1) * t152 - rSges(4,2) * t150) * t151;
t101 = t148 * rSges(3,3) + t149 * t164;
t100 = -t149 * rSges(3,3) + t148 * t164;
t92 = V_base(5) * rSges(2,3) - t132 * V_base(6) + t170;
t91 = t133 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t90 = t132 * V_base(4) - t133 * V_base(5) + t166;
t87 = rSges(4,1) * t117 - rSges(4,2) * t116 + rSges(4,3) * t175;
t85 = rSges(4,1) * t115 - rSges(4,2) * t114 + rSges(4,3) * t176;
t71 = t139 * t141 + (-t100 - t134) * V_base(6) + t170;
t70 = t101 * V_base(6) - t139 * t142 + t160;
t69 = t100 * t142 - t101 * t141 + t159;
t68 = t109 * t112 - t143 * t85 + t157;
t67 = -t109 * t113 + t143 * t87 + t156;
t66 = -t112 * t87 + t113 * t85 + t155;
t65 = qJD(4) * t116 + t112 * t172 - t143 * t181 + t157;
t64 = qJD(4) * t114 - t113 * t172 + t143 * t180 + t156;
t63 = qJD(4) * t150 * t151 - t112 * t180 + t113 * t181 + t155;
t1 = m(1) * (t123 ^ 2 + t124 ^ 2 + t125 ^ 2) / 0.2e1 + m(2) * (t90 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(3) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + t142 * (t158 * t148 + t154 * t149) / 0.2e1 + t141 * (t154 * t148 - t158 * t149) / 0.2e1 + m(4) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + m(5) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + ((t187 * t114 + t185 * t115 + t186 * t176) * t143 + (t194 * t114 + t190 * t115 + t192 * t176) * t113 + (t195 * t114 + t191 * t115 + t193 * t176) * t112) * t112 / 0.2e1 + ((t187 * t116 + t185 * t117 + t186 * t175) * t143 + (t194 * t116 + t190 * t117 + t192 * t175) * t113 + (t195 * t116 + t191 * t117 + t193 * t175) * t112) * t113 / 0.2e1 + ((-t193 * t112 - t192 * t113 - t186 * t143) * t153 + ((t187 * t150 + t185 * t152) * t143 + (t194 * t150 + t190 * t152) * t113 + (t195 * t150 + t191 * t152) * t112) * t151) * t143 / 0.2e1 + ((-t128 * t148 + t130 * t149 + Icges(1,4)) * V_base(5) + (-t148 * t129 + t149 * t131 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t149 * t128 + t148 * t130 + Icges(1,2)) * V_base(5) + (t129 * t149 + t131 * t148 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t151 * t99 + t153 * t97) * t142 + (t151 * t98 + t153 * t96) * t141 + (t153 * t137 + t151 * t138 + Icges(1,3) + Icges(2,3)) * V_base(6)) * V_base(6) / 0.2e1 + V_base(6) * V_base(4) * (Icges(2,5) * t149 - Icges(2,6) * t148 + Icges(1,5)) + V_base(6) * V_base(5) * (Icges(2,5) * t148 + Icges(2,6) * t149 + Icges(1,6));
T = t1;
