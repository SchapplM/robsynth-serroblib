% Calculate kinetic energy for
% S4RPPP1
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
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
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
% Datum: 2019-03-08 18:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPPP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RPPP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPP1_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPP1_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPP1_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:26:14
% EndTime: 2019-03-08 18:26:15
% DurationCPUTime: 1.40s
% Computational Cost: add. (683->184), mult. (1461->231), div. (0->0), fcn. (1481->6), ass. (0->87)
t200 = Icges(3,1) + Icges(4,2) + Icges(5,3);
t199 = Icges(3,4) + Icges(4,6) - Icges(5,6);
t198 = Icges(3,5) + Icges(5,5) - Icges(4,4);
t197 = Icges(3,2) + Icges(5,2) + Icges(4,3);
t196 = Icges(3,6) - Icges(4,5) - Icges(5,4);
t195 = Icges(3,3) + Icges(5,1) + Icges(4,1);
t194 = rSges(5,1) + pkin(3);
t193 = rSges(5,3) + qJ(4);
t152 = sin(pkin(6));
t155 = cos(pkin(4));
t157 = cos(qJ(1));
t154 = cos(pkin(6));
t156 = sin(qJ(1));
t174 = t156 * t154;
t126 = t152 * t157 + t155 * t174;
t178 = t152 * t156;
t127 = t154 * t157 - t155 * t178;
t153 = sin(pkin(4));
t177 = t153 * t156;
t192 = -t199 * t126 + t200 * t127 + t198 * t177;
t175 = t155 * t157;
t124 = -t154 * t175 + t178;
t125 = t152 * t175 + t174;
t176 = t153 * t157;
t191 = -t199 * t124 + t200 * t125 - t198 * t176;
t190 = t197 * t126 - t199 * t127 - t196 * t177;
t189 = t197 * t124 - t199 * t125 + t196 * t176;
t188 = -t196 * t126 + t198 * t127 + t195 * t177;
t187 = -t196 * t124 + t198 * t125 - t195 * t176;
t186 = t195 * t155 + (t198 * t152 + t196 * t154) * t153;
t185 = t196 * t155 + (t199 * t152 + t197 * t154) * t153;
t184 = t198 * t155 + (t200 * t152 + t199 * t154) * t153;
t182 = rSges(5,2) * t126 + t193 * t127 + t177 * t194;
t181 = t124 * rSges(5,2) + t193 * t125 - t176 * t194;
t180 = Icges(2,4) * t156;
t179 = qJ(2) * t155;
t101 = pkin(2) * t125 + qJ(3) * t124;
t131 = t156 * pkin(1) - qJ(2) * t176;
t173 = -t101 - t131;
t102 = pkin(2) * t127 + qJ(3) * t126;
t132 = pkin(1) * t157 + qJ(2) * t177;
t172 = -t102 - t132;
t171 = t194 * t155 + (-rSges(5,2) * t154 + t152 * t193) * t153;
t170 = qJD(2) * t153;
t169 = qJD(3) * t154;
t168 = V_base(5) * pkin(5) + V_base(1);
t165 = -pkin(5) - t179;
t164 = qJD(2) * t155 + V_base(4) * t131 + V_base(3);
t128 = (pkin(2) * t152 - qJ(3) * t154) * t153;
t163 = -t128 + t165;
t162 = V_base(4) * t101 + t164;
t161 = t156 * t170 + V_base(5) * t179 + t168;
t149 = V_base(6) + qJD(1);
t160 = t149 * t132 - t157 * t170 + V_base(2);
t159 = qJD(3) * t126 + V_base(5) * t128 + t161;
t158 = qJD(3) * t124 + t149 * t102 + t160;
t150 = Icges(2,4) * t157;
t143 = rSges(2,1) * t157 - t156 * rSges(2,2);
t142 = t156 * rSges(2,1) + rSges(2,2) * t157;
t141 = Icges(2,1) * t157 - t180;
t140 = Icges(2,1) * t156 + t150;
t139 = -Icges(2,2) * t156 + t150;
t138 = Icges(2,2) * t157 + t180;
t137 = Icges(2,5) * t157 - Icges(2,6) * t156;
t136 = Icges(2,5) * t156 + Icges(2,6) * t157;
t135 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t134 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t133 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t119 = rSges(4,1) * t155 + (-rSges(4,2) * t152 - rSges(4,3) * t154) * t153;
t117 = rSges(3,3) * t155 + (rSges(3,1) * t152 + rSges(3,2) * t154) * t153;
t105 = V_base(5) * rSges(2,3) - t142 * t149 + t168;
t104 = t143 * t149 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t103 = t142 * V_base(4) - t143 * V_base(5) + V_base(3);
t98 = rSges(3,1) * t127 - rSges(3,2) * t126 + rSges(3,3) * t177;
t97 = t125 * rSges(3,1) - t124 * rSges(3,2) - rSges(3,3) * t176;
t96 = -rSges(4,1) * t176 - t125 * rSges(4,2) + t124 * rSges(4,3);
t94 = rSges(4,1) * t177 - rSges(4,2) * t127 + rSges(4,3) * t126;
t74 = t117 * V_base(5) + (-t131 - t97) * t149 + t161;
t73 = t149 * t98 + (-t117 + t165) * V_base(4) + t160;
t72 = t97 * V_base(4) + (-t132 - t98) * V_base(5) + t164;
t71 = t119 * V_base(5) + (-t96 + t173) * t149 + t159;
t70 = t149 * t94 + (-t119 + t163) * V_base(4) + t158;
t69 = -t153 * t169 + t96 * V_base(4) + (-t94 + t172) * V_base(5) + t162;
t68 = qJD(4) * t127 + t171 * V_base(5) + (t173 - t181) * t149 + t159;
t67 = qJD(4) * t125 + t182 * t149 + (t163 - t171) * V_base(4) + t158;
t66 = t181 * V_base(4) + (qJD(4) * t152 - t169) * t153 + (t172 - t182) * V_base(5) + t162;
t1 = m(1) * (t133 ^ 2 + t134 ^ 2 + t135 ^ 2) / 0.2e1 + m(2) * (t103 ^ 2 + t104 ^ 2 + t105 ^ 2) / 0.2e1 + m(3) * (t72 ^ 2 + t73 ^ 2 + t74 ^ 2) / 0.2e1 + m(4) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + m(5) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + (Icges(2,3) * t149 + t136 * V_base(5) + t137 * V_base(4) + (t186 * t149 + t187 * V_base(5) + t188 * V_base(4)) * t155 + ((t191 * t152 - t189 * t154) * V_base(5) + (t192 * t152 - t190 * t154) * V_base(4) + (t184 * t152 + t185 * t154) * t149) * t153) * t149 / 0.2e1 + ((-t185 * t126 + t184 * t127 + t186 * t177 + t137) * t149 + (t189 * t126 + t191 * t127 - t156 * t138 + t157 * t140 + t187 * t177 + Icges(1,4)) * V_base(5) + (t190 * t126 + t192 * t127 - t156 * t139 + t157 * t141 + t188 * t177 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((-t185 * t124 + t184 * t125 - t186 * t176 + t136) * t149 + (t189 * t124 + t191 * t125 + t157 * t138 + t156 * t140 - t187 * t176 + Icges(1,2)) * V_base(5) + (t190 * t124 + t192 * t125 + t139 * t157 + t156 * t141 - t188 * t176 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
