% Calculate kinetic energy for
% S4RRRP7
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
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRRP7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP7_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RRRP7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP7_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP7_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP7_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:20:06
% EndTime: 2019-12-31 17:20:08
% DurationCPUTime: 1.58s
% Computational Cost: add. (661->199), mult. (1296->292), div. (0->0), fcn. (1250->6), ass. (0->97)
t199 = Icges(4,1) + Icges(5,1);
t198 = -Icges(4,4) + Icges(5,5);
t197 = Icges(5,4) + Icges(4,5);
t196 = Icges(4,2) + Icges(5,3);
t195 = -Icges(5,6) + Icges(4,6);
t194 = -Icges(4,3) - Icges(5,2);
t193 = rSges(5,1) + pkin(3);
t192 = rSges(5,3) + qJ(4);
t147 = sin(qJ(3));
t150 = cos(qJ(3));
t152 = cos(qJ(1));
t149 = sin(qJ(1));
t151 = cos(qJ(2));
t171 = t149 * t151;
t112 = t147 * t171 + t150 * t152;
t113 = -t147 * t152 + t150 * t171;
t148 = sin(qJ(2));
t173 = t148 * t149;
t191 = t196 * t112 + t198 * t113 - t195 * t173;
t170 = t151 * t152;
t114 = t147 * t170 - t149 * t150;
t115 = t149 * t147 + t150 * t170;
t172 = t148 * t152;
t190 = t196 * t114 + t198 * t115 - t195 * t172;
t189 = -t195 * t112 + t197 * t113 - t194 * t173;
t188 = -t195 * t114 + t197 * t115 - t194 * t172;
t187 = t198 * t112 + t199 * t113 + t197 * t173;
t186 = t198 * t114 + t199 * t115 + t197 * t172;
t185 = t195 * t151 + (t196 * t147 + t198 * t150) * t148;
t184 = t194 * t151 + (-t195 * t147 + t197 * t150) * t148;
t183 = -t197 * t151 + (t198 * t147 + t199 * t150) * t148;
t178 = rSges(5,2) * t173 + t192 * t112 + t193 * t113;
t177 = rSges(5,2) * t172 + t192 * t114 + t193 * t115;
t176 = Icges(2,4) * t149;
t175 = Icges(3,4) * t148;
t174 = Icges(3,4) * t151;
t169 = -rSges(5,2) * t151 + (t192 * t147 + t193 * t150) * t148;
t168 = qJD(3) * t148;
t167 = V_base(5) * pkin(4) + V_base(1);
t141 = qJD(2) * t149 + V_base(4);
t143 = V_base(6) + qJD(1);
t164 = pkin(2) * t151 + pkin(6) * t148;
t140 = -qJD(2) * t152 + V_base(5);
t163 = rSges(3,1) * t151 - rSges(3,2) * t148;
t162 = Icges(3,1) * t151 - t175;
t161 = -Icges(3,2) * t148 + t174;
t160 = Icges(3,5) * t151 - Icges(3,6) * t148;
t139 = pkin(1) * t152 + t149 * pkin(5);
t159 = -V_base(4) * pkin(4) + t143 * t139 + V_base(2);
t138 = t149 * pkin(1) - pkin(5) * t152;
t158 = V_base(4) * t138 - t139 * V_base(5) + V_base(3);
t157 = (Icges(3,5) * t148 + Icges(3,6) * t151) * t143 + (-Icges(3,3) * t152 + t149 * t160) * t140 + (Icges(3,3) * t149 + t152 * t160) * t141;
t118 = t164 * t149;
t137 = pkin(2) * t148 - pkin(6) * t151;
t156 = t140 * t137 + (-t118 - t138) * t143 + t167;
t119 = t164 * t152;
t155 = t143 * t119 - t137 * t141 + t159;
t154 = t141 * t118 - t119 * t140 + t158;
t100 = Icges(3,6) * t149 + t152 * t161;
t103 = -Icges(3,5) * t152 + t149 * t162;
t104 = Icges(3,5) * t149 + t152 * t162;
t127 = Icges(3,2) * t151 + t175;
t130 = Icges(3,1) * t148 + t174;
t99 = -Icges(3,6) * t152 + t149 * t161;
t153 = (-t100 * t148 + t104 * t151) * t141 + (t103 * t151 - t148 * t99) * t140 + (-t127 * t148 + t130 * t151) * t143;
t145 = Icges(2,4) * t152;
t136 = rSges(2,1) * t152 - t149 * rSges(2,2);
t135 = t149 * rSges(2,1) + rSges(2,2) * t152;
t134 = rSges(3,1) * t148 + rSges(3,2) * t151;
t133 = -qJD(3) * t151 + t143;
t132 = Icges(2,1) * t152 - t176;
t131 = Icges(2,1) * t149 + t145;
t129 = -Icges(2,2) * t149 + t145;
t128 = Icges(2,2) * t152 + t176;
t123 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t122 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t121 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t111 = t152 * t168 + t141;
t110 = t149 * t168 + t140;
t108 = t149 * rSges(3,3) + t152 * t163;
t107 = -rSges(3,3) * t152 + t149 * t163;
t106 = -rSges(4,3) * t151 + (rSges(4,1) * t150 - rSges(4,2) * t147) * t148;
t90 = V_base(5) * rSges(2,3) - t135 * t143 + t167;
t89 = t136 * t143 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t88 = t135 * V_base(4) - t136 * V_base(5) + V_base(3);
t85 = t115 * rSges(4,1) - t114 * rSges(4,2) + rSges(4,3) * t172;
t83 = rSges(4,1) * t113 - rSges(4,2) * t112 + rSges(4,3) * t173;
t69 = t134 * t140 + (-t107 - t138) * t143 + t167;
t68 = t108 * t143 - t134 * t141 + t159;
t67 = t107 * t141 - t108 * t140 + t158;
t66 = t106 * t110 - t133 * t83 + t156;
t65 = -t106 * t111 + t133 * t85 + t155;
t64 = -t110 * t85 + t111 * t83 + t154;
t63 = qJD(4) * t114 + t110 * t169 - t133 * t178 + t156;
t62 = qJD(4) * t112 - t111 * t169 + t133 * t177 + t155;
t61 = qJD(4) * t147 * t148 - t110 * t177 + t111 * t178 + t154;
t1 = m(1) * (t121 ^ 2 + t122 ^ 2 + t123 ^ 2) / 0.2e1 + m(2) * (t88 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + m(3) * (t67 ^ 2 + t68 ^ 2 + t69 ^ 2) / 0.2e1 + t141 * (t157 * t149 + t153 * t152) / 0.2e1 + t140 * (t153 * t149 - t157 * t152) / 0.2e1 + m(4) * (t64 ^ 2 + t65 ^ 2 + t66 ^ 2) / 0.2e1 + m(5) * (t61 ^ 2 + t62 ^ 2 + t63 ^ 2) / 0.2e1 + ((t185 * t112 + t183 * t113 + t184 * t173) * t133 + (t190 * t112 + t186 * t113 + t188 * t173) * t111 + (t191 * t112 + t187 * t113 + t189 * t173) * t110) * t110 / 0.2e1 + ((t185 * t114 + t183 * t115 + t184 * t172) * t133 + (t190 * t114 + t186 * t115 + t188 * t172) * t111 + (t191 * t114 + t187 * t115 + t189 * t172) * t110) * t111 / 0.2e1 + ((-t189 * t110 - t188 * t111 - t184 * t133) * t151 + ((t185 * t147 + t183 * t150) * t133 + (t190 * t147 + t186 * t150) * t111 + (t191 * t147 + t187 * t150) * t110) * t148) * t133 / 0.2e1 + ((t100 * t151 + t104 * t148) * t141 + (t103 * t148 + t151 * t99) * t140 + (t151 * t127 + t148 * t130 + Icges(2,3)) * t143) * t143 / 0.2e1 + ((-t149 * t128 + t131 * t152 + Icges(1,4)) * V_base(5) + (-t149 * t129 + t152 * t132 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t152 * t128 + t149 * t131 + Icges(1,2)) * V_base(5) + (t129 * t152 + t149 * t132 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t143 * (Icges(2,5) * t152 - Icges(2,6) * t149) + V_base(5) * t143 * (Icges(2,5) * t149 + Icges(2,6) * t152) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
