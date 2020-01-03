% Calculate kinetic energy for
% S5RRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRP3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP3_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRRP3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP3_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP3_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP3_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:04
% EndTime: 2019-12-31 21:49:05
% DurationCPUTime: 1.65s
% Computational Cost: add. (1070->198), mult. (751->255), div. (0->0), fcn. (529->8), ass. (0->106)
t226 = Icges(5,4) - Icges(6,5);
t225 = Icges(5,1) + Icges(6,1);
t224 = Icges(5,2) + Icges(6,3);
t158 = cos(qJ(4));
t223 = t226 * t158;
t156 = sin(qJ(4));
t222 = t226 * t156;
t221 = Icges(6,4) + Icges(5,5);
t220 = Icges(5,6) - Icges(6,6);
t219 = t224 * t156 - t223;
t218 = t225 * t158 - t222;
t217 = rSges(6,1) + pkin(4);
t216 = rSges(6,3) + qJ(5);
t155 = qJ(1) + qJ(2);
t151 = qJ(3) + t155;
t145 = sin(t151);
t146 = cos(t151);
t215 = t219 * t145 + t220 * t146;
t214 = -t220 * t145 + t219 * t146;
t213 = t218 * t145 - t221 * t146;
t212 = t221 * t145 + t218 * t146;
t211 = Icges(6,2) + Icges(5,3);
t210 = -t224 * t158 - t222;
t209 = t225 * t156 + t223;
t208 = -t220 * t156 + t221 * t158;
t207 = t216 * t156 + t217 * t158;
t119 = -qJD(4) * t146 + V_base(5);
t120 = qJD(4) * t145 + V_base(4);
t147 = V_base(6) + qJD(1);
t144 = qJD(2) + t147;
t141 = qJD(3) + t144;
t206 = (t210 * t156 + t209 * t158) * t141 + (t214 * t156 + t212 * t158) * t120 + (t215 * t156 + t213 * t158) * t119;
t205 = (t221 * t156 + t220 * t158) * t141 + (t211 * t145 + t208 * t146) * t120 + (t208 * t145 - t211 * t146) * t119;
t204 = -pkin(5) - pkin(6);
t157 = sin(qJ(1));
t200 = pkin(1) * t157;
t159 = cos(qJ(1));
t199 = pkin(1) * t159;
t148 = sin(t155);
t198 = pkin(2) * t148;
t149 = cos(t155);
t197 = pkin(2) * t149;
t196 = -t146 * rSges(6,2) + t207 * t145;
t195 = t145 * rSges(6,2) + t207 * t146;
t194 = Icges(2,4) * t157;
t193 = Icges(3,4) * t148;
t192 = Icges(4,4) * t145;
t187 = t217 * t156 - t216 * t158;
t186 = qJD(5) * t156;
t185 = -pkin(7) + t204;
t184 = t147 * t199 + V_base(2);
t183 = V_base(4) * t200 + V_base(3);
t182 = V_base(5) * pkin(5) + V_base(1);
t179 = t144 * t197 + t184;
t178 = V_base(4) * t198 + t183;
t177 = -t197 - t199;
t176 = rSges(5,1) * t158 - rSges(5,2) * t156;
t167 = V_base(5) * pkin(6) - t147 * t200 + t182;
t106 = pkin(3) * t146 + pkin(8) * t145;
t164 = t141 * t106 + t185 * V_base(4) + t179;
t163 = V_base(5) * pkin(7) - t144 * t198 + t167;
t105 = pkin(3) * t145 - pkin(8) * t146;
t162 = V_base(4) * t105 + (-t106 + t177) * V_base(5) + t178;
t150 = Icges(2,4) * t159;
t143 = Icges(3,4) * t149;
t140 = Icges(4,4) * t146;
t137 = rSges(2,1) * t159 - t157 * rSges(2,2);
t136 = t157 * rSges(2,1) + rSges(2,2) * t159;
t135 = rSges(5,1) * t156 + rSges(5,2) * t158;
t132 = Icges(2,1) * t159 - t194;
t131 = Icges(2,1) * t157 + t150;
t128 = -Icges(2,2) * t157 + t150;
t127 = Icges(2,2) * t159 + t194;
t118 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t117 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t116 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t114 = rSges(3,1) * t149 - rSges(3,2) * t148;
t113 = rSges(3,1) * t148 + rSges(3,2) * t149;
t112 = Icges(3,1) * t149 - t193;
t111 = Icges(3,1) * t148 + t143;
t110 = -Icges(3,2) * t148 + t143;
t109 = Icges(3,2) * t149 + t193;
t104 = rSges(4,1) * t146 - rSges(4,2) * t145;
t103 = rSges(4,1) * t145 + rSges(4,2) * t146;
t102 = Icges(4,1) * t146 - t192;
t101 = Icges(4,1) * t145 + t140;
t100 = -Icges(4,2) * t145 + t140;
t99 = Icges(4,2) * t146 + t192;
t92 = V_base(5) * rSges(2,3) - t136 * t147 + t182;
t91 = t147 * t137 + V_base(2) + (-pkin(5) - rSges(2,3)) * V_base(4);
t90 = t136 * V_base(4) - t137 * V_base(5) + V_base(3);
t89 = t145 * rSges(5,3) + t176 * t146;
t87 = -t146 * rSges(5,3) + t176 * t145;
t73 = V_base(5) * rSges(3,3) - t113 * t144 + t167;
t72 = t144 * t114 + (-rSges(3,3) + t204) * V_base(4) + t184;
t71 = V_base(4) * t113 + (-t114 - t199) * V_base(5) + t183;
t70 = V_base(5) * rSges(4,3) - t103 * t141 + t163;
t69 = t141 * t104 + (-rSges(4,3) + t185) * V_base(4) + t179;
t68 = V_base(4) * t103 + (-t104 + t177) * V_base(5) + t178;
t67 = t119 * t135 + (-t105 - t87) * t141 + t163;
t66 = -t120 * t135 + t141 * t89 + t164;
t65 = -t119 * t89 + t120 * t87 + t162;
t64 = t146 * t186 + t187 * t119 + (-t105 - t196) * t141 + t163;
t63 = -t187 * t120 + t195 * t141 + t145 * t186 + t164;
t62 = -qJD(5) * t158 - t195 * t119 + t196 * t120 + t162;
t1 = m(1) * (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + m(2) * (t90 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(3) * (t71 ^ 2 + t72 ^ 2 + t73 ^ 2) / 0.2e1 + m(4) * (t68 ^ 2 + t69 ^ 2 + t70 ^ 2) / 0.2e1 + m(5) * (t65 ^ 2 + t66 ^ 2 + t67 ^ 2) / 0.2e1 + m(6) * (t62 ^ 2 + t63 ^ 2 + t64 ^ 2) / 0.2e1 + (t145 * t206 - t205 * t146) * t119 / 0.2e1 + (t205 * t145 + t146 * t206) * t120 / 0.2e1 + ((t212 * t156 - t214 * t158) * t120 + (t213 * t156 - t215 * t158) * t119 + (t209 * t156 - t210 * t158 + Icges(4,3)) * t141) * t141 / 0.2e1 + ((t146 * t101 - t148 * t109 + t149 * t111 - t157 * t127 + t159 * t131 - t145 * t99 + Icges(1,4)) * V_base(5) + (-t145 * t100 + t146 * t102 - t148 * t110 + t149 * t112 - t157 * t128 + t159 * t132 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t145 * t101 + t149 * t109 + t148 * t111 + t159 * t127 + t157 * t131 + t146 * t99 + Icges(1,2)) * V_base(5) + (t146 * t100 + t145 * t102 + t149 * t110 + t148 * t112 + t159 * t128 + t157 * t132 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t141 * (Icges(4,5) * t146 - Icges(4,6) * t145) + V_base(5) * t141 * (Icges(4,5) * t145 + Icges(4,6) * t146) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t157 + Icges(2,6) * t159) * V_base(5) + (Icges(2,5) * t159 - Icges(2,6) * t157) * V_base(4) + Icges(2,3) * t147 / 0.2e1) * t147 + ((Icges(3,5) * t148 + Icges(3,6) * t149) * V_base(5) + (Icges(3,5) * t149 - Icges(3,6) * t148) * V_base(4) + Icges(3,3) * t144 / 0.2e1) * t144;
T = t1;
