% Calculate kinetic energy for
% S5RPPPR4
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
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
% Datum: 2019-12-31 17:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPPR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR4_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPPR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR4_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR4_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR4_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:45:04
% EndTime: 2019-12-31 17:45:05
% DurationCPUTime: 1.35s
% Computational Cost: add. (841->222), mult. (740->274), div. (0->0), fcn. (520->8), ass. (0->117)
t215 = Icges(3,4) + Icges(4,6);
t214 = Icges(3,1) + Icges(4,2);
t213 = -Icges(4,4) + Icges(3,5);
t212 = Icges(4,5) - Icges(3,6);
t211 = -Icges(3,2) - Icges(4,3);
t146 = qJ(1) + pkin(7);
t139 = cos(t146);
t210 = t215 * t139;
t137 = sin(t146);
t209 = t215 * t137;
t208 = t214 * t137 + t210;
t207 = t214 * t139 - t209;
t150 = sin(qJ(1));
t151 = cos(qJ(1));
t206 = Icges(2,5) * t150 + Icges(2,6) * t151 + t213 * t137 - t212 * t139;
t205 = Icges(2,5) * t151 - Icges(2,6) * t150 + t212 * t137 + t213 * t139;
t147 = sin(pkin(8));
t148 = cos(pkin(8));
t193 = Icges(5,4) * t147;
t164 = Icges(5,2) * t148 + t193;
t78 = Icges(5,6) * t137 - t164 * t139;
t192 = Icges(5,4) * t148;
t166 = Icges(5,1) * t147 + t192;
t80 = Icges(5,5) * t137 - t166 * t139;
t204 = t211 * t139 + t147 * t80 + t148 * t78 - t209;
t77 = Icges(5,6) * t139 + t164 * t137;
t79 = Icges(5,5) * t139 + t166 * t137;
t203 = t211 * t137 - t147 * t79 - t148 * t77 + t210;
t145 = pkin(8) + qJ(5);
t138 = cos(t145);
t136 = sin(t145);
t191 = Icges(6,4) * t136;
t101 = Icges(6,1) * t138 - t191;
t118 = qJD(5) * t137 + V_base(5);
t119 = qJD(5) * t139 + V_base(4);
t140 = V_base(6) + qJD(1);
t163 = Icges(6,2) * t138 + t191;
t68 = Icges(6,6) * t139 + t163 * t137;
t69 = Icges(6,6) * t137 - t163 * t139;
t190 = Icges(6,4) * t138;
t165 = Icges(6,1) * t136 + t190;
t70 = Icges(6,5) * t139 + t165 * t137;
t71 = Icges(6,5) * t137 - t165 * t139;
t96 = -Icges(6,2) * t136 + t190;
t202 = (t136 * t70 + t138 * t68) * t119 + (t136 * t71 + t138 * t69) * t118 + (t101 * t136 + t138 * t96) * t140;
t200 = pkin(1) * t150;
t199 = pkin(1) * t151;
t198 = pkin(4) * t147;
t197 = pkin(4) * t148;
t196 = -pkin(5) - qJ(2);
t195 = Icges(2,4) * t150;
t187 = qJ(4) * t137;
t186 = qJ(4) * t139;
t184 = -pkin(3) + t196;
t183 = t140 * t199 + V_base(2);
t182 = V_base(5) * pkin(5) + V_base(1);
t104 = pkin(2) * t137 - qJ(3) * t139;
t179 = -t104 - t200;
t108 = pkin(2) * t139 + qJ(3) * t137;
t178 = -t108 - t199;
t177 = V_base(5) * qJ(2) + t182;
t176 = V_base(4) * t200 + qJD(2) + V_base(3);
t175 = V_base(4) * t104 + t176;
t174 = qJD(3) * t137 + t177;
t173 = rSges(5,1) * t147 + rSges(5,2) * t148;
t172 = rSges(6,1) * t136 + rSges(6,2) * t138;
t162 = Icges(5,5) * t147 + Icges(5,6) * t148;
t161 = Icges(6,5) * t136 + Icges(6,6) * t138;
t116 = -Icges(5,2) * t147 + t192;
t117 = Icges(5,1) * t148 - t193;
t160 = t116 * t148 + t117 * t147;
t159 = V_base(4) * t187 + t175;
t158 = t179 - t187;
t157 = t178 - t186;
t156 = -qJD(3) * t139 + t140 * t108 + t183;
t155 = V_base(5) * pkin(3) + qJD(4) * t139 + t174;
t154 = (Icges(6,3) * t137 - t161 * t139) * t118 + (Icges(6,3) * t139 + t161 * t137) * t119 + (Icges(6,5) * t138 - Icges(6,6) * t136) * t140;
t153 = qJD(4) * t137 + t140 * t186 + t156;
t152 = (Icges(5,5) * t148 - Icges(5,6) * t147) * t140 + (Icges(5,3) * t139 + t162 * t137) * V_base(4) + (Icges(5,3) * t137 - t162 * t139) * V_base(5);
t142 = Icges(2,4) * t151;
t128 = rSges(2,1) * t151 - t150 * rSges(2,2);
t127 = t150 * rSges(2,1) + rSges(2,2) * t151;
t126 = Icges(2,1) * t151 - t195;
t125 = Icges(2,1) * t150 + t142;
t124 = -Icges(2,2) * t150 + t142;
t123 = Icges(2,2) * t151 + t195;
t120 = rSges(5,1) * t148 - rSges(5,2) * t147;
t114 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t113 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t112 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t110 = rSges(3,1) * t139 - rSges(3,2) * t137;
t109 = -rSges(4,2) * t139 + rSges(4,3) * t137;
t107 = rSges(6,1) * t138 - rSges(6,2) * t136;
t106 = rSges(3,1) * t137 + rSges(3,2) * t139;
t105 = -rSges(4,2) * t137 - rSges(4,3) * t139;
t86 = pkin(6) * t139 + t137 * t198;
t85 = pkin(6) * t137 - t139 * t198;
t84 = V_base(5) * rSges(2,3) - t127 * t140 + t182;
t83 = t128 * t140 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t82 = rSges(5,3) * t137 - t173 * t139;
t81 = rSges(5,3) * t139 + t173 * t137;
t74 = t127 * V_base(4) - t128 * V_base(5) + V_base(3);
t73 = rSges(6,3) * t137 - t172 * t139;
t72 = rSges(6,3) * t139 + t172 * t137;
t65 = V_base(5) * rSges(3,3) + (-t106 - t200) * t140 + t177;
t64 = t110 * t140 + (-rSges(3,3) + t196) * V_base(4) + t183;
t63 = V_base(4) * t106 + (-t110 - t199) * V_base(5) + t176;
t62 = V_base(5) * rSges(4,1) + (-t105 + t179) * t140 + t174;
t61 = t109 * t140 + (-rSges(4,1) + t196) * V_base(4) + t156;
t60 = V_base(4) * t105 + (-t109 + t178) * V_base(5) + t175;
t59 = t120 * V_base(5) + (t158 - t82) * t140 + t155;
t58 = t140 * t81 + (-t120 + t184) * V_base(4) + t153;
t57 = V_base(4) * t82 + (t157 - t81) * V_base(5) + t159;
t56 = V_base(5) * t197 + t107 * t118 + (t158 - t73 - t85) * t140 + t155;
t55 = -t107 * t119 + (t72 + t86) * t140 + (t184 - t197) * V_base(4) + t153;
t54 = -t118 * t72 + t119 * t73 + V_base(4) * t85 + (t157 - t86) * V_base(5) + t159;
t1 = m(1) * (t112 ^ 2 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + m(2) * (t74 ^ 2 + t83 ^ 2 + t84 ^ 2) / 0.2e1 + m(3) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + m(4) * (t60 ^ 2 + t61 ^ 2 + t62 ^ 2) / 0.2e1 + m(5) * (t57 ^ 2 + t58 ^ 2 + t59 ^ 2) / 0.2e1 + m(6) * (t54 ^ 2 + t55 ^ 2 + t56 ^ 2) / 0.2e1 + t119 * (t202 * t137 + t154 * t139) / 0.2e1 + t118 * (t154 * t137 - t202 * t139) / 0.2e1 + ((-t136 * t68 + t138 * t70) * t119 + (-t136 * t69 + t138 * t71) * t118 + (-t147 * t78 + t148 * t80 + t206) * V_base(5) + (-t147 * t77 + t148 * t79 + t205) * V_base(4) + (t101 * t138 - t116 * t147 + t117 * t148 - t136 * t96 + Icges(4,1) + Icges(2,3) + Icges(3,3)) * t140) * t140 / 0.2e1 + (t152 * t139 + (t160 * t137 + t205) * t140 + (-t150 * t123 + t125 * t151 + t204 * t137 + t208 * t139 + Icges(1,4)) * V_base(5) + (-t150 * t124 + t126 * t151 - t203 * t137 + t207 * t139 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t152 * t137 + (-t160 * t139 + t206) * t140 + (t123 * t151 + t150 * t125 + t208 * t137 - t204 * t139 + Icges(1,2)) * V_base(5) + (t124 * t151 + t150 * t126 + t207 * t137 + t203 * t139 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
