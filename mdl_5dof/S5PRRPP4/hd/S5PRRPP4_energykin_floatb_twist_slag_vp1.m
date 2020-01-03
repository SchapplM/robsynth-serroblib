% Calculate kinetic energy for
% S5PRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
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
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRPP4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP4_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRRPP4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP4_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPP4_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPP4_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:40:43
% EndTime: 2019-12-31 17:40:45
% DurationCPUTime: 1.73s
% Computational Cost: add. (907->186), mult. (922->238), div. (0->0), fcn. (702->6), ass. (0->99)
t243 = Icges(4,4) - Icges(6,4) - Icges(5,5);
t242 = Icges(4,1) + Icges(5,1) + Icges(6,1);
t241 = Icges(4,2) + Icges(6,2) + Icges(5,3);
t167 = cos(qJ(3));
t240 = t243 * t167;
t166 = sin(qJ(3));
t239 = t243 * t166;
t238 = Icges(5,4) + Icges(4,5) - Icges(6,5);
t237 = Icges(4,6) - Icges(5,6) + Icges(6,6);
t236 = t241 * t166 - t240;
t235 = t242 * t167 - t239;
t234 = rSges(6,1) + pkin(4);
t163 = pkin(7) + qJ(2);
t157 = sin(t163);
t158 = cos(t163);
t233 = t236 * t157 + t237 * t158;
t232 = -t237 * t157 + t236 * t158;
t231 = t235 * t157 - t238 * t158;
t230 = t238 * t157 + t235 * t158;
t229 = Icges(5,2) + Icges(4,3) + Icges(6,3);
t228 = -t241 * t167 - t239;
t227 = t242 * t166 + t240;
t226 = -t237 * t166 + t238 * t167;
t225 = rSges(6,3) + qJ(5);
t224 = rSges(6,2) * t166 + t234 * t167;
t135 = -qJD(3) * t158 + V_base(5);
t136 = qJD(3) * t157 + V_base(4);
t160 = V_base(6) + qJD(2);
t221 = (t228 * t166 + t227 * t167) * t160 + (t232 * t166 + t230 * t167) * t136 + (t233 * t166 + t231 * t167) * t135;
t220 = (t238 * t166 + t237 * t167) * t160 + (t229 * t157 + t226 * t158) * t136 + (t226 * t157 - t229 * t158) * t135;
t165 = cos(pkin(7));
t216 = pkin(1) * t165;
t214 = -pkin(5) - qJ(1);
t164 = sin(pkin(7));
t213 = Icges(2,4) * t164;
t212 = Icges(3,4) * t157;
t205 = t224 * t157 + t225 * t158;
t204 = -t225 * t157 + t224 * t158;
t189 = pkin(3) * t167 + qJ(4) * t166;
t110 = t189 * t157;
t124 = t157 * pkin(2) - t158 * pkin(6);
t203 = -t110 - t124;
t202 = qJD(4) * t166;
t195 = pkin(1) * V_base(6);
t201 = t165 * t195 + V_base(2);
t200 = V_base(5) * qJ(1) + V_base(1);
t196 = qJD(1) + V_base(3);
t194 = -t167 * rSges(6,2) + t234 * t166;
t193 = V_base(4) * t164 * pkin(1) + t196;
t192 = rSges(4,1) * t167 - rSges(4,2) * t166;
t191 = rSges(5,1) * t167 + rSges(5,3) * t166;
t176 = V_base(5) * pkin(5) - t164 * t195 + t200;
t125 = t158 * pkin(2) + t157 * pkin(6);
t175 = t160 * t125 + t214 * V_base(4) + t201;
t150 = t166 * pkin(3) - t167 * qJ(4);
t174 = t135 * t150 + t158 * t202 + t176;
t111 = t189 * t158;
t173 = t160 * t111 + t157 * t202 + t175;
t172 = V_base(4) * t124 + (-t125 - t216) * V_base(5) + t193;
t171 = -qJD(4) * t167 + t136 * t110 + t172;
t159 = Icges(2,4) * t165;
t156 = Icges(3,4) * t158;
t153 = t166 * rSges(4,1) + t167 * rSges(4,2);
t152 = t166 * rSges(5,1) - t167 * rSges(5,3);
t138 = t165 * rSges(2,1) - t164 * rSges(2,2);
t137 = t164 * rSges(2,1) + t165 * rSges(2,2);
t134 = Icges(2,1) * t165 - t213;
t133 = Icges(2,1) * t164 + t159;
t132 = -Icges(2,2) * t164 + t159;
t131 = Icges(2,2) * t165 + t213;
t128 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t127 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t126 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t123 = t158 * rSges(3,1) - t157 * rSges(3,2);
t122 = t157 * rSges(3,1) + t158 * rSges(3,2);
t121 = Icges(3,1) * t158 - t212;
t120 = Icges(3,1) * t157 + t156;
t119 = -Icges(3,2) * t157 + t156;
t118 = Icges(3,2) * t158 + t212;
t107 = V_base(5) * rSges(2,3) - V_base(6) * t137 + t200;
t106 = V_base(6) * t138 + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t105 = t157 * rSges(4,3) + t192 * t158;
t104 = t157 * rSges(5,2) + t191 * t158;
t102 = -t158 * rSges(4,3) + t192 * t157;
t101 = -t158 * rSges(5,2) + t191 * t157;
t80 = V_base(4) * t137 - V_base(5) * t138 + t196;
t79 = V_base(5) * rSges(3,3) - t160 * t122 + t176;
t78 = t160 * t123 + (-rSges(3,3) + t214) * V_base(4) + t201;
t77 = V_base(4) * t122 + (-t123 - t216) * V_base(5) + t193;
t76 = t135 * t153 + (-t102 - t124) * t160 + t176;
t75 = t160 * t105 - t136 * t153 + t175;
t74 = t136 * t102 - t135 * t105 + t172;
t73 = t135 * t152 + (-t101 + t203) * t160 + t174;
t72 = t160 * t104 + (-t150 - t152) * t136 + t173;
t71 = -qJD(5) * t157 + t194 * t135 + (t203 - t205) * t160 + t174;
t70 = qJD(5) * t158 + t204 * t160 + (-t150 - t194) * t136 + t173;
t69 = t136 * t101 + (-t104 - t111) * t135 + t171;
t68 = t205 * t136 + (-t111 - t204) * t135 + t171;
t1 = m(1) * (t126 ^ 2 + t127 ^ 2 + t128 ^ 2) / 0.2e1 + m(2) * (t106 ^ 2 + t107 ^ 2 + t80 ^ 2) / 0.2e1 + m(3) * (t77 ^ 2 + t78 ^ 2 + t79 ^ 2) / 0.2e1 + m(4) * (t74 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + m(5) * (t69 ^ 2 + t72 ^ 2 + t73 ^ 2) / 0.2e1 + m(6) * (t68 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + (t221 * t157 - t220 * t158) * t135 / 0.2e1 + (t220 * t157 + t221 * t158) * t136 / 0.2e1 + ((-t157 * t118 + t158 * t120 - t164 * t131 + t165 * t133 + Icges(1,4)) * V_base(5) + (-t157 * t119 + t158 * t121 - t164 * t132 + t165 * t134 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t158 * t118 + t157 * t120 + t165 * t131 + t164 * t133 + Icges(1,2)) * V_base(5) + (t158 * t119 + t157 * t121 + t165 * t132 + t164 * t134 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t230 * t166 - t232 * t167) * t136 + (t231 * t166 - t233 * t167) * t135 + (t227 * t166 - t228 * t167 + Icges(3,3)) * t160) * t160 / 0.2e1 + V_base(4) * t160 * (Icges(3,5) * t158 - Icges(3,6) * t157) + V_base(5) * t160 * (Icges(3,5) * t157 + Icges(3,6) * t158) + ((Icges(2,5) * t164 + Icges(2,6) * t165 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t165 - Icges(2,6) * t164 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T = t1;
