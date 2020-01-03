% Calculate kinetic energy for
% S5RPRPR10
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR10_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR10_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRPR10_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR10_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR10_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR10_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:25:52
% EndTime: 2019-12-31 18:25:54
% DurationCPUTime: 1.54s
% Computational Cost: add. (906->212), mult. (1000->261), div. (0->0), fcn. (962->8), ass. (0->107)
t218 = Icges(2,4) - Icges(3,5);
t217 = Icges(2,1) + Icges(3,1);
t216 = Icges(3,4) + Icges(2,5);
t215 = Icges(2,2) + Icges(3,3);
t214 = Icges(2,6) - Icges(3,6);
t200 = sin(qJ(1));
t213 = t218 * t200;
t165 = cos(qJ(1));
t212 = t218 * t165;
t209 = -t165 * t215 - t213;
t208 = t200 * t215 - t212;
t205 = t200 * t217 + t212;
t204 = t165 * t217 - t213;
t201 = cos(qJ(3));
t199 = sin(qJ(3));
t198 = t165 * pkin(2);
t197 = -pkin(6) - qJ(4);
t178 = t200 * t199;
t120 = -t165 * t201 - t178;
t196 = Icges(4,4) * t120;
t191 = qJ(3) + pkin(8);
t181 = sin(t191);
t182 = cos(t191);
t115 = -t165 * t182 - t181 * t200;
t195 = Icges(5,4) * t115;
t163 = sin(qJ(5));
t194 = Icges(6,4) * t163;
t164 = cos(qJ(5));
t193 = Icges(6,4) * t164;
t141 = pkin(1) * t200 - t165 * qJ(2);
t190 = V_base(4) * t141 + V_base(3);
t189 = V_base(5) * pkin(5) + V_base(1);
t186 = t200 * pkin(2);
t184 = t165 * t199;
t156 = V_base(6) + qJD(1);
t144 = t165 * pkin(1) + qJ(2) * t200;
t183 = -t144 - t198;
t180 = V_base(4) * t186 + t190;
t179 = qJD(2) * t200 + t189;
t155 = pkin(3) * t201 + pkin(2);
t113 = pkin(3) * t178 + (-pkin(2) + t155) * t165;
t177 = -t113 + t183;
t176 = -rSges(6,1) * t164 + rSges(6,2) * t163;
t175 = -Icges(6,1) * t164 + t194;
t174 = Icges(6,2) * t163 - t193;
t173 = -Icges(6,5) * t164 + Icges(6,6) * t163;
t172 = -qJD(2) * t165 + t156 * t144 + V_base(2);
t112 = -pkin(3) * t184 + t155 * t200 - t186;
t171 = V_base(4) * t112 - qJD(4) + t180;
t170 = V_base(4) * pkin(6) + t156 * t198 + t172;
t110 = -qJD(5) * t115 + V_base(5);
t116 = t165 * t181 - t182 * t200;
t111 = qJD(5) * t116 + V_base(4);
t152 = -qJD(3) + t156;
t169 = (-Icges(6,3) * t115 + t173 * t116) * t110 + (Icges(6,3) * t116 + t173 * t115) * t111 + (-Icges(6,5) * t163 - Icges(6,6) * t164) * t152;
t168 = V_base(4) * qJ(4) + t152 * t113 + t170;
t167 = (-t186 - t141) * t156 + t179;
t130 = -Icges(6,2) * t164 - t194;
t135 = -Icges(6,1) * t163 - t193;
t78 = -Icges(6,6) * t115 + t174 * t116;
t79 = Icges(6,6) * t116 + t174 * t115;
t80 = -Icges(6,5) * t115 + t175 * t116;
t81 = Icges(6,5) * t116 + t175 * t115;
t166 = (t163 * t79 - t164 * t81) * t111 + (t163 * t78 - t164 * t80) * t110 + (t130 * t163 - t135 * t164) * t152;
t146 = t165 * rSges(2,1) - rSges(2,2) * t200;
t145 = t165 * rSges(3,1) + rSges(3,3) * t200;
t143 = rSges(2,1) * t200 + t165 * rSges(2,2);
t142 = rSges(3,1) * t200 - t165 * rSges(3,3);
t140 = -rSges(6,1) * t163 - rSges(6,2) * t164;
t124 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t123 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t122 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t121 = t200 * t201 - t184;
t118 = Icges(4,4) * t121;
t114 = Icges(5,4) * t116;
t108 = V_base(5) * rSges(2,3) - t143 * t156 + t189;
t107 = t146 * t156 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t105 = t143 * V_base(4) - t146 * V_base(5) + V_base(3);
t104 = -rSges(4,1) * t120 + rSges(4,2) * t121;
t103 = rSges(4,1) * t121 + rSges(4,2) * t120;
t102 = -Icges(4,1) * t120 + t118;
t101 = Icges(4,1) * t121 + t196;
t100 = Icges(4,2) * t121 - t196;
t99 = Icges(4,2) * t120 + t118;
t96 = -pkin(4) * t115 + pkin(7) * t116;
t95 = -pkin(4) * t116 - pkin(7) * t115;
t94 = -rSges(5,1) * t115 - rSges(5,2) * t116;
t93 = -rSges(5,1) * t116 + rSges(5,2) * t115;
t92 = -Icges(5,1) * t115 - t114;
t91 = -Icges(5,1) * t116 + t195;
t90 = -Icges(5,2) * t116 - t195;
t89 = Icges(5,2) * t115 - t114;
t86 = V_base(5) * rSges(3,2) + (-t141 - t142) * t156 + t179;
t85 = t156 * t145 + (-rSges(3,2) - pkin(5)) * V_base(4) + t172;
t84 = t142 * V_base(4) + (-t144 - t145) * V_base(5) + t190;
t83 = rSges(6,3) * t116 + t176 * t115;
t82 = -rSges(6,3) * t115 + t176 * t116;
t75 = -t152 * t103 + (-pkin(6) - rSges(4,3)) * V_base(5) + t167;
t74 = t152 * t104 + (rSges(4,3) - pkin(5)) * V_base(4) + t170;
t73 = V_base(4) * t103 + (-t104 + t183) * V_base(5) + t180;
t72 = (-t112 - t93) * t152 + (-rSges(5,3) + t197) * V_base(5) + t167;
t71 = t152 * t94 + (rSges(5,3) - pkin(5)) * V_base(4) + t168;
t70 = V_base(4) * t93 + (t177 - t94) * V_base(5) + t171;
t69 = t110 * t140 + t197 * V_base(5) + (-t112 - t82 - t95) * t152 + t167;
t68 = -V_base(4) * pkin(5) - t111 * t140 + (t83 + t96) * t152 + t168;
t67 = -t110 * t83 + t111 * t82 + V_base(4) * t95 + (t177 - t96) * V_base(5) + t171;
t1 = m(1) * (t122 ^ 2 + t123 ^ 2 + t124 ^ 2) / 0.2e1 + m(2) * (t105 ^ 2 + t107 ^ 2 + t108 ^ 2) / 0.2e1 + m(3) * (t84 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + m(4) * (t73 ^ 2 + t74 ^ 2 + t75 ^ 2) / 0.2e1 + m(5) * (t70 ^ 2 + t71 ^ 2 + t72 ^ 2) / 0.2e1 + m(6) * (t67 ^ 2 + t68 ^ 2 + t69 ^ 2) / 0.2e1 + t111 * (t166 * t115 + t169 * t116) / 0.2e1 + t110 * (-t169 * t115 + t166 * t116) / 0.2e1 + ((-t163 * t81 - t164 * t79) * t111 + (-t163 * t80 - t164 * t78) * t110 + (-t164 * t130 - t163 * t135 + Icges(4,3) + Icges(5,3)) * t152) * t152 / 0.2e1 + ((-t101 * t120 - t115 * t91 - t116 * t89 + t121 * t99 + t205 * t165 + t209 * t200 + Icges(1,4)) * V_base(5) + (t100 * t121 - t102 * t120 - t115 * t92 - t116 * t90 + t204 * t165 + t208 * t200 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t101 * t121 + t115 * t89 - t116 * t91 + t120 * t99 - t209 * t165 + t205 * t200 + Icges(1,2)) * V_base(5) + (t100 * t120 + t102 * t121 + t115 * t90 - t116 * t92 - t208 * t165 + t204 * t200 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t152 * (-Icges(4,5) * t121 + Icges(5,5) * t116 - Icges(4,6) * t120 - Icges(5,6) * t115) + V_base(4) * t152 * (Icges(4,5) * t120 + Icges(5,5) * t115 - Icges(4,6) * t121 + Icges(5,6) * t116) + ((Icges(2,3) / 0.2e1 + Icges(3,2) / 0.2e1) * t156 + (-t214 * V_base(4) + t216 * V_base(5)) * t200 + (t214 * V_base(5) + t216 * V_base(4)) * t165) * t156 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
