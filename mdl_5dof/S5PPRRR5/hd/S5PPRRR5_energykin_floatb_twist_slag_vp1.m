% Calculate kinetic energy for
% S5PPRRR5
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
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
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
% Datum: 2019-12-31 17:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRRR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR5_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PPRRR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR5_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRR5_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRR5_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:35:36
% EndTime: 2019-12-31 17:35:37
% DurationCPUTime: 1.66s
% Computational Cost: add. (873->211), mult. (1000->265), div. (0->0), fcn. (962->8), ass. (0->109)
t218 = Icges(2,4) - Icges(3,5);
t217 = Icges(2,1) + Icges(3,1);
t216 = Icges(3,4) + Icges(2,5);
t215 = Icges(2,2) + Icges(3,3);
t214 = Icges(2,6) - Icges(3,6);
t200 = sin(pkin(8));
t213 = t218 * t200;
t163 = cos(pkin(8));
t212 = t218 * t163;
t211 = -t215 * t163 - t213;
t210 = t215 * t200 - t212;
t209 = t217 * t200 + t212;
t208 = t217 * t163 - t213;
t205 = -pkin(5) - pkin(6);
t203 = cos(qJ(3));
t202 = sin(qJ(3));
t201 = t163 * pkin(2);
t179 = t200 * t202;
t120 = -t163 * t203 - t179;
t199 = Icges(4,4) * t120;
t194 = qJ(3) + qJ(4);
t185 = sin(t194);
t186 = cos(t194);
t115 = -t163 * t186 - t185 * t200;
t198 = Icges(5,4) * t115;
t164 = sin(qJ(5));
t197 = Icges(6,4) * t164;
t165 = cos(qJ(5));
t196 = Icges(6,4) * t165;
t193 = V_base(5) * qJ(1) + V_base(1);
t189 = qJD(1) + V_base(3);
t188 = t163 * t202;
t187 = t200 * pkin(2);
t159 = V_base(6) - qJD(3);
t140 = t163 * pkin(1) + qJ(2) * t200;
t184 = -t140 - t201;
t182 = qJD(2) * t200 + t193;
t137 = pkin(1) * t200 - t163 * qJ(2);
t181 = V_base(4) * t137 + t189;
t155 = pkin(3) * t203 + pkin(2);
t113 = pkin(3) * t179 + (-pkin(2) + t155) * t163;
t180 = -t113 + t184;
t178 = V_base(4) * t187 + t181;
t177 = -rSges(6,1) * t165 + rSges(6,2) * t164;
t176 = -Icges(6,1) * t165 + t197;
t175 = Icges(6,2) * t164 - t196;
t174 = -Icges(6,5) * t165 + Icges(6,6) * t164;
t173 = -qJD(2) * t163 + V_base(6) * t140 + V_base(2);
t112 = -pkin(3) * t188 + t155 * t200 - t187;
t172 = V_base(4) * t112 + t178;
t171 = V_base(4) * pkin(5) + V_base(6) * t201 + t173;
t110 = -qJD(5) * t115 + V_base(5);
t116 = t163 * t185 - t186 * t200;
t111 = qJD(5) * t116 + V_base(4);
t154 = -qJD(4) + t159;
t170 = (-Icges(6,3) * t115 + t116 * t174) * t110 + (Icges(6,3) * t116 + t115 * t174) * t111 + (-Icges(6,5) * t164 - Icges(6,6) * t165) * t154;
t169 = V_base(4) * pkin(6) + t159 * t113 + t171;
t168 = (-t187 - t137) * V_base(6) + t182;
t167 = -t159 * t112 + t168;
t144 = -Icges(6,2) * t165 - t197;
t145 = -Icges(6,1) * t164 - t196;
t78 = -Icges(6,6) * t115 + t116 * t175;
t79 = Icges(6,6) * t116 + t115 * t175;
t80 = -Icges(6,5) * t115 + t116 * t176;
t81 = Icges(6,5) * t116 + t115 * t176;
t166 = (t164 * t79 - t165 * t81) * t111 + (t164 * t78 - t165 * t80) * t110 + (t144 * t164 - t145 * t165) * t154;
t146 = -t164 * rSges(6,1) - rSges(6,2) * t165;
t142 = t163 * rSges(2,1) - rSges(2,2) * t200;
t141 = t163 * rSges(3,1) + rSges(3,3) * t200;
t139 = rSges(2,1) * t200 + t163 * rSges(2,2);
t138 = rSges(3,1) * t200 - t163 * rSges(3,3);
t124 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t123 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t122 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t121 = t200 * t203 - t188;
t117 = Icges(4,4) * t121;
t114 = Icges(5,4) * t116;
t109 = V_base(5) * rSges(2,3) - t139 * V_base(6) + t193;
t108 = t142 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t105 = -rSges(4,1) * t120 + rSges(4,2) * t121;
t104 = rSges(4,1) * t121 + rSges(4,2) * t120;
t103 = t139 * V_base(4) - t142 * V_base(5) + t189;
t102 = -Icges(4,1) * t120 + t117;
t101 = Icges(4,1) * t121 + t199;
t100 = Icges(4,2) * t121 - t199;
t99 = Icges(4,2) * t120 + t117;
t96 = -pkin(4) * t115 + pkin(7) * t116;
t95 = -pkin(4) * t116 - pkin(7) * t115;
t94 = -rSges(5,1) * t115 - rSges(5,2) * t116;
t93 = -rSges(5,1) * t116 + rSges(5,2) * t115;
t92 = -Icges(5,1) * t115 - t114;
t91 = -Icges(5,1) * t116 + t198;
t90 = -Icges(5,2) * t116 - t198;
t89 = Icges(5,2) * t115 - t114;
t86 = V_base(5) * rSges(3,2) + (-t137 - t138) * V_base(6) + t182;
t85 = t141 * V_base(6) + (-rSges(3,2) - qJ(1)) * V_base(4) + t173;
t84 = t138 * V_base(4) + (-t140 - t141) * V_base(5) + t181;
t83 = t116 * rSges(6,3) + t115 * t177;
t82 = -t115 * rSges(6,3) + t116 * t177;
t75 = -t159 * t104 + (-pkin(5) - rSges(4,3)) * V_base(5) + t168;
t74 = t105 * t159 + (rSges(4,3) - qJ(1)) * V_base(4) + t171;
t73 = t104 * V_base(4) + (-t105 + t184) * V_base(5) + t178;
t72 = -t154 * t93 + (-rSges(5,3) + t205) * V_base(5) + t167;
t71 = t154 * t94 + (rSges(5,3) - qJ(1)) * V_base(4) + t169;
t70 = t93 * V_base(4) + (t180 - t94) * V_base(5) + t172;
t69 = t110 * t146 + t205 * V_base(5) + (-t82 - t95) * t154 + t167;
t68 = -V_base(4) * qJ(1) - t111 * t146 + (t83 + t96) * t154 + t169;
t67 = -t110 * t83 + t111 * t82 + t95 * V_base(4) + (t180 - t96) * V_base(5) + t172;
t1 = m(1) * (t122 ^ 2 + t123 ^ 2 + t124 ^ 2) / 0.2e1 + m(2) * (t103 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + m(3) * (t84 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + m(4) * (t73 ^ 2 + t74 ^ 2 + t75 ^ 2) / 0.2e1 + m(5) * (t70 ^ 2 + t71 ^ 2 + t72 ^ 2) / 0.2e1 + m(6) * (t67 ^ 2 + t68 ^ 2 + t69 ^ 2) / 0.2e1 + t111 * (t166 * t115 + t170 * t116) / 0.2e1 + t110 * (-t170 * t115 + t166 * t116) / 0.2e1 + ((-t164 * t81 - t165 * t79) * t111 + (-t164 * t80 - t165 * t78) * t110 + (-t165 * t144 - t164 * t145 + Icges(5,3)) * t154) * t154 / 0.2e1 + ((-t101 * t120 - t115 * t91 - t116 * t89 + t121 * t99 + t209 * t163 + t200 * t211 + Icges(1,4)) * V_base(5) + (t100 * t121 - t102 * t120 - t115 * t92 - t116 * t90 + t208 * t163 + t210 * t200 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t101 * t121 + t115 * t89 - t116 * t91 + t120 * t99 - t211 * t163 + t209 * t200 + Icges(1,2)) * V_base(5) + (t100 * t120 + t102 * t121 + t115 * t90 - t116 * t92 - t163 * t210 + t200 * t208 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t154 * (Icges(5,5) * t115 + Icges(5,6) * t116) + V_base(5) * t154 * (Icges(5,5) * t116 - Icges(5,6) * t115) + (Icges(1,6) * V_base(5) + Icges(1,5) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1 + Icges(3,2) / 0.2e1) * V_base(6) + (-t214 * V_base(4) + t216 * V_base(5)) * t200 + (t214 * V_base(5) + t216 * V_base(4)) * t163) * V_base(6) + ((-Icges(4,5) * t121 - Icges(4,6) * t120) * V_base(5) + (Icges(4,5) * t120 - Icges(4,6) * t121) * V_base(4) + Icges(4,3) * t159 / 0.2e1) * t159;
T = t1;
