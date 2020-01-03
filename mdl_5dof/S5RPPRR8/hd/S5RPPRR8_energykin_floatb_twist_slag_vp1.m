% Calculate kinetic energy for
% S5RPPRR8
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR8_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR8_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPRR8_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR8_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR8_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR8_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:01:01
% EndTime: 2019-12-31 18:01:03
% DurationCPUTime: 1.67s
% Computational Cost: add. (897->212), mult. (1000->264), div. (0->0), fcn. (962->8), ass. (0->108)
t217 = Icges(2,4) - Icges(3,5);
t216 = Icges(2,1) + Icges(3,1);
t215 = Icges(3,4) + Icges(2,5);
t214 = Icges(2,2) + Icges(3,3);
t213 = Icges(2,6) - Icges(3,6);
t202 = sin(qJ(1));
t212 = t217 * t202;
t165 = cos(qJ(1));
t211 = t217 * t165;
t210 = -t214 * t165 - t212;
t209 = t214 * t202 - t211;
t208 = t216 * t202 + t211;
t207 = t216 * t165 - t212;
t201 = t165 * pkin(2);
t200 = -pkin(6) - qJ(3);
t199 = cos(pkin(8));
t198 = sin(pkin(8));
t178 = t202 * t198;
t119 = -t165 * t199 - t178;
t197 = Icges(4,4) * t119;
t192 = pkin(8) + qJ(4);
t182 = sin(t192);
t183 = cos(t192);
t115 = -t165 * t183 - t182 * t202;
t196 = Icges(5,4) * t115;
t163 = sin(qJ(5));
t195 = Icges(6,4) * t163;
t164 = cos(qJ(5));
t194 = Icges(6,4) * t164;
t141 = pkin(1) * t202 - t165 * qJ(2);
t191 = V_base(4) * t141 + V_base(3);
t190 = V_base(5) * pkin(5) + V_base(1);
t187 = t202 * pkin(2);
t156 = V_base(6) + qJD(1);
t144 = t165 * pkin(1) + qJ(2) * t202;
t185 = -t144 - t201;
t184 = t165 * t198;
t181 = qJD(2) * t202 + t190;
t155 = pkin(3) * t199 + pkin(2);
t113 = pkin(3) * t178 + (-pkin(2) + t155) * t165;
t180 = -t113 + t185;
t179 = -t141 - t187;
t177 = V_base(4) * t187 - qJD(3) + t191;
t176 = -rSges(6,1) * t164 + rSges(6,2) * t163;
t175 = -Icges(6,1) * t164 + t195;
t174 = Icges(6,2) * t163 - t194;
t173 = -Icges(6,5) * t164 + Icges(6,6) * t163;
t172 = -qJD(2) * t165 + t156 * t144 + V_base(2);
t112 = -pkin(3) * t184 + t155 * t202 - t187;
t171 = V_base(4) * t112 + t177;
t170 = V_base(4) * qJ(3) + t156 * t201 + t172;
t110 = -qJD(5) * t115 + V_base(5);
t116 = t165 * t182 - t183 * t202;
t111 = qJD(5) * t116 + V_base(4);
t154 = -qJD(4) + t156;
t169 = (-Icges(6,3) * t115 + t116 * t173) * t110 + (Icges(6,3) * t116 + t115 * t173) * t111 + (-Icges(6,5) * t163 - Icges(6,6) * t164) * t154;
t168 = V_base(4) * pkin(6) + t156 * t113 + t170;
t167 = (-t112 + t179) * t156 + t181;
t130 = -Icges(6,2) * t164 - t195;
t135 = -Icges(6,1) * t163 - t194;
t78 = -Icges(6,6) * t115 + t116 * t174;
t79 = Icges(6,6) * t116 + t115 * t174;
t80 = -Icges(6,5) * t115 + t116 * t175;
t81 = Icges(6,5) * t116 + t115 * t175;
t166 = (t163 * t79 - t164 * t81) * t111 + (t163 * t78 - t164 * t80) * t110 + (t130 * t163 - t135 * t164) * t154;
t146 = t165 * rSges(2,1) - rSges(2,2) * t202;
t145 = t165 * rSges(3,1) + rSges(3,3) * t202;
t143 = rSges(2,1) * t202 + t165 * rSges(2,2);
t142 = rSges(3,1) * t202 - t165 * rSges(3,3);
t140 = -rSges(6,1) * t163 - rSges(6,2) * t164;
t124 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t123 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t122 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t120 = t199 * t202 - t184;
t117 = Icges(4,4) * t120;
t114 = Icges(5,4) * t116;
t108 = V_base(5) * rSges(2,3) - t143 * t156 + t190;
t107 = t146 * t156 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t105 = t143 * V_base(4) - t146 * V_base(5) + V_base(3);
t104 = -rSges(4,1) * t119 + rSges(4,2) * t120;
t103 = rSges(4,1) * t120 + rSges(4,2) * t119;
t102 = -Icges(4,1) * t119 + t117;
t101 = Icges(4,1) * t120 + t197;
t100 = Icges(4,2) * t120 - t197;
t99 = Icges(4,2) * t119 + t117;
t96 = -pkin(4) * t115 + pkin(7) * t116;
t95 = -pkin(4) * t116 - pkin(7) * t115;
t94 = -rSges(5,1) * t115 - rSges(5,2) * t116;
t93 = -rSges(5,1) * t116 + rSges(5,2) * t115;
t92 = -Icges(5,1) * t115 - t114;
t91 = -Icges(5,1) * t116 + t196;
t90 = -Icges(5,2) * t116 - t196;
t89 = Icges(5,2) * t115 - t114;
t86 = V_base(5) * rSges(3,2) + (-t141 - t142) * t156 + t181;
t85 = t156 * t145 + (-rSges(3,2) - pkin(5)) * V_base(4) + t172;
t84 = t142 * V_base(4) + (-t144 - t145) * V_base(5) + t191;
t83 = rSges(6,3) * t116 + t115 * t176;
t82 = -rSges(6,3) * t115 + t116 * t176;
t75 = (-qJ(3) - rSges(4,3)) * V_base(5) + (-t103 + t179) * t156 + t181;
t74 = t156 * t104 + (rSges(4,3) - pkin(5)) * V_base(4) + t170;
t73 = V_base(4) * t103 + (-t104 + t185) * V_base(5) + t177;
t72 = -t154 * t93 + (-rSges(5,3) + t200) * V_base(5) + t167;
t71 = t154 * t94 + (rSges(5,3) - pkin(5)) * V_base(4) + t168;
t70 = V_base(4) * t93 + (t180 - t94) * V_base(5) + t171;
t69 = t110 * t140 + t200 * V_base(5) + (-t82 - t95) * t154 + t167;
t68 = -V_base(4) * pkin(5) - t111 * t140 + (t83 + t96) * t154 + t168;
t67 = -t110 * t83 + t111 * t82 + V_base(4) * t95 + (t180 - t96) * V_base(5) + t171;
t1 = m(1) * (t122 ^ 2 + t123 ^ 2 + t124 ^ 2) / 0.2e1 + m(2) * (t105 ^ 2 + t107 ^ 2 + t108 ^ 2) / 0.2e1 + m(3) * (t84 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + m(4) * (t73 ^ 2 + t74 ^ 2 + t75 ^ 2) / 0.2e1 + m(5) * (t70 ^ 2 + t71 ^ 2 + t72 ^ 2) / 0.2e1 + m(6) * (t67 ^ 2 + t68 ^ 2 + t69 ^ 2) / 0.2e1 + t111 * (t166 * t115 + t169 * t116) / 0.2e1 + t110 * (-t169 * t115 + t166 * t116) / 0.2e1 + ((-t163 * t81 - t164 * t79) * t111 + (-t163 * t80 - t164 * t78) * t110 + (-t164 * t130 - t163 * t135 + Icges(5,3)) * t154) * t154 / 0.2e1 + ((-t101 * t119 - t115 * t91 - t116 * t89 + t120 * t99 + t208 * t165 + t202 * t210 + Icges(1,4)) * V_base(5) + (t100 * t120 - t102 * t119 - t115 * t92 - t116 * t90 + t207 * t165 + t209 * t202 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t101 * t120 + t115 * t89 - t116 * t91 + t119 * t99 - t210 * t165 + t208 * t202 + Icges(1,2)) * V_base(5) + (t100 * t119 + t102 * t120 + t115 * t90 - t116 * t92 - t165 * t209 + t202 * t207 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t154 * (Icges(5,5) * t115 + Icges(5,6) * t116) + V_base(5) * t154 * (Icges(5,5) * t116 - Icges(5,6) * t115) + ((Icges(4,5) * t119 - Icges(4,6) * t120) * V_base(4) + (-Icges(4,5) * t120 - Icges(4,6) * t119) * V_base(5) + (Icges(2,3) / 0.2e1 + Icges(3,2) / 0.2e1 + Icges(4,3) / 0.2e1) * t156 + (-t213 * V_base(4) + t215 * V_base(5)) * t202 + (t213 * V_base(5) + t215 * V_base(4)) * t165) * t156 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
