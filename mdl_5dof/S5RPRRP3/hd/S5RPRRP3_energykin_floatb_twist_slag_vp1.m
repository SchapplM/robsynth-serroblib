% Calculate kinetic energy for
% S5RPRRP3
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2020-01-03 11:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRRP3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP3_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP3_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP3_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:46:50
% EndTime: 2020-01-03 11:46:52
% DurationCPUTime: 2.58s
% Computational Cost: add. (1212->213), mult. (960->287), div. (0->0), fcn. (740->8), ass. (0->115)
t251 = Icges(5,4) + Icges(6,4);
t250 = Icges(5,1) + Icges(6,1);
t249 = Icges(5,2) + Icges(6,2);
t160 = qJ(3) + qJ(4);
t153 = cos(t160);
t248 = t251 * t153;
t152 = sin(t160);
t247 = t251 * t152;
t246 = Icges(5,5) + Icges(6,5);
t245 = -Icges(5,6) - Icges(6,6);
t244 = -t249 * t152 + t248;
t243 = t250 * t153 - t247;
t242 = rSges(6,1) + pkin(4);
t159 = qJ(1) + pkin(8);
t149 = sin(t159);
t150 = cos(t159);
t241 = t244 * t149 + t245 * t150;
t240 = t245 * t149 - t244 * t150;
t239 = t243 * t149 - t246 * t150;
t238 = -t246 * t149 - t243 * t150;
t237 = Icges(5,3) + Icges(6,3);
t236 = t249 * t153 + t247;
t235 = t250 * t152 + t248;
t234 = t245 * t152 + t246 * t153;
t233 = -rSges(6,3) - qJ(5);
t232 = -rSges(6,2) * t152 + t242 * t153;
t202 = -qJD(3) - qJD(4);
t105 = t149 * t202 + V_base(6);
t106 = t150 * t202 + V_base(5);
t151 = V_base(4) + qJD(1);
t231 = t105 * (t240 * t152 - t238 * t153) + t106 * (t241 * t152 - t239 * t153) + t151 * (t236 * t152 - t235 * t153);
t228 = (-t246 * t152 + t245 * t153) * t151 + (-t234 * t149 + t237 * t150) * t106 + (t237 * t149 + t234 * t150) * t105;
t131 = -qJD(3) * t149 + V_base(6);
t132 = -qJD(3) * t150 + V_base(5);
t163 = cos(qJ(3));
t161 = sin(qJ(3));
t210 = Icges(4,4) * t161;
t136 = Icges(4,2) * t163 + t210;
t209 = Icges(4,4) * t163;
t139 = Icges(4,1) * t161 + t209;
t183 = -Icges(4,2) * t161 + t209;
t96 = -Icges(4,6) * t150 + t149 * t183;
t97 = -Icges(4,6) * t149 - t150 * t183;
t186 = Icges(4,1) * t163 - t210;
t98 = -Icges(4,5) * t150 + t149 * t186;
t99 = -Icges(4,5) * t149 - t150 * t186;
t224 = (t136 * t161 - t139 * t163) * t151 + (t161 * t96 - t163 * t98) * t132 + (t161 * t97 - t163 * t99) * t131;
t220 = pkin(1) * t151;
t219 = pkin(3) * t161;
t218 = t163 * pkin(3);
t217 = -pkin(5) - qJ(2);
t215 = t232 * t149 + t150 * t233;
t214 = t149 * t233 - t232 * t150;
t118 = -t150 * pkin(2) - t149 * pkin(6);
t76 = -pkin(7) * t149 - t150 * t218;
t213 = -t118 - t76;
t164 = cos(qJ(1));
t212 = Icges(2,4) * t164;
t211 = Icges(3,4) * t150;
t162 = sin(qJ(1));
t201 = t162 * t220 + V_base(3);
t200 = V_base(6) * pkin(5) + V_base(2);
t197 = rSges(6,2) * t153 + t242 * t152;
t196 = V_base(6) * qJ(2) + t164 * t220 + t200;
t195 = rSges(4,1) * t163 - rSges(4,2) * t161;
t194 = rSges(5,1) * t153 - rSges(5,2) * t152;
t180 = Icges(4,5) * t163 - Icges(4,6) * t161;
t174 = t131 * t219 + t196;
t171 = -(-Icges(4,3) * t149 - t150 * t180) * t131 - (-Icges(4,3) * t150 + t149 * t180) * t132 - (Icges(4,5) * t161 + Icges(4,6) * t163) * t151;
t117 = t149 * pkin(2) - t150 * pkin(6);
t170 = t151 * t117 + t217 * V_base(5) + t201;
t169 = qJD(2) + V_base(1) + (-V_base(6) * t162 - t164 * V_base(5)) * pkin(1);
t75 = -pkin(7) * t150 + t149 * t218;
t168 = -t132 * t219 + t151 * t75 + t170;
t167 = -t117 * V_base(6) + V_base(5) * t118 + t169;
t166 = -t131 * t75 + t132 * t76 + t167;
t155 = Icges(2,4) * t162;
t147 = Icges(3,4) * t149;
t144 = -rSges(2,1) * t164 + rSges(2,2) * t162;
t143 = rSges(2,1) * t162 + rSges(2,2) * t164;
t142 = rSges(4,1) * t161 + rSges(4,2) * t163;
t141 = -Icges(2,1) * t164 + t155;
t140 = Icges(2,1) * t162 + t212;
t138 = Icges(2,2) * t162 - t212;
t137 = Icges(2,2) * t164 + t155;
t130 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t129 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t128 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t126 = rSges(5,1) * t152 + rSges(5,2) * t153;
t116 = -rSges(3,1) * t150 + rSges(3,2) * t149;
t115 = rSges(3,1) * t149 + rSges(3,2) * t150;
t114 = -Icges(3,1) * t150 + t147;
t113 = Icges(3,1) * t149 + t211;
t112 = Icges(3,2) * t149 - t211;
t111 = Icges(3,2) * t150 + t147;
t103 = -rSges(4,3) * t149 - t150 * t195;
t102 = -rSges(4,3) * t150 + t149 * t195;
t101 = V_base(6) * rSges(2,3) - t144 * t151 + t200;
t100 = t143 * t151 + V_base(3) + (-rSges(2,3) - pkin(5)) * V_base(5);
t93 = -t143 * V_base(6) + t144 * V_base(5) + V_base(1);
t92 = -rSges(5,3) * t149 - t150 * t194;
t90 = -rSges(5,3) * t150 + t149 * t194;
t72 = V_base(6) * rSges(3,3) - t116 * t151 + t196;
t71 = t115 * t151 + (-rSges(3,3) + t217) * V_base(5) + t201;
t68 = -t115 * V_base(6) + t116 * V_base(5) + t169;
t67 = t131 * t142 + (-t103 - t118) * t151 + t196;
t66 = t102 * t151 - t132 * t142 + t170;
t65 = -t102 * t131 + t103 * t132 + t167;
t64 = t105 * t126 + (-t92 + t213) * t151 + t174;
t63 = -t106 * t126 + t151 * t90 + t168;
t62 = -t105 * t90 + t106 * t92 + t166;
t61 = -qJD(5) * t150 + t197 * t105 + (t213 - t214) * t151 + t174;
t60 = -qJD(5) * t149 - t106 * t197 + t151 * t215 + t168;
t59 = -t105 * t215 + t106 * t214 + t166;
t1 = m(1) * (t128 ^ 2 + t129 ^ 2 + t130 ^ 2) / 0.2e1 + m(2) * (t100 ^ 2 + t101 ^ 2 + t93 ^ 2) / 0.2e1 + m(3) * (t68 ^ 2 + t71 ^ 2 + t72 ^ 2) / 0.2e1 + m(4) * (t65 ^ 2 + t66 ^ 2 + t67 ^ 2) / 0.2e1 + t132 * (-t224 * t149 + t171 * t150) / 0.2e1 + t131 * (t171 * t149 + t224 * t150) / 0.2e1 + m(5) * (t62 ^ 2 + t63 ^ 2 + t64 ^ 2) / 0.2e1 + m(6) * (t59 ^ 2 + t60 ^ 2 + t61 ^ 2) / 0.2e1 + (t228 * t149 + t231 * t150) * t105 / 0.2e1 + (-t231 * t149 + t228 * t150) * t106 / 0.2e1 + ((t112 * t150 + t114 * t149 + t138 * t164 + t141 * t162 + Icges(1,6)) * V_base(6) + (t150 * t111 + t149 * t113 + t164 * t137 + t162 * t140 + Icges(1,2)) * V_base(5)) * V_base(5) / 0.2e1 + ((t149 * t112 - t150 * t114 + t162 * t138 - t164 * t141 + Icges(1,3)) * V_base(6) + (t111 * t149 - t113 * t150 + t137 * t162 - t140 * t164 + Icges(1,6)) * V_base(5)) * V_base(6) / 0.2e1 + ((t161 * t98 + t163 * t96) * t132 + (t161 * t99 + t163 * t97) * t131 + (t239 * t152 + t241 * t153) * t106 + (t238 * t152 + t240 * t153) * t105 + (t163 * t136 + t161 * t139 + t235 * t152 + t236 * t153 + Icges(2,3) + Icges(3,3)) * t151) * t151 / 0.2e1 + t151 * V_base(6) * (-Icges(2,5) * t164 - Icges(3,5) * t150 + Icges(2,6) * t162 + Icges(3,6) * t149) + t151 * V_base(5) * (Icges(2,5) * t162 + Icges(3,5) * t149 + Icges(2,6) * t164 + Icges(3,6) * t150) + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + Icges(1,1) * V_base(4) / 0.2e1) * V_base(4);
T = t1;
