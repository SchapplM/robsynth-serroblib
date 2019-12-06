% Calculate kinetic energy for
% S5RRRRP2
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
% Datum: 2019-12-05 18:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRP2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRRP2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP2_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP2_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP2_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:47:29
% EndTime: 2019-12-05 18:47:31
% DurationCPUTime: 2.28s
% Computational Cost: add. (1244->212), mult. (960->292), div. (0->0), fcn. (740->8), ass. (0->116)
t252 = Icges(5,4) + Icges(6,4);
t251 = Icges(5,1) + Icges(6,1);
t250 = Icges(5,2) + Icges(6,2);
t160 = qJ(3) + qJ(4);
t154 = cos(t160);
t249 = t252 * t154;
t152 = sin(t160);
t248 = t252 * t152;
t247 = Icges(5,5) + Icges(6,5);
t246 = Icges(5,6) + Icges(6,6);
t245 = -t250 * t152 + t249;
t244 = t251 * t154 - t248;
t243 = rSges(6,1) + pkin(4);
t161 = qJ(1) + qJ(2);
t153 = sin(t161);
t155 = cos(t161);
t242 = -t245 * t153 + t246 * t155;
t241 = t246 * t153 + t245 * t155;
t240 = -t244 * t153 + t247 * t155;
t239 = t247 * t153 + t244 * t155;
t238 = Icges(5,3) + Icges(6,3);
t237 = t250 * t154 + t248;
t236 = t251 * t152 + t249;
t235 = -t246 * t152 + t247 * t154;
t234 = rSges(6,3) + qJ(5);
t233 = -rSges(6,2) * t152 + t243 * t154;
t131 = qJD(3) * t153 + V_base(6);
t105 = qJD(4) * t153 + t131;
t132 = qJD(3) * t155 + V_base(5);
t106 = qJD(4) * t155 + t132;
t151 = V_base(4) + qJD(1);
t149 = qJD(2) + t151;
t232 = t105 * (t241 * t152 - t239 * t154) + t106 * (t242 * t152 - t240 * t154) + t149 * (t237 * t152 - t236 * t154);
t231 = (t247 * t152 + t246 * t154) * t149 + (-t235 * t153 + t238 * t155) * t106 + (t238 * t153 + t235 * t155) * t105;
t164 = cos(qJ(3));
t162 = sin(qJ(3));
t210 = Icges(4,4) * t162;
t187 = Icges(4,1) * t164 - t210;
t100 = Icges(4,5) * t155 - t153 * t187;
t101 = Icges(4,5) * t153 + t155 * t187;
t136 = Icges(4,2) * t164 + t210;
t209 = Icges(4,4) * t164;
t139 = Icges(4,1) * t162 + t209;
t184 = -Icges(4,2) * t162 + t209;
t98 = Icges(4,6) * t155 - t153 * t184;
t99 = Icges(4,6) * t153 + t155 * t184;
t227 = (t101 * t164 - t162 * t99) * t131 + (t100 * t164 - t162 * t98) * t132 - (t136 * t162 - t139 * t164) * t149;
t226 = -pkin(5) - pkin(6);
t163 = sin(qJ(1));
t222 = pkin(1) * t163;
t165 = cos(qJ(1));
t221 = pkin(1) * t165;
t220 = pkin(3) * t162;
t219 = t164 * pkin(3);
t217 = -t233 * t153 + t155 * t234;
t216 = t153 * t234 + t233 * t155;
t126 = t155 * pkin(2) + t153 * pkin(7);
t76 = pkin(8) * t153 + t155 * t219;
t215 = -t126 - t76;
t214 = Icges(2,4) * t163;
t213 = Icges(2,4) * t165;
t212 = Icges(3,4) * t153;
t211 = Icges(3,4) * t155;
t202 = V_base(6) * pkin(5) + V_base(2);
t199 = rSges(6,2) * t154 + t243 * t152;
t198 = V_base(5) * t221 + V_base(6) * t222 + V_base(1);
t197 = rSges(4,1) * t164 - rSges(4,2) * t162;
t196 = rSges(5,1) * t154 - rSges(5,2) * t152;
t194 = -t151 * t222 + V_base(3);
t181 = Icges(4,5) * t164 - Icges(4,6) * t162;
t175 = V_base(6) * pkin(6) - t151 * t221 + t202;
t172 = (Icges(4,3) * t153 + t155 * t181) * t131 + (Icges(4,3) * t155 - t153 * t181) * t132 + (Icges(4,5) * t162 + Icges(4,6) * t164) * t149;
t171 = t131 * t220 + t175;
t125 = -t153 * pkin(2) + t155 * pkin(7);
t170 = -t125 * V_base(6) + V_base(5) * t126 + t198;
t169 = t149 * t125 + t226 * V_base(5) + t194;
t75 = pkin(8) * t155 - t153 * t219;
t168 = -t131 * t75 + t132 * t76 + t170;
t167 = -t132 * t220 + t149 * t75 + t169;
t144 = rSges(2,1) * t165 - rSges(2,2) * t163;
t143 = -rSges(2,1) * t163 - rSges(2,2) * t165;
t142 = rSges(4,1) * t162 + rSges(4,2) * t164;
t141 = Icges(2,1) * t165 - t214;
t140 = -Icges(2,1) * t163 - t213;
t138 = -Icges(2,2) * t163 + t213;
t137 = -Icges(2,2) * t165 - t214;
t130 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t129 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t128 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t124 = rSges(3,1) * t155 - rSges(3,2) * t153;
t123 = -rSges(3,1) * t153 - rSges(3,2) * t155;
t122 = rSges(5,1) * t152 + rSges(5,2) * t154;
t120 = Icges(3,1) * t155 - t212;
t119 = -Icges(3,1) * t153 - t211;
t116 = -Icges(3,2) * t153 + t211;
t115 = -Icges(3,2) * t155 - t212;
t103 = rSges(4,3) * t153 + t155 * t197;
t102 = rSges(4,3) * t155 - t153 * t197;
t95 = V_base(6) * rSges(2,3) - t144 * t151 + t202;
t94 = t143 * t151 + V_base(3) + (-rSges(2,3) - pkin(5)) * V_base(5);
t93 = -t143 * V_base(6) + t144 * V_base(5) + V_base(1);
t92 = rSges(5,3) * t153 + t155 * t196;
t90 = rSges(5,3) * t155 - t153 * t196;
t72 = V_base(6) * rSges(3,3) - t124 * t149 + t175;
t71 = t123 * t149 + (-rSges(3,3) + t226) * V_base(5) + t194;
t68 = -t123 * V_base(6) + t124 * V_base(5) + t198;
t67 = t131 * t142 + (-t103 - t126) * t149 + t175;
t66 = t102 * t149 - t132 * t142 + t169;
t65 = -t102 * t131 + t103 * t132 + t170;
t64 = t105 * t122 + (-t92 + t215) * t149 + t171;
t63 = -t106 * t122 + t149 * t90 + t167;
t62 = -t105 * t90 + t106 * t92 + t168;
t61 = qJD(5) * t155 + t199 * t105 + (t215 - t216) * t149 + t171;
t60 = qJD(5) * t153 - t106 * t199 + t149 * t217 + t167;
t59 = -t105 * t217 + t106 * t216 + t168;
t1 = m(1) * (t128 ^ 2 + t129 ^ 2 + t130 ^ 2) / 0.2e1 + m(2) * (t93 ^ 2 + t94 ^ 2 + t95 ^ 2) / 0.2e1 + m(3) * (t68 ^ 2 + t71 ^ 2 + t72 ^ 2) / 0.2e1 + m(4) * (t65 ^ 2 + t66 ^ 2 + t67 ^ 2) / 0.2e1 + t132 * (-t153 * t227 + t172 * t155) / 0.2e1 + t131 * (t172 * t153 + t155 * t227) / 0.2e1 + m(5) * (t62 ^ 2 + t63 ^ 2 + t64 ^ 2) / 0.2e1 + m(6) * (t59 ^ 2 + t60 ^ 2 + t61 ^ 2) / 0.2e1 + (t231 * t153 - t232 * t155) * t105 / 0.2e1 + (t232 * t153 + t231 * t155) * t106 / 0.2e1 + ((-t116 * t155 - t120 * t153 - t138 * t165 - t141 * t163 + Icges(1,6)) * V_base(6) + (-t155 * t115 - t153 * t119 - t165 * t137 - t163 * t140 + Icges(1,2)) * V_base(5)) * V_base(5) / 0.2e1 + ((-t153 * t116 + t155 * t120 - t163 * t138 + t165 * t141 + Icges(1,3)) * V_base(6) + (-t115 * t153 + t119 * t155 - t137 * t163 + t140 * t165 + Icges(1,6)) * V_base(5)) * V_base(6) / 0.2e1 + ((t100 * t162 + t164 * t98) * t132 + (t101 * t162 + t164 * t99) * t131 + (t240 * t152 + t242 * t154) * t106 + (t239 * t152 + t241 * t154) * t105 + (t164 * t136 + t162 * t139 + t236 * t152 + t237 * t154 + Icges(3,3)) * t149) * t149 / 0.2e1 + t149 * V_base(6) * (Icges(3,5) * t155 - Icges(3,6) * t153) + t149 * V_base(5) * (-Icges(3,5) * t153 - Icges(3,6) * t155) + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + Icges(1,1) * V_base(4) / 0.2e1) * V_base(4) + ((-Icges(2,5) * t163 - Icges(2,6) * t165) * V_base(5) + (Icges(2,5) * t165 - Icges(2,6) * t163) * V_base(6) + Icges(2,3) * t151 / 0.2e1) * t151;
T = t1;
