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
% Datum: 2020-01-03 12:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:11:01
% EndTime: 2020-01-03 12:11:03
% DurationCPUTime: 2.26s
% Computational Cost: add. (1244->213), mult. (960->292), div. (0->0), fcn. (740->8), ass. (0->116)
t250 = Icges(5,4) + Icges(6,4);
t249 = Icges(5,1) + Icges(6,1);
t248 = Icges(5,2) + Icges(6,2);
t160 = qJ(3) + qJ(4);
t153 = cos(t160);
t247 = t250 * t153;
t151 = sin(t160);
t246 = t250 * t151;
t245 = Icges(5,5) + Icges(6,5);
t244 = -Icges(5,6) - Icges(6,6);
t243 = -t248 * t151 + t247;
t242 = t249 * t153 - t246;
t241 = rSges(6,1) + pkin(4);
t161 = qJ(1) + qJ(2);
t152 = sin(t161);
t154 = cos(t161);
t240 = t243 * t152 + t244 * t154;
t239 = t244 * t152 - t243 * t154;
t238 = t242 * t152 - t245 * t154;
t237 = -t245 * t152 - t242 * t154;
t236 = Icges(5,3) + Icges(6,3);
t235 = t248 * t153 + t246;
t234 = t249 * t151 + t247;
t233 = t244 * t151 + t245 * t153;
t232 = -rSges(6,3) - qJ(5);
t231 = -rSges(6,2) * t151 + t153 * t241;
t203 = -qJD(3) - qJD(4);
t105 = t152 * t203 + V_base(6);
t106 = t154 * t203 + V_base(5);
t150 = V_base(4) + qJD(1);
t148 = qJD(2) + t150;
t230 = (t239 * t151 - t237 * t153) * t105 + (t240 * t151 - t238 * t153) * t106 + (t235 * t151 - t234 * t153) * t148;
t229 = (-t245 * t151 + t244 * t153) * t148 + (-t233 * t152 + t236 * t154) * t106 + (t236 * t152 + t233 * t154) * t105;
t164 = cos(qJ(3));
t162 = sin(qJ(3));
t211 = Icges(4,4) * t162;
t187 = Icges(4,1) * t164 - t211;
t100 = -Icges(4,5) * t154 + t152 * t187;
t101 = -Icges(4,5) * t152 - t154 * t187;
t131 = -qJD(3) * t152 + V_base(6);
t132 = -qJD(3) * t154 + V_base(5);
t136 = Icges(4,2) * t164 + t211;
t210 = Icges(4,4) * t164;
t139 = Icges(4,1) * t162 + t210;
t184 = -Icges(4,2) * t162 + t210;
t98 = -Icges(4,6) * t154 + t152 * t184;
t99 = -Icges(4,6) * t152 - t154 * t184;
t225 = (t101 * t164 - t162 * t99) * t131 + (t100 * t164 - t162 * t98) * t132 - (t136 * t162 - t139 * t164) * t148;
t224 = -pkin(5) - pkin(6);
t220 = pkin(1) * t150;
t219 = pkin(3) * t162;
t218 = t164 * pkin(3);
t216 = t231 * t152 + t232 * t154;
t215 = t232 * t152 - t231 * t154;
t126 = -t154 * pkin(2) - t152 * pkin(7);
t76 = -pkin(8) * t152 - t154 * t218;
t214 = -t126 - t76;
t165 = cos(qJ(1));
t213 = Icges(2,4) * t165;
t212 = Icges(3,4) * t154;
t163 = sin(qJ(1));
t202 = t163 * t220 + V_base(3);
t201 = V_base(6) * pkin(5) + V_base(2);
t198 = rSges(6,2) * t153 + t151 * t241;
t197 = V_base(6) * pkin(6) + t165 * t220 + t201;
t196 = rSges(4,1) * t164 - rSges(4,2) * t162;
t195 = rSges(5,1) * t153 - rSges(5,2) * t151;
t181 = Icges(4,5) * t164 - Icges(4,6) * t162;
t175 = t131 * t219 + t197;
t172 = -(-Icges(4,3) * t152 - t154 * t181) * t131 - (-Icges(4,3) * t154 + t152 * t181) * t132 - (Icges(4,5) * t162 + Icges(4,6) * t164) * t148;
t125 = t152 * pkin(2) - t154 * pkin(7);
t171 = t148 * t125 + t224 * V_base(5) + t202;
t170 = V_base(1) + (-V_base(6) * t163 - t165 * V_base(5)) * pkin(1);
t75 = -pkin(8) * t154 + t152 * t218;
t169 = -t132 * t219 + t148 * t75 + t171;
t168 = -t125 * V_base(6) + V_base(5) * t126 + t170;
t167 = -t131 * t75 + t132 * t76 + t168;
t155 = Icges(2,4) * t163;
t147 = Icges(3,4) * t152;
t144 = -rSges(2,1) * t165 + rSges(2,2) * t163;
t143 = rSges(2,1) * t163 + rSges(2,2) * t165;
t142 = rSges(4,1) * t162 + rSges(4,2) * t164;
t141 = -Icges(2,1) * t165 + t155;
t140 = Icges(2,1) * t163 + t213;
t138 = Icges(2,2) * t163 - t213;
t137 = Icges(2,2) * t165 + t155;
t130 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t129 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t128 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t124 = -rSges(3,1) * t154 + rSges(3,2) * t152;
t123 = rSges(3,1) * t152 + rSges(3,2) * t154;
t122 = rSges(5,1) * t151 + rSges(5,2) * t153;
t120 = -Icges(3,1) * t154 + t147;
t119 = Icges(3,1) * t152 + t212;
t116 = Icges(3,2) * t152 - t212;
t115 = Icges(3,2) * t154 + t147;
t103 = -rSges(4,3) * t152 - t154 * t196;
t102 = -rSges(4,3) * t154 + t152 * t196;
t95 = V_base(6) * rSges(2,3) - t144 * t150 + t201;
t94 = t143 * t150 + V_base(3) + (-rSges(2,3) - pkin(5)) * V_base(5);
t93 = -t143 * V_base(6) + t144 * V_base(5) + V_base(1);
t92 = -rSges(5,3) * t152 - t154 * t195;
t90 = -rSges(5,3) * t154 + t152 * t195;
t72 = V_base(6) * rSges(3,3) - t124 * t148 + t197;
t71 = t123 * t148 + (-rSges(3,3) + t224) * V_base(5) + t202;
t68 = -t123 * V_base(6) + t124 * V_base(5) + t170;
t67 = t131 * t142 + (-t103 - t126) * t148 + t197;
t66 = t102 * t148 - t132 * t142 + t171;
t65 = -t102 * t131 + t103 * t132 + t168;
t64 = t105 * t122 + (-t92 + t214) * t148 + t175;
t63 = -t106 * t122 + t148 * t90 + t169;
t62 = -t105 * t90 + t106 * t92 + t167;
t61 = -qJD(5) * t154 + t198 * t105 + (t214 - t215) * t148 + t175;
t60 = -qJD(5) * t152 - t106 * t198 + t148 * t216 + t169;
t59 = -t105 * t216 + t106 * t215 + t167;
t1 = m(1) * (t128 ^ 2 + t129 ^ 2 + t130 ^ 2) / 0.2e1 + m(2) * (t93 ^ 2 + t94 ^ 2 + t95 ^ 2) / 0.2e1 + m(3) * (t68 ^ 2 + t71 ^ 2 + t72 ^ 2) / 0.2e1 + m(4) * (t65 ^ 2 + t66 ^ 2 + t67 ^ 2) / 0.2e1 + t132 * (t225 * t152 + t172 * t154) / 0.2e1 + t131 * (t172 * t152 - t225 * t154) / 0.2e1 + m(5) * (t62 ^ 2 + t63 ^ 2 + t64 ^ 2) / 0.2e1 + m(6) * (t59 ^ 2 + t60 ^ 2 + t61 ^ 2) / 0.2e1 + (t229 * t152 + t230 * t154) * t105 / 0.2e1 + (-t230 * t152 + t229 * t154) * t106 / 0.2e1 + ((t116 * t154 + t120 * t152 + t138 * t165 + t141 * t163 + Icges(1,6)) * V_base(6) + (t154 * t115 + t152 * t119 + t165 * t137 + t163 * t140 + Icges(1,2)) * V_base(5)) * V_base(5) / 0.2e1 + ((t152 * t116 - t154 * t120 + t163 * t138 - t165 * t141 + Icges(1,3)) * V_base(6) + (t115 * t152 - t119 * t154 + t137 * t163 - t140 * t165 + Icges(1,6)) * V_base(5)) * V_base(6) / 0.2e1 + ((t100 * t162 + t164 * t98) * t132 + (t101 * t162 + t164 * t99) * t131 + (t238 * t151 + t240 * t153) * t106 + (t237 * t151 + t239 * t153) * t105 + (t164 * t136 + t162 * t139 + t234 * t151 + t235 * t153 + Icges(3,3)) * t148) * t148 / 0.2e1 + t148 * V_base(6) * (-Icges(3,5) * t154 + Icges(3,6) * t152) + t148 * V_base(5) * (Icges(3,5) * t152 + Icges(3,6) * t154) + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + Icges(1,1) * V_base(4) / 0.2e1) * V_base(4) + ((Icges(2,5) * t163 + Icges(2,6) * t165) * V_base(5) + (-Icges(2,5) * t165 + Icges(2,6) * t163) * V_base(6) + Icges(2,3) * t150 / 0.2e1) * t150;
T = t1;
