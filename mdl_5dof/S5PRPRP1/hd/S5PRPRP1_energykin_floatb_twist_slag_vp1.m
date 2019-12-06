% Calculate kinetic energy for
% S5PRPRP1
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP1_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP1_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRPRP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP1_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP1_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP1_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP1_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:28:15
% EndTime: 2019-12-05 15:28:17
% DurationCPUTime: 2.15s
% Computational Cost: add. (1087->218), mult. (911->282), div. (0->0), fcn. (691->8), ass. (0->119)
t247 = Icges(5,4) - Icges(6,5);
t246 = Icges(5,1) + Icges(6,1);
t245 = Icges(5,2) + Icges(6,3);
t161 = pkin(8) + qJ(4);
t155 = cos(t161);
t244 = t247 * t155;
t153 = sin(t161);
t243 = t247 * t153;
t242 = Icges(6,4) + Icges(5,5);
t241 = Icges(5,6) - Icges(6,6);
t240 = t245 * t153 - t244;
t239 = t246 * t155 - t243;
t238 = rSges(6,1) + pkin(4);
t237 = rSges(6,3) + qJ(5);
t162 = pkin(7) + qJ(2);
t154 = sin(t162);
t156 = cos(t162);
t236 = t240 * t154 + t241 * t156;
t235 = -t241 * t154 + t240 * t156;
t234 = t239 * t154 - t242 * t156;
t233 = t242 * t154 + t239 * t156;
t232 = Icges(6,2) + Icges(5,3);
t231 = -t245 * t155 - t243;
t230 = t246 * t153 + t244;
t229 = -t241 * t153 + t242 * t155;
t228 = t237 * t153 + t238 * t155;
t142 = -qJD(4) * t156 + V_base(5);
t143 = qJD(4) * t154 + V_base(4);
t158 = V_base(6) + qJD(2);
t225 = (t231 * t153 + t230 * t155) * t158 + (t235 * t153 + t233 * t155) * t143 + (t236 * t153 + t234 * t155) * t142;
t224 = (t242 * t153 + t241 * t155) * t158 + (t232 * t154 + t229 * t156) * t143 + (t229 * t154 - t232 * t156) * t142;
t166 = cos(pkin(7));
t220 = pkin(1) * t166;
t163 = sin(pkin(8));
t219 = pkin(3) * t163;
t165 = cos(pkin(8));
t218 = pkin(3) * t165;
t217 = -pkin(5) - qJ(1);
t216 = -rSges(6,2) * t156 + t228 * t154;
t215 = rSges(6,2) * t154 + t228 * t156;
t126 = t154 * pkin(2) - t156 * qJ(3);
t78 = -pkin(6) * t156 + t154 * t218;
t214 = -t126 - t78;
t164 = sin(pkin(7));
t213 = Icges(2,4) * t164;
t212 = Icges(3,4) * t154;
t211 = Icges(4,4) * t163;
t210 = Icges(4,4) * t165;
t204 = t238 * t153 - t237 * t155;
t203 = qJD(5) * t153;
t196 = pkin(1) * V_base(6);
t202 = t166 * t196 + V_base(2);
t201 = V_base(5) * qJ(1) + V_base(1);
t197 = qJD(1) + V_base(3);
t128 = t156 * pkin(2) + t154 * qJ(3);
t195 = -t128 - t220;
t194 = V_base(4) * t164 * pkin(1) + t197;
t193 = V_base(4) * t126 + t194;
t192 = rSges(4,1) * t165 - rSges(4,2) * t163;
t191 = rSges(5,1) * t155 - rSges(5,2) * t153;
t188 = Icges(4,1) * t165 - t211;
t185 = -Icges(4,2) * t163 + t210;
t182 = Icges(4,5) * t165 - Icges(4,6) * t163;
t179 = -qJD(3) * t156 + t158 * t128 + t202;
t176 = V_base(5) * pkin(5) - t164 * t196 + t201;
t175 = qJD(3) * t154 + t176;
t174 = (Icges(4,5) * t163 + Icges(4,6) * t165) * t158 + (-Icges(4,3) * t156 + t154 * t182) * V_base(5) + (Icges(4,3) * t154 + t156 * t182) * V_base(4);
t173 = V_base(5) * t219 + t175;
t79 = pkin(6) * t154 + t156 * t218;
t172 = V_base(4) * t78 + (t195 - t79) * V_base(5) + t193;
t171 = t158 * t79 + (t217 - t219) * V_base(4) + t179;
t100 = Icges(4,6) * t154 + t156 * t185;
t101 = -Icges(4,5) * t156 + t154 * t188;
t102 = Icges(4,5) * t154 + t156 * t188;
t136 = Icges(4,2) * t165 + t211;
t139 = Icges(4,1) * t163 + t210;
t99 = -Icges(4,6) * t156 + t154 * t185;
t168 = (-t100 * t163 + t102 * t165) * V_base(4) + (t101 * t165 - t163 * t99) * V_base(5) + (-t136 * t163 + t139 * t165) * t158;
t157 = Icges(2,4) * t166;
t151 = Icges(3,4) * t156;
t146 = rSges(2,1) * t166 - rSges(2,2) * t164;
t145 = rSges(2,1) * t164 + rSges(2,2) * t166;
t144 = rSges(4,1) * t163 + rSges(4,2) * t165;
t141 = Icges(2,1) * t166 - t213;
t140 = Icges(2,1) * t164 + t157;
t138 = -Icges(2,2) * t164 + t157;
t137 = Icges(2,2) * t166 + t213;
t132 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t131 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t130 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t129 = rSges(3,1) * t156 - rSges(3,2) * t154;
t127 = rSges(3,1) * t154 + rSges(3,2) * t156;
t125 = rSges(5,1) * t153 + rSges(5,2) * t155;
t122 = Icges(3,1) * t156 - t212;
t121 = Icges(3,1) * t154 + t151;
t118 = -Icges(3,2) * t154 + t151;
t117 = Icges(3,2) * t156 + t212;
t114 = Icges(3,5) * t156 - Icges(3,6) * t154;
t113 = Icges(3,5) * t154 + Icges(3,6) * t156;
t106 = V_base(5) * rSges(2,3) - t145 * V_base(6) + t201;
t105 = t146 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t104 = rSges(4,3) * t154 + t156 * t192;
t103 = -rSges(4,3) * t156 + t154 * t192;
t96 = rSges(5,3) * t154 + t156 * t191;
t94 = -rSges(5,3) * t156 + t154 * t191;
t80 = t145 * V_base(4) - t146 * V_base(5) + t197;
t75 = V_base(5) * rSges(3,3) - t127 * t158 + t176;
t74 = t129 * t158 + (-rSges(3,3) + t217) * V_base(4) + t202;
t73 = t127 * V_base(4) + (-t129 - t220) * V_base(5) + t194;
t72 = t144 * V_base(5) + (-t103 - t126) * t158 + t175;
t71 = t104 * t158 + (-t144 + t217) * V_base(4) + t179;
t70 = t103 * V_base(4) + (-t104 + t195) * V_base(5) + t193;
t69 = t125 * t142 + (-t94 + t214) * t158 + t173;
t68 = -t125 * t143 + t158 * t96 + t171;
t67 = -t142 * t96 + t143 * t94 + t172;
t66 = t156 * t203 + t204 * t142 + (t214 - t216) * t158 + t173;
t65 = -t143 * t204 + t154 * t203 + t158 * t215 + t171;
t64 = -qJD(5) * t155 - t142 * t215 + t143 * t216 + t172;
t1 = m(1) * (t130 ^ 2 + t131 ^ 2 + t132 ^ 2) / 0.2e1 + m(2) * (t105 ^ 2 + t106 ^ 2 + t80 ^ 2) / 0.2e1 + m(3) * (t73 ^ 2 + t74 ^ 2 + t75 ^ 2) / 0.2e1 + m(4) * (t70 ^ 2 + t71 ^ 2 + t72 ^ 2) / 0.2e1 + m(5) * (t67 ^ 2 + t68 ^ 2 + t69 ^ 2) / 0.2e1 + m(6) * (t64 ^ 2 + t65 ^ 2 + t66 ^ 2) / 0.2e1 + (t225 * t154 - t224 * t156) * t142 / 0.2e1 + (t224 * t154 + t225 * t156) * t143 / 0.2e1 + ((t101 * t163 + t165 * t99 + t113) * V_base(5) + (t100 * t165 + t102 * t163 + t114) * V_base(4) + (t233 * t153 - t235 * t155) * t143 + (t234 * t153 - t236 * t155) * t142 + (t165 * t136 + t163 * t139 + t230 * t153 - t231 * t155 + Icges(3,3)) * t158) * t158 / 0.2e1 + (t114 * t158 + t154 * t174 + t168 * t156 + (-t117 * t154 + t121 * t156 - t137 * t164 + t140 * t166 + Icges(1,4)) * V_base(5) + (-t154 * t118 + t156 * t122 - t164 * t138 + t166 * t141 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t113 * t158 + t154 * t168 - t174 * t156 + (t156 * t117 + t154 * t121 + t166 * t137 + t164 * t140 + Icges(1,2)) * V_base(5) + (t118 * t156 + t122 * t154 + t138 * t166 + t141 * t164 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t164 + Icges(2,6) * t166 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t166 - Icges(2,6) * t164 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T = t1;
