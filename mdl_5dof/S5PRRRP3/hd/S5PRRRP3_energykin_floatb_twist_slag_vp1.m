% Calculate kinetic energy for
% S5PRRRP3
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRP3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP3_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRRRP3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP3_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP3_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRP3_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:43:27
% EndTime: 2019-12-05 16:43:29
% DurationCPUTime: 2.27s
% Computational Cost: add. (1199->214), mult. (960->288), div. (0->0), fcn. (740->8), ass. (0->116)
t245 = Icges(5,4) + Icges(6,4);
t244 = Icges(5,1) + Icges(6,1);
t243 = Icges(5,2) + Icges(6,2);
t162 = qJ(3) + qJ(4);
t156 = cos(t162);
t242 = t245 * t156;
t155 = sin(t162);
t241 = t245 * t155;
t240 = Icges(5,5) + Icges(6,5);
t239 = Icges(5,6) + Icges(6,6);
t238 = -t243 * t155 + t242;
t237 = t244 * t156 - t241;
t236 = rSges(6,1) + pkin(4);
t161 = pkin(8) + qJ(2);
t151 = sin(t161);
t152 = cos(t161);
t235 = t238 * t151 - t239 * t152;
t234 = t239 * t151 + t238 * t152;
t233 = t237 * t151 - t240 * t152;
t232 = t240 * t151 + t237 * t152;
t231 = Icges(5,3) + Icges(6,3);
t230 = t243 * t156 + t241;
t229 = t244 * t155 + t242;
t228 = -t239 * t155 + t240 * t156;
t227 = rSges(6,3) + qJ(5);
t226 = -rSges(6,2) * t155 + t236 * t156;
t106 = V_base(5) + (-qJD(3) - qJD(4)) * t152;
t139 = qJD(3) * t151 + V_base(4);
t107 = qJD(4) * t151 + t139;
t154 = V_base(6) + qJD(2);
t223 = (-t230 * t155 + t229 * t156) * t154 + (-t234 * t155 + t232 * t156) * t107 + (-t235 * t155 + t233 * t156) * t106;
t222 = (t240 * t155 + t239 * t156) * t154 + (t231 * t151 + t228 * t152) * t107 + (t228 * t151 - t231 * t152) * t106;
t164 = cos(pkin(8));
t218 = pkin(1) * t164;
t165 = sin(qJ(3));
t217 = pkin(3) * t165;
t166 = cos(qJ(3));
t216 = t166 * pkin(3);
t215 = -pkin(5) - qJ(1);
t213 = t226 * t151 - t227 * t152;
t212 = t227 * t151 + t226 * t152;
t118 = t151 * pkin(2) - t152 * pkin(6);
t76 = -pkin(7) * t152 + t151 * t216;
t211 = -t118 - t76;
t163 = sin(pkin(8));
t210 = Icges(2,4) * t163;
t209 = Icges(3,4) * t151;
t208 = Icges(4,4) * t165;
t207 = Icges(4,4) * t166;
t194 = pkin(1) * V_base(6);
t200 = t164 * t194 + V_base(2);
t199 = V_base(5) * qJ(1) + V_base(1);
t195 = qJD(1) + V_base(3);
t193 = rSges(6,2) * t156 + t236 * t155;
t192 = V_base(4) * t163 * pkin(1) + t195;
t191 = rSges(4,1) * t166 - rSges(4,2) * t165;
t190 = rSges(5,1) * t156 - rSges(5,2) * t155;
t188 = Icges(4,1) * t166 - t208;
t185 = -Icges(4,2) * t165 + t207;
t182 = Icges(4,5) * t166 - Icges(4,6) * t165;
t138 = -qJD(3) * t152 + V_base(5);
t177 = t138 * (-Icges(4,3) * t152 + t151 * t182) + t139 * (Icges(4,3) * t151 + t152 * t182) + (Icges(4,5) * t165 + Icges(4,6) * t166) * t154;
t176 = V_base(5) * pkin(5) - t163 * t194 + t199;
t119 = t152 * pkin(2) + t151 * pkin(6);
t175 = t154 * t119 + t215 * V_base(4) + t200;
t174 = t138 * t217 + t176;
t173 = V_base(4) * t118 + (-t119 - t218) * V_base(5) + t192;
t77 = pkin(7) * t151 + t152 * t216;
t172 = -t139 * t217 + t154 * t77 + t175;
t171 = -t138 * t77 + t139 * t76 + t173;
t100 = Icges(4,5) * t151 + t152 * t188;
t143 = Icges(4,2) * t166 + t208;
t144 = Icges(4,1) * t165 + t207;
t97 = -Icges(4,6) * t152 + t151 * t185;
t98 = Icges(4,6) * t151 + t152 * t185;
t99 = -Icges(4,5) * t152 + t151 * t188;
t168 = (t100 * t166 - t165 * t98) * t139 + (-t165 * t97 + t166 * t99) * t138 + (-t143 * t165 + t144 * t166) * t154;
t153 = Icges(2,4) * t164;
t149 = Icges(3,4) * t152;
t145 = rSges(4,1) * t165 + rSges(4,2) * t166;
t141 = rSges(2,1) * t164 - rSges(2,2) * t163;
t140 = rSges(2,1) * t163 + rSges(2,2) * t164;
t137 = Icges(2,1) * t164 - t210;
t136 = Icges(2,1) * t163 + t153;
t135 = -Icges(2,2) * t163 + t153;
t134 = Icges(2,2) * t164 + t210;
t131 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t130 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t129 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t127 = rSges(5,1) * t155 + rSges(5,2) * t156;
t117 = rSges(3,1) * t152 - rSges(3,2) * t151;
t116 = rSges(3,1) * t151 + rSges(3,2) * t152;
t115 = Icges(3,1) * t152 - t209;
t114 = Icges(3,1) * t151 + t149;
t113 = -Icges(3,2) * t151 + t149;
t112 = Icges(3,2) * t152 + t209;
t104 = V_base(5) * rSges(2,3) - t140 * V_base(6) + t199;
t103 = t141 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t102 = rSges(4,3) * t151 + t152 * t191;
t101 = -rSges(4,3) * t152 + t151 * t191;
t94 = rSges(5,3) * t151 + t152 * t190;
t92 = -rSges(5,3) * t152 + t151 * t190;
t78 = t140 * V_base(4) - t141 * V_base(5) + t195;
t74 = V_base(5) * rSges(3,3) - t116 * t154 + t176;
t73 = t117 * t154 + (-rSges(3,3) + t215) * V_base(4) + t200;
t69 = t116 * V_base(4) + (-t117 - t218) * V_base(5) + t192;
t68 = t138 * t145 + (-t101 - t118) * t154 + t176;
t67 = t102 * t154 - t139 * t145 + t175;
t66 = t101 * t139 - t102 * t138 + t173;
t65 = t106 * t127 + (-t92 + t211) * t154 + t174;
t64 = -t107 * t127 + t154 * t94 + t172;
t63 = -t106 * t94 + t107 * t92 + t171;
t62 = qJD(5) * t151 + t193 * t106 + (t211 - t213) * t154 + t174;
t61 = -qJD(5) * t152 - t107 * t193 + t154 * t212 + t172;
t60 = -t106 * t212 + t107 * t213 + t171;
t1 = m(1) * (t129 ^ 2 + t130 ^ 2 + t131 ^ 2) / 0.2e1 + m(2) * (t103 ^ 2 + t104 ^ 2 + t78 ^ 2) / 0.2e1 + m(3) * (t69 ^ 2 + t73 ^ 2 + t74 ^ 2) / 0.2e1 + m(4) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + t139 * (t151 * t177 + t152 * t168) / 0.2e1 + t138 * (t151 * t168 - t152 * t177) / 0.2e1 + m(5) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + m(6) * (t60 ^ 2 + t61 ^ 2 + t62 ^ 2) / 0.2e1 + (t223 * t151 - t222 * t152) * t106 / 0.2e1 + (t222 * t151 + t223 * t152) * t107 / 0.2e1 + ((-t112 * t151 + t114 * t152 - t134 * t163 + t136 * t164 + Icges(1,4)) * V_base(5) + (-t151 * t113 + t152 * t115 - t163 * t135 + t164 * t137 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t152 * t112 + t151 * t114 + t164 * t134 + t163 * t136 + Icges(1,2)) * V_base(5) + (t113 * t152 + t115 * t151 + t135 * t164 + t137 * t163 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t100 * t165 + t166 * t98) * t139 + (t165 * t99 + t166 * t97) * t138 + (t232 * t155 + t234 * t156) * t107 + (t233 * t155 + t235 * t156) * t106 + (t166 * t143 + t165 * t144 + t229 * t155 + t230 * t156 + Icges(3,3)) * t154) * t154 / 0.2e1 + t154 * V_base(4) * (Icges(3,5) * t152 - Icges(3,6) * t151) + V_base(5) * t154 * (Icges(3,5) * t151 + Icges(3,6) * t152) + ((Icges(2,5) * t163 + Icges(2,6) * t164 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t164 - Icges(2,6) * t163 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T = t1;
