% Calculate kinetic energy for
% S5RPPRP6
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
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
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRP6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP6_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPRP6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP6_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP6_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP6_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:54:54
% EndTime: 2019-12-31 17:54:57
% DurationCPUTime: 2.27s
% Computational Cost: add. (716->211), mult. (913->260), div. (0->0), fcn. (695->6), ass. (0->117)
t261 = -Icges(5,4) + Icges(6,5);
t260 = Icges(5,1) + Icges(6,1);
t259 = Icges(5,2) + Icges(6,3);
t154 = pkin(7) + qJ(4);
t146 = cos(t154);
t258 = t261 * t146;
t145 = sin(t154);
t257 = t261 * t145;
t256 = Icges(6,4) + Icges(5,5);
t255 = Icges(5,6) - Icges(6,6);
t254 = -t259 * t146 + t257;
t253 = t260 * t145 - t258;
t252 = rSges(6,1) + pkin(4);
t251 = rSges(6,3) + qJ(5);
t250 = Icges(2,4) + Icges(3,6);
t158 = sin(qJ(1));
t159 = cos(qJ(1));
t249 = t254 * t158 - t255 * t159;
t248 = -t255 * t158 - t254 * t159;
t247 = t253 * t158 + t256 * t159;
t246 = t256 * t158 - t253 * t159;
t245 = Icges(2,1) + Icges(3,2);
t244 = -Icges(3,4) + Icges(2,5);
t243 = Icges(3,5) - Icges(2,6);
t242 = Icges(2,2) + Icges(3,3);
t241 = Icges(6,2) + Icges(5,3);
t240 = t259 * t145 + t258;
t239 = t260 * t146 + t257;
t238 = t256 * t145 + t255 * t146;
t237 = t250 * t159;
t236 = t252 * t145 - t251 * t146;
t235 = t250 * t158;
t141 = qJD(4) * t158 + V_base(5);
t142 = qJD(4) * t159 + V_base(4);
t147 = V_base(6) + qJD(1);
t234 = (t246 * t145 - t248 * t146) * t141 + (t247 * t145 - t249 * t146) * t142 + (t239 * t145 - t240 * t146) * t147;
t233 = t245 * t158 + t237;
t232 = t245 * t159 - t235;
t231 = t244 * t158 - t243 * t159;
t230 = t243 * t158 + t244 * t159;
t205 = qJ(3) * t159;
t229 = qJD(3) * t158 + t147 * t205;
t228 = (-t255 * t145 + t256 * t146) * t147 + (t238 * t158 + t241 * t159) * t142 + (t241 * t158 - t238 * t159) * t141;
t155 = sin(pkin(7));
t156 = cos(pkin(7));
t213 = Icges(4,4) * t155;
t176 = Icges(4,2) * t156 + t213;
t95 = Icges(4,6) * t158 - t176 * t159;
t212 = Icges(4,4) * t156;
t179 = Icges(4,1) * t155 + t212;
t97 = Icges(4,5) * t158 - t179 * t159;
t224 = t155 * t97 + t156 * t95 - t242 * t159 - t235;
t94 = Icges(4,6) * t159 + t176 * t158;
t96 = Icges(4,5) * t159 + t179 * t158;
t223 = t155 * t96 + t156 * t94 + t242 * t158 - t237;
t222 = -pkin(2) - pkin(5);
t218 = pkin(3) * t155;
t217 = pkin(3) * t156;
t216 = rSges(6,2) * t159 + t236 * t158;
t215 = t158 * rSges(6,2) - t236 * t159;
t204 = t158 * qJ(3);
t202 = t251 * t145 + t252 * t146;
t201 = qJD(2) * t159;
t200 = qJD(5) * t146;
t137 = pkin(1) * t159 + t158 * qJ(2);
t199 = t147 * t137 + V_base(2);
t134 = t158 * pkin(1) - qJ(2) * t159;
t198 = V_base(4) * t134 + V_base(3);
t197 = V_base(5) * pkin(5) + V_base(1);
t194 = V_base(4) * t204 + t198;
t193 = qJD(2) * t158 + t197;
t192 = -t134 - t204;
t191 = -t137 - t205;
t102 = pkin(6) * t158 - t159 * t218;
t190 = -t102 + t192;
t189 = rSges(4,1) * t155 + rSges(4,2) * t156;
t188 = rSges(5,1) * t145 + rSges(5,2) * t146;
t173 = Icges(4,5) * t155 + Icges(4,6) * t156;
t119 = -Icges(4,2) * t155 + t212;
t120 = Icges(4,1) * t156 - t213;
t168 = t119 * t156 + t120 * t155;
t167 = t199 - t201;
t166 = V_base(5) * pkin(2) + qJD(3) * t159 + t193;
t165 = V_base(5) * t217 + t166;
t162 = (Icges(4,5) * t156 - Icges(4,6) * t155) * t147 + (Icges(4,3) * t159 + t173 * t158) * V_base(4) + (Icges(4,3) * t158 - t173 * t159) * V_base(5);
t103 = pkin(6) * t159 + t158 * t218;
t161 = V_base(4) * t102 + (-t103 + t191) * V_base(5) + t194;
t160 = t147 * t103 + (-t217 + t222) * V_base(4) + t199 + t229;
t139 = rSges(2,1) * t159 - t158 * rSges(2,2);
t138 = -rSges(3,2) * t159 + t158 * rSges(3,3);
t136 = t158 * rSges(2,1) + rSges(2,2) * t159;
t135 = -t158 * rSges(3,2) - rSges(3,3) * t159;
t121 = rSges(4,1) * t156 - rSges(4,2) * t155;
t117 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t116 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t115 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t113 = rSges(5,1) * t146 - rSges(5,2) * t145;
t99 = t158 * rSges(4,3) - t189 * t159;
t98 = rSges(4,3) * t159 + t189 * t158;
t90 = t158 * rSges(5,3) - t188 * t159;
t88 = rSges(5,3) * t159 + t188 * t158;
t73 = V_base(5) * rSges(2,3) - t136 * t147 + t197;
t72 = t139 * t147 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t71 = t136 * V_base(4) - t139 * V_base(5) + V_base(3);
t70 = V_base(5) * rSges(3,1) + (-t134 - t135) * t147 + t193;
t69 = t147 * t138 + (-rSges(3,1) - pkin(5)) * V_base(4) + t167;
t68 = t135 * V_base(4) + (-t137 - t138) * V_base(5) + t198;
t67 = t121 * V_base(5) + (t192 - t99) * t147 + t166;
t66 = t147 * t98 + (-t121 + t222) * V_base(4) + t167 + t229;
t65 = V_base(4) * t99 + (t191 - t98) * V_base(5) + t194;
t64 = t113 * t141 + (t190 - t90) * t147 + t165;
t63 = -t142 * t113 + t147 * t88 + t160 - t201;
t62 = -t141 * t88 + t142 * t90 + t161;
t61 = -t158 * t200 + t202 * t141 + (t190 - t215) * t147 + t165;
t60 = (-qJD(2) + t200) * t159 + t216 * t147 - t202 * t142 + t160;
t59 = qJD(5) * t145 - t216 * t141 + t215 * t142 + t161;
t1 = m(1) * (t115 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + m(2) * (t71 ^ 2 + t72 ^ 2 + t73 ^ 2) / 0.2e1 + m(3) * (t68 ^ 2 + t69 ^ 2 + t70 ^ 2) / 0.2e1 + m(4) * (t65 ^ 2 + t66 ^ 2 + t67 ^ 2) / 0.2e1 + m(5) * (t62 ^ 2 + t63 ^ 2 + t64 ^ 2) / 0.2e1 + m(6) * (t59 ^ 2 + t60 ^ 2 + t61 ^ 2) / 0.2e1 + (t228 * t158 - t234 * t159) * t141 / 0.2e1 + (t234 * t158 + t228 * t159) * t142 / 0.2e1 + (t162 * t159 + (t168 * t158 + t230) * t147 + (t224 * t158 + t233 * t159 + Icges(1,4)) * V_base(5) + (t223 * t158 + t232 * t159 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t162 * t158 + (-t168 * t159 + t231) * t147 + (t233 * t158 - t224 * t159 + Icges(1,2)) * V_base(5) + (t232 * t158 - t223 * t159 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t249 * t145 + t247 * t146) * t142 + (t248 * t145 + t246 * t146) * t141 + (-t155 * t95 + t156 * t97 + t231) * V_base(5) + (-t155 * t94 + t156 * t96 + t230) * V_base(4) + (-t155 * t119 + t156 * t120 + t240 * t145 + t239 * t146 + Icges(3,1) + Icges(2,3)) * t147) * t147 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
