% Calculate kinetic energy for
% S5RPRPP4
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
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
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
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPP4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP4_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRPP4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP4_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP4_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPP4_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:14:14
% EndTime: 2019-12-31 18:14:16
% DurationCPUTime: 1.86s
% Computational Cost: add. (734->204), mult. (931->249), div. (0->0), fcn. (713->6), ass. (0->110)
t259 = -Icges(5,4) + Icges(6,5);
t258 = Icges(5,1) + Icges(6,1);
t257 = Icges(5,2) + Icges(6,3);
t154 = qJ(3) + pkin(7);
t146 = cos(t154);
t256 = t259 * t146;
t145 = sin(t154);
t255 = t259 * t145;
t254 = Icges(6,4) + Icges(5,5);
t253 = Icges(5,6) - Icges(6,6);
t252 = -t257 * t146 + t255;
t251 = t258 * t145 - t256;
t250 = rSges(6,1) + pkin(4);
t249 = rSges(6,3) + qJ(5);
t248 = Icges(2,4) + Icges(3,6);
t157 = sin(qJ(1));
t159 = cos(qJ(1));
t247 = t252 * t157 - t253 * t159;
t246 = -t253 * t157 - t252 * t159;
t245 = t251 * t157 + t254 * t159;
t244 = t254 * t157 - t251 * t159;
t243 = Icges(2,1) + Icges(3,2);
t242 = -Icges(3,4) + Icges(2,5);
t241 = Icges(3,5) - Icges(2,6);
t240 = Icges(2,2) + Icges(3,3);
t239 = t257 * t145 + t256;
t238 = t258 * t146 + t255;
t237 = Icges(6,2) + Icges(4,3) + Icges(5,3);
t156 = sin(qJ(3));
t158 = cos(qJ(3));
t236 = Icges(4,5) * t156 + Icges(4,6) * t158 + t254 * t145 + t253 * t146;
t235 = t248 * t159;
t234 = t250 * t145 - t249 * t146;
t233 = t248 * t157;
t208 = Icges(4,4) * t158;
t126 = -Icges(4,2) * t156 + t208;
t209 = Icges(4,4) * t156;
t131 = Icges(4,1) * t158 - t209;
t142 = qJD(3) * t157 + V_base(5);
t143 = qJD(3) * t159 + V_base(4);
t147 = V_base(6) + qJD(1);
t176 = Icges(4,2) * t158 + t209;
t94 = Icges(4,6) * t159 + t176 * t157;
t95 = Icges(4,6) * t157 - t176 * t159;
t179 = Icges(4,1) * t156 + t208;
t96 = Icges(4,5) * t159 + t179 * t157;
t97 = Icges(4,5) * t157 - t179 * t159;
t232 = t142 * (t244 * t145 - t246 * t146 + t156 * t97 + t158 * t95) + t143 * (t245 * t145 - t247 * t146 + t156 * t96 + t158 * t94) + t147 * (t126 * t158 + t131 * t156 + t238 * t145 - t239 * t146);
t217 = pkin(3) * t156;
t103 = qJ(4) * t159 + t157 * t217;
t231 = qJD(4) * t157 + t147 * t103;
t230 = -t240 * t159 - t233;
t229 = t240 * t157 - t235;
t228 = t243 * t157 + t235;
t227 = t243 * t159 - t233;
t224 = (Icges(4,5) * t158 - Icges(4,6) * t156 - t253 * t145 + t254 * t146) * t147 + (t236 * t157 + t237 * t159) * t143 + (t237 * t157 - t236 * t159) * t142;
t216 = pkin(3) * t158;
t215 = pkin(6) * t159;
t214 = t157 * pkin(6);
t212 = rSges(6,2) * t159 + t234 * t157;
t211 = t157 * rSges(6,2) - t234 * t159;
t201 = t249 * t145 + t250 * t146;
t200 = qJD(2) * t159;
t199 = qJD(5) * t146;
t138 = pkin(1) * t159 + t157 * qJ(2);
t198 = t147 * t138 + V_base(2);
t134 = t157 * pkin(1) - qJ(2) * t159;
t197 = V_base(4) * t134 + V_base(3);
t196 = V_base(5) * pkin(5) + V_base(1);
t193 = -t134 - t214;
t192 = qJD(2) * t157 + t196;
t102 = qJ(4) * t157 - t159 * t217;
t191 = -t102 + t193;
t190 = V_base(5) * pkin(2) + t192;
t189 = rSges(4,1) * t156 + rSges(4,2) * t158;
t188 = rSges(5,1) * t145 + rSges(5,2) * t146;
t167 = qJD(4) * t159 + t142 * t216 + t190;
t163 = t147 * t215 + (-pkin(2) - pkin(5)) * V_base(4) + t198;
t162 = V_base(4) * t214 + (-t138 - t215) * V_base(5) + t197;
t161 = t163 - t200;
t160 = t143 * t102 + t162;
t140 = rSges(2,1) * t159 - t157 * rSges(2,2);
t139 = -rSges(3,2) * t159 + t157 * rSges(3,3);
t137 = rSges(4,1) * t158 - rSges(4,2) * t156;
t136 = t157 * rSges(2,1) + rSges(2,2) * t159;
t135 = -t157 * rSges(3,2) - rSges(3,3) * t159;
t118 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t117 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t116 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t113 = rSges(5,1) * t146 - rSges(5,2) * t145;
t101 = t157 * rSges(4,3) - t189 * t159;
t100 = rSges(4,3) * t159 + t189 * t157;
t90 = t157 * rSges(5,3) - t188 * t159;
t88 = rSges(5,3) * t159 + t188 * t157;
t74 = V_base(5) * rSges(2,3) - t136 * t147 + t196;
t73 = t140 * t147 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t71 = t136 * V_base(4) - t140 * V_base(5) + V_base(3);
t70 = V_base(5) * rSges(3,1) + (-t134 - t135) * t147 + t192;
t69 = -t200 + t147 * t139 + (-rSges(3,1) - pkin(5)) * V_base(4) + t198;
t68 = t135 * V_base(4) + (-t138 - t139) * V_base(5) + t197;
t67 = t137 * t142 + (-t101 + t193) * t147 + t190;
t66 = t147 * t100 - t143 * t137 + t161;
t65 = -t142 * t100 + t143 * t101 + t162;
t64 = t113 * t142 + (t191 - t90) * t147 + t167;
t63 = t147 * t88 + (-t113 - t216) * t143 + t161 + t231;
t62 = t143 * t90 + (-t103 - t88) * t142 + t160;
t61 = -t157 * t199 + t201 * t142 + (t191 - t211) * t147 + t167;
t60 = (-qJD(2) + t199) * t159 + t212 * t147 + (-t201 - t216) * t143 + t163 + t231;
t59 = qJD(5) * t145 + t211 * t143 + (-t103 - t212) * t142 + t160;
t1 = m(1) * (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + m(2) * (t71 ^ 2 + t73 ^ 2 + t74 ^ 2) / 0.2e1 + m(3) * (t68 ^ 2 + t69 ^ 2 + t70 ^ 2) / 0.2e1 + m(4) * (t65 ^ 2 + t66 ^ 2 + t67 ^ 2) / 0.2e1 + m(5) * (t62 ^ 2 + t63 ^ 2 + t64 ^ 2) / 0.2e1 + m(6) * (t59 ^ 2 + t60 ^ 2 + t61 ^ 2) / 0.2e1 + (t224 * t157 - t232 * t159) * t142 / 0.2e1 + (t232 * t157 + t224 * t159) * t143 / 0.2e1 + ((t157 * t230 + t159 * t228 + Icges(1,4)) * V_base(5) + (t229 * t157 + t227 * t159 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t228 * t157 - t230 * t159 + Icges(1,2)) * V_base(5) + (t157 * t227 - t159 * t229 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t247 * t145 + t245 * t146 - t156 * t94 + t158 * t96) * t143 + (t246 * t145 + t244 * t146 - t156 * t95 + t158 * t97) * t142 + (-t156 * t126 + t158 * t131 + t239 * t145 + t238 * t146 + Icges(3,1) + Icges(2,3)) * t147) * t147 / 0.2e1 + t147 * V_base(5) * (t242 * t157 - t241 * t159) + t147 * V_base(4) * (t241 * t157 + t242 * t159) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
