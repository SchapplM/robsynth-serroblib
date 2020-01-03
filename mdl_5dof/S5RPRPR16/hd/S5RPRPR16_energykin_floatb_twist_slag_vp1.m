% Calculate kinetic energy for
% S5RPRPR16
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
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
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR16_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR16_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRPR16_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR16_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR16_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR16_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:38:40
% EndTime: 2019-12-31 18:38:42
% DurationCPUTime: 2.48s
% Computational Cost: add. (641->235), mult. (1154->322), div. (0->0), fcn. (992->6), ass. (0->117)
t251 = -Icges(4,4) - Icges(5,6);
t250 = Icges(4,1) + Icges(5,2);
t249 = Icges(4,2) + Icges(5,3);
t166 = cos(qJ(3));
t248 = t251 * t166;
t163 = sin(qJ(3));
t247 = t251 * t163;
t246 = -Icges(5,4) + Icges(4,5);
t245 = Icges(5,5) - Icges(4,6);
t244 = t249 * t166 - t247;
t243 = t250 * t163 - t248;
t242 = Icges(2,4) + Icges(3,6);
t164 = sin(qJ(1));
t167 = cos(qJ(1));
t241 = -t244 * t164 + t245 * t167;
t240 = t245 * t164 + t244 * t167;
t239 = t243 * t164 + t246 * t167;
t238 = t246 * t164 - t243 * t167;
t237 = Icges(2,1) + Icges(3,2);
t236 = Icges(5,1) + Icges(4,3);
t235 = -Icges(3,4) + Icges(2,5);
t234 = Icges(3,5) - Icges(2,6);
t233 = Icges(2,2) + Icges(3,3);
t232 = t249 * t163 + t248;
t231 = t250 * t166 + t247;
t230 = t246 * t163 - t245 * t166;
t229 = t242 * t167;
t228 = t242 * t164;
t150 = qJD(3) * t164 + V_base(5);
t151 = qJD(3) * t167 + V_base(4);
t154 = V_base(6) + qJD(1);
t227 = t150 * (-t238 * t163 + t240 * t166) + t151 * (-t239 * t163 + t241 * t166) - t154 * (t231 * t163 - t232 * t166);
t226 = -t233 * t167 - t228;
t225 = t233 * t164 - t229;
t224 = t237 * t164 + t229;
t223 = t237 * t167 - t228;
t220 = (t245 * t163 + t246 * t166) * t154 + (t230 * t164 + t236 * t167) * t151 + (t236 * t164 - t230 * t167) * t150;
t213 = pkin(6) * t164;
t212 = pkin(6) * t167;
t204 = t163 * t164;
t203 = t163 * t167;
t202 = t164 * t166;
t201 = t166 * t167;
t200 = qJD(4) * t166;
t199 = qJD(5) * t163;
t140 = t164 * pkin(1) - qJ(2) * t167;
t198 = V_base(4) * t140 + V_base(3);
t197 = V_base(5) * pkin(5) + V_base(1);
t194 = -t140 - t213;
t193 = qJD(2) * t164 + t197;
t188 = pkin(3) * t163 - qJ(4) * t166;
t114 = t188 * t167;
t192 = t114 + t194;
t191 = V_base(5) * pkin(2) + t193;
t190 = rSges(4,1) * t163 + rSges(4,2) * t166;
t189 = rSges(5,2) * t163 + rSges(5,3) * t166;
t146 = pkin(1) * t167 + t164 * qJ(2);
t175 = -qJD(2) * t167 + t154 * t146 + V_base(2);
t143 = pkin(3) * t166 + qJ(4) * t163;
t174 = t150 * t143 + t191;
t171 = V_base(4) * t213 + (-t146 - t212) * V_base(5) + t198;
t170 = t154 * t212 + (-pkin(2) - pkin(5)) * V_base(4) + t175;
t169 = qJD(4) * t163 - t151 * t114 + t171;
t113 = t188 * t164;
t168 = t154 * t113 + t167 * t200 + t170;
t165 = cos(qJ(5));
t162 = sin(qJ(5));
t148 = rSges(2,1) * t167 - t164 * rSges(2,2);
t147 = -rSges(3,2) * t167 + t164 * rSges(3,3);
t145 = rSges(4,1) * t166 - rSges(4,2) * t163;
t144 = -rSges(5,2) * t166 + rSges(5,3) * t163;
t142 = t164 * rSges(2,1) + rSges(2,2) * t167;
t141 = -t164 * rSges(3,2) - rSges(3,3) * t167;
t139 = qJD(5) * t166 + t154;
t120 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t119 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t118 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t117 = pkin(4) * t167 + pkin(7) * t204;
t116 = t164 * pkin(4) - pkin(7) * t203;
t111 = -t162 * t202 + t165 * t167;
t110 = -t162 * t167 - t165 * t202;
t109 = t162 * t201 + t164 * t165;
t108 = -t164 * t162 + t165 * t201;
t107 = t164 * t199 + t151;
t106 = -t167 * t199 + t150;
t104 = rSges(5,1) * t167 - t164 * t189;
t103 = t164 * rSges(5,1) + t167 * t189;
t102 = t164 * rSges(4,3) - t167 * t190;
t101 = rSges(4,3) * t167 + t164 * t190;
t100 = rSges(6,3) * t166 + (rSges(6,1) * t162 + rSges(6,2) * t165) * t163;
t91 = Icges(6,5) * t166 + (Icges(6,1) * t162 + Icges(6,4) * t165) * t163;
t88 = Icges(6,6) * t166 + (Icges(6,4) * t162 + Icges(6,2) * t165) * t163;
t85 = Icges(6,3) * t166 + (Icges(6,5) * t162 + Icges(6,6) * t165) * t163;
t82 = V_base(5) * rSges(2,3) - t142 * t154 + t197;
t81 = t148 * t154 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t80 = t142 * V_base(4) - t148 * V_base(5) + V_base(3);
t79 = rSges(6,1) * t111 + rSges(6,2) * t110 + rSges(6,3) * t204;
t78 = t109 * rSges(6,1) + t108 * rSges(6,2) - rSges(6,3) * t203;
t77 = Icges(6,1) * t111 + Icges(6,4) * t110 + Icges(6,5) * t204;
t76 = Icges(6,1) * t109 + Icges(6,4) * t108 - Icges(6,5) * t203;
t75 = Icges(6,4) * t111 + Icges(6,2) * t110 + Icges(6,6) * t204;
t74 = Icges(6,4) * t109 + Icges(6,2) * t108 - Icges(6,6) * t203;
t73 = Icges(6,5) * t111 + Icges(6,6) * t110 + Icges(6,3) * t204;
t72 = Icges(6,5) * t109 + Icges(6,6) * t108 - Icges(6,3) * t203;
t71 = V_base(5) * rSges(3,1) + (-t140 - t141) * t154 + t193;
t70 = t154 * t147 + (-rSges(3,1) - pkin(5)) * V_base(4) + t175;
t69 = t141 * V_base(4) + (-t146 - t147) * V_base(5) + t198;
t68 = t145 * t150 + (-t102 + t194) * t154 + t191;
t67 = t154 * t101 - t151 * t145 + t170;
t66 = -t150 * t101 + t151 * t102 + t171;
t65 = -t164 * t200 + t144 * t150 + (-t103 + t192) * t154 + t174;
t64 = t154 * t104 + (-t143 - t144) * t151 + t168;
t63 = t151 * t103 + (-t104 - t113) * t150 + t169;
t62 = t100 * t106 - t139 * t78 + (pkin(7) * t150 - qJD(4) * t164) * t166 + (-t116 + t192) * t154 + t174;
t61 = -t107 * t100 + t154 * t117 + t139 * t79 + (-pkin(7) * t166 - t143) * t151 + t168;
t60 = -t106 * t79 + t107 * t78 + t151 * t116 + (-t113 - t117) * t150 + t169;
t1 = m(1) * (t118 ^ 2 + t119 ^ 2 + t120 ^ 2) / 0.2e1 + m(2) * (t80 ^ 2 + t81 ^ 2 + t82 ^ 2) / 0.2e1 + m(3) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + m(4) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + m(5) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + m(6) * (t60 ^ 2 + t61 ^ 2 + t62 ^ 2) / 0.2e1 + t107 * ((t110 * t75 + t111 * t77 + t73 * t204) * t107 + (t110 * t74 + t111 * t76 + t204 * t72) * t106 + (t110 * t88 + t111 * t91 + t204 * t85) * t139) / 0.2e1 + t106 * ((t108 * t75 + t109 * t77 - t203 * t73) * t107 + (t108 * t74 + t109 * t76 - t72 * t203) * t106 + (t108 * t88 + t109 * t91 - t203 * t85) * t139) / 0.2e1 + t139 * ((t72 * t106 + t73 * t107 + t85 * t139) * t166 + ((t162 * t77 + t165 * t75) * t107 + (t162 * t76 + t165 * t74) * t106 + (t162 * t91 + t165 * t88) * t139) * t163) / 0.2e1 + (t220 * t164 + t227 * t167) * t150 / 0.2e1 + (-t227 * t164 + t220 * t167) * t151 / 0.2e1 + ((t164 * t226 + t224 * t167 + Icges(1,4)) * V_base(5) + (t164 * t225 + t167 * t223 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t224 * t164 - t167 * t226 + Icges(1,2)) * V_base(5) + (t164 * t223 - t167 * t225 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t241 * t163 + t239 * t166) * t151 + (t240 * t163 + t238 * t166) * t150 + (t232 * t163 + t231 * t166 + Icges(3,1) + Icges(2,3)) * t154) * t154 / 0.2e1 + t154 * V_base(5) * (t235 * t164 - t234 * t167) + t154 * V_base(4) * (t234 * t164 + t235 * t167) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
