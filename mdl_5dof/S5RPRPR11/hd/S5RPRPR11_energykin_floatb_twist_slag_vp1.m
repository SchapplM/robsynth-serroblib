% Calculate kinetic energy for
% S5RPRPR11
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR11_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR11_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRPR11_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR11_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR11_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR11_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:27:06
% EndTime: 2019-12-31 18:27:09
% DurationCPUTime: 2.80s
% Computational Cost: add. (1091->255), mult. (1320->349), div. (0->0), fcn. (1180->8), ass. (0->132)
t258 = Icges(4,4) - Icges(5,5);
t257 = Icges(4,1) + Icges(5,1);
t256 = Icges(4,2) + Icges(5,3);
t182 = pkin(8) + qJ(3);
t175 = sin(t182);
t255 = t258 * t175;
t176 = cos(t182);
t254 = t258 * t176;
t253 = Icges(5,4) + Icges(4,5);
t252 = Icges(4,6) - Icges(5,6);
t251 = t175 * t256 - t254;
t250 = t176 * t257 - t255;
t249 = Icges(5,2) + Icges(4,3);
t187 = sin(qJ(1));
t189 = cos(qJ(1));
t248 = t187 * t251 + t189 * t252;
t247 = -t187 * t252 + t189 * t251;
t246 = t187 * t250 - t189 * t253;
t245 = t187 * t253 + t189 * t250;
t244 = -t176 * t256 - t255;
t243 = t175 * t257 + t254;
t242 = -t175 * t252 + t176 * t253;
t171 = -qJD(3) * t189 + V_base(5);
t172 = qJD(3) * t187 + V_base(4);
t177 = V_base(6) + qJD(1);
t241 = (t175 * t244 + t176 * t243) * t177 + (t175 * t247 + t176 * t245) * t172 + (t175 * t248 + t176 * t246) * t171;
t240 = (t175 * t253 + t176 * t252) * t177 + (t187 * t249 + t189 * t242) * t172 + (t187 * t242 - t189 * t249) * t171;
t183 = sin(pkin(8));
t236 = pkin(2) * t183;
t235 = pkin(4) * t175;
t234 = pkin(4) * t176;
t184 = cos(pkin(8));
t233 = pkin(2) * t184;
t232 = Icges(2,4) * t187;
t231 = Icges(3,4) * t183;
t230 = Icges(3,4) * t184;
t101 = -pkin(6) * t189 + t187 * t233;
t167 = t187 * pkin(1) - qJ(2) * t189;
t224 = -t101 - t167;
t223 = qJD(4) * t175;
t222 = t167 * V_base(4) + V_base(3);
t221 = V_base(5) * pkin(5) + V_base(1);
t212 = pkin(3) * t176 + qJ(4) * t175;
t133 = t212 * t187;
t218 = -t133 + t224;
t217 = qJD(2) * t187 + t221;
t216 = t236 * V_base(5) + t217;
t215 = rSges(3,1) * t184 - rSges(3,2) * t183;
t214 = rSges(4,1) * t176 - rSges(4,2) * t175;
t213 = rSges(5,1) * t176 + rSges(5,3) * t175;
t211 = Icges(3,1) * t184 - t231;
t208 = -Icges(3,2) * t183 + t230;
t205 = Icges(3,5) * t184 - Icges(3,6) * t183;
t186 = sin(qJ(5));
t188 = cos(qJ(5));
t136 = t175 * t188 - t176 * t186;
t202 = t175 * t186 + t176 * t188;
t169 = pkin(1) * t189 + qJ(2) * t187;
t201 = -qJD(2) * t189 + t169 * t177 + V_base(2);
t146 = pkin(3) * t175 - qJ(4) * t176;
t200 = t146 * t171 + t189 * t223 + t216;
t102 = pkin(6) * t187 + t189 * t233;
t197 = V_base(4) * t101 + (-t102 - t169) * V_base(5) + t222;
t196 = (-Icges(3,3) * t189 + t187 * t205) * V_base(5) + (Icges(3,3) * t187 + t189 * t205) * V_base(4) + (Icges(3,5) * t183 + Icges(3,6) * t184) * t177;
t195 = -qJD(4) * t176 + t133 * t172 + t197;
t194 = t177 * t102 + (-pkin(5) - t236) * V_base(4) + t201;
t134 = t212 * t189;
t193 = t134 * t177 + t187 * t223 + t194;
t127 = -Icges(3,6) * t189 + t187 * t208;
t128 = Icges(3,6) * t187 + t189 * t208;
t129 = -Icges(3,5) * t189 + t187 * t211;
t130 = Icges(3,5) * t187 + t189 * t211;
t156 = Icges(3,2) * t184 + t231;
t157 = Icges(3,1) * t183 + t230;
t190 = (-t128 * t183 + t130 * t184) * V_base(4) + (-t127 * t183 + t129 * t184) * V_base(5) + (-t156 * t183 + t157 * t184) * t177;
t180 = Icges(2,4) * t189;
t170 = rSges(2,1) * t189 - rSges(2,2) * t187;
t168 = rSges(2,1) * t187 + rSges(2,2) * t189;
t164 = Icges(2,1) * t189 - t232;
t163 = Icges(2,1) * t187 + t180;
t162 = -Icges(2,2) * t187 + t180;
t161 = Icges(2,2) * t189 + t232;
t160 = Icges(2,5) * t189 - Icges(2,6) * t187;
t159 = Icges(2,5) * t187 + Icges(2,6) * t189;
t158 = rSges(3,1) * t183 + rSges(3,2) * t184;
t154 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t153 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t152 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t150 = -qJD(5) * t187 + t172;
t149 = V_base(5) + (-qJD(3) + qJD(5)) * t189;
t148 = rSges(4,1) * t175 + rSges(4,2) * t176;
t147 = rSges(5,1) * t175 - rSges(5,3) * t176;
t145 = -pkin(7) * t187 + t189 * t234;
t144 = pkin(7) * t189 + t187 * t234;
t132 = rSges(3,3) * t187 + t189 * t215;
t131 = -rSges(3,3) * t189 + t187 * t215;
t124 = t202 * t189;
t123 = t136 * t189;
t122 = t202 * t187;
t121 = t136 * t187;
t119 = rSges(4,3) * t187 + t189 * t214;
t118 = rSges(5,2) * t187 + t189 * t213;
t117 = -rSges(4,3) * t189 + t187 * t214;
t116 = -rSges(5,2) * t189 + t187 * t213;
t100 = V_base(5) * rSges(2,3) - t168 * t177 + t221;
t99 = t170 * t177 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t97 = t168 * V_base(4) - t170 * V_base(5) + V_base(3);
t94 = rSges(6,1) * t136 - rSges(6,2) * t202;
t93 = Icges(6,1) * t136 - Icges(6,4) * t202;
t92 = Icges(6,4) * t136 - Icges(6,2) * t202;
t91 = Icges(6,5) * t136 - Icges(6,6) * t202;
t90 = rSges(6,1) * t124 + rSges(6,2) * t123 - rSges(6,3) * t187;
t89 = rSges(6,1) * t122 + rSges(6,2) * t121 + rSges(6,3) * t189;
t88 = Icges(6,1) * t124 + Icges(6,4) * t123 - Icges(6,5) * t187;
t87 = Icges(6,1) * t122 + Icges(6,4) * t121 + Icges(6,5) * t189;
t86 = Icges(6,4) * t124 + Icges(6,2) * t123 - Icges(6,6) * t187;
t85 = Icges(6,4) * t122 + Icges(6,2) * t121 + Icges(6,6) * t189;
t84 = Icges(6,5) * t124 + Icges(6,6) * t123 - Icges(6,3) * t187;
t83 = Icges(6,5) * t122 + Icges(6,6) * t121 + Icges(6,3) * t189;
t82 = t158 * V_base(5) + (-t131 - t167) * t177 + t217;
t81 = t177 * t132 + (-pkin(5) - t158) * V_base(4) + t201;
t80 = t131 * V_base(4) + (-t132 - t169) * V_base(5) + t222;
t79 = t148 * t171 + (-t117 + t224) * t177 + t216;
t78 = t119 * t177 - t148 * t172 + t194;
t77 = t117 * t172 - t119 * t171 + t197;
t76 = t147 * t171 + (-t116 + t218) * t177 + t200;
t75 = t177 * t118 + (-t146 - t147) * t172 + t193;
t74 = t116 * t172 + (-t118 - t134) * t171 + t195;
t73 = t171 * t235 + t149 * t94 + (-t144 - t89 + t218) * t177 + t200;
t72 = -t150 * t94 + (t145 + t90) * t177 + (-t146 - t235) * t172 + t193;
t71 = t144 * t172 - t149 * t90 + t150 * t89 + (-t134 - t145) * t171 + t195;
t1 = m(1) * (t152 ^ 2 + t153 ^ 2 + t154 ^ 2) / 0.2e1 + m(2) * (t100 ^ 2 + t97 ^ 2 + t99 ^ 2) / 0.2e1 + m(3) * (t80 ^ 2 + t81 ^ 2 + t82 ^ 2) / 0.2e1 + m(4) * (t77 ^ 2 + t78 ^ 2 + t79 ^ 2) / 0.2e1 + m(5) * (t74 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + m(6) * (t71 ^ 2 + t72 ^ 2 + t73 ^ 2) / 0.2e1 + t150 * ((t123 * t86 + t124 * t88 - t187 * t84) * t150 + (t123 * t85 + t124 * t87 - t187 * t83) * t149 + (t123 * t92 + t124 * t93 - t187 * t91) * t177) / 0.2e1 + t149 * ((t121 * t86 + t122 * t88 + t189 * t84) * t150 + (t121 * t85 + t122 * t87 + t189 * t83) * t149 + (t121 * t92 + t122 * t93 + t189 * t91) * t177) / 0.2e1 + (t241 * t187 - t240 * t189) * t171 / 0.2e1 + (t240 * t187 + t241 * t189) * t172 / 0.2e1 + (t160 * t177 + t196 * t187 + t190 * t189 + (-t161 * t187 + t163 * t189 + Icges(1,4)) * V_base(5) + (-t187 * t162 + t164 * t189 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t159 * t177 + t190 * t187 - t196 * t189 + (t161 * t189 + t187 * t163 + Icges(1,2)) * V_base(5) + (t162 * t189 + t164 * t187 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t136 * t88 - t202 * t86) * t150 + (t136 * t87 - t202 * t85) * t149 + (t127 * t184 + t129 * t183 + t159) * V_base(5) + (t128 * t184 + t130 * t183 + t160) * V_base(4) + (t175 * t245 - t176 * t247) * t172 + (t175 * t246 - t176 * t248) * t171 + (t136 * t93 + t156 * t184 + t157 * t183 + t243 * t175 - t244 * t176 - t202 * t92 + Icges(2,3)) * t177) * t177 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
