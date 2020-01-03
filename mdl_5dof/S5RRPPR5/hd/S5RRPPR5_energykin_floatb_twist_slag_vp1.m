% Calculate kinetic energy for
% S5RRPPR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPPR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR5_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPPR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR5_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR5_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR5_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:28:49
% EndTime: 2019-12-31 19:28:51
% DurationCPUTime: 2.54s
% Computational Cost: add. (1113->248), mult. (1342->335), div. (0->0), fcn. (1202->8), ass. (0->128)
t259 = Icges(4,4) - Icges(5,5);
t258 = Icges(4,1) + Icges(5,1);
t257 = Icges(4,2) + Icges(5,3);
t182 = qJ(2) + pkin(8);
t176 = cos(t182);
t256 = t259 * t176;
t175 = sin(t182);
t255 = t259 * t175;
t254 = Icges(5,4) + Icges(4,5);
t253 = Icges(4,6) - Icges(5,6);
t252 = t257 * t175 - t256;
t251 = t258 * t176 - t255;
t186 = sin(qJ(1));
t189 = cos(qJ(1));
t250 = t252 * t186 + t253 * t189;
t249 = -t253 * t186 + t252 * t189;
t248 = t251 * t186 - t254 * t189;
t247 = t254 * t186 + t251 * t189;
t246 = -t257 * t176 - t255;
t245 = t258 * t175 + t256;
t244 = Icges(5,2) + Icges(3,3) + Icges(4,3);
t185 = sin(qJ(2));
t188 = cos(qJ(2));
t243 = Icges(3,5) * t188 - Icges(3,6) * t185 - t253 * t175 + t254 * t176;
t230 = Icges(3,4) * t188;
t209 = -Icges(3,2) * t185 + t230;
t127 = -Icges(3,6) * t189 + t209 * t186;
t128 = Icges(3,6) * t186 + t209 * t189;
t231 = Icges(3,4) * t185;
t212 = Icges(3,1) * t188 - t231;
t129 = -Icges(3,5) * t189 + t212 * t186;
t130 = Icges(3,5) * t186 + t212 * t189;
t159 = Icges(3,2) * t188 + t231;
t162 = Icges(3,1) * t185 + t230;
t172 = -qJD(2) * t189 + V_base(5);
t173 = qJD(2) * t186 + V_base(4);
t177 = V_base(6) + qJD(1);
t242 = (-t159 * t185 + t162 * t188 + t246 * t175 + t245 * t176) * t177 + (-t128 * t185 + t130 * t188 + t249 * t175 + t247 * t176) * t173 + (-t127 * t185 + t129 * t188 + t250 * t175 + t248 * t176) * t172;
t241 = (Icges(3,5) * t185 + Icges(3,6) * t188 + t254 * t175 + t253 * t176) * t177 + (t244 * t186 + t243 * t189) * t173 + (t243 * t186 - t244 * t189) * t172;
t237 = pkin(2) * t185;
t236 = pkin(4) * t175;
t235 = pkin(4) * t176;
t234 = pkin(2) * t188;
t232 = Icges(2,4) * t186;
t102 = -qJ(3) * t189 + t234 * t186;
t170 = t186 * pkin(1) - pkin(6) * t189;
t225 = -t102 - t170;
t103 = qJ(3) * t186 + t234 * t189;
t213 = pkin(3) * t176 + qJ(4) * t175;
t132 = t213 * t189;
t224 = -t103 - t132;
t223 = qJD(4) * t175;
t222 = V_base(5) * pkin(5) + V_base(1);
t131 = t213 * t186;
t219 = -t131 + t225;
t146 = pkin(3) * t175 - qJ(4) * t176;
t218 = -t146 - t237;
t217 = qJD(3) * t186 + t172 * t237 + t222;
t216 = rSges(3,1) * t188 - rSges(3,2) * t185;
t215 = rSges(4,1) * t176 - rSges(4,2) * t175;
t214 = rSges(5,1) * t176 + rSges(5,3) * t175;
t184 = sin(qJ(5));
t187 = cos(qJ(5));
t136 = t175 * t187 - t176 * t184;
t203 = t175 * t184 + t176 * t187;
t171 = pkin(1) * t189 + t186 * pkin(6);
t202 = -V_base(4) * pkin(5) + t177 * t171 + V_base(2);
t201 = V_base(4) * t170 - t171 * V_base(5) + V_base(3);
t200 = t172 * t146 + t189 * t223 + t217;
t199 = t173 * t102 + t201;
t195 = -qJD(3) * t189 + t177 * t103 + t202;
t194 = -qJD(4) * t176 + t173 * t131 + t199;
t193 = t177 * t132 + t186 * t223 + t195;
t180 = Icges(2,4) * t189;
t169 = rSges(2,1) * t189 - t186 * rSges(2,2);
t168 = t186 * rSges(2,1) + rSges(2,2) * t189;
t167 = rSges(3,1) * t185 + rSges(3,2) * t188;
t164 = Icges(2,1) * t189 - t232;
t163 = Icges(2,1) * t186 + t180;
t161 = -Icges(2,2) * t186 + t180;
t160 = Icges(2,2) * t189 + t232;
t155 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t154 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t153 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t150 = -qJD(5) * t186 + t173;
t149 = V_base(5) + (-qJD(2) + qJD(5)) * t189;
t148 = rSges(4,1) * t175 + rSges(4,2) * t176;
t147 = rSges(5,1) * t175 - rSges(5,3) * t176;
t145 = -t186 * pkin(7) + t189 * t235;
t144 = pkin(7) * t189 + t186 * t235;
t134 = t186 * rSges(3,3) + t216 * t189;
t133 = -rSges(3,3) * t189 + t216 * t186;
t124 = t203 * t189;
t123 = t136 * t189;
t122 = t203 * t186;
t121 = t136 * t186;
t119 = t186 * rSges(4,3) + t215 * t189;
t118 = t186 * rSges(5,2) + t214 * t189;
t117 = -rSges(4,3) * t189 + t215 * t186;
t116 = -rSges(5,2) * t189 + t214 * t186;
t100 = V_base(5) * rSges(2,3) - t168 * t177 + t222;
t99 = t169 * t177 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t97 = t168 * V_base(4) - t169 * V_base(5) + V_base(3);
t94 = rSges(6,1) * t136 - rSges(6,2) * t203;
t93 = Icges(6,1) * t136 - Icges(6,4) * t203;
t92 = Icges(6,4) * t136 - Icges(6,2) * t203;
t91 = Icges(6,5) * t136 - Icges(6,6) * t203;
t90 = rSges(6,1) * t124 + rSges(6,2) * t123 - rSges(6,3) * t186;
t89 = t122 * rSges(6,1) + t121 * rSges(6,2) + rSges(6,3) * t189;
t88 = Icges(6,1) * t124 + Icges(6,4) * t123 - Icges(6,5) * t186;
t87 = Icges(6,1) * t122 + Icges(6,4) * t121 + Icges(6,5) * t189;
t86 = Icges(6,4) * t124 + Icges(6,2) * t123 - Icges(6,6) * t186;
t85 = Icges(6,4) * t122 + Icges(6,2) * t121 + Icges(6,6) * t189;
t84 = Icges(6,5) * t124 + Icges(6,6) * t123 - Icges(6,3) * t186;
t83 = Icges(6,5) * t122 + Icges(6,6) * t121 + Icges(6,3) * t189;
t82 = t167 * t172 + (-t133 - t170) * t177 + t222;
t81 = t134 * t177 - t167 * t173 + t202;
t80 = t133 * t173 - t134 * t172 + t201;
t79 = t148 * t172 + (-t117 + t225) * t177 + t217;
t78 = t177 * t119 + (-t148 - t237) * t173 + t195;
t77 = t117 * t173 + (-t103 - t119) * t172 + t199;
t76 = t147 * t172 + (-t116 + t219) * t177 + t200;
t75 = t177 * t118 + (-t147 + t218) * t173 + t193;
t74 = t116 * t173 + (-t118 + t224) * t172 + t194;
t73 = t172 * t236 + t149 * t94 + (-t144 - t89 + t219) * t177 + t200;
t72 = -t150 * t94 + (t145 + t90) * t177 + (t218 - t236) * t173 + t193;
t71 = t144 * t173 - t149 * t90 + t150 * t89 + (-t145 + t224) * t172 + t194;
t1 = m(1) * (t153 ^ 2 + t154 ^ 2 + t155 ^ 2) / 0.2e1 + m(2) * (t100 ^ 2 + t97 ^ 2 + t99 ^ 2) / 0.2e1 + m(3) * (t80 ^ 2 + t81 ^ 2 + t82 ^ 2) / 0.2e1 + m(4) * (t77 ^ 2 + t78 ^ 2 + t79 ^ 2) / 0.2e1 + m(5) * (t74 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + m(6) * (t71 ^ 2 + t72 ^ 2 + t73 ^ 2) / 0.2e1 + t150 * ((t123 * t86 + t124 * t88 - t186 * t84) * t150 + (t123 * t85 + t124 * t87 - t186 * t83) * t149 + (t123 * t92 + t124 * t93 - t186 * t91) * t177) / 0.2e1 + t149 * ((t121 * t86 + t122 * t88 + t189 * t84) * t150 + (t121 * t85 + t122 * t87 + t189 * t83) * t149 + (t121 * t92 + t122 * t93 + t189 * t91) * t177) / 0.2e1 + ((-t186 * t160 + t163 * t189 + Icges(1,4)) * V_base(5) + (-t186 * t161 + t164 * t189 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t160 * t189 + t186 * t163 + Icges(1,2)) * V_base(5) + (t161 * t189 + t186 * t164 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (t242 * t186 - t241 * t189) * t172 / 0.2e1 + (t241 * t186 + t242 * t189) * t173 / 0.2e1 + ((t136 * t88 - t203 * t86) * t150 + (t136 * t87 - t203 * t85) * t149 + (t128 * t188 + t130 * t185 + t247 * t175 - t249 * t176) * t173 + (t127 * t188 + t129 * t185 + t248 * t175 - t250 * t176) * t172 + (t136 * t93 + t159 * t188 + t162 * t185 + t245 * t175 - t246 * t176 - t203 * t92 + Icges(2,3)) * t177) * t177 / 0.2e1 + V_base(4) * t177 * (Icges(2,5) * t189 - Icges(2,6) * t186) + V_base(5) * t177 * (Icges(2,5) * t186 + Icges(2,6) * t189) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
