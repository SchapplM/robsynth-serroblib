% Calculate kinetic energy for
% S5RRRRP4
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
% Datum: 2019-12-31 21:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRP4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP4_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP4_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRRP4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP4_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP4_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP4_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP4_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:50:37
% EndTime: 2019-12-31 21:50:38
% DurationCPUTime: 1.94s
% Computational Cost: add. (1216->215), mult. (953->293), div. (0->0), fcn. (733->8), ass. (0->117)
t244 = Icges(5,4) - Icges(6,5);
t243 = Icges(5,1) + Icges(6,1);
t242 = Icges(5,2) + Icges(6,3);
t164 = qJ(3) + qJ(4);
t159 = cos(t164);
t241 = t244 * t159;
t157 = sin(t164);
t240 = t244 * t157;
t239 = Icges(6,4) + Icges(5,5);
t238 = Icges(5,6) - Icges(6,6);
t237 = t242 * t157 - t241;
t236 = t243 * t159 - t240;
t235 = rSges(6,1) + pkin(4);
t234 = rSges(6,3) + qJ(5);
t165 = qJ(1) + qJ(2);
t158 = sin(t165);
t160 = cos(t165);
t233 = t237 * t158 + t238 * t160;
t232 = -t238 * t158 + t237 * t160;
t231 = t236 * t158 - t239 * t160;
t230 = t239 * t158 + t236 * t160;
t229 = Icges(6,2) + Icges(5,3);
t228 = -t242 * t159 - t240;
t227 = t243 * t157 + t241;
t226 = -t238 * t157 + t239 * t159;
t225 = t234 * t157 + t235 * t159;
t110 = V_base(5) + (-qJD(3) - qJD(4)) * t160;
t137 = qJD(3) * t158 + V_base(4);
t111 = qJD(4) * t158 + t137;
t156 = V_base(6) + qJD(1);
t154 = qJD(2) + t156;
t224 = (t228 * t157 + t227 * t159) * t154 + (t232 * t157 + t230 * t159) * t111 + (t233 * t157 + t231 * t159) * t110;
t223 = (t239 * t157 + t238 * t159) * t154 + (t229 * t158 + t226 * t160) * t111 + (t226 * t158 - t229 * t160) * t110;
t222 = -pkin(5) - pkin(6);
t167 = sin(qJ(1));
t218 = pkin(1) * t167;
t169 = cos(qJ(1));
t217 = pkin(1) * t169;
t166 = sin(qJ(3));
t216 = pkin(3) * t166;
t168 = cos(qJ(3));
t215 = pkin(3) * t168;
t213 = -rSges(6,2) * t160 + t158 * t225;
t212 = rSges(6,2) * t158 + t160 * t225;
t131 = t158 * pkin(2) - t160 * pkin(7);
t78 = -pkin(8) * t160 + t215 * t158;
t211 = -t131 - t78;
t210 = Icges(2,4) * t167;
t209 = Icges(3,4) * t158;
t208 = Icges(4,4) * t166;
t207 = Icges(4,4) * t168;
t202 = t235 * t157 - t234 * t159;
t201 = qJD(5) * t157;
t200 = t156 * t217 + V_base(2);
t199 = V_base(4) * t218 + V_base(3);
t198 = V_base(5) * pkin(5) + V_base(1);
t195 = rSges(4,1) * t168 - rSges(4,2) * t166;
t194 = rSges(5,1) * t159 - rSges(5,2) * t157;
t191 = Icges(4,1) * t168 - t208;
t188 = -Icges(4,2) * t166 + t207;
t185 = Icges(4,5) * t168 - Icges(4,6) * t166;
t182 = V_base(5) * pkin(6) - t156 * t218 + t198;
t136 = -qJD(3) * t160 + V_base(5);
t179 = (Icges(4,3) * t158 + t185 * t160) * t137 + (-Icges(4,3) * t160 + t185 * t158) * t136 + (Icges(4,5) * t166 + Icges(4,6) * t168) * t154;
t178 = t136 * t216 + t182;
t132 = t160 * pkin(2) + t158 * pkin(7);
t177 = t154 * t132 + t222 * V_base(4) + t200;
t176 = V_base(4) * t131 + (-t132 - t217) * V_base(5) + t199;
t79 = pkin(8) * t158 + t215 * t160;
t175 = -t137 * t216 + t154 * t79 + t177;
t174 = -t136 * t79 + t137 * t78 + t176;
t101 = -Icges(4,6) * t160 + t188 * t158;
t102 = Icges(4,6) * t158 + t188 * t160;
t103 = -Icges(4,5) * t160 + t191 * t158;
t104 = Icges(4,5) * t158 + t191 * t160;
t141 = Icges(4,2) * t168 + t208;
t144 = Icges(4,1) * t166 + t207;
t171 = (-t102 * t166 + t104 * t168) * t137 + (-t101 * t166 + t103 * t168) * t136 + (-t141 * t166 + t144 * t168) * t154;
t161 = Icges(2,4) * t169;
t153 = Icges(3,4) * t160;
t149 = rSges(2,1) * t169 - rSges(2,2) * t167;
t148 = rSges(2,1) * t167 + rSges(2,2) * t169;
t147 = rSges(4,1) * t166 + rSges(4,2) * t168;
t146 = Icges(2,1) * t169 - t210;
t145 = Icges(2,1) * t167 + t161;
t143 = -Icges(2,2) * t167 + t161;
t142 = Icges(2,2) * t169 + t210;
t135 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t134 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t133 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t130 = rSges(3,1) * t160 - rSges(3,2) * t158;
t129 = rSges(3,1) * t158 + rSges(3,2) * t160;
t128 = rSges(5,1) * t157 + rSges(5,2) * t159;
t125 = Icges(3,1) * t160 - t209;
t124 = Icges(3,1) * t158 + t153;
t121 = -Icges(3,2) * t158 + t153;
t120 = Icges(3,2) * t160 + t209;
t106 = rSges(4,3) * t158 + t195 * t160;
t105 = -rSges(4,3) * t160 + t195 * t158;
t98 = V_base(5) * rSges(2,3) - t148 * t156 + t198;
t97 = t149 * t156 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t96 = t148 * V_base(4) - t149 * V_base(5) + V_base(3);
t95 = rSges(5,3) * t158 + t194 * t160;
t93 = -rSges(5,3) * t160 + t194 * t158;
t75 = V_base(5) * rSges(3,3) - t129 * t154 + t182;
t74 = t130 * t154 + (-rSges(3,3) + t222) * V_base(4) + t200;
t73 = t129 * V_base(4) + (-t130 - t217) * V_base(5) + t199;
t72 = t136 * t147 + (-t105 - t131) * t154 + t182;
t71 = t106 * t154 - t137 * t147 + t177;
t70 = t105 * t137 - t106 * t136 + t176;
t69 = t110 * t128 + (-t93 + t211) * t154 + t178;
t68 = -t111 * t128 + t154 * t95 + t175;
t67 = -t110 * t95 + t111 * t93 + t174;
t66 = t160 * t201 + t202 * t110 + (t211 - t213) * t154 + t178;
t65 = -t202 * t111 + t212 * t154 + t158 * t201 + t175;
t64 = -qJD(5) * t159 - t212 * t110 + t213 * t111 + t174;
t1 = m(1) * (t133 ^ 2 + t134 ^ 2 + t135 ^ 2) / 0.2e1 + m(2) * (t96 ^ 2 + t97 ^ 2 + t98 ^ 2) / 0.2e1 + m(3) * (t73 ^ 2 + t74 ^ 2 + t75 ^ 2) / 0.2e1 + m(4) * (t70 ^ 2 + t71 ^ 2 + t72 ^ 2) / 0.2e1 + t137 * (t179 * t158 + t171 * t160) / 0.2e1 + t136 * (t171 * t158 - t179 * t160) / 0.2e1 + m(5) * (t67 ^ 2 + t68 ^ 2 + t69 ^ 2) / 0.2e1 + m(6) * (t64 ^ 2 + t65 ^ 2 + t66 ^ 2) / 0.2e1 + (t224 * t158 - t223 * t160) * t110 / 0.2e1 + (t223 * t158 + t224 * t160) * t111 / 0.2e1 + ((-t120 * t158 + t124 * t160 - t142 * t167 + t145 * t169 + Icges(1,4)) * V_base(5) + (-t158 * t121 + t160 * t125 - t167 * t143 + t169 * t146 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t160 * t120 + t158 * t124 + t169 * t142 + t167 * t145 + Icges(1,2)) * V_base(5) + (t121 * t160 + t125 * t158 + t143 * t169 + t146 * t167 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t102 * t168 + t104 * t166) * t137 + (t101 * t168 + t103 * t166) * t136 + (t230 * t157 - t232 * t159) * t111 + (t231 * t157 - t233 * t159) * t110 + (t168 * t141 + t166 * t144 + t227 * t157 - t228 * t159 + Icges(3,3)) * t154) * t154 / 0.2e1 + t154 * V_base(4) * (Icges(3,5) * t160 - Icges(3,6) * t158) + V_base(5) * t154 * (Icges(3,5) * t158 + Icges(3,6) * t160) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t167 + Icges(2,6) * t169) * V_base(5) + (Icges(2,5) * t169 - Icges(2,6) * t167) * V_base(4) + Icges(2,3) * t156 / 0.2e1) * t156;
T = t1;
