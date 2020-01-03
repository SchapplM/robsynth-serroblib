% Calculate kinetic energy for
% S5RRPRP3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRP3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP3_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPRP3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP3_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP3_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP3_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:50:52
% EndTime: 2019-12-31 19:50:54
% DurationCPUTime: 2.22s
% Computational Cost: add. (1132->218), mult. (911->285), div. (0->0), fcn. (691->8), ass. (0->119)
t244 = Icges(5,4) - Icges(6,5);
t243 = Icges(5,1) + Icges(6,1);
t242 = Icges(5,2) + Icges(6,3);
t162 = pkin(8) + qJ(4);
t155 = cos(t162);
t241 = t244 * t155;
t154 = sin(t162);
t240 = t244 * t154;
t239 = Icges(6,4) + Icges(5,5);
t238 = Icges(5,6) - Icges(6,6);
t237 = t242 * t154 - t241;
t236 = t243 * t155 - t240;
t235 = rSges(6,1) + pkin(4);
t234 = rSges(6,3) + qJ(5);
t163 = qJ(1) + qJ(2);
t157 = sin(t163);
t158 = cos(t163);
t233 = t237 * t157 + t238 * t158;
t232 = -t238 * t157 + t237 * t158;
t231 = t236 * t157 - t239 * t158;
t230 = t239 * t157 + t236 * t158;
t229 = Icges(6,2) + Icges(5,3);
t228 = -t242 * t155 - t240;
t227 = t243 * t154 + t241;
t226 = -t238 * t154 + t239 * t155;
t225 = t234 * t154 + t235 * t155;
t137 = -qJD(4) * t158 + V_base(5);
t138 = qJD(4) * t157 + V_base(4);
t156 = V_base(6) + qJD(1);
t152 = qJD(2) + t156;
t224 = (t228 * t154 + t227 * t155) * t152 + (t232 * t154 + t230 * t155) * t138 + (t233 * t154 + t231 * t155) * t137;
t223 = (t239 * t154 + t238 * t155) * t152 + (t229 * t157 + t226 * t158) * t138 + (t226 * t157 - t229 * t158) * t137;
t222 = -pkin(5) - pkin(6);
t167 = sin(qJ(1));
t218 = pkin(1) * t167;
t168 = cos(qJ(1));
t217 = pkin(1) * t168;
t164 = sin(pkin(8));
t216 = pkin(3) * t164;
t165 = cos(pkin(8));
t215 = pkin(3) * t165;
t214 = -rSges(6,2) * t158 + t225 * t157;
t213 = rSges(6,2) * t157 + t225 * t158;
t126 = pkin(2) * t157 - qJ(3) * t158;
t78 = -pkin(7) * t158 + t215 * t157;
t212 = -t126 - t78;
t211 = Icges(2,4) * t167;
t210 = Icges(3,4) * t157;
t209 = Icges(4,4) * t164;
t208 = Icges(4,4) * t165;
t202 = t235 * t154 - t234 * t155;
t201 = qJD(5) * t154;
t200 = t156 * t217 + V_base(2);
t199 = V_base(4) * t218 + V_base(3);
t198 = V_base(5) * pkin(5) + V_base(1);
t128 = pkin(2) * t158 + qJ(3) * t157;
t195 = -t128 - t217;
t194 = V_base(4) * t126 + t199;
t193 = rSges(4,1) * t165 - rSges(4,2) * t164;
t192 = rSges(5,1) * t155 - rSges(5,2) * t154;
t189 = Icges(4,1) * t165 - t209;
t186 = -Icges(4,2) * t164 + t208;
t183 = Icges(4,5) * t165 - Icges(4,6) * t164;
t180 = -qJD(3) * t158 + t152 * t128 + t200;
t179 = V_base(5) * pkin(6) - t156 * t218 + t198;
t176 = qJD(3) * t157 + t179;
t175 = V_base(5) * t216 + t176;
t174 = (Icges(4,3) * t157 + t183 * t158) * V_base(4) + (Icges(4,5) * t164 + Icges(4,6) * t165) * t152 + (-Icges(4,3) * t158 + t183 * t157) * V_base(5);
t79 = pkin(7) * t157 + t215 * t158;
t173 = V_base(4) * t78 + (t195 - t79) * V_base(5) + t194;
t172 = t152 * t79 + (-t216 + t222) * V_base(4) + t180;
t101 = -Icges(4,6) * t158 + t186 * t157;
t102 = Icges(4,6) * t157 + t186 * t158;
t103 = -Icges(4,5) * t158 + t189 * t157;
t104 = Icges(4,5) * t157 + t189 * t158;
t134 = Icges(4,2) * t165 + t209;
t135 = Icges(4,1) * t164 + t208;
t169 = (-t102 * t164 + t104 * t165) * V_base(4) + (-t101 * t164 + t103 * t165) * V_base(5) + (-t134 * t164 + t135 * t165) * t152;
t159 = Icges(2,4) * t168;
t151 = Icges(3,4) * t158;
t146 = rSges(2,1) * t168 - t167 * rSges(2,2);
t145 = t167 * rSges(2,1) + rSges(2,2) * t168;
t144 = Icges(2,1) * t168 - t211;
t143 = Icges(2,1) * t167 + t159;
t142 = -Icges(2,2) * t167 + t159;
t141 = Icges(2,2) * t168 + t211;
t136 = rSges(4,1) * t164 + rSges(4,2) * t165;
t132 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t131 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t130 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t129 = rSges(3,1) * t158 - rSges(3,2) * t157;
t127 = rSges(3,1) * t157 + rSges(3,2) * t158;
t125 = Icges(3,1) * t158 - t210;
t124 = Icges(3,1) * t157 + t151;
t123 = -Icges(3,2) * t157 + t151;
t122 = Icges(3,2) * t158 + t210;
t121 = Icges(3,5) * t158 - Icges(3,6) * t157;
t120 = Icges(3,5) * t157 + Icges(3,6) * t158;
t119 = rSges(5,1) * t154 + rSges(5,2) * t155;
t106 = rSges(4,3) * t157 + t193 * t158;
t105 = -rSges(4,3) * t158 + t193 * t157;
t98 = V_base(5) * rSges(2,3) - t145 * t156 + t198;
t97 = t146 * t156 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t96 = t145 * V_base(4) - t146 * V_base(5) + V_base(3);
t95 = rSges(5,3) * t157 + t192 * t158;
t93 = -rSges(5,3) * t158 + t192 * t157;
t75 = V_base(5) * rSges(3,3) - t127 * t152 + t179;
t74 = t129 * t152 + (-rSges(3,3) + t222) * V_base(4) + t200;
t73 = V_base(4) * t127 + (-t129 - t217) * V_base(5) + t199;
t72 = t136 * V_base(5) + (-t105 - t126) * t152 + t176;
t71 = t106 * t152 + (-t136 + t222) * V_base(4) + t180;
t70 = V_base(4) * t105 + (-t106 + t195) * V_base(5) + t194;
t69 = t119 * t137 + (-t93 + t212) * t152 + t175;
t68 = -t119 * t138 + t152 * t95 + t172;
t67 = -t137 * t95 + t138 * t93 + t173;
t66 = t158 * t201 + t202 * t137 + (t212 - t214) * t152 + t175;
t65 = -t202 * t138 + t213 * t152 + t157 * t201 + t172;
t64 = -qJD(5) * t155 - t213 * t137 + t214 * t138 + t173;
t1 = m(1) * (t130 ^ 2 + t131 ^ 2 + t132 ^ 2) / 0.2e1 + m(2) * (t96 ^ 2 + t97 ^ 2 + t98 ^ 2) / 0.2e1 + m(3) * (t73 ^ 2 + t74 ^ 2 + t75 ^ 2) / 0.2e1 + m(4) * (t70 ^ 2 + t71 ^ 2 + t72 ^ 2) / 0.2e1 + m(5) * (t67 ^ 2 + t68 ^ 2 + t69 ^ 2) / 0.2e1 + m(6) * (t64 ^ 2 + t65 ^ 2 + t66 ^ 2) / 0.2e1 + (t224 * t157 - t223 * t158) * t137 / 0.2e1 + (t223 * t157 + t224 * t158) * t138 / 0.2e1 + ((t101 * t165 + t103 * t164 + t120) * V_base(5) + (t102 * t165 + t104 * t164 + t121) * V_base(4) + (t230 * t154 - t232 * t155) * t138 + (t231 * t154 - t233 * t155) * t137 + (t165 * t134 + t164 * t135 + t227 * t154 - t228 * t155 + Icges(3,3)) * t152) * t152 / 0.2e1 + (t121 * t152 + t174 * t157 + t169 * t158 + (-t122 * t157 + t124 * t158 - t167 * t141 + t143 * t168 + Icges(1,4)) * V_base(5) + (-t157 * t123 + t158 * t125 - t167 * t142 + t168 * t144 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t120 * t152 + t169 * t157 - t174 * t158 + (t158 * t122 + t157 * t124 + t168 * t141 + t167 * t143 + Icges(1,2)) * V_base(5) + (t123 * t158 + t125 * t157 + t142 * t168 + t167 * t144 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t167 + Icges(2,6) * t168) * V_base(5) + (Icges(2,5) * t168 - Icges(2,6) * t167) * V_base(4) + Icges(2,3) * t156 / 0.2e1) * t156;
T = t1;
