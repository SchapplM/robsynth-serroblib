% Calculate kinetic energy for
% S5RPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR7_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRPR7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR7_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR7_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR7_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:18:50
% EndTime: 2019-12-31 18:18:52
% DurationCPUTime: 1.92s
% Computational Cost: add. (1361->264), mult. (1148->363), div. (0->0), fcn. (984->10), ass. (0->132)
t241 = Icges(4,3) + Icges(5,3);
t179 = qJ(3) + pkin(9);
t171 = sin(t179);
t173 = cos(t179);
t183 = sin(qJ(3));
t186 = cos(qJ(3));
t240 = Icges(4,5) * t186 + Icges(5,5) * t173 - Icges(4,6) * t183 - Icges(5,6) * t171;
t180 = qJ(1) + pkin(8);
t172 = sin(t180);
t174 = cos(t180);
t221 = Icges(5,4) * t173;
t199 = -Icges(5,2) * t171 + t221;
t100 = -Icges(5,6) * t174 + t172 * t199;
t101 = Icges(5,6) * t172 + t174 * t199;
t222 = Icges(5,4) * t171;
t201 = Icges(5,1) * t173 - t222;
t102 = -Icges(5,5) * t174 + t172 * t201;
t103 = Icges(5,5) * t172 + t174 * t201;
t223 = Icges(4,4) * t186;
t200 = -Icges(4,2) * t183 + t223;
t111 = -Icges(4,6) * t174 + t172 * t200;
t112 = Icges(4,6) * t172 + t174 * t200;
t224 = Icges(4,4) * t183;
t202 = Icges(4,1) * t186 - t224;
t114 = -Icges(4,5) * t174 + t172 * t202;
t115 = Icges(4,5) * t172 + t174 * t202;
t135 = Icges(5,2) * t173 + t222;
t138 = Icges(5,1) * t171 + t221;
t151 = -qJD(3) * t174 + V_base(5);
t152 = qJD(3) * t172 + V_base(4);
t156 = Icges(4,2) * t186 + t224;
t159 = Icges(4,1) * t183 + t223;
t175 = V_base(6) + qJD(1);
t237 = (-t135 * t171 + t138 * t173 - t156 * t183 + t159 * t186) * t175 + (-t101 * t171 + t103 * t173 - t112 * t183 + t115 * t186) * t152 + (-t100 * t171 + t102 * t173 - t111 * t183 + t114 * t186) * t151;
t236 = (Icges(4,5) * t183 + Icges(5,5) * t171 + Icges(4,6) * t186 + Icges(5,6) * t173) * t175 + (t241 * t172 + t240 * t174) * t152 + (t240 * t172 - t241 * t174) * t151;
t184 = sin(qJ(1));
t232 = pkin(1) * t184;
t187 = cos(qJ(1));
t231 = pkin(1) * t187;
t230 = pkin(3) * t183;
t229 = pkin(3) * t186;
t228 = -pkin(5) - qJ(2);
t226 = Icges(2,4) * t184;
t225 = Icges(3,4) * t172;
t220 = t171 * t172;
t219 = t171 * t174;
t182 = sin(qJ(5));
t218 = t172 * t182;
t185 = cos(qJ(5));
t217 = t172 * t185;
t216 = t174 * t182;
t215 = t174 * t185;
t214 = qJD(5) * t171;
t213 = t175 * t231 + V_base(2);
t212 = V_base(5) * pkin(5) + V_base(1);
t145 = pkin(2) * t172 - pkin(6) * t174;
t209 = -t145 - t232;
t208 = V_base(5) * qJ(2) + t212;
t207 = V_base(4) * t232 + qJD(2) + V_base(3);
t96 = -qJ(4) * t174 + t172 * t229;
t206 = t209 - t96;
t205 = pkin(4) * t173 + pkin(7) * t171;
t204 = rSges(4,1) * t186 - rSges(4,2) * t183;
t203 = rSges(5,1) * t173 - rSges(5,2) * t171;
t196 = qJD(4) * t172 + t151 * t230 + t208;
t146 = pkin(2) * t174 + pkin(6) * t172;
t193 = t175 * t146 + t228 * V_base(4) + t213;
t192 = V_base(4) * t145 + (-t146 - t231) * V_base(5) + t207;
t191 = t152 * t96 + t192;
t97 = qJ(4) * t172 + t174 * t229;
t190 = -qJD(4) * t174 + t175 * t97 + t193;
t177 = Icges(2,4) * t187;
t169 = Icges(3,4) * t174;
t164 = rSges(2,1) * t187 - t184 * rSges(2,2);
t163 = t184 * rSges(2,1) + rSges(2,2) * t187;
t162 = rSges(4,1) * t183 + rSges(4,2) * t186;
t161 = Icges(2,1) * t187 - t226;
t160 = Icges(2,1) * t184 + t177;
t158 = -Icges(2,2) * t184 + t177;
t157 = Icges(2,2) * t187 + t226;
t150 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t149 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t148 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t147 = -qJD(5) * t173 + t175;
t144 = pkin(4) * t171 - pkin(7) * t173;
t143 = rSges(3,1) * t174 - rSges(3,2) * t172;
t142 = rSges(3,1) * t172 + rSges(3,2) * t174;
t141 = rSges(5,1) * t171 + rSges(5,2) * t173;
t140 = Icges(3,1) * t174 - t225;
t139 = Icges(3,1) * t172 + t169;
t137 = -Icges(3,2) * t172 + t169;
t136 = Icges(3,2) * t174 + t225;
t128 = t173 * t215 + t218;
t127 = -t173 * t216 + t217;
t126 = t173 * t217 - t216;
t125 = -t173 * t218 - t215;
t124 = t174 * t214 + t152;
t123 = t172 * t214 + t151;
t122 = t205 * t174;
t121 = t205 * t172;
t120 = rSges(4,3) * t172 + t174 * t204;
t119 = -rSges(4,3) * t174 + t172 * t204;
t118 = -rSges(6,3) * t173 + (rSges(6,1) * t185 - rSges(6,2) * t182) * t171;
t117 = V_base(5) * rSges(2,3) - t163 * t175 + t212;
t116 = t164 * t175 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t113 = -Icges(6,5) * t173 + (Icges(6,1) * t185 - Icges(6,4) * t182) * t171;
t110 = -Icges(6,6) * t173 + (Icges(6,4) * t185 - Icges(6,2) * t182) * t171;
t107 = -Icges(6,3) * t173 + (Icges(6,5) * t185 - Icges(6,6) * t182) * t171;
t106 = t163 * V_base(4) - t164 * V_base(5) + V_base(3);
t105 = rSges(5,3) * t172 + t174 * t203;
t104 = -rSges(5,3) * t174 + t172 * t203;
t94 = V_base(5) * rSges(3,3) + (-t142 - t232) * t175 + t208;
t93 = t143 * t175 + (-rSges(3,3) + t228) * V_base(4) + t213;
t91 = V_base(4) * t142 + (-t143 - t231) * V_base(5) + t207;
t90 = rSges(6,1) * t128 + rSges(6,2) * t127 + rSges(6,3) * t219;
t89 = rSges(6,1) * t126 + rSges(6,2) * t125 + rSges(6,3) * t220;
t88 = Icges(6,1) * t128 + Icges(6,4) * t127 + Icges(6,5) * t219;
t87 = Icges(6,1) * t126 + Icges(6,4) * t125 + Icges(6,5) * t220;
t86 = Icges(6,4) * t128 + Icges(6,2) * t127 + Icges(6,6) * t219;
t85 = Icges(6,4) * t126 + Icges(6,2) * t125 + Icges(6,6) * t220;
t84 = Icges(6,5) * t128 + Icges(6,6) * t127 + Icges(6,3) * t219;
t83 = Icges(6,5) * t126 + Icges(6,6) * t125 + Icges(6,3) * t220;
t82 = t151 * t162 + (-t119 + t209) * t175 + t208;
t81 = t120 * t175 - t152 * t162 + t193;
t80 = t152 * t119 - t151 * t120 + t192;
t79 = t141 * t151 + (-t104 + t206) * t175 + t196;
t78 = t105 * t175 + (-t141 - t230) * t152 + t190;
t77 = t152 * t104 + (-t105 - t97) * t151 + t191;
t76 = t118 * t123 + t144 * t151 - t147 * t89 + (-t121 + t206) * t175 + t196;
t75 = -t118 * t124 + t122 * t175 + t147 * t90 + (-t144 - t230) * t152 + t190;
t74 = t152 * t121 - t123 * t90 + t124 * t89 + (-t122 - t97) * t151 + t191;
t1 = m(1) * (t148 ^ 2 + t149 ^ 2 + t150 ^ 2) / 0.2e1 + m(2) * (t106 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + m(3) * (t91 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + m(4) * (t80 ^ 2 + t81 ^ 2 + t82 ^ 2) / 0.2e1 + m(5) * (t77 ^ 2 + t78 ^ 2 + t79 ^ 2) / 0.2e1 + m(6) * (t74 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + t124 * ((t127 * t86 + t128 * t88 + t84 * t219) * t124 + (t127 * t85 + t128 * t87 + t219 * t83) * t123 + (t107 * t219 + t110 * t127 + t113 * t128) * t147) / 0.2e1 + t123 * ((t125 * t86 + t126 * t88 + t220 * t84) * t124 + (t125 * t85 + t126 * t87 + t83 * t220) * t123 + (t107 * t220 + t110 * t125 + t113 * t126) * t147) / 0.2e1 + t147 * ((-t107 * t147 - t83 * t123 - t84 * t124) * t173 + ((-t182 * t86 + t185 * t88) * t124 + (-t182 * t85 + t185 * t87) * t123 + (-t110 * t182 + t113 * t185) * t147) * t171) / 0.2e1 + (t237 * t172 - t236 * t174) * t151 / 0.2e1 + (t236 * t172 + t237 * t174) * t152 / 0.2e1 + ((-t136 * t172 + t139 * t174 - t184 * t157 + t160 * t187 + Icges(1,4)) * V_base(5) + (-t137 * t172 + t140 * t174 - t184 * t158 + t161 * t187 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t136 * t174 + t139 * t172 + t157 * t187 + t184 * t160 + Icges(1,2)) * V_base(5) + (t137 * t174 + t140 * t172 + t158 * t187 + t184 * t161 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t101 * t173 + t103 * t171 + t112 * t186 + t115 * t183) * t152 + (t100 * t173 + t102 * t171 + t111 * t186 + t114 * t183) * t151 + (t135 * t173 + t138 * t171 + t156 * t186 + t159 * t183 + Icges(2,3) + Icges(3,3)) * t175) * t175 / 0.2e1 + t175 * V_base(5) * (Icges(2,5) * t184 + Icges(3,5) * t172 + Icges(2,6) * t187 + Icges(3,6) * t174) + t175 * V_base(4) * (Icges(2,5) * t187 + Icges(3,5) * t174 - Icges(2,6) * t184 - Icges(3,6) * t172) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
