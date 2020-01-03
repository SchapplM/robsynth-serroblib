% Calculate kinetic energy for
% S5RPRPR9
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
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR9_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR9_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRPR9_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR9_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR9_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR9_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:23:38
% EndTime: 2019-12-31 18:23:40
% DurationCPUTime: 2.48s
% Computational Cost: add. (1083->245), mult. (1144->336), div. (0->0), fcn. (980->8), ass. (0->122)
t249 = Icges(4,4) + Icges(5,6);
t248 = Icges(4,1) + Icges(5,2);
t247 = Icges(4,2) + Icges(5,3);
t182 = cos(qJ(3));
t246 = t249 * t182;
t179 = sin(qJ(3));
t245 = t249 * t179;
t244 = Icges(5,4) - Icges(4,5);
t243 = Icges(5,5) - Icges(4,6);
t242 = t179 * t247 - t246;
t241 = t182 * t248 - t245;
t177 = qJ(1) + pkin(8);
t171 = sin(t177);
t172 = cos(t177);
t240 = t243 * t171 + t242 * t172;
t239 = t242 * t171 - t243 * t172;
t238 = Icges(5,1) + Icges(4,3);
t237 = t241 * t171 + t244 * t172;
t236 = -t244 * t171 + t241 * t172;
t235 = -t182 * t247 - t245;
t234 = t179 * t248 + t246;
t233 = t243 * t179 - t244 * t182;
t145 = -qJD(3) * t172 + V_base(5);
t146 = qJD(3) * t171 + V_base(4);
t173 = V_base(6) + qJD(1);
t230 = (t235 * t179 + t234 * t182) * t173 + (t240 * t179 + t236 * t182) * t146 + (t239 * t179 + t237 * t182) * t145;
t229 = (-t244 * t179 - t243 * t182) * t173 + (t238 * t171 + t233 * t172) * t146 + (t233 * t171 - t238 * t172) * t145;
t180 = sin(qJ(1));
t225 = pkin(1) * t180;
t183 = cos(qJ(1));
t224 = pkin(1) * t183;
t223 = pkin(7) * t179;
t222 = -pkin(5) - qJ(2);
t221 = Icges(2,4) * t180;
t220 = Icges(3,4) * t171;
t215 = t171 * t182;
t214 = t172 * t182;
t178 = sin(qJ(5));
t213 = t178 * t179;
t181 = cos(qJ(5));
t212 = t179 * t181;
t211 = qJD(4) * t179;
t210 = qJD(5) * t182;
t209 = t173 * t224 + V_base(2);
t208 = V_base(5) * pkin(5) + V_base(1);
t140 = pkin(2) * t171 - pkin(6) * t172;
t205 = -t140 - t225;
t204 = V_base(5) * qJ(2) + t208;
t203 = V_base(4) * t225 + qJD(2) + V_base(3);
t199 = pkin(3) * t182 + qJ(4) * t179;
t125 = t199 * t171;
t202 = -t125 + t205;
t201 = rSges(4,1) * t182 - rSges(4,2) * t179;
t200 = -rSges(5,2) * t182 + rSges(5,3) * t179;
t162 = pkin(3) * t179 - qJ(4) * t182;
t192 = t145 * t162 + t172 * t211 + t204;
t141 = pkin(2) * t172 + pkin(6) * t171;
t189 = t173 * t141 + t222 * V_base(4) + t209;
t126 = t199 * t172;
t188 = t173 * t126 + t171 * t211 + t189;
t187 = V_base(4) * t140 + (-t141 - t224) * V_base(5) + t203;
t186 = -qJD(4) * t182 + t146 * t125 + t187;
t175 = Icges(2,4) * t183;
t170 = Icges(3,4) * t172;
t166 = rSges(2,1) * t183 - t180 * rSges(2,2);
t165 = t180 * rSges(2,1) + rSges(2,2) * t183;
t164 = rSges(4,1) * t179 + rSges(4,2) * t182;
t163 = -rSges(5,2) * t179 - rSges(5,3) * t182;
t161 = qJD(5) * t179 + t173;
t158 = Icges(2,1) * t183 - t221;
t157 = Icges(2,1) * t180 + t175;
t155 = -Icges(2,2) * t180 + t175;
t154 = Icges(2,2) * t183 + t221;
t144 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t143 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t142 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t139 = rSges(3,1) * t172 - rSges(3,2) * t171;
t138 = rSges(3,1) * t171 + rSges(3,2) * t172;
t137 = Icges(3,1) * t172 - t220;
t136 = Icges(3,1) * t171 + t170;
t135 = -Icges(3,2) * t171 + t170;
t134 = Icges(3,2) * t172 + t220;
t131 = -pkin(4) * t172 + pkin(7) * t215;
t130 = pkin(4) * t171 + pkin(7) * t214;
t127 = rSges(6,3) * t179 + (-rSges(6,1) * t178 - rSges(6,2) * t181) * t182;
t124 = Icges(6,5) * t179 + (-Icges(6,1) * t178 - Icges(6,4) * t181) * t182;
t123 = Icges(6,6) * t179 + (-Icges(6,4) * t178 - Icges(6,2) * t181) * t182;
t122 = Icges(6,3) * t179 + (-Icges(6,5) * t178 - Icges(6,6) * t181) * t182;
t121 = t171 * t213 - t172 * t181;
t120 = t171 * t212 + t172 * t178;
t119 = t171 * t181 + t172 * t213;
t118 = -t171 * t178 + t172 * t212;
t117 = t172 * t210 + t146;
t116 = t171 * t210 + t145;
t113 = -rSges(5,1) * t172 + t200 * t171;
t112 = rSges(5,1) * t171 + t200 * t172;
t111 = rSges(4,3) * t171 + t201 * t172;
t110 = -rSges(4,3) * t172 + t201 * t171;
t109 = V_base(5) * rSges(2,3) - t165 * t173 + t208;
t108 = t166 * t173 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t95 = t165 * V_base(4) - t166 * V_base(5) + V_base(3);
t93 = V_base(5) * rSges(3,3) + (-t138 - t225) * t173 + t204;
t92 = t139 * t173 + (-rSges(3,3) + t222) * V_base(4) + t209;
t91 = V_base(4) * t138 + (-t139 - t224) * V_base(5) + t203;
t90 = rSges(6,1) * t121 + rSges(6,2) * t120 + rSges(6,3) * t215;
t89 = rSges(6,1) * t119 + rSges(6,2) * t118 + rSges(6,3) * t214;
t88 = Icges(6,1) * t121 + Icges(6,4) * t120 + Icges(6,5) * t215;
t87 = Icges(6,1) * t119 + Icges(6,4) * t118 + Icges(6,5) * t214;
t86 = Icges(6,4) * t121 + Icges(6,2) * t120 + Icges(6,6) * t215;
t85 = Icges(6,4) * t119 + Icges(6,2) * t118 + Icges(6,6) * t214;
t84 = Icges(6,5) * t121 + Icges(6,6) * t120 + Icges(6,3) * t215;
t83 = Icges(6,5) * t119 + Icges(6,6) * t118 + Icges(6,3) * t214;
t82 = t145 * t164 + (-t110 + t205) * t173 + t204;
t81 = t111 * t173 - t146 * t164 + t189;
t80 = t146 * t110 - t145 * t111 + t187;
t79 = t145 * t163 + (-t113 + t202) * t173 + t192;
t78 = t112 * t173 + (-t162 - t163) * t146 + t188;
t77 = t146 * t113 + (-t112 - t126) * t145 + t186;
t76 = t145 * t223 + t116 * t127 - t161 * t90 + (-t131 + t202) * t173 + t192;
t75 = -t117 * t127 + t130 * t173 + t161 * t89 + (-t162 - t223) * t146 + t188;
t74 = -t116 * t89 + t117 * t90 + t146 * t131 + (-t126 - t130) * t145 + t186;
t1 = m(1) * (t142 ^ 2 + t143 ^ 2 + t144 ^ 2) / 0.2e1 + m(2) * (t108 ^ 2 + t109 ^ 2 + t95 ^ 2) / 0.2e1 + m(3) * (t91 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + m(4) * (t80 ^ 2 + t81 ^ 2 + t82 ^ 2) / 0.2e1 + m(5) * (t77 ^ 2 + t78 ^ 2 + t79 ^ 2) / 0.2e1 + m(6) * (t74 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + t117 * ((t118 * t85 + t119 * t87 + t83 * t214) * t117 + (t118 * t86 + t119 * t88 + t84 * t214) * t116 + (t118 * t123 + t119 * t124 + t122 * t214) * t161) / 0.2e1 + t116 * ((t120 * t85 + t121 * t87 + t83 * t215) * t117 + (t120 * t86 + t121 * t88 + t84 * t215) * t116 + (t120 * t123 + t121 * t124 + t122 * t215) * t161) / 0.2e1 + t161 * ((t84 * t116 + t83 * t117 + t122 * t161) * t179 + ((-t178 * t87 - t181 * t85) * t117 + (-t178 * t88 - t181 * t86) * t116 + (-t123 * t181 - t124 * t178) * t161) * t182) / 0.2e1 + (t230 * t171 - t229 * t172) * t145 / 0.2e1 + (t229 * t171 + t230 * t172) * t146 / 0.2e1 + ((-t134 * t171 + t136 * t172 - t180 * t154 + t157 * t183 + Icges(1,4)) * V_base(5) + (-t171 * t135 + t172 * t137 - t180 * t155 + t183 * t158 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t172 * t134 + t171 * t136 + t183 * t154 + t180 * t157 + Icges(1,2)) * V_base(5) + (t135 * t172 + t137 * t171 + t155 * t183 + t180 * t158 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t236 * t179 - t240 * t182) * t146 + (t237 * t179 - t239 * t182) * t145 + (t234 * t179 - t235 * t182 + Icges(2,3) + Icges(3,3)) * t173) * t173 / 0.2e1 + t173 * V_base(5) * (Icges(2,5) * t180 + Icges(3,5) * t171 + Icges(2,6) * t183 + Icges(3,6) * t172) + t173 * V_base(4) * (Icges(2,5) * t183 + Icges(3,5) * t172 - Icges(2,6) * t180 - Icges(3,6) * t171) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
