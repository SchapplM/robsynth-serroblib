% Calculate kinetic energy for
% S5RPPPR3
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3]';
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
% Datum: 2019-12-31 17:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPPR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR3_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPPR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR3_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR3_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR3_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:43:48
% EndTime: 2019-12-31 17:43:51
% DurationCPUTime: 2.54s
% Computational Cost: add. (986->239), mult. (1098->305), div. (0->0), fcn. (956->8), ass. (0->123)
t247 = Icges(4,4) - Icges(5,5);
t246 = Icges(4,1) + Icges(5,1);
t245 = Icges(4,2) + Icges(5,3);
t177 = sin(pkin(8));
t244 = t247 * t177;
t178 = cos(pkin(8));
t243 = t247 * t178;
t242 = Icges(5,4) + Icges(4,5);
t241 = Icges(4,6) - Icges(5,6);
t240 = t245 * t177 - t243;
t239 = t246 * t178 - t244;
t176 = qJ(1) + pkin(7);
t170 = sin(t176);
t171 = cos(t176);
t238 = -t240 * t170 - t241 * t171;
t237 = Icges(5,2) + Icges(4,3);
t236 = -t241 * t170 + t240 * t171;
t235 = t239 * t170 - t242 * t171;
t234 = t242 * t170 + t239 * t171;
t233 = -t245 * t178 - t244;
t232 = t246 * t177 + t243;
t231 = -t241 * t177 + t242 * t178;
t180 = sin(qJ(1));
t182 = cos(qJ(1));
t230 = Icges(2,5) * t180 + Icges(3,5) * t170 + Icges(2,6) * t182 + Icges(3,6) * t171;
t229 = Icges(2,5) * t182 + Icges(3,5) * t171 - Icges(2,6) * t180 - Icges(3,6) * t170;
t172 = V_base(6) + qJD(1);
t228 = (-t238 * t177 + t235 * t178) * V_base(5) + (t236 * t177 + t234 * t178) * V_base(4) + (t233 * t177 + t232 * t178) * t172;
t227 = (t231 * t170 - t237 * t171) * V_base(5) + (t237 * t170 + t231 * t171) * V_base(4) + (t242 * t177 + t241 * t178) * t172;
t225 = pkin(1) * t180;
t224 = pkin(1) * t182;
t223 = pkin(4) * t177;
t222 = pkin(4) * t178;
t221 = -pkin(5) - qJ(2);
t220 = Icges(2,4) * t180;
t219 = Icges(3,4) * t170;
t214 = qJD(4) * t177;
t213 = t172 * t224 + V_base(2);
t212 = V_base(5) * pkin(5) + V_base(1);
t153 = pkin(3) * t177 - qJ(4) * t178;
t209 = -t153 + t221;
t135 = pkin(2) * t170 - qJ(3) * t171;
t208 = -t135 - t225;
t137 = pkin(2) * t171 + qJ(3) * t170;
t207 = -t137 - t224;
t206 = V_base(5) * qJ(2) + t212;
t205 = V_base(4) * t225 + qJD(2) + V_base(3);
t198 = pkin(3) * t178 + qJ(4) * t177;
t123 = t198 * t170;
t204 = -t123 + t208;
t124 = t198 * t171;
t203 = -t124 + t207;
t202 = qJD(3) * t170 + t206;
t201 = V_base(4) * t135 + t205;
t200 = rSges(4,1) * t178 - rSges(4,2) * t177;
t199 = rSges(5,1) * t178 + rSges(5,3) * t177;
t179 = sin(qJ(5));
t181 = cos(qJ(5));
t141 = t177 * t181 - t178 * t179;
t191 = t177 * t179 + t178 * t181;
t190 = -qJD(3) * t171 + t172 * t137 + t213;
t189 = V_base(5) * t153 + t171 * t214 + t202;
t188 = t172 * t124 + t170 * t214 + t190;
t187 = -qJD(4) * t178 + V_base(4) * t123 + t201;
t174 = Icges(2,4) * t182;
t169 = Icges(3,4) * t171;
t165 = rSges(2,1) * t182 - t180 * rSges(2,2);
t164 = t180 * rSges(2,1) + rSges(2,2) * t182;
t163 = Icges(2,1) * t182 - t220;
t162 = Icges(2,1) * t180 + t174;
t161 = -Icges(2,2) * t180 + t174;
t160 = Icges(2,2) * t182 + t220;
t155 = rSges(4,1) * t177 + rSges(4,2) * t178;
t154 = rSges(5,1) * t177 - rSges(5,3) * t178;
t152 = -qJD(5) * t170 + V_base(4);
t151 = qJD(5) * t171 + V_base(5);
t144 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t143 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t142 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t138 = rSges(3,1) * t171 - rSges(3,2) * t170;
t136 = rSges(3,1) * t170 + rSges(3,2) * t171;
t134 = Icges(3,1) * t171 - t219;
t133 = Icges(3,1) * t170 + t169;
t132 = -Icges(3,2) * t170 + t169;
t131 = Icges(3,2) * t171 + t219;
t128 = -pkin(6) * t170 + t171 * t222;
t127 = pkin(6) * t171 + t170 * t222;
t122 = t191 * t171;
t121 = t141 * t171;
t120 = t191 * t170;
t119 = t141 * t170;
t116 = V_base(5) * rSges(2,3) - t164 * t172 + t212;
t115 = t165 * t172 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t114 = rSges(4,3) * t170 + t171 * t200;
t113 = rSges(5,2) * t170 + t171 * t199;
t112 = -rSges(4,3) * t171 + t170 * t200;
t111 = -rSges(5,2) * t171 + t170 * t199;
t98 = t164 * V_base(4) - t165 * V_base(5) + V_base(3);
t97 = rSges(6,1) * t141 - rSges(6,2) * t191;
t96 = Icges(6,1) * t141 - Icges(6,4) * t191;
t95 = Icges(6,4) * t141 - Icges(6,2) * t191;
t94 = Icges(6,5) * t141 - Icges(6,6) * t191;
t93 = V_base(5) * rSges(3,3) + (-t136 - t225) * t172 + t206;
t92 = t138 * t172 + (-rSges(3,3) + t221) * V_base(4) + t213;
t91 = V_base(4) * t136 + (-t138 - t224) * V_base(5) + t205;
t90 = rSges(6,1) * t122 + rSges(6,2) * t121 - rSges(6,3) * t170;
t89 = rSges(6,1) * t120 + rSges(6,2) * t119 + rSges(6,3) * t171;
t88 = Icges(6,1) * t122 + Icges(6,4) * t121 - Icges(6,5) * t170;
t87 = Icges(6,1) * t120 + Icges(6,4) * t119 + Icges(6,5) * t171;
t86 = Icges(6,4) * t122 + Icges(6,2) * t121 - Icges(6,6) * t170;
t85 = Icges(6,4) * t120 + Icges(6,2) * t119 + Icges(6,6) * t171;
t84 = Icges(6,5) * t122 + Icges(6,6) * t121 - Icges(6,3) * t170;
t83 = Icges(6,5) * t120 + Icges(6,6) * t119 + Icges(6,3) * t171;
t82 = t155 * V_base(5) + (-t112 + t208) * t172 + t202;
t81 = t114 * t172 + (-t155 + t221) * V_base(4) + t190;
t80 = V_base(4) * t112 + (-t114 + t207) * V_base(5) + t201;
t79 = t154 * V_base(5) + (-t111 + t204) * t172 + t189;
t78 = t113 * t172 + (-t154 + t209) * V_base(4) + t188;
t77 = V_base(4) * t111 + (-t113 + t203) * V_base(5) + t187;
t76 = V_base(5) * t223 + t151 * t97 + (-t127 + t204 - t89) * t172 + t189;
t75 = -t152 * t97 + (t128 + t90) * t172 + (t209 - t223) * V_base(4) + t188;
t74 = V_base(4) * t127 - t151 * t90 + t152 * t89 + (-t128 + t203) * V_base(5) + t187;
t1 = m(1) * (t142 ^ 2 + t143 ^ 2 + t144 ^ 2) / 0.2e1 + m(2) * (t115 ^ 2 + t116 ^ 2 + t98 ^ 2) / 0.2e1 + m(3) * (t91 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + m(4) * (t80 ^ 2 + t81 ^ 2 + t82 ^ 2) / 0.2e1 + m(5) * (t77 ^ 2 + t78 ^ 2 + t79 ^ 2) / 0.2e1 + m(6) * (t74 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + t152 * ((t121 * t86 + t122 * t88 - t170 * t84) * t152 + (t121 * t85 + t122 * t87 - t170 * t83) * t151 + (t121 * t95 + t122 * t96 - t170 * t94) * t172) / 0.2e1 + t151 * ((t119 * t86 + t120 * t88 + t171 * t84) * t152 + (t119 * t85 + t120 * t87 + t171 * t83) * t151 + (t119 * t95 + t120 * t96 + t171 * t94) * t172) / 0.2e1 + ((t141 * t88 - t191 * t86) * t152 + (t141 * t87 - t191 * t85) * t151 + (t235 * t177 + t238 * t178 + t230) * V_base(5) + (t234 * t177 - t236 * t178 + t229) * V_base(4) + (t141 * t96 + t232 * t177 - t233 * t178 - t191 * t95 + Icges(2,3) + Icges(3,3)) * t172) * t172 / 0.2e1 + (t229 * t172 + t228 * t171 + t227 * t170 + (-t131 * t170 + t133 * t171 - t180 * t160 + t162 * t182 + Icges(1,4)) * V_base(5) + (-t132 * t170 + t134 * t171 - t180 * t161 + t163 * t182 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t230 * t172 - t227 * t171 + t228 * t170 + (t131 * t171 + t133 * t170 + t160 * t182 + t180 * t162 + Icges(1,2)) * V_base(5) + (t132 * t171 + t134 * t170 + t161 * t182 + t180 * t163 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
