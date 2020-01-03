% Calculate kinetic energy for
% S5PRPPR4
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
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
% Datum: 2019-12-31 17:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPPR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR4_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR4_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRPPR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR4_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR4_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR4_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPPR4_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:36:43
% EndTime: 2019-12-31 17:36:45
% DurationCPUTime: 2.38s
% Computational Cost: add. (973->238), mult. (1098->308), div. (0->0), fcn. (956->8), ass. (0->123)
t248 = Icges(4,4) - Icges(5,5);
t247 = Icges(4,1) + Icges(5,1);
t246 = Icges(4,2) + Icges(5,3);
t177 = sin(pkin(8));
t245 = t248 * t177;
t179 = cos(pkin(8));
t244 = t248 * t179;
t243 = Icges(5,4) + Icges(4,5);
t242 = Icges(4,6) - Icges(5,6);
t241 = t246 * t177 - t244;
t240 = t247 * t179 - t245;
t176 = pkin(7) + qJ(2);
t170 = sin(t176);
t171 = cos(t176);
t239 = -t241 * t170 - t242 * t171;
t238 = Icges(5,2) + Icges(4,3);
t237 = -t242 * t170 + t241 * t171;
t236 = t240 * t170 - t243 * t171;
t235 = t243 * t170 + t240 * t171;
t234 = -t246 * t179 - t245;
t233 = t247 * t177 + t244;
t232 = -t242 * t177 + t243 * t179;
t173 = V_base(6) + qJD(2);
t229 = (-t239 * t177 + t236 * t179) * V_base(5) + (t237 * t177 + t235 * t179) * V_base(4) + (t234 * t177 + t233 * t179) * t173;
t228 = (t232 * t170 - t238 * t171) * V_base(5) + (t238 * t170 + t232 * t171) * V_base(4) + (t243 * t177 + t242 * t179) * t173;
t180 = cos(pkin(7));
t226 = pkin(1) * t180;
t225 = pkin(4) * t177;
t224 = pkin(4) * t179;
t223 = -pkin(5) - qJ(1);
t178 = sin(pkin(7));
t222 = Icges(2,4) * t178;
t221 = Icges(3,4) * t170;
t200 = pkin(3) * t179 + qJ(4) * t177;
t123 = t200 * t170;
t135 = pkin(2) * t170 - qJ(3) * t171;
t216 = -t123 - t135;
t215 = qJD(4) * t177;
t207 = pkin(1) * V_base(6);
t214 = t180 * t207 + V_base(2);
t213 = V_base(5) * qJ(1) + V_base(1);
t209 = qJD(1) + V_base(3);
t159 = pkin(3) * t177 - qJ(4) * t179;
t208 = -t159 + t223;
t137 = pkin(2) * t171 + qJ(3) * t170;
t206 = -t137 - t226;
t205 = V_base(4) * t178 * pkin(1) + t209;
t124 = t200 * t171;
t204 = -t124 + t206;
t203 = V_base(4) * t135 + t205;
t202 = rSges(4,1) * t179 - rSges(4,2) * t177;
t201 = rSges(5,1) * t179 + rSges(5,3) * t177;
t181 = sin(qJ(5));
t182 = cos(qJ(5));
t141 = t177 * t182 - t179 * t181;
t193 = t177 * t181 + t179 * t182;
t192 = -qJD(3) * t171 + t173 * t137 + t214;
t191 = V_base(5) * pkin(5) - t178 * t207 + t213;
t190 = t173 * t124 + t170 * t215 + t192;
t189 = -qJD(4) * t179 + V_base(4) * t123 + t203;
t188 = qJD(3) * t170 + t191;
t185 = V_base(5) * t159 + t171 * t215 + t188;
t172 = Icges(2,4) * t180;
t169 = Icges(3,4) * t171;
t163 = rSges(2,1) * t180 - rSges(2,2) * t178;
t162 = rSges(2,1) * t178 + rSges(2,2) * t180;
t161 = rSges(4,1) * t177 + rSges(4,2) * t179;
t160 = rSges(5,1) * t177 - rSges(5,3) * t179;
t158 = -qJD(5) * t170 + V_base(4);
t157 = qJD(5) * t171 + V_base(5);
t156 = Icges(2,1) * t180 - t222;
t155 = Icges(2,1) * t178 + t172;
t152 = -Icges(2,2) * t178 + t172;
t151 = Icges(2,2) * t180 + t222;
t144 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t143 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t142 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t138 = rSges(3,1) * t171 - rSges(3,2) * t170;
t136 = rSges(3,1) * t170 + rSges(3,2) * t171;
t134 = Icges(3,1) * t171 - t221;
t133 = Icges(3,1) * t170 + t169;
t132 = -Icges(3,2) * t170 + t169;
t131 = Icges(3,2) * t171 + t221;
t130 = Icges(3,5) * t171 - Icges(3,6) * t170;
t129 = Icges(3,5) * t170 + Icges(3,6) * t171;
t128 = -pkin(6) * t170 + t171 * t224;
t127 = pkin(6) * t171 + t170 * t224;
t122 = t193 * t171;
t121 = t141 * t171;
t120 = t193 * t170;
t119 = t141 * t170;
t117 = V_base(5) * rSges(2,3) - t162 * V_base(6) + t213;
t116 = t163 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t114 = rSges(4,3) * t170 + t171 * t202;
t113 = rSges(5,2) * t170 + t171 * t201;
t112 = -rSges(4,3) * t171 + t170 * t202;
t111 = -rSges(5,2) * t171 + t170 * t201;
t98 = rSges(6,1) * t141 - rSges(6,2) * t193;
t97 = t162 * V_base(4) - t163 * V_base(5) + t209;
t96 = Icges(6,1) * t141 - Icges(6,4) * t193;
t95 = Icges(6,4) * t141 - Icges(6,2) * t193;
t94 = Icges(6,5) * t141 - Icges(6,6) * t193;
t93 = V_base(5) * rSges(3,3) - t136 * t173 + t191;
t92 = t138 * t173 + (-rSges(3,3) + t223) * V_base(4) + t214;
t91 = t136 * V_base(4) + (-t138 - t226) * V_base(5) + t205;
t90 = rSges(6,1) * t122 + rSges(6,2) * t121 - rSges(6,3) * t170;
t89 = rSges(6,1) * t120 + rSges(6,2) * t119 + rSges(6,3) * t171;
t88 = Icges(6,1) * t122 + Icges(6,4) * t121 - Icges(6,5) * t170;
t87 = Icges(6,1) * t120 + Icges(6,4) * t119 + Icges(6,5) * t171;
t86 = Icges(6,4) * t122 + Icges(6,2) * t121 - Icges(6,6) * t170;
t85 = Icges(6,4) * t120 + Icges(6,2) * t119 + Icges(6,6) * t171;
t84 = Icges(6,5) * t122 + Icges(6,6) * t121 - Icges(6,3) * t170;
t83 = Icges(6,5) * t120 + Icges(6,6) * t119 + Icges(6,3) * t171;
t82 = t161 * V_base(5) + (-t112 - t135) * t173 + t188;
t81 = t114 * t173 + (-t161 + t223) * V_base(4) + t192;
t80 = t112 * V_base(4) + (-t114 + t206) * V_base(5) + t203;
t79 = t160 * V_base(5) + (-t111 + t216) * t173 + t185;
t78 = t113 * t173 + (-t160 + t208) * V_base(4) + t190;
t77 = t111 * V_base(4) + (-t113 + t204) * V_base(5) + t189;
t76 = V_base(5) * t225 + t157 * t98 + (-t127 - t89 + t216) * t173 + t185;
t75 = -t158 * t98 + (t128 + t90) * t173 + (t208 - t225) * V_base(4) + t190;
t74 = t127 * V_base(4) - t157 * t90 + t158 * t89 + (-t128 + t204) * V_base(5) + t189;
t1 = m(1) * (t142 ^ 2 + t143 ^ 2 + t144 ^ 2) / 0.2e1 + m(2) * (t116 ^ 2 + t117 ^ 2 + t97 ^ 2) / 0.2e1 + m(3) * (t91 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + m(4) * (t80 ^ 2 + t81 ^ 2 + t82 ^ 2) / 0.2e1 + m(5) * (t77 ^ 2 + t78 ^ 2 + t79 ^ 2) / 0.2e1 + m(6) * (t74 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + t158 * ((t121 * t86 + t122 * t88 - t170 * t84) * t158 + (t121 * t85 + t122 * t87 - t170 * t83) * t157 + (t121 * t95 + t122 * t96 - t170 * t94) * t173) / 0.2e1 + t157 * ((t119 * t86 + t120 * t88 + t171 * t84) * t158 + (t119 * t85 + t120 * t87 + t171 * t83) * t157 + (t119 * t95 + t120 * t96 + t171 * t94) * t173) / 0.2e1 + ((t141 * t88 - t193 * t86) * t158 + (t141 * t87 - t193 * t85) * t157 + (t236 * t177 + t239 * t179 + t129) * V_base(5) + (t235 * t177 - t237 * t179 + t130) * V_base(4) + (t141 * t96 + t233 * t177 - t234 * t179 - t193 * t95 + Icges(3,3)) * t173) * t173 / 0.2e1 + (t130 * t173 + t229 * t171 + t228 * t170 + (-t131 * t170 + t133 * t171 - t151 * t178 + t155 * t180 + Icges(1,4)) * V_base(5) + (-t132 * t170 + t134 * t171 - t152 * t178 + t156 * t180 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t129 * t173 - t228 * t171 + t229 * t170 + (t131 * t171 + t133 * t170 + t151 * t180 + t155 * t178 + Icges(1,2)) * V_base(5) + (t132 * t171 + t134 * t170 + t152 * t180 + t156 * t178 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t178 + Icges(2,6) * t180 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t180 - Icges(2,6) * t178 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T = t1;
