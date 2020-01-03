% Calculate kinetic energy for
% S5RPRPR15
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR15_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR15_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRPR15_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR15_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR15_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR15_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:36:34
% EndTime: 2019-12-31 18:36:36
% DurationCPUTime: 2.59s
% Computational Cost: add. (886->282), mult. (1388->399), div. (0->0), fcn. (1282->8), ass. (0->133)
t246 = Icges(2,4) + Icges(3,6);
t245 = Icges(2,1) + Icges(3,2);
t244 = -Icges(3,4) + Icges(2,5);
t243 = Icges(3,5) - Icges(2,6);
t242 = Icges(2,2) + Icges(3,3);
t186 = cos(qJ(1));
t241 = t246 * t186;
t184 = sin(qJ(1));
t240 = t246 * t184;
t239 = -t242 * t186 - t240;
t238 = t242 * t184 - t241;
t237 = t245 * t184 + t241;
t236 = t245 * t186 - t240;
t183 = sin(qJ(3));
t185 = cos(qJ(3));
t181 = cos(pkin(8));
t226 = pkin(4) * t181;
t233 = -pkin(7) * t185 + t183 * t226;
t224 = Icges(4,4) * t183;
t198 = Icges(4,2) * t185 + t224;
t119 = Icges(4,6) * t186 + t184 * t198;
t120 = Icges(4,6) * t184 - t186 * t198;
t223 = Icges(4,4) * t185;
t199 = Icges(4,1) * t183 + t223;
t121 = Icges(4,5) * t186 + t184 * t199;
t122 = Icges(4,5) * t184 - t186 * t199;
t146 = -Icges(4,2) * t183 + t223;
t151 = Icges(4,1) * t185 - t224;
t164 = qJD(3) * t184 + V_base(5);
t165 = qJD(3) * t186 + V_base(4);
t171 = V_base(6) + qJD(1);
t232 = (t119 * t185 + t121 * t183) * t165 + (t120 * t185 + t122 * t183) * t164 + (t146 * t185 + t151 * t183) * t171;
t228 = pkin(6) * t184;
t227 = pkin(6) * t186;
t180 = sin(pkin(8));
t220 = t180 * t186;
t219 = t183 * t186;
t179 = pkin(8) + qJ(5);
t169 = sin(t179);
t218 = t184 * t169;
t170 = cos(t179);
t217 = t184 * t170;
t216 = t184 * t180;
t215 = t184 * t181;
t214 = t184 * t185;
t213 = t185 * t186;
t211 = qJD(4) * t185;
t210 = qJD(5) * t185;
t155 = t184 * pkin(1) - qJ(2) * t186;
t209 = V_base(4) * t155 + V_base(3);
t208 = V_base(5) * pkin(5) + V_base(1);
t205 = -t155 - t228;
t204 = qJD(2) * t184 + t208;
t200 = pkin(3) * t183 - qJ(4) * t185;
t134 = t200 * t186;
t203 = t134 + t205;
t202 = V_base(5) * pkin(2) + t204;
t201 = rSges(4,1) * t183 + rSges(4,2) * t185;
t197 = Icges(4,5) * t183 + Icges(4,6) * t185;
t160 = pkin(1) * t186 + t184 * qJ(2);
t193 = -qJD(2) * t186 + t171 * t160 + V_base(2);
t192 = (Icges(4,3) * t186 + t184 * t197) * t165 + (Icges(4,3) * t184 - t186 * t197) * t164 + (Icges(4,5) * t185 - Icges(4,6) * t183) * t171;
t158 = pkin(3) * t185 + qJ(4) * t183;
t191 = t164 * t158 - t184 * t211 + t202;
t190 = V_base(4) * t228 + (-t160 - t227) * V_base(5) + t209;
t189 = t171 * t227 + (-pkin(2) - pkin(5)) * V_base(4) + t193;
t188 = qJD(4) * t183 - t165 * t134 + t190;
t133 = t200 * t184;
t187 = t171 * t133 + t186 * t211 + t189;
t162 = rSges(2,1) * t186 - t184 * rSges(2,2);
t161 = -rSges(3,2) * t186 + t184 * rSges(3,3);
t159 = rSges(4,1) * t185 - rSges(4,2) * t183;
t157 = t184 * rSges(2,1) + rSges(2,2) * t186;
t156 = -t184 * rSges(3,2) - rSges(3,3) * t186;
t154 = qJD(5) * t183 + t171;
t138 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t137 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t136 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t131 = -t184 * t210 + t165;
t130 = t186 * t210 + t164;
t129 = -t181 * t219 + t216;
t128 = t180 * t219 + t215;
t127 = t183 * t215 + t220;
t126 = t181 * t186 - t183 * t216;
t124 = t184 * rSges(4,3) - t186 * t201;
t123 = rSges(4,3) * t186 + t184 * t201;
t116 = -t170 * t219 + t218;
t115 = t169 * t219 + t217;
t114 = t169 * t186 + t183 * t217;
t113 = t170 * t186 - t183 * t218;
t112 = rSges(5,3) * t183 + (rSges(5,1) * t181 - rSges(5,2) * t180) * t185;
t110 = Icges(5,5) * t183 + (Icges(5,1) * t181 - Icges(5,4) * t180) * t185;
t109 = Icges(5,6) * t183 + (Icges(5,4) * t181 - Icges(5,2) * t180) * t185;
t108 = Icges(5,3) * t183 + (Icges(5,5) * t181 - Icges(5,6) * t180) * t185;
t106 = rSges(6,3) * t183 + (rSges(6,1) * t170 - rSges(6,2) * t169) * t185;
t105 = Icges(6,5) * t183 + (Icges(6,1) * t170 - Icges(6,4) * t169) * t185;
t104 = Icges(6,6) * t183 + (Icges(6,4) * t170 - Icges(6,2) * t169) * t185;
t103 = Icges(6,3) * t183 + (Icges(6,5) * t170 - Icges(6,6) * t169) * t185;
t102 = pkin(7) * t183 + t185 * t226;
t101 = V_base(5) * rSges(2,3) - t157 * t171 + t208;
t100 = t162 * t171 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t99 = t157 * V_base(4) - t162 * V_base(5) + V_base(3);
t98 = pkin(4) * t216 - t186 * t233;
t97 = pkin(4) * t220 + t184 * t233;
t96 = t129 * rSges(5,1) + t128 * rSges(5,2) + rSges(5,3) * t213;
t95 = rSges(5,1) * t127 + rSges(5,2) * t126 - rSges(5,3) * t214;
t94 = Icges(5,1) * t129 + Icges(5,4) * t128 + Icges(5,5) * t213;
t93 = Icges(5,1) * t127 + Icges(5,4) * t126 - Icges(5,5) * t214;
t92 = Icges(5,4) * t129 + Icges(5,2) * t128 + Icges(5,6) * t213;
t91 = Icges(5,4) * t127 + Icges(5,2) * t126 - Icges(5,6) * t214;
t90 = Icges(5,5) * t129 + Icges(5,6) * t128 + Icges(5,3) * t213;
t89 = Icges(5,5) * t127 + Icges(5,6) * t126 - Icges(5,3) * t214;
t88 = V_base(5) * rSges(3,1) + (-t155 - t156) * t171 + t204;
t87 = t171 * t161 + (-rSges(3,1) - pkin(5)) * V_base(4) + t193;
t86 = t116 * rSges(6,1) + t115 * rSges(6,2) + rSges(6,3) * t213;
t85 = rSges(6,1) * t114 + rSges(6,2) * t113 - rSges(6,3) * t214;
t84 = Icges(6,1) * t116 + Icges(6,4) * t115 + Icges(6,5) * t213;
t83 = Icges(6,1) * t114 + Icges(6,4) * t113 - Icges(6,5) * t214;
t82 = Icges(6,4) * t116 + Icges(6,2) * t115 + Icges(6,6) * t213;
t81 = Icges(6,4) * t114 + Icges(6,2) * t113 - Icges(6,6) * t214;
t80 = Icges(6,5) * t116 + Icges(6,6) * t115 + Icges(6,3) * t213;
t79 = Icges(6,5) * t114 + Icges(6,6) * t113 - Icges(6,3) * t214;
t78 = t156 * V_base(4) + (-t160 - t161) * V_base(5) + t209;
t77 = t159 * t164 + (-t124 + t205) * t171 + t202;
t76 = t171 * t123 - t165 * t159 + t189;
t75 = -t164 * t123 + t165 * t124 + t190;
t74 = t112 * t164 + (t203 - t96) * t171 + t191;
t73 = t171 * t95 + (-t112 - t158) * t165 + t187;
t72 = t165 * t96 + (-t133 - t95) * t164 + t188;
t71 = t102 * t164 + t106 * t130 - t154 * t86 + (t203 - t98) * t171 + t191;
t70 = -t131 * t106 + t154 * t85 + t171 * t97 + (-t102 - t158) * t165 + t187;
t69 = -t130 * t85 + t131 * t86 + t165 * t98 + (-t133 - t97) * t164 + t188;
t1 = m(1) * (t136 ^ 2 + t137 ^ 2 + t138 ^ 2) / 0.2e1 + m(2) * (t100 ^ 2 + t101 ^ 2 + t99 ^ 2) / 0.2e1 + m(3) * (t78 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(4) * (t75 ^ 2 + t76 ^ 2 + t77 ^ 2) / 0.2e1 + m(5) * (t72 ^ 2 + t73 ^ 2 + t74 ^ 2) / 0.2e1 + m(6) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + t131 * ((t113 * t81 + t114 * t83 - t79 * t214) * t131 + (t113 * t82 + t114 * t84 - t214 * t80) * t130 + (-t103 * t214 + t104 * t113 + t105 * t114) * t154) / 0.2e1 + t130 * ((t115 * t81 + t116 * t83 + t213 * t79) * t131 + (t115 * t82 + t116 * t84 + t80 * t213) * t130 + (t103 * t213 + t115 * t104 + t116 * t105) * t154) / 0.2e1 + t154 * ((t103 * t154 + t80 * t130 + t79 * t131) * t183 + ((-t169 * t81 + t170 * t83) * t131 + (-t169 * t82 + t170 * t84) * t130 + (-t104 * t169 + t105 * t170) * t154) * t185) / 0.2e1 + (t192 * t184 - t232 * t186 + (t128 * t91 + t129 * t93 + t213 * t89) * t165 + (t128 * t92 + t129 * t94 + t213 * t90) * t164 + (t108 * t213 + t128 * t109 + t129 * t110) * t171) * t164 / 0.2e1 + (t192 * t186 + t232 * t184 + (t126 * t91 + t127 * t93 - t214 * t89) * t165 + (t126 * t92 + t127 * t94 - t214 * t90) * t164 + (-t108 * t214 + t109 * t126 + t110 * t127) * t171) * t165 / 0.2e1 + ((t184 * t239 + t237 * t186 + Icges(1,4)) * V_base(5) + (t238 * t184 + t236 * t186 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t237 * t184 - t239 * t186 + Icges(1,2)) * V_base(5) + (t184 * t236 - t186 * t238 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t119 * t183 + t121 * t185) * t165 + (-t120 * t183 + t122 * t185) * t164 + (t90 * t164 + t89 * t165) * t183 + ((-t180 * t91 + t181 * t93) * t165 + (-t180 * t92 + t181 * t94) * t164) * t185 + (Icges(2,3) + Icges(3,1) + (-t109 * t180 + t110 * t181 + t151) * t185 + (-t146 + t108) * t183) * t171) * t171 / 0.2e1 + t171 * V_base(5) * (t244 * t184 - t243 * t186) + t171 * V_base(4) * (t243 * t184 + t244 * t186) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
