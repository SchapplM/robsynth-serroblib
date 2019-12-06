% Calculate kinetic energy for
% S5RRRPR1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-05 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR1_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRPR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR1_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR1_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR1_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:37:44
% EndTime: 2019-12-05 18:37:46
% DurationCPUTime: 1.79s
% Computational Cost: add. (1335->259), mult. (1214->367), div. (0->0), fcn. (996->10), ass. (0->139)
t254 = Icges(4,3) + Icges(5,3);
t186 = qJ(2) + qJ(3);
t174 = pkin(9) + t186;
t169 = sin(t174);
t170 = cos(t174);
t179 = sin(t186);
t180 = cos(t186);
t253 = Icges(4,5) * t180 + Icges(5,5) * t170 - Icges(4,6) * t179 - Icges(5,6) * t169;
t188 = sin(qJ(1));
t190 = cos(qJ(1));
t236 = Icges(5,4) * t170;
t212 = -Icges(5,2) * t169 + t236;
t100 = -Icges(5,6) * t190 + t188 * t212;
t101 = Icges(5,6) * t188 + t190 * t212;
t237 = Icges(5,4) * t169;
t216 = Icges(5,1) * t170 - t237;
t102 = -Icges(5,5) * t190 + t188 * t216;
t103 = Icges(5,5) * t188 + t190 * t216;
t238 = Icges(4,4) * t180;
t213 = -Icges(4,2) * t179 + t238;
t112 = -Icges(4,6) * t190 + t188 * t213;
t113 = Icges(4,6) * t188 + t190 * t213;
t239 = Icges(4,4) * t179;
t217 = Icges(4,1) * t180 - t239;
t114 = -Icges(4,5) * t190 + t188 * t217;
t115 = Icges(4,5) * t188 + t190 * t217;
t135 = Icges(5,2) * t170 + t237;
t136 = Icges(5,1) * t169 + t236;
t140 = Icges(4,2) * t180 + t239;
t141 = Icges(4,1) * t179 + t238;
t228 = -qJD(2) - qJD(3);
t143 = t190 * t228 + V_base(5);
t166 = qJD(2) * t188 + V_base(4);
t144 = qJD(3) * t188 + t166;
t175 = V_base(6) + qJD(1);
t252 = (-t135 * t169 + t136 * t170 - t140 * t179 + t141 * t180) * t175 + (-t101 * t169 + t103 * t170 - t113 * t179 + t115 * t180) * t144 + (-t100 * t169 + t102 * t170 - t112 * t179 + t114 * t180) * t143;
t251 = (Icges(4,5) * t179 + Icges(5,5) * t169 + Icges(4,6) * t180 + Icges(5,6) * t170) * t175 + (t254 * t188 + t253 * t190) * t144 + (t253 * t188 - t254 * t190) * t143;
t187 = sin(qJ(2));
t247 = pkin(2) * t187;
t246 = pkin(3) * t179;
t245 = pkin(4) * t169;
t189 = cos(qJ(2));
t244 = t189 * pkin(2);
t242 = Icges(2,4) * t188;
t241 = Icges(3,4) * t187;
t240 = Icges(3,4) * t189;
t172 = qJ(5) + t174;
t167 = sin(t172);
t235 = Icges(6,4) * t167;
t168 = cos(t172);
t234 = Icges(6,4) * t168;
t108 = -pkin(7) * t190 + t188 * t244;
t163 = t188 * pkin(1) - t190 * pkin(6);
t233 = -t108 - t163;
t232 = pkin(4) * t170;
t231 = pkin(3) * t180;
t227 = V_base(5) * pkin(5) + V_base(1);
t85 = -qJ(4) * t190 + t188 * t231;
t224 = -t85 + t233;
t165 = -qJD(2) * t190 + V_base(5);
t223 = t165 * t247 + t227;
t222 = rSges(3,1) * t189 - rSges(3,2) * t187;
t221 = rSges(4,1) * t180 - rSges(4,2) * t179;
t220 = rSges(5,1) * t170 - rSges(5,2) * t169;
t219 = rSges(6,1) * t168 - rSges(6,2) * t167;
t218 = Icges(3,1) * t189 - t241;
t215 = Icges(6,1) * t168 - t235;
t214 = -Icges(3,2) * t187 + t240;
t211 = -Icges(6,2) * t167 + t234;
t210 = Icges(3,5) * t189 - Icges(3,6) * t187;
t207 = Icges(6,5) * t168 - Icges(6,6) * t167;
t206 = qJD(4) * t188 + t143 * t246 + t223;
t164 = t190 * pkin(1) + t188 * pkin(6);
t205 = -V_base(4) * pkin(5) + t175 * t164 + V_base(2);
t204 = V_base(4) * t163 - t164 * V_base(5) + V_base(3);
t132 = V_base(5) + (-qJD(5) + t228) * t190;
t133 = qJD(5) * t188 + t144;
t203 = (Icges(6,5) * t167 + Icges(6,6) * t168) * t175 + (-Icges(6,3) * t190 + t188 * t207) * t132 + (Icges(6,3) * t188 + t190 * t207) * t133;
t200 = (-Icges(3,3) * t190 + t188 * t210) * t165 + (Icges(3,3) * t188 + t190 * t210) * t166 + (Icges(3,5) * t187 + Icges(3,6) * t189) * t175;
t109 = pkin(7) * t188 + t190 * t244;
t199 = t166 * t108 - t109 * t165 + t204;
t198 = t175 * t109 - t166 * t247 + t205;
t197 = t144 * t85 + t199;
t86 = qJ(4) * t188 + t190 * t231;
t196 = -qJD(4) * t190 + t175 * t86 + t198;
t126 = Icges(6,2) * t168 + t235;
t127 = Icges(6,1) * t167 + t234;
t91 = -Icges(6,6) * t190 + t188 * t211;
t92 = Icges(6,6) * t188 + t190 * t211;
t93 = -Icges(6,5) * t190 + t188 * t215;
t94 = Icges(6,5) * t188 + t190 * t215;
t195 = (-t167 * t92 + t168 * t94) * t133 + (-t167 * t91 + t168 * t93) * t132 + (-t126 * t167 + t127 * t168) * t175;
t121 = -Icges(3,6) * t190 + t188 * t214;
t122 = Icges(3,6) * t188 + t190 * t214;
t123 = -Icges(3,5) * t190 + t188 * t218;
t124 = Icges(3,5) * t188 + t190 * t218;
t154 = Icges(3,2) * t189 + t241;
t157 = Icges(3,1) * t187 + t240;
t192 = (-t122 * t187 + t124 * t189) * t166 + (-t121 * t187 + t123 * t189) * t165 + (-t154 * t187 + t157 * t189) * t175;
t182 = Icges(2,4) * t190;
t162 = rSges(2,1) * t190 - rSges(2,2) * t188;
t161 = rSges(2,1) * t188 + rSges(2,2) * t190;
t160 = rSges(3,1) * t187 + rSges(3,2) * t189;
t159 = Icges(2,1) * t190 - t242;
t158 = Icges(2,1) * t188 + t182;
t156 = -Icges(2,2) * t188 + t182;
t155 = Icges(2,2) * t190 + t242;
t150 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t149 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t148 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t142 = rSges(4,1) * t179 + rSges(4,2) * t180;
t137 = rSges(5,1) * t169 + rSges(5,2) * t170;
t130 = rSges(6,1) * t167 + rSges(6,2) * t168;
t129 = rSges(3,3) * t188 + t190 * t222;
t128 = -rSges(3,3) * t190 + t188 * t222;
t117 = rSges(4,3) * t188 + t190 * t221;
t116 = -rSges(4,3) * t190 + t188 * t221;
t107 = rSges(5,3) * t188 + t190 * t220;
t106 = -rSges(5,3) * t190 + t188 * t220;
t105 = V_base(5) * rSges(2,3) - t161 * t175 + t227;
t104 = t162 * t175 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t97 = t161 * V_base(4) - t162 * V_base(5) + V_base(3);
t96 = rSges(6,3) * t188 + t190 * t219;
t95 = -rSges(6,3) * t190 + t188 * t219;
t82 = pkin(8) * t188 + t190 * t232;
t81 = -pkin(8) * t190 + t188 * t232;
t80 = t160 * t165 + (-t128 - t163) * t175 + t227;
t79 = t129 * t175 - t160 * t166 + t205;
t78 = t128 * t166 - t129 * t165 + t204;
t77 = t142 * t143 + (-t116 + t233) * t175 + t223;
t76 = t117 * t175 - t142 * t144 + t198;
t75 = t116 * t144 - t117 * t143 + t199;
t74 = t137 * t143 + (-t106 + t224) * t175 + t206;
t73 = t107 * t175 + (-t137 - t246) * t144 + t196;
t72 = t106 * t144 + (-t107 - t86) * t143 + t197;
t71 = t143 * t245 + t130 * t132 + (-t81 - t95 + t224) * t175 + t206;
t70 = -t130 * t133 + (t82 + t96) * t175 + (-t245 - t246) * t144 + t196;
t69 = -t132 * t96 + t133 * t95 + t144 * t81 + (-t82 - t86) * t143 + t197;
t1 = m(1) * (t148 ^ 2 + t149 ^ 2 + t150 ^ 2) / 0.2e1 + m(2) * (t104 ^ 2 + t105 ^ 2 + t97 ^ 2) / 0.2e1 + m(3) * (t78 ^ 2 + t79 ^ 2 + t80 ^ 2) / 0.2e1 + t166 * (t188 * t200 + t190 * t192) / 0.2e1 + t165 * (t188 * t192 - t190 * t200) / 0.2e1 + m(4) * (t75 ^ 2 + t76 ^ 2 + t77 ^ 2) / 0.2e1 + m(5) * (t72 ^ 2 + t73 ^ 2 + t74 ^ 2) / 0.2e1 + m(6) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + t133 * (t188 * t203 + t190 * t195) / 0.2e1 + t132 * (t188 * t195 - t190 * t203) / 0.2e1 + (t252 * t188 - t251 * t190) * t143 / 0.2e1 + (t251 * t188 + t252 * t190) * t144 / 0.2e1 + ((-t155 * t188 + t158 * t190 + Icges(1,4)) * V_base(5) + (-t156 * t188 + t159 * t190 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t155 * t190 + t158 * t188 + Icges(1,2)) * V_base(5) + (t156 * t190 + t159 * t188 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t122 * t189 + t124 * t187) * t166 + (t121 * t189 + t123 * t187) * t165 + (t167 * t94 + t168 * t92) * t133 + (t167 * t93 + t168 * t91) * t132 + (t101 * t170 + t103 * t169 + t113 * t180 + t115 * t179) * t144 + (t100 * t170 + t102 * t169 + t112 * t180 + t114 * t179) * t143 + (t126 * t168 + t127 * t167 + t135 * t170 + t136 * t169 + t140 * t180 + t141 * t179 + t154 * t189 + t157 * t187 + Icges(2,3)) * t175) * t175 / 0.2e1 + V_base(4) * t175 * (Icges(2,5) * t190 - Icges(2,6) * t188) + V_base(5) * t175 * (Icges(2,5) * t188 + Icges(2,6) * t190) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
