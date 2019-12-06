% Calculate kinetic energy for
% S5RRPRR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR2_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPRR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR2_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR2_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR2_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:27:00
% EndTime: 2019-12-05 18:27:02
% DurationCPUTime: 1.92s
% Computational Cost: add. (1317->259), mult. (1196->367), div. (0->0), fcn. (978->10), ass. (0->139)
t254 = Icges(3,3) + Icges(4,3);
t186 = qJ(2) + pkin(9);
t174 = sin(t186);
t175 = cos(t186);
t188 = sin(qJ(2));
t190 = cos(qJ(2));
t253 = Icges(3,5) * t190 + Icges(4,5) * t175 - Icges(3,6) * t188 - Icges(4,6) * t174;
t189 = sin(qJ(1));
t191 = cos(qJ(1));
t238 = Icges(4,4) * t175;
t213 = -Icges(4,2) * t174 + t238;
t112 = -Icges(4,6) * t191 + t213 * t189;
t113 = Icges(4,6) * t189 + t213 * t191;
t239 = Icges(4,4) * t174;
t217 = Icges(4,1) * t175 - t239;
t114 = -Icges(4,5) * t191 + t217 * t189;
t115 = Icges(4,5) * t189 + t217 * t191;
t240 = Icges(3,4) * t190;
t214 = -Icges(3,2) * t188 + t240;
t121 = -Icges(3,6) * t191 + t214 * t189;
t122 = Icges(3,6) * t189 + t214 * t191;
t241 = Icges(3,4) * t188;
t218 = Icges(3,1) * t190 - t241;
t123 = -Icges(3,5) * t191 + t218 * t189;
t124 = Icges(3,5) * t189 + t218 * t191;
t140 = Icges(4,2) * t175 + t239;
t141 = Icges(4,1) * t174 + t238;
t154 = Icges(3,2) * t190 + t241;
t157 = Icges(3,1) * t188 + t240;
t165 = -qJD(2) * t191 + V_base(5);
t166 = qJD(2) * t189 + V_base(4);
t177 = V_base(6) + qJD(1);
t252 = (-t140 * t174 + t141 * t175 - t154 * t188 + t157 * t190) * t177 + (-t113 * t174 + t115 * t175 - t122 * t188 + t124 * t190) * t166 + (-t112 * t174 + t114 * t175 - t121 * t188 + t123 * t190) * t165;
t251 = (Icges(3,5) * t188 + Icges(4,5) * t174 + Icges(3,6) * t190 + Icges(4,6) * t175) * t177 + (t254 * t189 + t253 * t191) * t166 + (t253 * t189 - t254 * t191) * t165;
t247 = pkin(2) * t188;
t246 = pkin(3) * t174;
t176 = qJ(4) + t186;
t170 = sin(t176);
t245 = pkin(4) * t170;
t244 = t190 * pkin(2);
t242 = Icges(2,4) * t189;
t237 = Icges(5,4) * t170;
t171 = cos(t176);
t236 = Icges(5,4) * t171;
t172 = qJ(5) + t176;
t167 = sin(t172);
t235 = Icges(6,4) * t167;
t168 = cos(t172);
t234 = Icges(6,4) * t168;
t108 = -qJ(3) * t191 + t244 * t189;
t163 = t189 * pkin(1) - pkin(6) * t191;
t233 = -t108 - t163;
t232 = pkin(4) * t171;
t231 = pkin(3) * t175;
t228 = -qJD(2) - qJD(4);
t227 = V_base(5) * pkin(5) + V_base(1);
t85 = -pkin(7) * t191 + t231 * t189;
t224 = -t85 + t233;
t144 = qJD(4) * t189 + t166;
t223 = qJD(3) * t189 + t165 * t247 + t227;
t222 = rSges(3,1) * t190 - rSges(3,2) * t188;
t221 = rSges(4,1) * t175 - rSges(4,2) * t174;
t220 = rSges(5,1) * t171 - rSges(5,2) * t170;
t219 = rSges(6,1) * t168 - rSges(6,2) * t167;
t216 = Icges(5,1) * t171 - t237;
t215 = Icges(6,1) * t168 - t235;
t212 = -Icges(5,2) * t170 + t236;
t211 = -Icges(6,2) * t167 + t234;
t208 = Icges(5,5) * t171 - Icges(5,6) * t170;
t207 = Icges(6,5) * t168 - Icges(6,6) * t167;
t206 = t165 * t246 + t223;
t164 = pkin(1) * t191 + t189 * pkin(6);
t205 = -V_base(4) * pkin(5) + t177 * t164 + V_base(2);
t204 = V_base(4) * t163 - t164 * V_base(5) + V_base(3);
t203 = t166 * t108 + t204;
t131 = V_base(5) + (-qJD(5) + t228) * t191;
t132 = qJD(5) * t189 + t144;
t202 = (Icges(6,5) * t167 + Icges(6,6) * t168) * t177 + (-Icges(6,3) * t191 + t207 * t189) * t131 + (Icges(6,3) * t189 + t207 * t191) * t132;
t143 = t228 * t191 + V_base(5);
t201 = (Icges(5,5) * t170 + Icges(5,6) * t171) * t177 + (-Icges(5,3) * t191 + t208 * t189) * t143 + (Icges(5,3) * t189 + t208 * t191) * t144;
t109 = qJ(3) * t189 + t244 * t191;
t198 = -qJD(3) * t191 + t177 * t109 + t205;
t86 = pkin(7) * t189 + t231 * t191;
t197 = t166 * t85 + (-t109 - t86) * t165 + t203;
t196 = t177 * t86 + (-t246 - t247) * t166 + t198;
t126 = Icges(6,2) * t168 + t235;
t127 = Icges(6,1) * t167 + t234;
t91 = -Icges(6,6) * t191 + t211 * t189;
t92 = Icges(6,6) * t189 + t211 * t191;
t93 = -Icges(6,5) * t191 + t215 * t189;
t94 = Icges(6,5) * t189 + t215 * t191;
t195 = (-t167 * t92 + t168 * t94) * t132 + (-t167 * t91 + t168 * t93) * t131 + (-t126 * t167 + t127 * t168) * t177;
t100 = -Icges(5,6) * t191 + t212 * t189;
t101 = Icges(5,6) * t189 + t212 * t191;
t102 = -Icges(5,5) * t191 + t216 * t189;
t103 = Icges(5,5) * t189 + t216 * t191;
t134 = Icges(5,2) * t171 + t237;
t135 = Icges(5,1) * t170 + t236;
t194 = (-t101 * t170 + t103 * t171) * t144 + (-t100 * t170 + t102 * t171) * t143 + (-t134 * t170 + t135 * t171) * t177;
t182 = Icges(2,4) * t191;
t162 = rSges(2,1) * t191 - t189 * rSges(2,2);
t161 = t189 * rSges(2,1) + rSges(2,2) * t191;
t160 = rSges(3,1) * t188 + rSges(3,2) * t190;
t159 = Icges(2,1) * t191 - t242;
t158 = Icges(2,1) * t189 + t182;
t156 = -Icges(2,2) * t189 + t182;
t155 = Icges(2,2) * t191 + t242;
t150 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t149 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t148 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t142 = rSges(4,1) * t174 + rSges(4,2) * t175;
t136 = rSges(5,1) * t170 + rSges(5,2) * t171;
t130 = rSges(6,1) * t167 + rSges(6,2) * t168;
t129 = t189 * rSges(3,3) + t222 * t191;
t128 = -rSges(3,3) * t191 + t222 * t189;
t117 = t189 * rSges(4,3) + t221 * t191;
t116 = -rSges(4,3) * t191 + t221 * t189;
t107 = t189 * rSges(5,3) + t220 * t191;
t106 = -rSges(5,3) * t191 + t220 * t189;
t105 = V_base(5) * rSges(2,3) - t161 * t177 + t227;
t104 = t162 * t177 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t97 = t161 * V_base(4) - t162 * V_base(5) + V_base(3);
t96 = t189 * rSges(6,3) + t219 * t191;
t95 = -rSges(6,3) * t191 + t219 * t189;
t82 = pkin(8) * t189 + t232 * t191;
t81 = -pkin(8) * t191 + t232 * t189;
t80 = t160 * t165 + (-t128 - t163) * t177 + t227;
t79 = t129 * t177 - t160 * t166 + t205;
t78 = t128 * t166 - t129 * t165 + t204;
t77 = t142 * t165 + (-t116 + t233) * t177 + t223;
t76 = t177 * t117 + (-t142 - t247) * t166 + t198;
t75 = t116 * t166 + (-t109 - t117) * t165 + t203;
t74 = t136 * t143 + (-t106 + t224) * t177 + t206;
t73 = t177 * t107 - t144 * t136 + t196;
t72 = t106 * t144 - t107 * t143 + t197;
t71 = t143 * t245 + t130 * t131 + (-t81 - t95 + t224) * t177 + t206;
t70 = -t144 * t245 - t132 * t130 + (t82 + t96) * t177 + t196;
t69 = -t131 * t96 + t132 * t95 - t143 * t82 + t144 * t81 + t197;
t1 = m(1) * (t148 ^ 2 + t149 ^ 2 + t150 ^ 2) / 0.2e1 + m(2) * (t104 ^ 2 + t105 ^ 2 + t97 ^ 2) / 0.2e1 + m(3) * (t78 ^ 2 + t79 ^ 2 + t80 ^ 2) / 0.2e1 + m(4) * (t75 ^ 2 + t76 ^ 2 + t77 ^ 2) / 0.2e1 + m(5) * (t72 ^ 2 + t73 ^ 2 + t74 ^ 2) / 0.2e1 + t144 * (t201 * t189 + t194 * t191) / 0.2e1 + t143 * (t194 * t189 - t201 * t191) / 0.2e1 + m(6) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + t132 * (t202 * t189 + t195 * t191) / 0.2e1 + t131 * (t195 * t189 - t202 * t191) / 0.2e1 + (t252 * t189 - t251 * t191) * t165 / 0.2e1 + (t251 * t189 + t252 * t191) * t166 / 0.2e1 + ((-t189 * t155 + t158 * t191 + Icges(1,4)) * V_base(5) + (-t189 * t156 + t159 * t191 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t155 * t191 + t189 * t158 + Icges(1,2)) * V_base(5) + (t156 * t191 + t189 * t159 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t101 * t171 + t103 * t170) * t144 + (t100 * t171 + t102 * t170) * t143 + (t167 * t94 + t168 * t92) * t132 + (t167 * t93 + t168 * t91) * t131 + (t113 * t175 + t115 * t174 + t122 * t190 + t124 * t188) * t166 + (t112 * t175 + t114 * t174 + t121 * t190 + t123 * t188) * t165 + (t126 * t168 + t127 * t167 + t134 * t171 + t135 * t170 + t140 * t175 + t141 * t174 + t154 * t190 + t157 * t188 + Icges(2,3)) * t177) * t177 / 0.2e1 + t177 * V_base(4) * (Icges(2,5) * t191 - Icges(2,6) * t189) + V_base(5) * t177 * (Icges(2,5) * t189 + Icges(2,6) * t191) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
