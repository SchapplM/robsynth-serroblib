% Calculate kinetic energy for
% S5RRRPR8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPR8_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR8_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRPR8_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR8_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR8_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR8_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:19:02
% EndTime: 2019-12-31 21:19:05
% DurationCPUTime: 2.97s
% Computational Cost: add. (1146->261), mult. (1388->378), div. (0->0), fcn. (1226->8), ass. (0->136)
t267 = Icges(4,4) + Icges(5,6);
t266 = Icges(4,1) + Icges(5,2);
t265 = -Icges(4,2) - Icges(5,3);
t187 = qJ(2) + qJ(3);
t184 = cos(t187);
t264 = t267 * t184;
t183 = sin(t187);
t263 = t267 * t183;
t262 = Icges(5,4) - Icges(4,5);
t261 = Icges(5,5) - Icges(4,6);
t260 = t265 * t183 + t264;
t259 = t266 * t184 - t263;
t258 = Icges(5,1) + Icges(4,3);
t190 = sin(qJ(1));
t193 = cos(qJ(1));
t257 = t260 * t190 + t261 * t193;
t256 = -t261 * t190 + t260 * t193;
t255 = t259 * t190 + t262 * t193;
t254 = -t262 * t190 + t259 * t193;
t253 = t265 * t184 - t263;
t252 = t266 * t183 + t264;
t251 = t261 * t183 - t262 * t184;
t153 = V_base(5) + (-qJD(2) - qJD(3)) * t193;
t178 = qJD(2) * t190 + V_base(4);
t154 = qJD(3) * t190 + t178;
t180 = V_base(6) + qJD(1);
t250 = (t253 * t183 + t252 * t184) * t180 + (-t256 * t183 + t254 * t184) * t154 + (-t257 * t183 + t255 * t184) * t153;
t249 = (-t262 * t183 - t261 * t184) * t180 + (t258 * t190 + t251 * t193) * t154 + (t251 * t190 - t258 * t193) * t153;
t189 = sin(qJ(2));
t245 = pkin(2) * t189;
t244 = pkin(8) * t183;
t192 = cos(qJ(2));
t243 = pkin(2) * t192;
t241 = Icges(2,4) * t190;
t240 = Icges(3,4) * t189;
t239 = Icges(3,4) * t192;
t234 = t184 * t190;
t233 = t184 * t193;
t188 = sin(qJ(5));
t232 = t188 * t190;
t231 = t188 * t193;
t191 = cos(qJ(5));
t230 = t190 * t191;
t229 = t191 * t193;
t106 = -pkin(7) * t193 + t190 * t243;
t175 = t190 * pkin(1) - t193 * pkin(6);
t228 = -t106 - t175;
t227 = qJD(4) * t183;
t226 = qJD(5) * t184;
t225 = V_base(5) * pkin(5) + V_base(1);
t217 = pkin(3) * t184 + qJ(4) * t183;
t135 = t217 * t190;
t222 = -t135 + t228;
t177 = -qJD(2) * t193 + V_base(5);
t221 = t177 * t245 + t225;
t220 = rSges(3,1) * t192 - rSges(3,2) * t189;
t219 = rSges(4,1) * t184 - rSges(4,2) * t183;
t218 = -rSges(5,2) * t184 + rSges(5,3) * t183;
t216 = Icges(3,1) * t192 - t240;
t214 = -Icges(3,2) * t189 + t239;
t211 = Icges(3,5) * t192 - Icges(3,6) * t189;
t150 = pkin(3) * t183 - qJ(4) * t184;
t207 = t153 * t150 + t193 * t227 + t221;
t176 = t193 * pkin(1) + t190 * pkin(6);
t206 = -V_base(4) * pkin(5) + t180 * t176 + V_base(2);
t205 = V_base(4) * t175 - t176 * V_base(5) + V_base(3);
t202 = (-Icges(3,3) * t193 + t190 * t211) * t177 + (Icges(3,3) * t190 + t193 * t211) * t178 + (Icges(3,5) * t189 + Icges(3,6) * t192) * t180;
t107 = pkin(7) * t190 + t193 * t243;
t201 = t178 * t106 - t107 * t177 + t205;
t200 = t180 * t107 - t178 * t245 + t206;
t136 = t217 * t193;
t199 = t180 * t136 + t190 * t227 + t200;
t198 = -qJD(4) * t184 + t154 * t135 + t201;
t129 = -Icges(3,6) * t193 + t190 * t214;
t130 = Icges(3,6) * t190 + t193 * t214;
t131 = -Icges(3,5) * t193 + t190 * t216;
t132 = Icges(3,5) * t190 + t193 * t216;
t164 = Icges(3,2) * t192 + t240;
t167 = Icges(3,1) * t189 + t239;
t195 = (-t130 * t189 + t132 * t192) * t178 + (-t129 * t189 + t131 * t192) * t177 + (-t164 * t189 + t167 * t192) * t180;
t185 = Icges(2,4) * t193;
t172 = rSges(2,1) * t193 - rSges(2,2) * t190;
t171 = rSges(2,1) * t190 + rSges(2,2) * t193;
t170 = rSges(3,1) * t189 + rSges(3,2) * t192;
t169 = Icges(2,1) * t193 - t241;
t168 = Icges(2,1) * t190 + t185;
t166 = -Icges(2,2) * t190 + t185;
t165 = Icges(2,2) * t193 + t241;
t160 = qJD(5) * t183 + t180;
t159 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t158 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t157 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t152 = rSges(4,1) * t183 + rSges(4,2) * t184;
t151 = -rSges(5,2) * t183 - rSges(5,3) * t184;
t143 = -pkin(4) * t193 + pkin(8) * t234;
t142 = pkin(4) * t190 + pkin(8) * t233;
t140 = t183 * t232 - t229;
t139 = t183 * t230 + t231;
t138 = t183 * t231 + t230;
t137 = t183 * t229 - t232;
t134 = rSges(3,3) * t190 + t193 * t220;
t133 = -rSges(3,3) * t193 + t190 * t220;
t126 = t193 * t226 + t154;
t125 = t190 * t226 + t153;
t124 = -rSges(5,1) * t193 + t190 * t218;
t123 = rSges(5,1) * t190 + t193 * t218;
t122 = rSges(4,3) * t190 + t193 * t219;
t121 = -rSges(4,3) * t193 + t190 * t219;
t105 = rSges(6,3) * t183 + (-rSges(6,1) * t188 - rSges(6,2) * t191) * t184;
t104 = Icges(6,5) * t183 + (-Icges(6,1) * t188 - Icges(6,4) * t191) * t184;
t103 = Icges(6,6) * t183 + (-Icges(6,4) * t188 - Icges(6,2) * t191) * t184;
t102 = Icges(6,3) * t183 + (-Icges(6,5) * t188 - Icges(6,6) * t191) * t184;
t101 = V_base(5) * rSges(2,3) - t171 * t180 + t225;
t100 = t172 * t180 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t98 = t171 * V_base(4) - t172 * V_base(5) + V_base(3);
t94 = rSges(6,1) * t140 + rSges(6,2) * t139 + rSges(6,3) * t234;
t93 = rSges(6,1) * t138 + rSges(6,2) * t137 + rSges(6,3) * t233;
t92 = Icges(6,1) * t140 + Icges(6,4) * t139 + Icges(6,5) * t234;
t91 = Icges(6,1) * t138 + Icges(6,4) * t137 + Icges(6,5) * t233;
t90 = Icges(6,4) * t140 + Icges(6,2) * t139 + Icges(6,6) * t234;
t89 = Icges(6,4) * t138 + Icges(6,2) * t137 + Icges(6,6) * t233;
t88 = Icges(6,5) * t140 + Icges(6,6) * t139 + Icges(6,3) * t234;
t87 = Icges(6,5) * t138 + Icges(6,6) * t137 + Icges(6,3) * t233;
t86 = t170 * t177 + (-t133 - t175) * t180 + t225;
t85 = t134 * t180 - t170 * t178 + t206;
t84 = t133 * t178 - t134 * t177 + t205;
t83 = t152 * t153 + (-t121 + t228) * t180 + t221;
t82 = t122 * t180 - t152 * t154 + t200;
t81 = t121 * t154 - t122 * t153 + t201;
t80 = t151 * t153 + (-t124 + t222) * t180 + t207;
t79 = t123 * t180 + (-t150 - t151) * t154 + t199;
t78 = t124 * t154 + (-t123 - t136) * t153 + t198;
t77 = t153 * t244 + t105 * t125 - t160 * t94 + (-t143 + t222) * t180 + t207;
t76 = -t105 * t126 + t142 * t180 + t160 * t93 + (-t150 - t244) * t154 + t199;
t75 = -t125 * t93 + t126 * t94 + t143 * t154 + (-t136 - t142) * t153 + t198;
t1 = m(1) * (t157 ^ 2 + t158 ^ 2 + t159 ^ 2) / 0.2e1 + m(2) * (t100 ^ 2 + t101 ^ 2 + t98 ^ 2) / 0.2e1 + m(3) * (t84 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + t178 * (t202 * t190 + t195 * t193) / 0.2e1 + t177 * (t195 * t190 - t202 * t193) / 0.2e1 + m(4) * (t81 ^ 2 + t82 ^ 2 + t83 ^ 2) / 0.2e1 + m(5) * (t78 ^ 2 + t79 ^ 2 + t80 ^ 2) / 0.2e1 + m(6) * (t75 ^ 2 + t76 ^ 2 + t77 ^ 2) / 0.2e1 + t126 * ((t137 * t89 + t138 * t91 + t87 * t233) * t126 + (t137 * t90 + t138 * t92 + t233 * t88) * t125 + (t102 * t233 + t103 * t137 + t104 * t138) * t160) / 0.2e1 + t125 * ((t139 * t89 + t140 * t91 + t234 * t87) * t126 + (t139 * t90 + t140 * t92 + t88 * t234) * t125 + (t102 * t234 + t103 * t139 + t104 * t140) * t160) / 0.2e1 + t160 * ((t102 * t160 + t88 * t125 + t87 * t126) * t183 + ((-t188 * t91 - t191 * t89) * t126 + (-t188 * t92 - t191 * t90) * t125 + (-t103 * t191 - t104 * t188) * t160) * t184) / 0.2e1 + (t190 * t250 - t249 * t193) * t153 / 0.2e1 + (t249 * t190 + t193 * t250) * t154 / 0.2e1 + ((-t165 * t190 + t168 * t193 + Icges(1,4)) * V_base(5) + (-t166 * t190 + t169 * t193 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t165 * t193 + t168 * t190 + Icges(1,2)) * V_base(5) + (t166 * t193 + t169 * t190 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t130 * t192 + t132 * t189) * t178 + (t129 * t192 + t131 * t189) * t177 + (t254 * t183 + t256 * t184) * t154 + (t255 * t183 + t257 * t184) * t153 + (t164 * t192 + t167 * t189 + t252 * t183 - t253 * t184 + Icges(2,3)) * t180) * t180 / 0.2e1 + t180 * V_base(4) * (Icges(2,5) * t193 - Icges(2,6) * t190) + V_base(5) * t180 * (Icges(2,5) * t190 + Icges(2,6) * t193) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
