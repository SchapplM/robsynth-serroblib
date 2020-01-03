% Calculate kinetic energy for
% S5RPRPR13
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
% Datum: 2019-12-31 18:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR13_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR13_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRPR13_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR13_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR13_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR13_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:31:52
% EndTime: 2019-12-31 18:31:55
% DurationCPUTime: 2.99s
% Computational Cost: add. (1082->264), mult. (1324->370), div. (0->0), fcn. (1162->8), ass. (0->137)
t266 = Icges(4,4) + Icges(5,6);
t265 = Icges(4,1) + Icges(5,2);
t264 = -Icges(4,2) - Icges(5,3);
t185 = pkin(8) + qJ(3);
t179 = cos(t185);
t263 = t266 * t179;
t178 = sin(t185);
t262 = t266 * t178;
t261 = Icges(5,4) - Icges(4,5);
t260 = Icges(5,5) - Icges(4,6);
t259 = t264 * t178 + t263;
t258 = t265 * t179 - t262;
t257 = Icges(5,1) + Icges(4,3);
t190 = sin(qJ(1));
t192 = cos(qJ(1));
t256 = t259 * t190 + t260 * t192;
t255 = -t260 * t190 + t259 * t192;
t254 = t258 * t190 + t261 * t192;
t253 = -t261 * t190 + t258 * t192;
t252 = t264 * t179 - t262;
t251 = t265 * t178 + t263;
t250 = t260 * t178 - t261 * t179;
t174 = -qJD(3) * t192 + V_base(5);
t175 = qJD(3) * t190 + V_base(4);
t180 = V_base(6) + qJD(1);
t249 = (t252 * t178 + t251 * t179) * t180 + (-t255 * t178 + t253 * t179) * t175 + (-t256 * t178 + t254 * t179) * t174;
t248 = (-t261 * t178 - t260 * t179) * t180 + (t257 * t190 + t250 * t192) * t175 + (t250 * t190 - t257 * t192) * t174;
t186 = sin(pkin(8));
t244 = pkin(2) * t186;
t243 = pkin(7) * t178;
t187 = cos(pkin(8));
t242 = pkin(2) * t187;
t241 = Icges(2,4) * t190;
t240 = Icges(3,4) * t186;
t239 = Icges(3,4) * t187;
t234 = t179 * t190;
t233 = t179 * t192;
t189 = sin(qJ(5));
t232 = t189 * t192;
t231 = t190 * t189;
t191 = cos(qJ(5));
t230 = t190 * t191;
t229 = t191 * t192;
t105 = -pkin(6) * t192 + t242 * t190;
t170 = t190 * pkin(1) - qJ(2) * t192;
t227 = -t105 - t170;
t226 = qJD(4) * t178;
t225 = qJD(5) * t179;
t224 = V_base(4) * t170 + V_base(3);
t223 = V_base(5) * pkin(5) + V_base(1);
t214 = pkin(3) * t179 + qJ(4) * t178;
t133 = t214 * t190;
t220 = -t133 + t227;
t219 = qJD(2) * t190 + t223;
t218 = V_base(5) * t244 + t219;
t217 = rSges(3,1) * t187 - rSges(3,2) * t186;
t216 = rSges(4,1) * t179 - rSges(4,2) * t178;
t215 = -rSges(5,2) * t179 + rSges(5,3) * t178;
t213 = Icges(3,1) * t187 - t240;
t211 = -Icges(3,2) * t186 + t239;
t208 = Icges(3,5) * t187 - Icges(3,6) * t186;
t172 = pkin(1) * t192 + t190 * qJ(2);
t204 = -qJD(2) * t192 + t180 * t172 + V_base(2);
t148 = pkin(3) * t178 - qJ(4) * t179;
t203 = t174 * t148 + t192 * t226 + t218;
t106 = pkin(6) * t190 + t242 * t192;
t200 = V_base(4) * t105 + (-t106 - t172) * V_base(5) + t224;
t199 = (-Icges(3,3) * t192 + t208 * t190) * V_base(5) + (Icges(3,3) * t190 + t208 * t192) * V_base(4) + (Icges(3,5) * t186 + Icges(3,6) * t187) * t180;
t198 = -qJD(4) * t179 + t175 * t133 + t200;
t197 = t180 * t106 + (-pkin(5) - t244) * V_base(4) + t204;
t134 = t214 * t192;
t196 = t180 * t134 + t190 * t226 + t197;
t127 = -Icges(3,6) * t192 + t211 * t190;
t128 = Icges(3,6) * t190 + t211 * t192;
t129 = -Icges(3,5) * t192 + t213 * t190;
t130 = Icges(3,5) * t190 + t213 * t192;
t159 = Icges(3,2) * t187 + t240;
t160 = Icges(3,1) * t186 + t239;
t193 = (-t128 * t186 + t130 * t187) * V_base(4) + (-t127 * t186 + t129 * t187) * V_base(5) + (-t159 * t186 + t160 * t187) * t180;
t183 = Icges(2,4) * t192;
t173 = rSges(2,1) * t192 - t190 * rSges(2,2);
t171 = t190 * rSges(2,1) + rSges(2,2) * t192;
t167 = Icges(2,1) * t192 - t241;
t166 = Icges(2,1) * t190 + t183;
t165 = -Icges(2,2) * t190 + t183;
t164 = Icges(2,2) * t192 + t241;
t163 = Icges(2,5) * t192 - Icges(2,6) * t190;
t162 = Icges(2,5) * t190 + Icges(2,6) * t192;
t161 = rSges(3,1) * t186 + rSges(3,2) * t187;
t157 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t156 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t155 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t154 = qJD(5) * t178 + t180;
t152 = -pkin(4) * t192 + pkin(7) * t234;
t151 = t190 * pkin(4) + pkin(7) * t233;
t150 = rSges(4,1) * t178 + rSges(4,2) * t179;
t149 = -rSges(5,2) * t178 - rSges(5,3) * t179;
t140 = t178 * t231 - t229;
t139 = t178 * t230 + t232;
t138 = t178 * t232 + t230;
t137 = t178 * t229 - t231;
t136 = t192 * t225 + t175;
t135 = t190 * t225 + t174;
t132 = t190 * rSges(3,3) + t217 * t192;
t131 = -rSges(3,3) * t192 + t217 * t190;
t123 = -rSges(5,1) * t192 + t215 * t190;
t122 = t190 * rSges(5,1) + t215 * t192;
t121 = t190 * rSges(4,3) + t216 * t192;
t120 = -rSges(4,3) * t192 + t216 * t190;
t104 = rSges(6,3) * t178 + (-rSges(6,1) * t189 - rSges(6,2) * t191) * t179;
t103 = V_base(5) * rSges(2,3) - t171 * t180 + t223;
t102 = t173 * t180 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t101 = Icges(6,5) * t178 + (-Icges(6,1) * t189 - Icges(6,4) * t191) * t179;
t100 = Icges(6,6) * t178 + (-Icges(6,4) * t189 - Icges(6,2) * t191) * t179;
t99 = Icges(6,3) * t178 + (-Icges(6,5) * t189 - Icges(6,6) * t191) * t179;
t97 = t171 * V_base(4) - t173 * V_base(5) + V_base(3);
t94 = rSges(6,1) * t140 + rSges(6,2) * t139 + rSges(6,3) * t234;
t93 = t138 * rSges(6,1) + t137 * rSges(6,2) + rSges(6,3) * t233;
t92 = Icges(6,1) * t140 + Icges(6,4) * t139 + Icges(6,5) * t234;
t91 = Icges(6,1) * t138 + Icges(6,4) * t137 + Icges(6,5) * t233;
t90 = Icges(6,4) * t140 + Icges(6,2) * t139 + Icges(6,6) * t234;
t89 = Icges(6,4) * t138 + Icges(6,2) * t137 + Icges(6,6) * t233;
t88 = Icges(6,5) * t140 + Icges(6,6) * t139 + Icges(6,3) * t234;
t87 = Icges(6,5) * t138 + Icges(6,6) * t137 + Icges(6,3) * t233;
t86 = t161 * V_base(5) + (-t131 - t170) * t180 + t219;
t85 = t180 * t132 + (-pkin(5) - t161) * V_base(4) + t204;
t84 = t131 * V_base(4) + (-t132 - t172) * V_base(5) + t224;
t83 = t150 * t174 + (-t120 + t227) * t180 + t218;
t82 = t180 * t121 - t175 * t150 + t197;
t81 = t120 * t175 - t121 * t174 + t200;
t80 = t149 * t174 + (-t123 + t220) * t180 + t203;
t79 = t180 * t122 + (-t148 - t149) * t175 + t196;
t78 = t123 * t175 + (-t122 - t134) * t174 + t198;
t77 = t174 * t243 + t104 * t135 - t154 * t94 + (-t152 + t220) * t180 + t203;
t76 = -t136 * t104 + t180 * t151 + t154 * t93 + (-t148 - t243) * t175 + t196;
t75 = -t135 * t93 + t136 * t94 + t152 * t175 + (-t134 - t151) * t174 + t198;
t1 = m(1) * (t155 ^ 2 + t156 ^ 2 + t157 ^ 2) / 0.2e1 + m(2) * (t102 ^ 2 + t103 ^ 2 + t97 ^ 2) / 0.2e1 + m(3) * (t84 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + m(4) * (t81 ^ 2 + t82 ^ 2 + t83 ^ 2) / 0.2e1 + m(5) * (t78 ^ 2 + t79 ^ 2 + t80 ^ 2) / 0.2e1 + m(6) * (t75 ^ 2 + t76 ^ 2 + t77 ^ 2) / 0.2e1 + t136 * ((t137 * t89 + t138 * t91 + t87 * t233) * t136 + (t137 * t90 + t138 * t92 + t88 * t233) * t135 + (t137 * t100 + t138 * t101 + t99 * t233) * t154) / 0.2e1 + t135 * ((t139 * t89 + t140 * t91 + t87 * t234) * t136 + (t139 * t90 + t140 * t92 + t88 * t234) * t135 + (t100 * t139 + t101 * t140 + t99 * t234) * t154) / 0.2e1 + t154 * ((t88 * t135 + t87 * t136 + t99 * t154) * t178 + ((-t189 * t91 - t191 * t89) * t136 + (-t189 * t92 - t191 * t90) * t135 + (-t100 * t191 - t101 * t189) * t154) * t179) / 0.2e1 + (t249 * t190 - t248 * t192) * t174 / 0.2e1 + (t248 * t190 + t249 * t192) * t175 / 0.2e1 + (t163 * t180 + t199 * t190 + t193 * t192 + (-t190 * t164 + t166 * t192 + Icges(1,4)) * V_base(5) + (-t190 * t165 + t192 * t167 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t162 * t180 + t193 * t190 - t199 * t192 + (t192 * t164 + t190 * t166 + Icges(1,2)) * V_base(5) + (t165 * t192 + t190 * t167 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t127 * t187 + t129 * t186 + t162) * V_base(5) + (t128 * t187 + t130 * t186 + t163) * V_base(4) + (t253 * t178 + t255 * t179) * t175 + (t254 * t178 + t256 * t179) * t174 + (t187 * t159 + t186 * t160 + t251 * t178 - t252 * t179 + Icges(2,3)) * t180) * t180 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
