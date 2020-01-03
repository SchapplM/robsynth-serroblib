% Calculate kinetic energy for
% S5RRPPR7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPPR7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR7_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPPR7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR7_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR7_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR7_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:34:51
% EndTime: 2019-12-31 19:34:54
% DurationCPUTime: 2.82s
% Computational Cost: add. (1104->257), mult. (1346->356), div. (0->0), fcn. (1184->8), ass. (0->133)
t267 = Icges(4,4) + Icges(5,6);
t266 = Icges(4,1) + Icges(5,2);
t265 = -Icges(4,2) - Icges(5,3);
t185 = qJ(2) + pkin(8);
t179 = cos(t185);
t264 = t267 * t179;
t178 = sin(t185);
t263 = t267 * t178;
t262 = Icges(5,4) - Icges(4,5);
t261 = Icges(5,5) - Icges(4,6);
t260 = t265 * t178 + t264;
t259 = t266 * t179 - t263;
t189 = sin(qJ(1));
t192 = cos(qJ(1));
t258 = t260 * t189 + t261 * t192;
t257 = -t261 * t189 + t260 * t192;
t256 = t259 * t189 + t262 * t192;
t255 = -t262 * t189 + t259 * t192;
t254 = t265 * t179 - t263;
t253 = t266 * t178 + t264;
t252 = Icges(5,1) + Icges(3,3) + Icges(4,3);
t188 = sin(qJ(2));
t191 = cos(qJ(2));
t251 = Icges(3,5) * t191 - Icges(3,6) * t188 + t261 * t178 - t262 * t179;
t239 = Icges(3,4) * t191;
t212 = -Icges(3,2) * t188 + t239;
t127 = -Icges(3,6) * t192 + t189 * t212;
t128 = Icges(3,6) * t189 + t192 * t212;
t240 = Icges(3,4) * t188;
t214 = Icges(3,1) * t191 - t240;
t129 = -Icges(3,5) * t192 + t189 * t214;
t130 = Icges(3,5) * t189 + t192 * t214;
t162 = Icges(3,2) * t191 + t240;
t165 = Icges(3,1) * t188 + t239;
t175 = -qJD(2) * t192 + V_base(5);
t176 = qJD(2) * t189 + V_base(4);
t180 = V_base(6) + qJD(1);
t250 = (-t162 * t188 + t165 * t191 + t254 * t178 + t253 * t179) * t180 + (-t128 * t188 + t130 * t191 - t257 * t178 + t255 * t179) * t176 + (-t127 * t188 + t129 * t191 - t258 * t178 + t256 * t179) * t175;
t249 = (Icges(3,5) * t188 + Icges(3,6) * t191 - t262 * t178 - t261 * t179) * t180 + (t252 * t189 + t251 * t192) * t176 + (t251 * t189 - t252 * t192) * t175;
t245 = pkin(2) * t188;
t244 = pkin(7) * t178;
t243 = pkin(2) * t191;
t241 = Icges(2,4) * t189;
t234 = t179 * t189;
t233 = t179 * t192;
t187 = sin(qJ(5));
t232 = t187 * t192;
t231 = t189 * t187;
t190 = cos(qJ(5));
t230 = t189 * t190;
t229 = t190 * t192;
t106 = -qJ(3) * t192 + t189 * t243;
t173 = t189 * pkin(1) - pkin(6) * t192;
t228 = -t106 - t173;
t107 = qJ(3) * t189 + t192 * t243;
t215 = pkin(3) * t179 + qJ(4) * t178;
t132 = t215 * t192;
t227 = -t107 - t132;
t226 = qJD(4) * t178;
t225 = qJD(5) * t179;
t224 = V_base(5) * pkin(5) + V_base(1);
t131 = t215 * t189;
t221 = -t131 + t228;
t148 = pkin(3) * t178 - qJ(4) * t179;
t220 = -t148 - t245;
t219 = qJD(3) * t189 + t175 * t245 + t224;
t218 = rSges(3,1) * t191 - rSges(3,2) * t188;
t217 = rSges(4,1) * t179 - rSges(4,2) * t178;
t216 = -rSges(5,2) * t179 + rSges(5,3) * t178;
t174 = pkin(1) * t192 + t189 * pkin(6);
t205 = -V_base(4) * pkin(5) + t180 * t174 + V_base(2);
t204 = V_base(4) * t173 - t174 * V_base(5) + V_base(3);
t203 = t175 * t148 + t192 * t226 + t219;
t202 = t176 * t106 + t204;
t198 = -qJD(3) * t192 + t180 * t107 + t205;
t197 = -qJD(4) * t179 + t176 * t131 + t202;
t196 = t180 * t132 + t189 * t226 + t198;
t183 = Icges(2,4) * t192;
t172 = rSges(2,1) * t192 - t189 * rSges(2,2);
t171 = t189 * rSges(2,1) + rSges(2,2) * t192;
t170 = rSges(3,1) * t188 + rSges(3,2) * t191;
t167 = Icges(2,1) * t192 - t241;
t166 = Icges(2,1) * t189 + t183;
t164 = -Icges(2,2) * t189 + t183;
t163 = Icges(2,2) * t192 + t241;
t158 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t157 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t156 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t155 = qJD(5) * t178 + t180;
t152 = -pkin(4) * t192 + pkin(7) * t234;
t151 = t189 * pkin(4) + pkin(7) * t233;
t150 = rSges(4,1) * t178 + rSges(4,2) * t179;
t149 = -rSges(5,2) * t178 - rSges(5,3) * t179;
t140 = t178 * t231 - t229;
t139 = t178 * t230 + t232;
t138 = t178 * t232 + t230;
t137 = t178 * t229 - t231;
t136 = t192 * t225 + t176;
t135 = t189 * t225 + t175;
t134 = t189 * rSges(3,3) + t192 * t218;
t133 = -rSges(3,3) * t192 + t189 * t218;
t123 = -rSges(5,1) * t192 + t189 * t216;
t122 = t189 * rSges(5,1) + t192 * t216;
t121 = t189 * rSges(4,3) + t192 * t217;
t120 = -rSges(4,3) * t192 + t189 * t217;
t104 = rSges(6,3) * t178 + (-rSges(6,1) * t187 - rSges(6,2) * t190) * t179;
t103 = V_base(5) * rSges(2,3) - t171 * t180 + t224;
t102 = t172 * t180 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t101 = Icges(6,5) * t178 + (-Icges(6,1) * t187 - Icges(6,4) * t190) * t179;
t100 = Icges(6,6) * t178 + (-Icges(6,4) * t187 - Icges(6,2) * t190) * t179;
t99 = Icges(6,3) * t178 + (-Icges(6,5) * t187 - Icges(6,6) * t190) * t179;
t97 = t171 * V_base(4) - t172 * V_base(5) + V_base(3);
t94 = rSges(6,1) * t140 + rSges(6,2) * t139 + rSges(6,3) * t234;
t93 = t138 * rSges(6,1) + t137 * rSges(6,2) + rSges(6,3) * t233;
t92 = Icges(6,1) * t140 + Icges(6,4) * t139 + Icges(6,5) * t234;
t91 = Icges(6,1) * t138 + Icges(6,4) * t137 + Icges(6,5) * t233;
t90 = Icges(6,4) * t140 + Icges(6,2) * t139 + Icges(6,6) * t234;
t89 = Icges(6,4) * t138 + Icges(6,2) * t137 + Icges(6,6) * t233;
t88 = Icges(6,5) * t140 + Icges(6,6) * t139 + Icges(6,3) * t234;
t87 = Icges(6,5) * t138 + Icges(6,6) * t137 + Icges(6,3) * t233;
t86 = t170 * t175 + (-t133 - t173) * t180 + t224;
t85 = t134 * t180 - t170 * t176 + t205;
t84 = t133 * t176 - t134 * t175 + t204;
t83 = t150 * t175 + (-t120 + t228) * t180 + t219;
t82 = t180 * t121 + (-t150 - t245) * t176 + t198;
t81 = t120 * t176 + (-t107 - t121) * t175 + t202;
t80 = t149 * t175 + (-t123 + t221) * t180 + t203;
t79 = t180 * t122 + (-t149 + t220) * t176 + t196;
t78 = t123 * t176 + (-t122 + t227) * t175 + t197;
t77 = t175 * t244 + t104 * t135 - t155 * t94 + (-t152 + t221) * t180 + t203;
t76 = -t136 * t104 + t180 * t151 + t155 * t93 + (t220 - t244) * t176 + t196;
t75 = -t135 * t93 + t136 * t94 + t152 * t176 + (-t151 + t227) * t175 + t197;
t1 = m(1) * (t156 ^ 2 + t157 ^ 2 + t158 ^ 2) / 0.2e1 + m(2) * (t102 ^ 2 + t103 ^ 2 + t97 ^ 2) / 0.2e1 + m(3) * (t84 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + m(4) * (t81 ^ 2 + t82 ^ 2 + t83 ^ 2) / 0.2e1 + m(5) * (t78 ^ 2 + t79 ^ 2 + t80 ^ 2) / 0.2e1 + m(6) * (t75 ^ 2 + t76 ^ 2 + t77 ^ 2) / 0.2e1 + t136 * ((t137 * t89 + t138 * t91 + t87 * t233) * t136 + (t137 * t90 + t138 * t92 + t233 * t88) * t135 + (t137 * t100 + t138 * t101 + t233 * t99) * t155) / 0.2e1 + t135 * ((t139 * t89 + t140 * t91 + t234 * t87) * t136 + (t139 * t90 + t140 * t92 + t88 * t234) * t135 + (t100 * t139 + t101 * t140 + t234 * t99) * t155) / 0.2e1 + t155 * ((t88 * t135 + t87 * t136 + t99 * t155) * t178 + ((-t187 * t91 - t190 * t89) * t136 + (-t187 * t92 - t190 * t90) * t135 + (-t100 * t190 - t101 * t187) * t155) * t179) / 0.2e1 + ((-t189 * t163 + t166 * t192 + Icges(1,4)) * V_base(5) + (-t189 * t164 + t167 * t192 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t163 * t192 + t189 * t166 + Icges(1,2)) * V_base(5) + (t164 * t192 + t189 * t167 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (t250 * t189 - t249 * t192) * t175 / 0.2e1 + (t249 * t189 + t250 * t192) * t176 / 0.2e1 + ((t128 * t191 + t130 * t188 + t255 * t178 + t257 * t179) * t176 + (t127 * t191 + t129 * t188 + t256 * t178 + t258 * t179) * t175 + (t191 * t162 + t188 * t165 + t253 * t178 - t254 * t179 + Icges(2,3)) * t180) * t180 / 0.2e1 + V_base(4) * t180 * (Icges(2,5) * t192 - Icges(2,6) * t189) + V_base(5) * t180 * (Icges(2,5) * t189 + Icges(2,6) * t192) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
