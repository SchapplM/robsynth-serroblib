% Calculate kinetic energy for
% S5PPRPR2
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRPR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR2_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PPRPR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR2_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR2_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRPR2_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:02:54
% EndTime: 2019-12-05 15:02:58
% DurationCPUTime: 3.40s
% Computational Cost: add. (1037->264), mult. (1324->365), div. (0->0), fcn. (1162->8), ass. (0->137)
t268 = Icges(4,4) + Icges(5,6);
t267 = Icges(4,1) + Icges(5,2);
t266 = -Icges(4,2) - Icges(5,3);
t184 = pkin(8) + qJ(3);
t181 = cos(t184);
t265 = t268 * t181;
t180 = sin(t184);
t264 = t268 * t180;
t263 = Icges(5,4) - Icges(4,5);
t262 = Icges(5,5) - Icges(4,6);
t261 = t180 * t266 + t265;
t260 = t181 * t267 - t264;
t259 = Icges(5,1) + Icges(4,3);
t186 = sin(pkin(7));
t188 = cos(pkin(7));
t258 = t186 * t261 + t188 * t262;
t257 = -t186 * t262 + t188 * t261;
t256 = t186 * t260 + t188 * t263;
t255 = -t186 * t263 + t188 * t260;
t254 = t181 * t266 - t264;
t253 = t180 * t267 + t265;
t252 = t180 * t262 - t181 * t263;
t251 = Icges(2,5) * t188 - Icges(2,6) * t186 + Icges(1,5);
t250 = Icges(2,5) * t186 + Icges(2,6) * t188 + Icges(1,6);
t174 = -qJD(3) * t188 + V_base(5);
t175 = qJD(3) * t186 + V_base(4);
t249 = (t180 * t254 + t181 * t253) * V_base(6) + (-t180 * t257 + t181 * t255) * t175 + (-t180 * t258 + t181 * t256) * t174;
t248 = (-t180 * t263 - t181 * t262) * V_base(6) + (t186 * t259 + t188 * t252) * t175 + (t186 * t252 - t188 * t259) * t174;
t185 = sin(pkin(8));
t245 = pkin(2) * t185;
t244 = pkin(6) * t180;
t187 = cos(pkin(8));
t243 = pkin(2) * t187;
t168 = pkin(1) * t186 - qJ(2) * t188;
t99 = -pkin(5) * t188 + t186 * t243;
t242 = -t168 - t99;
t241 = Icges(2,4) * t186;
t240 = Icges(3,4) * t185;
t239 = Icges(3,4) * t187;
t234 = t181 * t186;
t233 = t181 * t188;
t190 = sin(qJ(5));
t232 = t186 * t190;
t191 = cos(qJ(5));
t231 = t186 * t191;
t230 = t188 * t190;
t229 = t188 * t191;
t227 = qJD(4) * t180;
t226 = qJD(5) * t181;
t225 = V_base(5) * qJ(1) + V_base(1);
t221 = qJD(1) + V_base(3);
t213 = pkin(3) * t181 + qJ(4) * t180;
t133 = t213 * t186;
t220 = -t133 + t242;
t219 = qJD(2) * t186 + t225;
t218 = t168 * V_base(4) + t221;
t217 = t245 * V_base(5) + t219;
t216 = rSges(3,1) * t187 - rSges(3,2) * t185;
t215 = rSges(4,1) * t181 - rSges(4,2) * t180;
t214 = -rSges(5,2) * t181 + rSges(5,3) * t180;
t212 = Icges(3,1) * t187 - t240;
t210 = -Icges(3,2) * t185 + t239;
t207 = Icges(3,5) * t187 - Icges(3,6) * t185;
t170 = pkin(1) * t188 + qJ(2) * t186;
t203 = -qJD(2) * t188 + t170 * V_base(6) + V_base(2);
t149 = pkin(3) * t180 - qJ(4) * t181;
t202 = t149 * t174 + t188 * t227 + t217;
t100 = pkin(5) * t186 + t188 * t243;
t199 = V_base(4) * t99 + (-t100 - t170) * V_base(5) + t218;
t198 = (-Icges(3,3) * t188 + t186 * t207) * V_base(5) + (Icges(3,3) * t186 + t188 * t207) * V_base(4) + (Icges(3,5) * t185 + Icges(3,6) * t187) * V_base(6);
t197 = V_base(6) * t100 + (-qJ(1) - t245) * V_base(4) + t203;
t196 = -qJD(4) * t181 + t133 * t175 + t199;
t134 = t213 * t188;
t195 = t134 * V_base(6) + t186 * t227 + t197;
t127 = -Icges(3,6) * t188 + t186 * t210;
t128 = Icges(3,6) * t186 + t188 * t210;
t129 = -Icges(3,5) * t188 + t186 * t212;
t130 = Icges(3,5) * t186 + t188 * t212;
t160 = Icges(3,2) * t187 + t240;
t163 = Icges(3,1) * t185 + t239;
t192 = (-t128 * t185 + t130 * t187) * V_base(4) + (-t127 * t185 + t129 * t187) * V_base(5) + (-t160 * t185 + t163 * t187) * V_base(6);
t182 = Icges(2,4) * t188;
t171 = rSges(2,1) * t188 - rSges(2,2) * t186;
t169 = rSges(2,1) * t186 + rSges(2,2) * t188;
t167 = rSges(3,1) * t185 + rSges(3,2) * t187;
t166 = qJD(5) * t180 + V_base(6);
t165 = Icges(2,1) * t188 - t241;
t164 = Icges(2,1) * t186 + t182;
t162 = -Icges(2,2) * t186 + t182;
t161 = Icges(2,2) * t188 + t241;
t156 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t155 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t154 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t151 = rSges(4,1) * t180 + rSges(4,2) * t181;
t150 = -rSges(5,2) * t180 - rSges(5,3) * t181;
t142 = -pkin(4) * t188 + pkin(6) * t234;
t141 = pkin(4) * t186 + pkin(6) * t233;
t140 = t180 * t232 - t229;
t139 = t180 * t231 + t230;
t138 = t180 * t230 + t231;
t137 = t180 * t229 - t232;
t136 = t188 * t226 + t175;
t135 = t186 * t226 + t174;
t132 = rSges(3,3) * t186 + t188 * t216;
t131 = -rSges(3,3) * t188 + t186 * t216;
t122 = -rSges(5,1) * t188 + t186 * t214;
t121 = rSges(5,1) * t186 + t188 * t214;
t120 = rSges(4,3) * t186 + t188 * t215;
t119 = -rSges(4,3) * t188 + t186 * t215;
t118 = V_base(5) * rSges(2,3) - t169 * V_base(6) + t225;
t117 = t171 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t104 = t180 * rSges(6,3) + (-rSges(6,1) * t190 - rSges(6,2) * t191) * t181;
t103 = Icges(6,5) * t180 + (-Icges(6,1) * t190 - Icges(6,4) * t191) * t181;
t102 = Icges(6,6) * t180 + (-Icges(6,4) * t190 - Icges(6,2) * t191) * t181;
t101 = Icges(6,3) * t180 + (-Icges(6,5) * t190 - Icges(6,6) * t191) * t181;
t95 = t169 * V_base(4) - t171 * V_base(5) + t221;
t94 = rSges(6,1) * t140 + rSges(6,2) * t139 + rSges(6,3) * t234;
t93 = rSges(6,1) * t138 + rSges(6,2) * t137 + rSges(6,3) * t233;
t92 = Icges(6,1) * t140 + Icges(6,4) * t139 + Icges(6,5) * t234;
t91 = Icges(6,1) * t138 + Icges(6,4) * t137 + Icges(6,5) * t233;
t90 = Icges(6,4) * t140 + Icges(6,2) * t139 + Icges(6,6) * t234;
t89 = Icges(6,4) * t138 + Icges(6,2) * t137 + Icges(6,6) * t233;
t88 = Icges(6,5) * t140 + Icges(6,6) * t139 + Icges(6,3) * t234;
t87 = Icges(6,5) * t138 + Icges(6,6) * t137 + Icges(6,3) * t233;
t86 = t167 * V_base(5) + (-t131 - t168) * V_base(6) + t219;
t85 = t132 * V_base(6) + (-qJ(1) - t167) * V_base(4) + t203;
t84 = t131 * V_base(4) + (-t132 - t170) * V_base(5) + t218;
t83 = t151 * t174 + (-t119 + t242) * V_base(6) + t217;
t82 = t120 * V_base(6) - t151 * t175 + t197;
t81 = t119 * t175 - t120 * t174 + t199;
t80 = t150 * t174 + (-t122 + t220) * V_base(6) + t202;
t79 = t121 * V_base(6) + (-t149 - t150) * t175 + t195;
t78 = t122 * t175 + (-t121 - t134) * t174 + t196;
t77 = t174 * t244 + t104 * t135 - t166 * t94 + (-t142 + t220) * V_base(6) + t202;
t76 = -t104 * t136 + t141 * V_base(6) + t166 * t93 + (-t149 - t244) * t175 + t195;
t75 = -t135 * t93 + t136 * t94 + t142 * t175 + (-t134 - t141) * t174 + t196;
t1 = m(1) * (t154 ^ 2 + t155 ^ 2 + t156 ^ 2) / 0.2e1 + m(2) * (t117 ^ 2 + t118 ^ 2 + t95 ^ 2) / 0.2e1 + m(3) * (t84 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + m(4) * (t81 ^ 2 + t82 ^ 2 + t83 ^ 2) / 0.2e1 + m(5) * (t78 ^ 2 + t79 ^ 2 + t80 ^ 2) / 0.2e1 + m(6) * (t75 ^ 2 + t76 ^ 2 + t77 ^ 2) / 0.2e1 + t136 * ((t137 * t89 + t138 * t91 + t87 * t233) * t136 + (t137 * t90 + t138 * t92 + t233 * t88) * t135 + (t101 * t233 + t102 * t137 + t103 * t138) * t166) / 0.2e1 + t135 * ((t139 * t89 + t140 * t91 + t234 * t87) * t136 + (t139 * t90 + t140 * t92 + t88 * t234) * t135 + (t101 * t234 + t102 * t139 + t103 * t140) * t166) / 0.2e1 + t166 * ((t101 * t166 + t135 * t88 + t136 * t87) * t180 + ((-t190 * t91 - t191 * t89) * t136 + (-t190 * t92 - t191 * t90) * t135 + (-t102 * t191 - t103 * t190) * t166) * t181) / 0.2e1 + (t249 * t186 - t248 * t188) * t174 / 0.2e1 + (t248 * t186 + t249 * t188) * t175 / 0.2e1 + (t198 * t186 + t192 * t188 + t251 * V_base(6) + (-t161 * t186 + t164 * t188 + Icges(1,4)) * V_base(5) + (-t186 * t162 + t188 * t165 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t192 * t186 - t198 * t188 + t250 * V_base(6) + (t188 * t161 + t186 * t164 + Icges(1,2)) * V_base(5) + (t162 * t188 + t165 * t186 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t180 * t255 + t181 * t257) * t175 + (t180 * t256 + t181 * t258) * t174 + (t127 * t187 + t129 * t185 + t250) * V_base(5) + (t128 * t187 + t130 * t185 + t251) * V_base(4) + (t187 * t160 + t185 * t163 + t253 * t180 - t254 * t181 + Icges(1,3) + Icges(2,3)) * V_base(6)) * V_base(6) / 0.2e1;
T = t1;
