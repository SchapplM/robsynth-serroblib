% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRPR2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR2_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR2_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR2_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR2_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:17:18
% EndTime: 2019-12-05 16:17:27
% DurationCPUTime: 4.21s
% Computational Cost: add. (10771->347), mult. (8931->483), div. (0->0), fcn. (8565->8), ass. (0->188)
t169 = sin(pkin(9));
t170 = cos(pkin(9));
t171 = sin(qJ(5));
t172 = cos(qJ(5));
t167 = pkin(8) + qJ(2);
t166 = qJ(3) + t167;
t162 = sin(t166);
t163 = cos(t166);
t240 = t170 * t171;
t113 = t162 * t240 + t163 * t172;
t239 = t170 * t172;
t114 = t162 * t239 - t163 * t171;
t245 = t162 * t169;
t68 = Icges(6,5) * t114 - Icges(6,6) * t113 + Icges(6,3) * t245;
t107 = Icges(6,4) * t114;
t71 = -Icges(6,2) * t113 + Icges(6,6) * t245 + t107;
t106 = Icges(6,4) * t113;
t75 = -Icges(6,1) * t114 - Icges(6,5) * t245 + t106;
t30 = -(t171 * t71 + t172 * t75) * t169 - t170 * t68;
t190 = (-rSges(6,1) * t171 - rSges(6,2) * t172) * t169;
t127 = qJD(5) * t190;
t230 = qJD(5) * t169;
t279 = t127 * t230;
t278 = -t113 * t71 - t114 * t75;
t115 = t162 * t172 - t163 * t240;
t116 = t162 * t171 + t163 * t239;
t277 = t115 * t71 - t116 * t75;
t243 = t163 * t169;
t274 = t68 * t243;
t164 = sin(t167);
t254 = pkin(2) * qJD(2);
t226 = t164 * t254;
t130 = rSges(4,1) * t162 + rSges(4,2) * t163;
t168 = qJD(2) + qJD(3);
t247 = t130 * t168;
t103 = -t226 - t247;
t152 = -qJD(5) * t170 + t168;
t123 = -rSges(6,3) * t170 + (rSges(6,1) * t172 - rSges(6,2) * t171) * t169;
t220 = t123 * t230;
t204 = rSges(6,1) * t114 - rSges(6,2) * t113;
t77 = rSges(6,3) * t245 + t204;
t273 = -t152 * t77 + t162 * t220;
t259 = pkin(4) * t170;
t186 = -t259 - pkin(3) + (-rSges(6,3) - pkin(7)) * t169;
t206 = pkin(7) * t169 + t259;
t180 = -t168 * t206 + t220;
t165 = cos(t167);
t260 = pkin(2) * qJD(2) ^ 2;
t228 = t165 * t260;
t231 = qJD(4) * t168;
t211 = t163 * t231 - t228;
t84 = -qJD(5) * t114 + t115 * t168;
t85 = -qJD(5) * t113 + t116 * t168;
t207 = rSges(6,1) * t85 + rSges(6,2) * t84;
t241 = t168 * t169;
t222 = t163 * t241;
t49 = rSges(6,3) * t222 + t207;
t131 = t163 * pkin(3) + t162 * qJ(4);
t154 = qJD(4) * t163;
t98 = t131 * t168 - t154;
t22 = t162 * t279 - t152 * t49 + (t163 * t180 - t98) * t168 + t211;
t251 = t22 * t162;
t117 = t206 * t162;
t156 = t163 * qJ(4);
t129 = pkin(3) * t162 - t156;
t153 = qJD(4) * t162;
t210 = t153 - t226;
t34 = (-t117 - t129) * t168 + t210 + t273;
t225 = t165 * t254;
t242 = t163 * t170;
t238 = pkin(4) * t242 + pkin(7) * t243 + t131;
t79 = t116 * rSges(6,1) + t115 * rSges(6,2) + rSges(6,3) * t243;
t265 = t152 * t79 - t163 * t220 + t168 * t238 - t154;
t35 = t225 + t265;
t272 = (t34 * t186 * t163 + (-t34 * qJ(4) + t35 * (-rSges(6,3) * t169 - pkin(3) - t206)) * t162) * t168 + t186 * t251;
t26 = t274 + t277;
t271 = t26 - t274;
t255 = rSges(5,1) * t170;
t227 = t162 * t255;
t233 = rSges(5,2) * t245 + t163 * rSges(5,3);
t101 = t227 - t233;
t223 = t162 * t241;
t244 = t163 * t168;
t235 = rSges(5,2) * t223 + rSges(5,3) * t244;
t270 = t168 * t101 + t235;
t268 = -rSges(5,1) * t242 - t162 * rSges(5,3);
t184 = -rSges(5,2) * t243 + t131 - t268;
t269 = t168 * t184;
t122 = t168 * t129;
t143 = qJ(4) * t244;
t267 = t143 + t122;
t234 = t143 + t153;
t266 = t234 - t153 + t122;
t70 = Icges(6,5) * t116 + Icges(6,6) * t115 + Icges(6,3) * t243;
t250 = Icges(6,4) * t116;
t73 = Icges(6,2) * t115 + Icges(6,6) * t243 + t250;
t108 = Icges(6,4) * t115;
t76 = Icges(6,1) * t116 + Icges(6,5) * t243 + t108;
t25 = -t113 * t73 + t114 * t76 + t70 * t245;
t82 = -qJD(5) * t116 + t113 * t168;
t83 = qJD(5) * t115 - t114 * t168;
t256 = t83 * rSges(6,1) + t82 * rSges(6,2);
t264 = t168 * t117 + t256 - t273;
t263 = t162 * (-Icges(6,2) * t114 - t106 - t75) + t163 * (-Icges(6,2) * t116 + t108 + t76);
t262 = -t163 / 0.2e1;
t261 = pkin(2) * t164;
t119 = -Icges(6,3) * t170 + (Icges(6,5) * t172 - Icges(6,6) * t171) * t169;
t248 = Icges(6,4) * t172;
t120 = -Icges(6,6) * t170 + (-Icges(6,2) * t171 + t248) * t169;
t249 = Icges(6,4) * t171;
t121 = -Icges(6,5) * t170 + (Icges(6,1) * t172 - t249) * t169;
t50 = -t113 * t120 + t114 * t121 + t119 * t245;
t253 = t152 * t50;
t24 = t245 * t68 + t278;
t252 = t162 * t24;
t246 = t162 * t168;
t135 = (-Icges(6,1) * t171 - t248) * t169;
t237 = -t120 + t135;
t134 = (-Icges(6,2) * t172 - t249) * t169;
t236 = t121 + t134;
t27 = t115 * t73 + t116 * t76 + t70 * t243;
t229 = t164 * t260;
t218 = -pkin(3) - t255;
t217 = -t230 / 0.2e1;
t216 = t230 / 0.2e1;
t132 = t163 * rSges(4,1) - rSges(4,2) * t162;
t215 = t162 * t217;
t214 = t162 * t216;
t213 = t163 * t217;
t212 = t163 * t216;
t112 = rSges(4,1) * t244 - rSges(4,2) * t246;
t209 = -t154 + t225;
t203 = t34 * t162 - t163 * t35;
t202 = -t162 * t79 + t163 * t77;
t201 = t162 * (-Icges(6,5) * t113 - Icges(6,6) * t114) + t163 * (Icges(6,5) * t115 - Icges(6,6) * t116);
t200 = t162 * t231 + t168 * (-pkin(3) * t246 + t234) - t229;
t199 = t168 * t215;
t198 = t168 * t212;
t195 = t154 - t207;
t194 = t156 - t204;
t191 = t169 * t201;
t189 = (t163 * t25 + t252) * t169;
t188 = (t162 * t26 + t163 * t27) * t169;
t133 = (-Icges(6,5) * t171 - Icges(6,6) * t172) * t169;
t185 = t79 + t238;
t183 = (Icges(6,1) * t115 - t250 - t73) * t163 + (-Icges(6,1) * t113 - t107 - t71) * t162;
t181 = t162 * t218 + t156 + t233;
t51 = t115 * t120 + t116 * t121 + t119 * t243;
t41 = t51 * t152;
t10 = qJD(5) * t188 + t41;
t43 = Icges(6,5) * t85 + Icges(6,6) * t84 + Icges(6,3) * t222;
t45 = Icges(6,4) * t85 + Icges(6,2) * t84 + Icges(6,6) * t222;
t47 = Icges(6,1) * t85 + Icges(6,4) * t84 + Icges(6,5) * t222;
t13 = -t170 * t43 + (-t171 * t45 + t172 * t47 + (t171 * t75 - t172 * t71) * qJD(5)) * t169;
t42 = Icges(6,5) * t83 + Icges(6,6) * t82 - Icges(6,3) * t223;
t44 = Icges(6,4) * t83 + Icges(6,2) * t82 - Icges(6,6) * t223;
t46 = Icges(6,1) * t83 + Icges(6,4) * t82 - Icges(6,5) * t223;
t14 = -t170 * t42 + (-t171 * t44 + t172 * t46 + (-t171 * t76 - t172 * t73) * qJD(5)) * t169;
t124 = qJD(5) * t133;
t125 = qJD(5) * t134;
t126 = qJD(5) * t135;
t19 = t115 * t125 + t116 * t126 + t120 * t82 + t121 * t83 + (-t119 * t246 + t124 * t163) * t169;
t20 = -t113 * t125 + t114 * t126 + t120 * t84 + t121 * t85 + (t119 * t244 + t124 * t162) * t169;
t31 = -t170 * t70 + (-t171 * t73 + t172 * t76) * t169;
t40 = -t124 * t170 + (-t125 * t171 + t126 * t172 + (-t120 * t172 - t121 * t171) * qJD(5)) * t169;
t37 = t40 * t152;
t9 = qJD(5) * t189 + t253;
t179 = (t41 + ((t24 + t27 - t278) * t163 + t271 * t162) * t230) * t215 + t37 + (t14 + t19) * t212 + (t31 + t51) * t199 + (t30 + t50) * t198 + (t13 + t20 + t10) * t214 + (-t253 + ((-t25 + t271 - t277) * t163 - t252) * t230 + t9) * t213;
t178 = ((t168 * t24 - t113 * t44 + t114 * t46 + t73 * t84 + t76 * t85 + (t162 * t42 + t244 * t70) * t169) * t163 + (-t168 * t25 - t113 * t45 + t114 * t47 + t71 * t84 - t75 * t85 + (t162 * t43 + t244 * t68) * t169) * t162) * t169;
t177 = ((t168 * t26 + t115 * t44 + t116 * t46 + t73 * t82 + t76 * t83 + (t163 * t42 - t246 * t70) * t169) * t163 + (-t168 * t27 + t115 * t45 + t116 * t47 + t71 * t82 - t75 * t83 + (t163 * t43 - t246 * t68) * t169) * t162) * t169;
t176 = ((t168 * t30 + t14) * t163 + (-t168 * t31 + t13) * t162) * t169;
t62 = (-t101 - t129) * t168 + t210;
t63 = t209 + t269;
t175 = (t62 * t218 * t163 + (t62 * (-rSges(5,3) - qJ(4)) + t63 * t218) * t162) * t168;
t161 = pkin(2) * t165;
t138 = rSges(5,2) * t222;
t104 = t132 * t168 + t225;
t95 = -t112 * t168 - t228;
t94 = -t168 * t247 - t229;
t93 = rSges(6,1) * t115 - rSges(6,2) * t116;
t92 = -rSges(6,1) * t113 - rSges(6,2) * t114;
t53 = (t268 * t168 + t138 - t98) * t168 + t211;
t52 = t168 * (-t168 * t227 + t235) + t200;
t48 = -rSges(6,3) * t223 + t256;
t36 = t202 * t230 + qJD(1);
t21 = t152 * t48 - t163 * t279 + t180 * t246 + t200;
t16 = ((-t168 * t79 + t49) * t163 + (-t168 * t77 - t48) * t162) * t230;
t1 = [m(6) * t16; m(4) * (t95 * (-t130 - t261) + t94 * (t132 + t161) + (-t112 - t225 + t104) * t103) + t179 + (t22 * (t194 - t261) + t34 * (t195 - t225) + t21 * (t161 + t185) + (t34 + t264 + t267) * t35 + t272) * m(6) + (t53 * (t181 - t261) + t62 * (t138 - t209) + t52 * (t161 + t184) + t175 + (t62 + t267 + t270) * t63) * m(5); t179 + (t21 * t185 + t22 * t194 + (t264 + t266) * t35 + (t195 + t265) * t34 + t272) * m(6) + (t53 * t181 + t52 * t184 + t175 + (t266 + t270) * t63 + (t138 + t269) * t62) * m(5) + (-t103 * t112 - t104 * t247 - t130 * t95 + t132 * t94 - (-t103 * t132 - t104 * t130) * t168) * m(4); 0.2e1 * (t251 / 0.2e1 + t21 * t262) * m(6) + 0.2e1 * (t53 * t162 / 0.2e1 + t52 * t262) * m(5); -t10 * t223 / 0.2e1 + (-t51 * t170 + t188) * t199 + (-t19 * t170 + t177) * t212 + (qJD(5) * t178 + t152 * t20) * t245 / 0.2e1 + (-t50 * t170 + t189) * t198 + (-t170 * t20 + t178) * t214 - t170 * (qJD(5) * t176 + t37) / 0.2e1 + t152 * (-t40 * t170 + t176) / 0.2e1 + ((t115 * t236 + t116 * t237 + t133 * t243) * t152 + (t263 * t115 + t183 * t116 + t163 * t191) * t230) * t213 + ((-t113 * t236 + t114 * t237 + t133 * t245) * t152 + (-t113 * t263 + t114 * t183 + t162 * t191) * t230) * t215 - t152 * (-t170 * t133 * t152 + ((-t171 * t236 + t172 * t237) * t152 + ((-t171 * t263 + t172 * t183) * t169 - t201 * t170) * qJD(5)) * t169) / 0.2e1 + (qJD(5) * t177 + t152 * t19 + t168 * t9) * t243 / 0.2e1 + ((-t21 * t79 + t22 * t77 + t34 * t49 - t35 * t48) * t170 + (t16 * t202 + t36 * (-t162 * t48 + t163 * t49 - t244 * t79 - t246 * t77) + t203 * t127 + ((t168 * t34 - t21) * t163 + (t168 * t35 + t22) * t162) * t123) * t169 - (-t34 * t92 + t35 * t93) * t152 - (t36 * (-t162 * t93 + t163 * t92) + t203 * t190) * t230) * m(6);];
tauc = t1(:);
