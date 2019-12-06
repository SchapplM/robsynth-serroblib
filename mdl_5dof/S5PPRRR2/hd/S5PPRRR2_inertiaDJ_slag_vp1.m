% Calculate time derivative of joint inertia matrix for
% S5PPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRR2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR2_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR2_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRR2_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRR2_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:14:11
% EndTime: 2019-12-05 15:14:22
% DurationCPUTime: 5.18s
% Computational Cost: add. (17378->460), mult. (18743->730), div. (0->0), fcn. (17954->8), ass. (0->248)
t189 = sin(pkin(8));
t190 = cos(pkin(8));
t299 = t190 ^ 2;
t306 = t189 ^ 2 + t299;
t186 = pkin(9) + qJ(3);
t182 = sin(t186);
t183 = cos(t186);
t305 = qJD(3) * (t182 * rSges(4,1) + t183 * rSges(4,2));
t192 = cos(qJ(4));
t294 = t192 * pkin(4);
t146 = -pkin(7) * t183 + t294 * t182;
t302 = 2 * m(5);
t301 = 2 * m(6);
t298 = -t183 / 0.2e1;
t297 = t189 / 0.2e1;
t296 = -t190 / 0.2e1;
t295 = t190 / 0.2e1;
t188 = qJ(4) + qJ(5);
t184 = sin(t188);
t271 = t189 * t184;
t245 = t183 * t271;
t185 = cos(t188);
t266 = t190 * t185;
t169 = -t245 - t266;
t270 = t189 * t185;
t244 = t183 * t270;
t267 = t190 * t184;
t170 = t244 - t267;
t278 = t182 * t189;
t110 = Icges(6,5) * t170 + Icges(6,6) * t169 + Icges(6,3) * t278;
t112 = Icges(6,4) * t170 + Icges(6,2) * t169 + Icges(6,6) * t278;
t114 = Icges(6,1) * t170 + Icges(6,4) * t169 + Icges(6,5) * t278;
t243 = t183 * t267;
t171 = -t243 + t270;
t242 = t183 * t266;
t172 = t242 + t271;
t277 = t182 * t190;
t49 = t110 * t277 + t171 * t112 + t172 * t114;
t292 = t189 * t49;
t264 = t190 * t192;
t191 = sin(qJ(4));
t269 = t189 * t191;
t173 = -t183 * t269 - t264;
t265 = t190 * t191;
t268 = t189 * t192;
t174 = t183 * t268 - t265;
t125 = Icges(5,5) * t174 + Icges(5,6) * t173 + Icges(5,3) * t278;
t127 = Icges(5,4) * t174 + Icges(5,2) * t173 + Icges(5,6) * t278;
t129 = Icges(5,1) * t174 + Icges(5,4) * t173 + Icges(5,5) * t278;
t175 = -t183 * t265 + t268;
t176 = t183 * t264 + t269;
t58 = t125 * t277 + t175 * t127 + t176 * t129;
t291 = t189 * t58;
t111 = Icges(6,5) * t172 + Icges(6,6) * t171 + Icges(6,3) * t277;
t113 = Icges(6,4) * t172 + Icges(6,2) * t171 + Icges(6,6) * t277;
t115 = Icges(6,1) * t172 + Icges(6,4) * t171 + Icges(6,5) * t277;
t48 = t111 * t278 + t169 * t113 + t170 * t115;
t290 = t190 * t48;
t126 = Icges(5,5) * t176 + Icges(5,6) * t175 + Icges(5,3) * t277;
t128 = Icges(5,4) * t176 + Icges(5,2) * t175 + Icges(5,6) * t277;
t130 = Icges(5,1) * t176 + Icges(5,4) * t175 + Icges(5,5) * t277;
t57 = t126 * t278 + t173 * t128 + t174 * t130;
t289 = t190 * t57;
t250 = qJD(4) * t191;
t248 = pkin(4) * t250;
t194 = -t146 * qJD(3) - t183 * t248;
t249 = qJD(4) * t192;
t247 = pkin(4) * t249;
t187 = qJD(4) + qJD(5);
t252 = qJD(3) * t190;
t203 = t182 * t252 - t187 * t189;
t138 = t203 * t184 - t187 * t242;
t139 = -t203 * t185 - t187 * t243;
t236 = t183 * t252;
t83 = t139 * rSges(6,1) + t138 * rSges(6,2) + rSges(6,3) * t236;
t288 = t189 * t247 + t194 * t190 + t83;
t116 = t170 * rSges(6,1) + t169 * rSges(6,2) + rSges(6,3) * t278;
t253 = qJD(3) * t189;
t202 = t182 * t253 + t187 * t190;
t136 = t202 * t184 - t187 * t244;
t137 = -t202 * t185 - t187 * t245;
t237 = t183 * t253;
t82 = t137 * rSges(6,1) + t136 * rSges(6,2) + rSges(6,3) * t237;
t287 = t116 * t236 + t82 * t277;
t284 = Icges(5,4) * t191;
t283 = Icges(5,4) * t192;
t282 = Icges(6,4) * t184;
t281 = Icges(6,4) * t185;
t279 = t182 * t187;
t219 = Icges(6,5) * t185 - Icges(6,6) * t184;
t276 = t183 * ((-Icges(6,5) * t184 - Icges(6,6) * t185) * t279 + (Icges(6,3) * t182 + t219 * t183) * qJD(3));
t220 = Icges(5,5) * t192 - Icges(5,6) * t191;
t251 = qJD(4) * t182;
t275 = t183 * ((-Icges(5,5) * t191 - Icges(5,6) * t192) * t251 + (Icges(5,3) * t182 + t220 * t183) * qJD(3));
t147 = -Icges(6,3) * t183 + t219 * t182;
t274 = t183 * t147;
t151 = -Icges(5,3) * t183 + t220 * t182;
t273 = t183 * t151;
t272 = t183 * t189;
t229 = rSges(6,1) * t185 - rSges(6,2) * t184;
t109 = (-rSges(6,1) * t184 - rSges(6,2) * t185) * t279 + (rSges(6,3) * t182 + t229 * t183) * qJD(3);
t195 = pkin(7) * t182 + t294 * t183;
t133 = t195 * qJD(3) - t182 * t248;
t263 = -t109 - t133;
t122 = -pkin(4) * t265 + t195 * t189;
t262 = t116 + t122;
t117 = t172 * rSges(6,1) + t171 * rSges(6,2) + rSges(6,3) * t277;
t123 = pkin(4) * t269 + t195 * t190;
t261 = t117 + t123;
t230 = rSges(5,1) * t192 - rSges(5,2) * t191;
t124 = (-rSges(5,1) * t191 - rSges(5,2) * t192) * t251 + (rSges(5,3) * t182 + t230 * t183) * qJD(3);
t232 = pkin(3) * t183 + pkin(6) * t182;
t178 = t232 * qJD(3);
t260 = -t124 - t178;
t150 = -t183 * rSges(6,3) + t229 * t182;
t80 = t183 * t116 + t150 * t278;
t259 = -t146 - t150;
t180 = t182 * pkin(3) - t183 * pkin(6);
t258 = t306 * qJD(3) * t180;
t156 = -t183 * rSges(5,3) + t230 * t182;
t257 = -t156 - t180;
t256 = t306 * t232;
t255 = qJD(3) * t182;
t254 = qJD(3) * t183;
t246 = t109 * t278 + t150 * t237 + t183 * t82;
t241 = -t178 + t263;
t240 = -t180 + t259;
t239 = t191 * t255;
t238 = t192 * t255;
t234 = t261 * t183;
t233 = t259 * t190;
t218 = -t112 * t184 + t114 * t185;
t54 = -t183 * t110 + t218 * t182;
t217 = -t113 * t184 + t115 * t185;
t55 = -t183 * t111 + t217 * t182;
t228 = t54 * t189 + t55 * t190;
t216 = -t127 * t191 + t129 * t192;
t60 = -t183 * t125 + t216 * t182;
t215 = -t128 * t191 + t130 * t192;
t61 = -t183 * t126 + t215 * t182;
t227 = t60 * t189 + t61 * t190;
t225 = Icges(5,1) * t192 - t284;
t224 = Icges(6,1) * t185 - t282;
t222 = -Icges(5,2) * t191 + t283;
t221 = -Icges(6,2) * t184 + t281;
t131 = t174 * rSges(5,1) + t173 * rSges(5,2) + rSges(5,3) * t278;
t132 = t176 * rSges(5,1) + t175 * rSges(5,2) + rSges(5,3) * t277;
t214 = t131 * t190 - t132 * t189;
t148 = -Icges(6,6) * t183 + t221 * t182;
t149 = -Icges(6,5) * t183 + t224 * t182;
t213 = t148 * t184 - t149 * t185;
t152 = -Icges(5,6) * t183 + t222 * t182;
t153 = -Icges(5,5) * t183 + t225 * t182;
t212 = t152 * t191 - t153 * t192;
t74 = Icges(6,5) * t137 + Icges(6,6) * t136 + Icges(6,3) * t237;
t208 = t110 * t254 + t182 * t74;
t75 = Icges(6,5) * t139 + Icges(6,6) * t138 + Icges(6,3) * t236;
t207 = t111 * t254 + t182 * t75;
t142 = -t174 * qJD(4) + t189 * t239;
t143 = t173 * qJD(4) - t189 * t238;
t89 = Icges(5,5) * t143 + Icges(5,6) * t142 + Icges(5,3) * t237;
t206 = t125 * t254 + t182 * t89;
t144 = -t176 * qJD(4) + t190 * t239;
t145 = t175 * qJD(4) - t190 * t238;
t90 = Icges(5,5) * t145 + Icges(5,6) * t144 + Icges(5,3) * t236;
t205 = t126 * t254 + t182 * t90;
t106 = (-Icges(6,2) * t185 - t282) * t279 + (Icges(6,6) * t182 + t221 * t183) * qJD(3);
t107 = (-Icges(6,1) * t184 - t281) * t279 + (Icges(6,5) * t182 + t224 * t183) * qJD(3);
t76 = Icges(6,4) * t137 + Icges(6,2) * t136 + Icges(6,6) * t237;
t78 = Icges(6,1) * t137 + Icges(6,4) * t136 + Icges(6,5) * t237;
t17 = (t218 * qJD(3) - t74) * t183 + (qJD(3) * t110 + (-t112 * t187 + t78) * t185 + (-t114 * t187 - t76) * t184) * t182;
t77 = Icges(6,4) * t139 + Icges(6,2) * t138 + Icges(6,6) * t236;
t79 = Icges(6,1) * t139 + Icges(6,4) * t138 + Icges(6,5) * t236;
t18 = (t217 * qJD(3) - t75) * t183 + (qJD(3) * t111 + (-t113 * t187 + t79) * t185 + (-t115 * t187 - t77) * t184) * t182;
t47 = t110 * t278 + t169 * t112 + t170 * t114;
t50 = t111 * t277 + t171 * t113 + t172 * t115;
t22 = t136 * t112 + t137 * t114 + t169 * t76 + t170 * t78 + t208 * t189;
t23 = t136 * t113 + t137 * t115 + t169 * t77 + t170 * t79 + t207 * t189;
t63 = t147 * t278 + t169 * t148 + t170 * t149;
t6 = -(t169 * t106 + t170 * t107 + t136 * t148 + t137 * t149) * t183 + (t23 * t190 + (t22 - t276) * t189) * t182 + (t63 * t182 + (t290 + (t47 - t274) * t189) * t183) * qJD(3);
t64 = t147 * t277 + t171 * t148 + t172 * t149;
t24 = t138 * t112 + t139 * t114 + t171 * t76 + t172 * t78 + t208 * t190;
t25 = t138 * t113 + t139 * t115 + t171 * t77 + t172 * t79 + t207 * t190;
t7 = -(t171 * t106 + t172 * t107 + t138 * t148 + t139 * t149) * t183 + (t24 * t189 + (t25 - t276) * t190) * t182 + (t64 * t182 + (t292 + (t50 - t274) * t190) * t183) * qJD(3);
t71 = -t213 * t182 - t274;
t201 = -t183 * ((t276 + (t213 * t183 + t228) * qJD(3)) * t183 + (t18 * t190 + t17 * t189 - (qJD(3) * t147 + (-t148 * t187 + t107) * t185 + (-t149 * t187 - t106) * t184) * t183 + t71 * qJD(3)) * t182) + t7 * t277 + t6 * t278 + (-t63 * t183 + (t189 * t47 + t290) * t182) * t237 + (-t64 * t183 + (t190 * t50 + t292) * t182) * t236 + (t228 * t182 - t71 * t183) * t255;
t13 = t23 * t189 - t22 * t190;
t14 = t25 * t189 - t24 * t190;
t200 = t6 * t296 + t7 * t297 + (-t17 * t190 + t18 * t189) * t298 + t13 * t278 / 0.2e1 + t14 * t277 / 0.2e1 + (t55 * t189 - t54 * t190) * t255 / 0.2e1 + (t189 * (t48 * t189 - t47 * t190) + t190 * (t50 * t189 - t49 * t190)) * t254 / 0.2e1;
t196 = qJD(3) * (-Icges(4,5) * t182 - Icges(4,6) * t183);
t163 = t189 * t196;
t135 = t257 * t190;
t134 = t257 * t189;
t121 = (-Icges(5,1) * t191 - t283) * t251 + (Icges(5,5) * t182 + t225 * t183) * qJD(3);
t120 = (-Icges(5,2) * t192 - t284) * t251 + (Icges(5,6) * t182 + t222 * t183) * qJD(3);
t118 = t306 * t305;
t104 = t117 * t255;
t103 = t116 * t277;
t99 = t194 * t189 - t190 * t247;
t98 = t260 * t190;
t97 = t260 * t189;
t96 = t145 * rSges(5,1) + t144 * rSges(5,2) + rSges(5,3) * t236;
t95 = t143 * rSges(5,1) + t142 * rSges(5,2) + rSges(5,3) * t237;
t94 = Icges(5,1) * t145 + Icges(5,4) * t144 + Icges(5,5) * t236;
t93 = Icges(5,1) * t143 + Icges(5,4) * t142 + Icges(5,5) * t237;
t92 = Icges(5,4) * t145 + Icges(5,2) * t144 + Icges(5,6) * t236;
t91 = Icges(5,4) * t143 + Icges(5,2) * t142 + Icges(5,6) * t237;
t88 = t240 * t190;
t87 = t240 * t189;
t86 = -t183 * t132 - t156 * t277;
t85 = t183 * t131 + t156 * t278;
t84 = -t212 * t182 - t273;
t81 = -t183 * t117 - t150 * t277;
t70 = t214 * t182;
t69 = t151 * t277 + t175 * t152 + t176 * t153;
t68 = t151 * t278 + t173 * t152 + t174 * t153;
t67 = -t117 * t278 + t103;
t66 = t241 * t190;
t65 = t241 * t189;
t62 = t189 * t131 + t190 * t132 + t256;
t59 = t126 * t277 + t175 * t128 + t176 * t130;
t56 = t125 * t278 + t173 * t127 + t174 * t129;
t53 = t182 * t233 - t234;
t52 = t183 * t122 + t146 * t278 + t80;
t51 = t189 * t95 + t190 * t96 - t258;
t46 = -t124 * t277 - t183 * t96 + (-t156 * t183 * t190 + t132 * t182) * qJD(3);
t45 = t124 * t278 + t183 * t95 + (-t131 * t182 + t156 * t272) * qJD(3);
t44 = t103 + (t122 * t190 - t261 * t189) * t182;
t43 = t262 * t189 + t261 * t190 + t256;
t42 = -t109 * t277 + t104 + (-t150 * t252 - t83) * t183;
t41 = -t116 * t255 + t246;
t40 = (-t189 * t96 + t190 * t95) * t182 + t214 * t254;
t39 = (-t117 * t254 - t182 * t83) * t189 + t287;
t38 = t288 * t190 + (t82 + t99) * t189 - t258;
t34 = t104 + (qJD(3) * t123 + t263 * t190) * t182 + (qJD(3) * t233 - t288) * t183;
t33 = t133 * t278 + t183 * t99 + (t146 * t272 - t262 * t182) * qJD(3) + t246;
t32 = t144 * t128 + t145 * t130 + t175 * t92 + t176 * t94 + t205 * t190;
t31 = t144 * t127 + t145 * t129 + t175 * t91 + t176 * t93 + t206 * t190;
t30 = t142 * t128 + t143 * t130 + t173 * t92 + t174 * t94 + t205 * t189;
t29 = t142 * t127 + t143 * t129 + t173 * t91 + t174 * t93 + t206 * t189;
t27 = (t215 * qJD(3) - t90) * t183 + (qJD(3) * t126 - t191 * t92 + t192 * t94 + (-t128 * t192 - t130 * t191) * qJD(4)) * t182;
t26 = (t216 * qJD(3) - t89) * t183 + (qJD(3) * t125 - t191 * t91 + t192 * t93 + (-t127 * t192 - t129 * t191) * qJD(4)) * t182;
t19 = (t122 * t254 + t182 * t99) * t190 + (-qJD(3) * t234 - t288 * t182) * t189 + t287;
t16 = t32 * t189 - t31 * t190;
t15 = t30 * t189 - t29 * t190;
t9 = -(t175 * t120 + t176 * t121 + t144 * t152 + t145 * t153) * t183 + (t31 * t189 + (t32 - t275) * t190) * t182 + (t69 * t182 + (t291 + (t59 - t273) * t190) * t183) * qJD(3);
t8 = -(t173 * t120 + t174 * t121 + t142 * t152 + t143 * t153) * t183 + (t30 * t190 + (t29 - t275) * t189) * t182 + (t68 * t182 + (t289 + (t56 - t273) * t189) * t183) * qJD(3);
t1 = [0; 0; 0; -m(4) * t118 + m(5) * t51 + m(6) * t38; m(5) * (t98 * t189 - t97 * t190) + m(6) * (t66 * t189 - t65 * t190); (t38 * t43 + t65 * t87 + t66 * t88) * t301 + (t134 * t97 + t135 * t98 + t62 * t51) * t302 + 0.2e1 * m(4) * (-t118 + t305) * t306 * (rSges(4,1) * t183 - rSges(4,2) * t182) + (-t299 * t163 - t13 - t15) * t190 + (t14 + t16 + (-t189 * t163 + t306 * t196) * t190) * t189; m(5) * t40 + m(6) * t19; m(5) * (t45 * t189 - t46 * t190) + m(6) * (t33 * t189 - t34 * t190); (t27 * t189 - t26 * t190) * t298 + t9 * t297 + t8 * t296 + (t15 * t297 + t16 * t295) * t182 + m(6) * (t19 * t43 + t33 * t88 + t34 * t87 + t38 * t44 + t52 * t66 + t53 * t65) + m(5) * (t46 * t134 + t45 * t135 + t40 * t62 + t70 * t51 + t85 * t98 + t86 * t97) + (t182 * (t61 * t189 - t60 * t190) / 0.2e1 + ((t57 * t189 - t56 * t190) * t297 + (t59 * t189 - t58 * t190) * t295) * t183) * qJD(3) + t200; (t19 * t44 + t33 * t52 + t34 * t53) * t301 + (t227 * t182 - t84 * t183) * t255 - t183 * ((t275 + (t212 * t183 + t227) * qJD(3)) * t183 + (t27 * t190 + t26 * t189 - (qJD(3) * t151 - t120 * t191 + t121 * t192 - t152 * t249 - t153 * t250) * t183 + t84 * qJD(3)) * t182) + (-t68 * t183 + (t189 * t56 + t289) * t182) * t237 + t8 * t278 + (t40 * t70 + t45 * t85 + t46 * t86) * t302 + (-t69 * t183 + (t190 * t59 + t291) * t182) * t236 + t9 * t277 + t201; m(6) * t39; m(6) * (t41 * t189 - t42 * t190); m(6) * (t38 * t67 + t39 * t43 + t41 * t88 + t42 * t87 + t65 * t81 + t66 * t80) + t200; m(6) * (t19 * t67 + t33 * t80 + t34 * t81 + t39 * t44 + t41 * t52 + t42 * t53) + t201; (t39 * t67 + t41 * t80 + t42 * t81) * t301 + t201;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
