% Calculate vector of inverse dynamics joint torques for
% S5PPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPPRR2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR2_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPPRR2_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR2_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPPRR2_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPPRR2_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:26
% EndTime: 2019-12-05 14:59:36
% DurationCPUTime: 7.59s
% Computational Cost: add. (13098->630), mult. (35900->906), div. (0->0), fcn. (42924->10), ass. (0->217)
t218 = sin(qJ(5));
t219 = cos(qJ(5));
t293 = -rSges(6,1) * t219 + rSges(6,2) * t218;
t292 = 2 * qJD(4);
t291 = 2 * qJDD(4);
t215 = cos(pkin(9));
t216 = cos(pkin(8));
t213 = sin(pkin(8));
t279 = cos(qJ(4));
t253 = t213 * t279;
t278 = sin(qJ(4));
t197 = t215 * t253 - t216 * t278;
t212 = sin(pkin(9));
t217 = cos(pkin(7));
t214 = sin(pkin(7));
t262 = t214 * t216;
t193 = -t212 * t217 + t215 * t262;
t252 = t213 * t278;
t179 = t193 * t279 + t214 * t252;
t168 = t179 * qJD(4);
t192 = t212 * t262 + t215 * t217;
t233 = -t193 * t278 + t214 * t253;
t95 = qJD(5) * t168 + qJDD(4) * t192 - qJDD(5) * t233;
t290 = t95 / 0.2e1;
t261 = t216 * t217;
t195 = t212 * t214 + t215 * t261;
t181 = t195 * t279 + t217 * t252;
t170 = t181 * qJD(4);
t194 = t212 * t261 - t214 * t215;
t232 = -t195 * t278 + t217 * t253;
t96 = qJD(5) * t170 + qJDD(4) * t194 - qJDD(5) * t232;
t289 = t96 / 0.2e1;
t147 = qJD(4) * t192 - qJD(5) * t233;
t288 = -t147 / 0.2e1;
t287 = t147 / 0.2e1;
t148 = qJD(4) * t194 - qJD(5) * t232;
t286 = -t148 / 0.2e1;
t285 = t148 / 0.2e1;
t188 = t197 * qJD(4);
t196 = t215 * t252 + t216 * t279;
t263 = t212 * t213;
t149 = qJD(5) * t188 + qJDD(4) * t263 + qJDD(5) * t196;
t284 = t149 / 0.2e1;
t184 = qJD(4) * t263 + qJD(5) * t196;
t283 = -t184 / 0.2e1;
t282 = t184 / 0.2e1;
t167 = t233 * qJD(4);
t117 = pkin(4) * t167 + pkin(6) * t168;
t144 = t179 * t219 + t192 * t218;
t91 = -qJD(5) * t144 - t167 * t218;
t143 = -t179 * t218 + t192 * t219;
t92 = qJD(5) * t143 + t167 * t219;
t55 = rSges(6,1) * t92 + rSges(6,2) * t91 + rSges(6,3) * t168;
t274 = t117 + t55;
t169 = t232 * qJD(4);
t118 = pkin(4) * t169 + pkin(6) * t170;
t146 = t181 * t219 + t194 * t218;
t93 = -qJD(5) * t146 - t169 * t218;
t145 = -t181 * t218 + t194 * t219;
t94 = qJD(5) * t145 + t169 * t219;
t56 = rSges(6,1) * t94 + rSges(6,2) * t93 + rSges(6,3) * t170;
t273 = t118 + t56;
t128 = pkin(4) * t179 - pkin(6) * t233;
t69 = rSges(6,1) * t144 + rSges(6,2) * t143 - rSges(6,3) * t233;
t272 = t128 + t69;
t130 = pkin(4) * t181 - pkin(6) * t232;
t70 = rSges(6,1) * t146 + rSges(6,2) * t145 - rSges(6,3) * t232;
t271 = t130 + t70;
t187 = t196 * qJD(4);
t158 = -pkin(4) * t187 + pkin(6) * t188;
t183 = t197 * t219 + t218 * t263;
t141 = -qJD(5) * t183 + t187 * t218;
t182 = -t197 * t218 + t219 * t263;
t142 = qJD(5) * t182 - t187 * t219;
t74 = rSges(6,1) * t142 + rSges(6,2) * t141 + rSges(6,3) * t188;
t270 = t158 + t74;
t269 = Icges(5,4) * t179;
t268 = Icges(5,4) * t181;
t267 = Icges(5,4) * t197;
t266 = Icges(6,4) * t144;
t265 = Icges(6,4) * t146;
t264 = Icges(6,4) * t183;
t108 = rSges(6,1) * t183 + rSges(6,2) * t182 + rSges(6,3) * t196;
t164 = pkin(4) * t197 + pkin(6) * t196;
t260 = t108 + t164;
t258 = qJD(3) * t213;
t259 = qJD(2) * t214 + t217 * t258;
t256 = qJDD(3) * t213;
t198 = qJDD(2) * t214 + t217 * t256;
t257 = -m(4) - m(5) - m(6);
t254 = -m(3) + t257;
t127 = pkin(4) * t233 + pkin(6) * t179;
t129 = pkin(4) * t232 + pkin(6) * t181;
t163 = -t196 * pkin(4) + pkin(6) * t197;
t249 = -qJD(2) * t217 + t214 * t258;
t248 = -qJD(3) * t216 + qJD(1);
t199 = -qJDD(2) * t217 + t214 * t256;
t208 = -qJDD(3) * t216 + qJDD(1);
t246 = -Icges(6,1) * t219 + Icges(6,4) * t218;
t245 = -Icges(6,4) * t219 + Icges(6,2) * t218;
t244 = -Icges(6,5) * t219 + Icges(6,6) * t218;
t103 = rSges(5,1) * t179 + rSges(5,2) * t233 + rSges(5,3) * t192;
t104 = rSges(5,1) * t181 + rSges(5,2) * t232 + rSges(5,3) * t194;
t243 = t103 * t194 - t104 * t192;
t115 = rSges(5,1) * t167 - rSges(5,2) * t168;
t116 = rSges(5,1) * t169 - rSges(5,2) * t170;
t242 = t115 * t194 - t116 * t192;
t241 = t128 * t194 - t130 * t192;
t89 = t179 * rSges(6,3) - t293 * t233;
t90 = t181 * rSges(6,3) - t232 * t293;
t138 = t197 * rSges(6,3) + t196 * t293;
t153 = rSges(5,1) * t197 - rSges(5,2) * t196 + rSges(5,3) * t263;
t240 = -t103 * t263 + t153 * t192;
t239 = t104 * t263 - t153 * t194;
t157 = -rSges(5,1) * t187 - rSges(5,2) * t188;
t238 = -t115 * t263 + t157 * t192;
t237 = t116 * t263 - t157 * t194;
t236 = -t128 * t263 + t164 * t192;
t235 = t130 * t263 - t164 * t194;
t234 = (Icges(6,5) * t182 - Icges(6,6) * t183) * t184 + t147 * (Icges(6,5) * t143 - Icges(6,6) * t144) + t148 * (Icges(6,5) * t145 - Icges(6,6) * t146);
t231 = (Icges(5,5) * t233 - Icges(5,6) * t179) * t192 + (Icges(5,5) * t232 - Icges(5,6) * t181) * t194 + (-Icges(5,5) * t196 - Icges(5,6) * t197) * t263;
t106 = Icges(6,2) * t182 + Icges(6,6) * t196 + t264;
t65 = Icges(6,2) * t143 - Icges(6,6) * t233 + t266;
t66 = Icges(6,2) * t145 - Icges(6,6) * t232 + t265;
t224 = (Icges(6,1) * t145 - t265 - t66) * t148 + (Icges(6,1) * t143 - t266 - t65) * t147 + (Icges(6,1) * t182 - t106 - t264) * t184;
t177 = Icges(6,4) * t182;
t107 = Icges(6,1) * t183 + Icges(6,5) * t196 + t177;
t139 = Icges(6,4) * t143;
t140 = Icges(6,4) * t145;
t67 = Icges(6,1) * t144 - Icges(6,5) * t233 + t139;
t68 = Icges(6,1) * t146 - Icges(6,5) * t232 + t140;
t223 = (-Icges(6,2) * t146 + t140 + t68) * t148 + (-Icges(6,2) * t144 + t139 + t67) * t147 + (-Icges(6,2) * t183 + t107 + t177) * t184;
t100 = Icges(5,2) * t232 + Icges(5,6) * t194 + t268;
t151 = -Icges(5,2) * t196 + Icges(5,6) * t263 + t267;
t99 = Icges(5,2) * t233 + Icges(5,6) * t192 + t269;
t222 = (Icges(5,1) * t232 - t100 - t268) * t194 + (Icges(5,1) * t233 - t269 - t99) * t192 + (-Icges(5,1) * t196 - t151 - t267) * t263;
t171 = Icges(5,4) * t233;
t101 = Icges(5,1) * t179 + Icges(5,5) * t192 + t171;
t172 = Icges(5,4) * t232;
t102 = Icges(5,1) * t181 + Icges(5,5) * t194 + t172;
t189 = Icges(5,4) * t196;
t152 = Icges(5,1) * t197 + Icges(5,5) * t263 - t189;
t221 = (Icges(5,2) * t181 - t102 - t172) * t194 + (Icges(5,2) * t179 - t101 - t171) * t192 + (Icges(5,2) * t197 - t152 + t189) * t263;
t220 = (Icges(6,3) * t181 + t218 * t66 - t219 * t68 - t232 * t244) * t148 + (Icges(6,3) * t179 + t218 * t65 - t219 * t67 - t233 * t244) * t147 + (Icges(6,3) * t197 + t106 * t218 - t107 * t219 + t196 * t244) * t184;
t162 = -rSges(5,1) * t196 - rSges(5,2) * t197;
t156 = -Icges(5,1) * t187 - Icges(5,4) * t188;
t155 = -Icges(5,4) * t187 - Icges(5,2) * t188;
t154 = -Icges(5,5) * t187 - Icges(5,6) * t188;
t150 = Icges(5,5) * t197 - Icges(5,6) * t196 + Icges(5,3) * t263;
t137 = Icges(6,5) * t197 + t196 * t246;
t136 = Icges(6,6) * t197 + t196 * t245;
t134 = rSges(6,1) * t182 - rSges(6,2) * t183;
t126 = rSges(5,1) * t232 - rSges(5,2) * t181;
t125 = rSges(5,1) * t233 - rSges(5,2) * t179;
t114 = Icges(5,1) * t169 - Icges(5,4) * t170;
t113 = Icges(5,1) * t167 - Icges(5,4) * t168;
t112 = Icges(5,4) * t169 - Icges(5,2) * t170;
t111 = Icges(5,4) * t167 - Icges(5,2) * t168;
t110 = Icges(5,5) * t169 - Icges(5,6) * t170;
t109 = Icges(5,5) * t167 - Icges(5,6) * t168;
t105 = Icges(6,5) * t183 + Icges(6,6) * t182 + Icges(6,3) * t196;
t98 = Icges(5,5) * t181 + Icges(5,6) * t232 + Icges(5,3) * t194;
t97 = Icges(5,5) * t179 + Icges(5,6) * t233 + Icges(5,3) * t192;
t88 = Icges(6,5) * t181 - t232 * t246;
t87 = Icges(6,5) * t179 - t233 * t246;
t86 = Icges(6,6) * t181 - t232 * t245;
t85 = Icges(6,6) * t179 - t233 * t245;
t82 = rSges(6,1) * t145 - rSges(6,2) * t146;
t81 = rSges(6,1) * t143 - rSges(6,2) * t144;
t73 = Icges(6,1) * t142 + Icges(6,4) * t141 + Icges(6,5) * t188;
t72 = Icges(6,4) * t142 + Icges(6,2) * t141 + Icges(6,6) * t188;
t71 = Icges(6,5) * t142 + Icges(6,6) * t141 + Icges(6,3) * t188;
t64 = Icges(6,5) * t146 + Icges(6,6) * t145 - Icges(6,3) * t232;
t63 = Icges(6,5) * t144 + Icges(6,6) * t143 - Icges(6,3) * t233;
t62 = qJD(4) * t239 + t249;
t61 = qJD(4) * t240 + t259;
t59 = qJD(4) * t243 + t248;
t54 = Icges(6,1) * t94 + Icges(6,4) * t93 + Icges(6,5) * t170;
t53 = Icges(6,1) * t92 + Icges(6,4) * t91 + Icges(6,5) * t168;
t52 = Icges(6,4) * t94 + Icges(6,2) * t93 + Icges(6,6) * t170;
t51 = Icges(6,4) * t92 + Icges(6,2) * t91 + Icges(6,6) * t168;
t50 = Icges(6,5) * t94 + Icges(6,6) * t93 + Icges(6,3) * t170;
t49 = Icges(6,5) * t92 + Icges(6,6) * t91 + Icges(6,3) * t168;
t45 = qJD(4) * t237 + qJDD(4) * t239 + t199;
t44 = qJD(4) * t238 + qJDD(4) * t240 + t198;
t43 = t105 * t196 + t106 * t182 + t107 * t183;
t36 = -t105 * t232 + t106 * t145 + t107 * t146;
t35 = -t105 * t233 + t106 * t143 + t107 * t144;
t34 = qJD(4) * t242 + qJDD(4) * t243 + t208;
t33 = qJD(4) * t235 - t108 * t148 + t184 * t70 + t249;
t32 = qJD(4) * t236 + t108 * t147 - t184 * t69 + t259;
t29 = t182 * t66 + t183 * t68 + t196 * t64;
t28 = t182 * t65 + t183 * t67 + t196 * t63;
t23 = qJD(4) * t241 - t147 * t70 + t148 * t69 + t248;
t22 = t145 * t66 + t146 * t68 - t232 * t64;
t21 = t145 * t65 + t146 * t67 - t232 * t63;
t20 = t143 * t66 + t144 * t68 - t233 * t64;
t19 = t143 * t65 + t144 * t67 - t233 * t63;
t18 = t105 * t188 + t106 * t141 + t107 * t142 + t182 * t72 + t183 * t73 + t196 * t71;
t17 = t105 * t170 + t106 * t93 + t107 * t94 + t145 * t72 + t146 * t73 - t232 * t71;
t16 = t105 * t168 + t106 * t91 + t107 * t92 + t143 * t72 + t144 * t73 - t233 * t71;
t15 = -t108 * t96 - t148 * t74 + t149 * t70 + t184 * t56 + t235 * qJDD(4) + (t118 * t263 - t158 * t194) * qJD(4) + t199;
t14 = t108 * t95 + t147 * t74 - t149 * t69 - t184 * t55 + t236 * qJDD(4) + (-t117 * t263 + t158 * t192) * qJD(4) + t198;
t13 = t141 * t66 + t142 * t68 + t182 * t52 + t183 * t54 + t188 * t64 + t196 * t50;
t12 = t141 * t65 + t142 * t67 + t182 * t51 + t183 * t53 + t188 * t63 + t196 * t49;
t11 = -t147 * t56 + t148 * t55 + t69 * t96 - t70 * t95 + t241 * qJDD(4) + (t117 * t194 - t118 * t192) * qJD(4) + t208;
t10 = t145 * t52 + t146 * t54 + t170 * t64 - t232 * t50 + t66 * t93 + t68 * t94;
t9 = t145 * t51 + t146 * t53 + t170 * t63 - t232 * t49 + t65 * t93 + t67 * t94;
t8 = t143 * t52 + t144 * t54 + t168 * t64 - t233 * t50 + t66 * t91 + t68 * t92;
t7 = t143 * t51 + t144 * t53 + t168 * t63 - t233 * t49 + t65 * t91 + t67 * t92;
t6 = t147 * t28 + t148 * t29 + t184 * t43;
t5 = t147 * t21 + t22 * t148 + t184 * t36;
t4 = t147 * t19 + t148 * t20 + t184 * t35;
t3 = t12 * t147 + t13 * t148 + t149 * t43 + t18 * t184 + t28 * t95 + t29 * t96;
t2 = t10 * t148 + t147 * t9 + t149 * t36 + t17 * t184 + t21 * t95 + t22 * t96;
t1 = t7 * t147 + t148 * t8 + t149 * t35 + t16 * t184 + t19 * t95 + t20 * t96;
t24 = [(m(2) + m(3)) * qJDD(1) + m(4) * t208 + m(5) * t34 + m(6) * t11 + (-m(2) + t254) * g(3); t254 * (g(1) * t214 - g(2) * t217) + m(4) * (t198 * t214 - t199 * t217) + m(5) * (t214 * t44 - t217 * t45) + m(6) * (t14 * t214 - t15 * t217) + m(3) * (t214 ^ 2 + t217 ^ 2) * qJDD(2); t257 * (-g(3) * t216 + (g(1) * t217 + g(2) * t214) * t213) + m(4) * (-t208 * t216 + (t198 * t217 + t199 * t214) * t213) + m(5) * (-t216 * t34 + (t214 * t45 + t217 * t44) * t213) + m(6) * (-t11 * t216 + (t14 * t217 + t15 * t214) * t213); (t12 * t192 + t13 * t194 + t18 * t263) * t282 + ((t182 * t86 + t183 * t88 + t197 * t64) * t148 + (t182 * t85 + t183 * t87 + t197 * t63) * t147 + (t197 * t105 + t182 * t136 + t183 * t137) * t184 + (t179 * t28 + t181 * t29 + t197 * t43) * qJD(5) + t220 * t196) * t283 + (t192 * t28 + t194 * t29 + t263 * t43) * t284 + (t10 * t194 + t17 * t263 + t192 * t9) * t285 + ((t145 * t86 + t146 * t88 + t181 * t64) * t148 + (t145 * t85 + t146 * t87 + t181 * t63) * t147 + (t181 * t105 + t145 * t136 + t146 * t137) * t184 + (t179 * t21 + t181 * t22 + t197 * t36) * qJD(5) - t220 * t232) * t286 + (t16 * t263 + t192 * t7 + t194 * t8) * t287 + ((t143 * t86 + t144 * t88 + t179 * t64) * t148 + (t143 * t85 + t144 * t87 + t179 * t63) * t147 + (t179 * t105 + t143 * t136 + t144 * t137) * t184 + (t179 * t19 + t181 * t20 + t197 * t35) * qJD(5) - t220 * t233) * t288 + (t192 * t21 + t194 * t22 + t263 * t36) * t289 + (t19 * t192 + t194 * t20 + t263 * t35) * t290 + ((t192 * (t101 * t179 + t192 * t97 + t233 * t99) + t194 * (t100 * t233 + t102 * t179 + t192 * t98) + t263 * (t150 * t192 + t151 * t233 + t152 * t179)) * t291 + (t192 * (t101 * t167 + t109 * t192 + t111 * t233 + t113 * t179 - t168 * t99) + t194 * (-t100 * t168 + t102 * t167 + t110 * t192 + t112 * t233 + t114 * t179) + t263 * (-t151 * t168 + t152 * t167 + t154 * t192 + t155 * t233 + t156 * t179)) * t292 + t1) * t192 / 0.2e1 + ((t192 * (t101 * t181 + t194 * t97 + t232 * t99) + t194 * (t100 * t232 + t102 * t181 + t194 * t98) + t263 * (t150 * t194 + t151 * t232 + t152 * t181)) * t291 + (t192 * (t101 * t169 + t109 * t194 + t111 * t232 + t113 * t181 - t170 * t99) + t194 * (-t100 * t170 + t102 * t169 + t110 * t194 + t112 * t232 + t114 * t181) + t263 * (-t151 * t170 + t152 * t169 + t154 * t194 + t155 * t232 + t156 * t181)) * t292 + t2) * t194 / 0.2e1 - (t179 * t4 + t181 * t5 + t197 * t6) * qJD(5) / 0.2e1 - (t194 * (t181 * t222 + t194 * t231 - t221 * t232) + t192 * (t179 * t222 + t192 * t231 - t221 * t233) + (t196 * t221 + t197 * t222 + t231 * t263) * t263) * (qJD(4) ^ 2) / 0.2e1 + ((t192 * (t101 * t197 - t196 * t99 + t263 * t97) + t194 * (-t100 * t196 + t102 * t197 + t263 * t98) + t263 * (t150 * t263 - t151 * t196 + t152 * t197)) * t291 + t3 + (t192 * (-t101 * t187 + t109 * t263 - t111 * t196 + t113 * t197 - t188 * t99) + t194 * (-t100 * t188 - t102 * t187 + t110 * t263 - t112 * t196 + t114 * t197) + t263 * (-t151 * t188 - t152 * t187 + t154 * t263 - t155 * t196 + t156 * t197)) * t292) * t263 / 0.2e1 + (-g(1) * (t129 + t90) - g(2) * (t127 + t89) - g(3) * (t138 + t163) - t32 * (t138 * t147 - t184 * t89) - t33 * (-t138 * t148 + t184 * t90) - t23 * (-t147 * t90 + t148 * t89) - (t32 * (t108 * t179 - t197 * t69) + t33 * (-t108 * t181 + t197 * t70) + t23 * (-t179 * t70 + t181 * t69)) * qJD(5) - (t32 * (-t127 * t263 + t163 * t192) + t33 * (t129 * t263 - t163 * t194) + t23 * (t127 * t194 - t129 * t192)) * qJD(4) + (-t14 * t272 + t15 * t271 + t273 * t33 - t274 * t32) * t263 + (t11 * t272 - t15 * t260 + t23 * t274 - t270 * t33) * t194 + (-t11 * t271 + t14 * t260 - t23 * t273 + t270 * t32) * t192) * m(6) + (t237 * t62 + t238 * t61 + t239 * t45 + t240 * t44 + t242 * t59 + t243 * t34 - (t61 * (-t125 * t263 + t162 * t192) + t62 * (t126 * t263 - t162 * t194) + t59 * (t125 * t194 - t126 * t192)) * qJD(4) - g(1) * t126 - g(2) * t125 - g(3) * t162) * m(5); t170 * t5 / 0.2e1 - t232 * t2 / 0.2e1 + (t196 * t36 - t21 * t233 - t22 * t232) * t289 + (-t10 * t232 + t168 * t21 + t17 * t196 + t170 * t22 + t188 * t36 - t233 * t9) * t285 + t168 * t4 / 0.2e1 - t233 * t1 / 0.2e1 + (-t19 * t233 + t196 * t35 - t20 * t232) * t290 + (t16 * t196 + t19 * t168 + t20 * t170 + t35 * t188 - t232 * t8 - t233 * t7) * t287 + t188 * t6 / 0.2e1 + t196 * t3 / 0.2e1 + (t196 * t43 - t232 * t29 - t233 * t28) * t284 + (-t12 * t233 - t13 * t232 + t28 * t168 + t29 * t170 + t18 * t196 + t43 * t188) * t282 + (t145 * t223 + t146 * t224 - t232 * t234) * t286 + (t143 * t223 + t144 * t224 - t233 * t234) * t288 + (t182 * t223 + t183 * t224 + t196 * t234) * t283 + (t14 * (-t108 * t233 - t196 * t69) + t15 * (t108 * t232 + t196 * t70) + t11 * (-t232 * t69 + t233 * t70) - g(1) * t82 - g(2) * t81 - g(3) * t134 + (-t108 * t170 + t134 * t148 - t184 * t82 + t188 * t70 + t196 * t56 + t232 * t74) * t33 + (t108 * t168 - t134 * t147 + t184 * t81 - t188 * t69 - t196 * t55 - t233 * t74) * t32 + (t147 * t82 - t148 * t81 - t168 * t70 + t170 * t69 - t232 * t55 + t233 * t56) * t23) * m(6);];
tau = t24;
