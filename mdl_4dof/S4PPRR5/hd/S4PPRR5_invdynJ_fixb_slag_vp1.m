% Calculate vector of inverse dynamics joint torques for
% S4PPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PPRR5_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR5_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR5_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR5_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR5_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR5_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRR5_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:19:41
% EndTime: 2019-12-31 16:19:51
% DurationCPUTime: 7.86s
% Computational Cost: add. (3634->434), mult. (10043->659), div. (0->0), fcn. (9714->6), ass. (0->210)
t159 = cos(qJ(3));
t291 = Icges(4,1) - Icges(4,2);
t293 = t159 * t291;
t292 = -2 * Icges(4,4);
t154 = sin(pkin(6));
t157 = sin(qJ(3));
t186 = Icges(4,5) * t159 - Icges(4,6) * t157;
t118 = t186 * t154;
t290 = t157 * t292 + t293;
t102 = qJD(3) * t118;
t155 = cos(pkin(6));
t276 = t155 ^ 2;
t289 = t276 * t102;
t277 = t154 ^ 2;
t119 = t186 * t155;
t283 = qJD(3) * t119;
t288 = t277 * t283;
t287 = t276 + t277;
t286 = t154 * t155;
t285 = 0.2e1 * qJD(3);
t284 = 2 * qJDD(3);
t153 = t159 * pkin(5);
t142 = -pkin(3) * t157 + t153;
t227 = qJD(3) * t154;
t226 = qJD(3) * t155;
t225 = qJD(3) * t157;
t224 = qJD(3) * t159;
t156 = sin(qJ(4));
t158 = cos(qJ(4));
t184 = Icges(5,5) * t158 - Icges(5,6) * t156;
t95 = Icges(5,3) * t157 + t159 * t184;
t240 = Icges(5,4) * t158;
t187 = -Icges(5,2) * t156 + t240;
t97 = Icges(5,6) * t157 + t159 * t187;
t241 = Icges(5,4) * t156;
t190 = Icges(5,1) * t158 - t241;
t99 = Icges(5,5) * t157 + t159 * t190;
t222 = qJD(4) * t159;
t132 = t155 * t222 + t227;
t133 = -t154 * t222 + t226;
t233 = t156 * t157;
t116 = t154 * t158 + t155 * t233;
t234 = t155 * t159;
t232 = t157 * t158;
t117 = -t154 * t156 + t155 * t232;
t242 = Icges(5,4) * t117;
t51 = Icges(5,2) * t116 + Icges(5,6) * t234 - t242;
t113 = Icges(5,4) * t116;
t53 = -Icges(5,1) * t117 + Icges(5,5) * t234 + t113;
t197 = t156 * t51 - t158 * t53;
t114 = -t154 * t233 + t155 * t158;
t235 = t154 * t159;
t115 = t154 * t232 + t155 * t156;
t243 = Icges(5,4) * t115;
t50 = Icges(5,2) * t114 - Icges(5,6) * t235 + t243;
t112 = Icges(5,4) * t114;
t52 = Icges(5,1) * t115 - Icges(5,5) * t235 + t112;
t198 = t156 * t50 - t158 * t52;
t160 = t132 * (-t155 * t95 + t197) + t133 * (t95 * t154 + t198);
t280 = t118 + (t159 ^ 2 * t292 + (-t290 - t293) * t157) * t155;
t236 = t154 * t157;
t279 = -t157 * (-Icges(4,6) * t155 + t154 * t290) - t159 * (0.2e1 * Icges(4,4) * t235 + Icges(4,5) * t155 + t236 * t291);
t130 = (-Icges(5,1) * t156 - t240) * t159;
t223 = qJD(4) * t157;
t278 = t132 * (-Icges(5,1) * t116 - t242 + t51) + t133 * (-Icges(5,1) * t114 + t243 + t50) - t223 * (t130 - t97);
t129 = (-Icges(5,2) * t158 - t241) * t159;
t161 = t132 * (Icges(5,2) * t117 + t113 + t53) + t133 * (-Icges(5,2) * t115 + t112 + t52) + t223 * (t129 + t99);
t220 = qJD(3) * qJD(4);
t165 = -qJDD(4) * t159 + t157 * t220;
t82 = qJDD(3) * t154 - t155 * t165;
t275 = t82 / 0.2e1;
t83 = qJDD(3) * t155 + t154 * t165;
t274 = t83 / 0.2e1;
t273 = -t132 / 0.2e1;
t272 = t132 / 0.2e1;
t271 = -t133 / 0.2e1;
t270 = t133 / 0.2e1;
t136 = qJDD(4) * t157 + t159 * t220;
t269 = t136 / 0.2e1;
t265 = g(2) * t155;
t260 = rSges(5,1) * t158;
t259 = rSges(5,2) * t156;
t258 = t157 * rSges(5,3);
t257 = t157 * t95;
t152 = t159 * rSges(5,3);
t49 = -Icges(5,5) * t117 + Icges(5,6) * t116 + Icges(5,3) * t234;
t19 = t114 * t51 + t115 * t53 - t235 * t49;
t256 = t19 * t155;
t48 = Icges(5,5) * t115 + Icges(5,6) * t114 - Icges(5,3) * t235;
t20 = t116 * t50 - t117 * t52 + t234 * t48;
t255 = t20 * t154;
t32 = t114 * t97 + t115 * t99 - t235 * t95;
t254 = t32 * t159;
t33 = t116 * t97 - t117 * t99 + t234 * t95;
t253 = t33 * t159;
t135 = t142 * qJD(3);
t131 = (-rSges(5,1) * t156 - rSges(5,2) * t158) * t159;
t207 = -t259 + t260;
t59 = qJD(4) * t131 + (-t157 * t207 + t152) * qJD(3);
t246 = t135 + t59;
t101 = t159 * t207 + t258;
t143 = pkin(3) * t159 + t157 * pkin(5);
t231 = t101 + t143;
t230 = pkin(3) * t235 + pkin(5) * t236;
t229 = rSges(5,2) * t233 + t152;
t228 = qJD(3) * t143;
t221 = -m(3) - m(4) - m(5);
t219 = t159 * t260;
t218 = t159 * t259;
t110 = t143 * t227;
t217 = t154 * t225;
t216 = t155 * t225;
t215 = t156 * t224;
t214 = t158 * t224;
t213 = -pkin(3) - t260;
t212 = -t226 / 0.2e1;
t211 = -t223 / 0.2e1;
t210 = t223 / 0.2e1;
t141 = rSges(4,1) * t159 - t157 * rSges(4,2);
t208 = t157 * rSges(4,1) + rSges(4,2) * t159;
t206 = -t49 * t132 - t48 * t133;
t18 = t114 * t50 + t115 * t52 - t235 * t48;
t205 = t18 * t154 - t256;
t21 = t116 * t51 - t117 * t53 + t234 * t49;
t204 = -t21 * t155 + t255;
t23 = t157 * t48 - t159 * t198;
t24 = t157 * t49 - t159 * t197;
t203 = t23 * t154 - t24 * t155;
t54 = t115 * rSges(5,1) + t114 * rSges(5,2) - rSges(5,3) * t235;
t55 = -t117 * rSges(5,1) + t116 * rSges(5,2) + rSges(5,3) * t234;
t202 = t154 * t55 + t155 * t54;
t150 = qJDD(2) * t154;
t134 = t208 * qJD(3);
t179 = -qJD(3) * t134 + qJDD(3) * t141;
t201 = t154 * (t154 * t179 + t150) - t276 * (-qJDD(2) - t179);
t151 = qJD(2) * t154;
t200 = (t141 * t227 + t151) * t154 - (-qJD(3) * t141 - qJD(2)) * t276;
t199 = t287 * t208;
t196 = t156 * t97 - t158 * t99;
t185 = Icges(4,5) * t157 + Icges(4,6) * t159;
t124 = t141 * t154;
t125 = t141 * t155;
t181 = -t124 * t227 - t125 * t226;
t126 = t142 * t154;
t127 = t142 * t155;
t180 = t126 * t154 + t127 * t155;
t178 = qJD(3) * t135 + qJDD(3) * t143;
t94 = Icges(5,3) * t159 - t157 * t184;
t175 = t196 + t94;
t70 = -qJD(4) * t115 - t154 * t215;
t71 = qJD(4) * t114 + t154 * t214;
t39 = Icges(5,5) * t71 + Icges(5,6) * t70 + Icges(5,3) * t217;
t174 = -t159 * t39 + t225 * t48;
t72 = qJD(4) * t117 + t155 * t215;
t73 = qJD(4) * t116 - t155 * t214;
t40 = Icges(5,5) * t73 + Icges(5,6) * t72 - Icges(5,3) * t216;
t173 = -t159 * t40 + t225 * t49;
t128 = (-Icges(5,5) * t156 - Icges(5,6) * t158) * t159;
t56 = qJD(3) * t94 + qJD(4) * t128;
t172 = -t159 * t56 + t225 * t95;
t22 = qJD(3) * t180 - t132 * t54 + t133 * t55 + qJD(1);
t170 = t22 * t202;
t80 = rSges(5,3) * t236 + (-t218 + t219) * t154;
t164 = t175 * t157;
t163 = t128 * t223 + t132 * (Icges(5,5) * t116 + Icges(5,6) * t117) + t133 * (Icges(5,5) * t114 - Icges(5,6) * t115);
t98 = Icges(5,5) * t159 - t157 * t190;
t96 = Icges(5,6) * t159 - t157 * t187;
t138 = t155 * t218;
t111 = t155 * t228;
t100 = -rSges(5,1) * t232 + t229;
t85 = Icges(4,3) * t154 - t155 * t185;
t84 = Icges(4,3) * t155 + t154 * t185;
t81 = t138 + (-t219 - t258) * t155;
t79 = t99 * t155;
t78 = t99 * t154;
t77 = t97 * t155;
t76 = t97 * t154;
t69 = rSges(5,1) * t116 + rSges(5,2) * t117;
t68 = rSges(5,1) * t114 - rSges(5,2) * t115;
t58 = qJD(3) * t98 + qJD(4) * t130;
t57 = qJD(3) * t96 + qJD(4) * t129;
t47 = -qJD(3) * t199 + qJD(1);
t46 = rSges(5,1) * t73 + rSges(5,2) * t72 - rSges(5,3) * t216;
t45 = rSges(5,1) * t71 + rSges(5,2) * t70 + rSges(5,3) * t217;
t44 = Icges(5,1) * t73 + Icges(5,4) * t72 - Icges(5,5) * t216;
t43 = Icges(5,1) * t71 + Icges(5,4) * t70 + Icges(5,5) * t217;
t42 = Icges(5,4) * t73 + Icges(5,2) * t72 - Icges(5,6) * t216;
t41 = Icges(5,4) * t71 + Icges(5,2) * t70 + Icges(5,6) * t217;
t38 = -t159 * t196 + t257;
t31 = qJD(3) * t181 - qJDD(3) * t199 + qJDD(1);
t30 = t54 * t223 - t101 * t133 + (-qJD(2) - t228) * t155;
t29 = t101 * t132 - t223 * t55 + t110 + t151;
t17 = t45 * t223 - t101 * t83 - t133 * t59 + t136 * t54 + (-qJDD(2) - t178) * t155;
t16 = t101 * t82 + t132 * t59 - t136 * t55 + t154 * t178 - t223 * t46 + t150;
t15 = (qJD(3) * t196 + t56) * t157 + (qJD(3) * t95 - t156 * t57 + t158 * t58 + (-t156 * t99 - t158 * t97) * qJD(4)) * t159;
t14 = t116 * t57 - t117 * t58 - t155 * t172 + t72 * t97 + t73 * t99;
t13 = t114 * t57 + t115 * t58 + t154 * t172 + t70 * t97 + t71 * t99;
t12 = -t132 * t45 + t133 * t46 - t54 * t82 + t55 * t83 + qJDD(1) + t180 * qJDD(3) + (-t110 * t154 - t111 * t155) * qJD(3);
t11 = t132 * t24 + t133 * t23 + t223 * t38;
t10 = t116 * t42 - t117 * t44 - t155 * t173 + t72 * t51 + t73 * t53;
t9 = t116 * t41 - t117 * t43 - t155 * t174 + t72 * t50 + t73 * t52;
t8 = t114 * t42 + t115 * t44 + t154 * t173 + t70 * t51 + t71 * t53;
t7 = t114 * t41 + t115 * t43 + t154 * t174 + t70 * t50 + t71 * t52;
t6 = t132 * t21 + t133 * t20 + t223 * t33;
t5 = t132 * t19 + t133 * t18 + t223 * t32;
t4 = (qJD(3) * t197 + t40) * t157 + (qJD(3) * t49 - t156 * t42 + t158 * t44 + (-t156 * t53 - t158 * t51) * qJD(4)) * t159;
t3 = (qJD(3) * t198 + t39) * t157 + (qJD(3) * t48 - t156 * t41 + t158 * t43 + (-t156 * t52 - t158 * t50) * qJD(4)) * t159;
t2 = t132 * t10 + t133 * t9 + t136 * t33 + t14 * t223 + t20 * t83 + t82 * t21;
t1 = t13 * t223 + t132 * t8 + t133 * t7 + t136 * t32 + t83 * t18 + t19 * t82;
t25 = [(m(2) + m(3)) * qJDD(1) + m(4) * t31 + m(5) * t12 + (-m(2) + t221) * g(3); t221 * (g(1) * t154 - t265) + m(4) * t201 + m(5) * (t154 * t16 - t155 * t17) + m(3) * t287 * qJDD(2); (t289 + (t280 * t154 + (-t119 - t279) * t155) * t227) * t212 - (-t288 + (t279 * t155 + (t118 - t280) * t154) * t226) * t227 / 0.2e1 + (t154 * t19 + t155 * t18) * t274 + (t154 * t21 + t155 * t20) * t275 + (t154 * t24 + t155 * t23) * t269 + (t154 * t8 + t155 * t7) * t270 + ((t114 * t76 + t115 * t78) * t133 + (-t114 * t77 - t115 * t79) * t132 + (t254 + (t114 * t96 + t115 * t98 - t256) * t157) * qJD(4) + (((t18 + t257) * qJD(4) - t206) * t157 + (-t175 * t223 - t160) * t159) * t154) * t271 + (t10 * t154 + t155 * t9) * t272 + ((t116 * t76 - t117 * t78) * t133 + (-t116 * t77 + t117 * t79) * t132 + (t253 + (t116 * t96 - t117 * t98 + t255) * t157) * qJD(4) + (((-t21 - t257) * qJD(4) + t206) * t157 + (qJD(4) * t164 + t160) * t159) * t155) * t273 - t11 * t222 / 0.2e1 + (t154 * t5 + ((-t156 * t76 + t158 * t78 + t48) * t133 + (t156 * t77 - t158 * t79 + t49) * t132 + t38 * qJD(4)) * t159 + ((t164 + (-t156 * t96 + t158 * t98 + t95) * t159 + t203) * qJD(4) + t160) * t157) * t211 + (t154 * t4 + (t3 + t6) * t155) * t210 + ((t277 * t85 + t286 * t84) * t284 + (t102 * t286 - t288) * t285 + t2) * t154 / 0.2e1 + ((t276 * t84 + t286 * t85) * t284 + (-t283 * t286 + t289) * t285 + t1) * t155 / 0.2e1 + (-t29 * (t132 * t100 + t142 * t227) - t30 * (-t133 * t100 - t142 * t226) - t22 * (-t132 * t80 + t133 * t81 - t227 * t230 - t228 * t276) - ((-t29 * t55 + t30 * t54) * t159 + (t29 * (-t101 * t155 - t81) + t30 * (-t101 * t154 + t80) + t170) * t157) * qJD(4) - g(1) * (t80 + t230) - g(2) * t138 - g(3) * (t157 * t213 + t153 + t229) - (t213 * t159 + (-rSges(5,3) - pkin(5)) * t157) * t265 + (-t17 * t231 - t30 * t246 + t12 * (t127 + t55) + t22 * (-t111 + t46)) * t155 + (t16 * t231 + t29 * t246 + t12 * (t126 - t54) + t22 * (-t110 - t45)) * t154) * m(5) + (-t134 * t200 + t141 * t201 + t181 * t47 - t199 * t31 - (t47 * (-t124 * t154 - t125 * t155) - t200 * t208) * qJD(3) - g(1) * t124 + g(2) * t125 + g(3) * t208) * m(4); t5 * t217 / 0.2e1 - t1 * t235 / 0.2e1 + (t32 * t157 - t159 * t205) * t274 + (t13 * t157 + (-t154 * t7 + t155 * t8) * t159 + (t157 * t205 + t254) * qJD(3)) * t270 + t157 * t6 * t212 + t2 * t234 / 0.2e1 + (t33 * t157 - t159 * t204) * t275 + (t14 * t157 + (t10 * t155 - t154 * t9) * t159 + (t157 * t204 + t253) * qJD(3)) * t272 + t11 * t224 / 0.2e1 + t157 * (t4 * t132 + t3 * t133 + t38 * t136 + t15 * t223 + t23 * t83 + t24 * t82) / 0.2e1 + (t38 * t157 - t159 * t203) * t269 + (t15 * t157 + (-t154 * t3 + t155 * t4) * t159 + (t157 * t203 + t38 * t159) * qJD(3)) * t210 + (t114 * t161 - t115 * t278 - t163 * t235) * t271 + (t161 * t116 + t117 * t278 + t163 * t234) * t273 + (t163 * t157 + (-t156 * t161 - t158 * t278) * t159) * t211 + ((-t16 * t55 + t17 * t54 - t29 * t46 + t30 * t45 + (t170 + (-t154 * t30 - t155 * t29) * t101) * qJD(3)) * t157 + (t29 * (-qJD(3) * t55 + t155 * t59) + t30 * (qJD(3) * t54 + t154 * t59) - t12 * t202 + t22 * (-t154 * t46 - t155 * t45) + (t154 * t17 + t155 * t16) * t101) * t159 - t29 * (t131 * t132 - t69 * t223) - t30 * (-t131 * t133 + t68 * t223) - t22 * (-t132 * t68 + t133 * t69) - g(1) * t68 - g(2) * t69 - g(3) * t131) * m(5);];
tau = t25;
