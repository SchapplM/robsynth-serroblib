% Calculate vector of inverse dynamics joint torques for
% S4PPRR4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-31 16:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PPRR4_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR4_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR4_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR4_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_invdynJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR4_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR4_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRR4_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:18:26
% EndTime: 2019-12-31 16:18:39
% DurationCPUTime: 8.15s
% Computational Cost: add. (6864->420), mult. (10043->650), div. (0->0), fcn. (9714->6), ass. (0->209)
t159 = cos(pkin(6));
t278 = t159 ^ 2;
t158 = sin(pkin(6));
t279 = t158 ^ 2;
t287 = t278 + t279;
t157 = pkin(7) + qJ(3);
t155 = sin(t157);
t156 = cos(t157);
t196 = -Icges(4,5) * t155 - Icges(4,6) * t156;
t291 = t158 * t159;
t160 = sin(qJ(4));
t265 = rSges(5,2) * t160;
t161 = cos(qJ(4));
t266 = rSges(5,1) * t161;
t216 = -t265 + t266;
t91 = t155 * rSges(5,3) + t156 * t216;
t290 = 2 * qJD(3);
t289 = 2 * qJDD(3);
t236 = qJD(3) * t155;
t235 = qJD(3) * t156;
t281 = t156 * pkin(3) + t155 * pkin(5);
t195 = Icges(5,5) * t161 - Icges(5,6) * t160;
t84 = -Icges(5,3) * t156 + t155 * t195;
t251 = Icges(5,4) * t161;
t198 = -Icges(5,2) * t160 + t251;
t86 = -Icges(5,6) * t156 + t155 * t198;
t252 = Icges(5,4) * t160;
t201 = Icges(5,1) * t161 - t252;
t88 = -Icges(5,5) * t156 + t155 * t201;
t232 = qJD(4) * t155;
t234 = qJD(3) * t158;
t135 = t159 * t232 + t234;
t233 = qJD(3) * t159;
t136 = t158 * t232 - t233;
t241 = t159 * t160;
t242 = t158 * t161;
t130 = -t156 * t241 + t242;
t246 = t155 * t159;
t240 = t159 * t161;
t243 = t158 * t160;
t131 = t156 * t240 + t243;
t253 = Icges(5,4) * t131;
t55 = Icges(5,2) * t130 + Icges(5,6) * t246 + t253;
t119 = Icges(5,4) * t130;
t57 = Icges(5,1) * t131 + Icges(5,5) * t246 + t119;
t205 = -t160 * t55 + t161 * t57;
t128 = -t156 * t243 - t240;
t247 = t155 * t158;
t129 = t156 * t242 - t241;
t254 = Icges(5,4) * t129;
t54 = Icges(5,2) * t128 + Icges(5,6) * t247 + t254;
t118 = Icges(5,4) * t128;
t56 = Icges(5,1) * t129 + Icges(5,5) * t247 + t118;
t206 = -t160 * t54 + t161 * t56;
t280 = -(-t159 * t84 - t205) * t135 - (-t158 * t84 - t206) * t136;
t125 = (-Icges(5,2) * t161 - t252) * t155;
t231 = qJD(4) * t156;
t164 = t135 * (-Icges(5,2) * t131 + t119 + t57) + t136 * (-Icges(5,2) * t129 + t118 + t56) - t231 * (t125 + t88);
t229 = qJD(3) * qJD(4);
t176 = qJDD(4) * t155 + t156 * t229;
t82 = qJDD(3) * t158 + t159 * t176;
t277 = t82 / 0.2e1;
t83 = -qJDD(3) * t159 + t158 * t176;
t276 = t83 / 0.2e1;
t134 = -qJDD(4) * t156 + t155 * t229;
t275 = t134 / 0.2e1;
t274 = -t135 / 0.2e1;
t273 = t135 / 0.2e1;
t272 = -t136 / 0.2e1;
t271 = t136 / 0.2e1;
t264 = t156 * t84;
t52 = Icges(5,5) * t129 + Icges(5,6) * t128 + Icges(5,3) * t247;
t20 = t130 * t54 + t131 * t56 + t246 * t52;
t263 = t158 * t20;
t53 = Icges(5,5) * t131 + Icges(5,6) * t130 + Icges(5,3) * t246;
t19 = t128 * t55 + t129 * t57 + t247 * t53;
t262 = t19 * t159;
t29 = t128 * t86 + t129 * t88 + t247 * t84;
t261 = t29 * t155;
t30 = t130 * t86 + t131 * t88 + t246 * t84;
t260 = t30 * t155;
t133 = t281 * qJD(3);
t127 = (-rSges(5,1) * t160 - rSges(5,2) * t161) * t155;
t51 = qJD(3) * t91 + qJD(4) * t127;
t258 = -t133 - t51;
t139 = pkin(3) * t155 - pkin(5) * t156;
t90 = -t156 * rSges(5,3) + t155 * t216;
t257 = -t139 - t90;
t245 = t156 * t158;
t244 = t156 * t159;
t226 = t155 * t265;
t239 = rSges(5,3) * t245 + t158 * t226;
t238 = rSges(5,3) * t244 + t159 * t226;
t237 = qJD(2) * t159;
t230 = -m(3) - m(4) - m(5);
t228 = qJDD(2) * t159;
t227 = t155 * t266;
t225 = t160 * t236;
t224 = t161 * t236;
t223 = t156 * t234;
t222 = t156 * t233;
t220 = t233 / 0.2e1;
t219 = -t231 / 0.2e1;
t218 = t231 / 0.2e1;
t138 = rSges(4,1) * t156 - rSges(4,2) * t155;
t137 = rSges(4,1) * t155 + rSges(4,2) * t156;
t215 = t53 * t135 + t52 * t136;
t153 = qJDD(2) * t158;
t190 = -qJD(3) * t133 - qJDD(3) * t139;
t78 = -qJD(4) * t129 + t158 * t225;
t79 = qJD(4) * t128 - t158 * t224;
t45 = rSges(5,1) * t79 + rSges(5,2) * t78 + rSges(5,3) * t223;
t58 = rSges(5,1) * t129 + rSges(5,2) * t128 + rSges(5,3) * t247;
t16 = -t134 * t58 + t136 * t51 + t159 * t190 + t231 * t45 + t83 * t90 + t153;
t80 = -qJD(4) * t131 + t159 * t225;
t81 = qJD(4) * t130 - t159 * t224;
t46 = rSges(5,1) * t81 + rSges(5,2) * t80 + rSges(5,3) * t222;
t59 = rSges(5,1) * t131 + rSges(5,2) * t130 + rSges(5,3) * t246;
t17 = t134 * t59 - t135 * t51 + t158 * t190 - t231 * t46 - t82 * t90 - t228;
t212 = t158 * t16 - t159 * t17;
t18 = t128 * t54 + t129 * t56 + t247 * t52;
t211 = t18 * t158 + t262;
t21 = t130 * t55 + t131 * t57 + t246 * t53;
t210 = t159 * t21 + t263;
t23 = t155 * t206 - t156 * t52;
t24 = t155 * t205 - t156 * t53;
t209 = t23 * t158 + t24 * t159;
t208 = -t158 * t59 + t159 * t58;
t154 = qJD(2) * t158;
t207 = -t158 * (-t137 * t234 - t237) - t159 * (-t137 * t233 + t154);
t204 = -t160 * t86 + t161 * t88;
t197 = Icges(4,5) * t156 - Icges(4,6) * t155;
t194 = t287 * t138;
t193 = t287 * qJD(3) * t137;
t122 = t281 * t158;
t123 = t281 * t159;
t192 = t122 * t158 + t123 * t159;
t132 = t138 * qJD(3);
t191 = -qJD(3) * t132 - qJDD(3) * t137;
t85 = Icges(5,3) * t155 + t156 * t195;
t187 = -t204 + t85;
t39 = Icges(5,5) * t79 + Icges(5,6) * t78 + Icges(5,3) * t223;
t186 = t155 * t39 + t235 * t52;
t40 = Icges(5,5) * t81 + Icges(5,6) * t80 + Icges(5,3) * t222;
t185 = t155 * t40 + t235 * t53;
t124 = (-Icges(5,5) * t160 - Icges(5,6) * t161) * t155;
t48 = qJD(3) * t85 + qJD(4) * t124;
t184 = t155 * t48 + t235 * t84;
t183 = qJD(3) * t139;
t22 = qJD(3) * t192 + t135 * t58 - t136 * t59 + qJD(1);
t182 = t22 * t208;
t126 = (-Icges(5,1) * t160 - t251) * t155;
t173 = qJD(3) * t196;
t172 = t124 * t231 - t135 * (Icges(5,5) * t130 - Icges(5,6) * t131) - t136 * (Icges(5,5) * t128 - Icges(5,6) * t129);
t89 = Icges(5,5) * t155 + t156 * t201;
t87 = Icges(5,6) * t155 + t156 * t198;
t168 = t155 * t172;
t165 = (Icges(5,1) * t130 - t253 - t55) * t135 + (Icges(5,1) * t128 - t254 - t54) * t136 - (t126 - t86) * t231;
t163 = t287 * t196;
t162 = (-t187 * t231 - t280) * t155;
t147 = pkin(5) * t244;
t146 = pkin(5) * t245;
t121 = t137 * t159;
t120 = t137 * t158;
t113 = t196 * t159;
t112 = t196 * t158;
t111 = t159 * t183;
t110 = t158 * t183;
t103 = t159 * t173;
t102 = t158 * t173;
t93 = Icges(4,3) * t158 + t159 * t197;
t92 = -Icges(4,3) * t159 + t158 * t197;
t77 = -t159 * t227 + t238;
t76 = -t158 * t227 + t239;
t75 = t88 * t159;
t74 = t88 * t158;
t73 = t86 * t159;
t72 = t86 * t158;
t69 = rSges(5,1) * t130 - rSges(5,2) * t131;
t68 = rSges(5,1) * t128 - rSges(5,2) * t129;
t61 = t158 * t191 - t228;
t60 = t159 * t191 + t153;
t50 = qJD(3) * t89 + qJD(4) * t126;
t49 = qJD(3) * t87 + qJD(4) * t125;
t47 = qJD(3) * t194 + qJD(1);
t44 = Icges(5,1) * t81 + Icges(5,4) * t80 + Icges(5,5) * t222;
t43 = Icges(5,1) * t79 + Icges(5,4) * t78 + Icges(5,5) * t223;
t42 = Icges(5,4) * t81 + Icges(5,2) * t80 + Icges(5,6) * t222;
t41 = Icges(5,4) * t79 + Icges(5,2) * t78 + Icges(5,6) * t223;
t34 = t155 * t204 - t264;
t33 = -qJD(3) * t193 + qJDD(3) * t194 + qJDD(1);
t32 = -t135 * t90 - t139 * t234 - t231 * t59 - t237;
t31 = t136 * t90 - t139 * t233 + t231 * t58 + t154;
t15 = (qJD(3) * t204 - t48) * t156 + (qJD(3) * t84 - t160 * t49 + t161 * t50 + (-t160 * t88 - t161 * t86) * qJD(4)) * t155;
t14 = t130 * t49 + t131 * t50 + t159 * t184 + t80 * t86 + t81 * t88;
t13 = t128 * t49 + t129 * t50 + t158 * t184 + t78 * t86 + t79 * t88;
t12 = t135 * t45 - t136 * t46 + t58 * t82 - t59 * t83 + qJDD(1) + t192 * qJDD(3) + (-t110 * t158 - t111 * t159) * qJD(3);
t11 = t135 * t24 + t136 * t23 - t231 * t34;
t10 = t130 * t42 + t131 * t44 + t159 * t185 + t55 * t80 + t57 * t81;
t9 = t130 * t41 + t131 * t43 + t159 * t186 + t54 * t80 + t56 * t81;
t8 = t128 * t42 + t129 * t44 + t158 * t185 + t55 * t78 + t57 * t79;
t7 = t128 * t41 + t129 * t43 + t158 * t186 + t54 * t78 + t56 * t79;
t6 = t135 * t21 + t136 * t20 - t231 * t30;
t5 = t135 * t19 + t136 * t18 - t231 * t29;
t4 = (qJD(3) * t205 - t40) * t156 + (qJD(3) * t53 - t160 * t42 + t161 * t44 + (-t160 * t57 - t161 * t55) * qJD(4)) * t155;
t3 = (qJD(3) * t206 - t39) * t156 + (qJD(3) * t52 - t160 * t41 + t161 * t43 + (-t160 * t56 - t161 * t54) * qJD(4)) * t155;
t2 = t135 * t10 + t134 * t30 + t136 * t9 - t14 * t231 + t20 * t83 + t82 * t21;
t1 = -t13 * t231 + t134 * t29 + t135 * t8 + t136 * t7 + t83 * t18 + t19 * t82;
t25 = [(m(2) + m(3)) * qJDD(1) + m(4) * t33 + m(5) * t12 + (-m(2) + t230) * g(3); t230 * (g(1) * t158 - g(2) * t159) + m(4) * (t158 * t60 - t159 * t61) + m(5) * t212 + m(3) * t287 * qJDD(2); (t158 * t24 - t159 * t23) * t275 + (t158 * t19 - t159 * t18) * t276 + (t158 * t21 - t159 * t20) * t277 + (t158 * t8 - t159 * t7) * t271 + ((-t128 * t73 - t129 * t75) * t135 + (-t128 * t72 - t129 * t74) * t136 + (t261 + (-t128 * t87 - t129 * t89 + t262) * t156) * qJD(4) + (((t18 - t264) * qJD(4) + t215) * t156 + t162) * t158) * t272 + (t10 * t158 - t159 * t9) * t273 + ((-t130 * t73 - t131 * t75) * t135 + (-t130 * t72 - t131 * t74) * t136 + (t260 + (-t130 * t87 - t131 * t89 + t263) * t156) * qJD(4) + (((t21 - t264) * qJD(4) + t215) * t156 + t162) * t159) * t274 - t11 * t232 / 0.2e1 + (((t160 * t73 - t161 * t75 + t53) * t135 + (t160 * t72 - t161 * t74 + t52) * t136 + t34 * qJD(4)) * t155 + ((t187 * t156 + (t160 * t87 - t161 * t89 - t84) * t155 + t209) * qJD(4) + t280) * t156) * t218 + (t112 * qJD(3) * t278 + (-t159 * t113 + t163) * t234) * t220 - (t113 * qJD(3) * t279 + (-t158 * t112 + t163) * t233) * t234 / 0.2e1 + ((t279 * t93 - t291 * t92) * t289 + (-t102 * t291 + t103 * t279) * t290 + t2) * t158 / 0.2e1 - ((t278 * t92 - t291 * t93) * t289 + (t102 * t278 - t103 * t291) * t290 + t1) * t159 / 0.2e1 + ((-t3 + t6) * t159 + (t4 + t5) * t158) * t219 + (-g(1) * (t147 + t238) - g(2) * (t146 + t239) - g(3) * (t91 + t281) - (g(1) * t159 + g(2) * t158) * t155 * (-pkin(3) - t266) + (t16 * t257 + t31 * t258 + t12 * (t123 + t59) + t22 * (-t111 + t46)) * t159 + (t17 * t257 + t32 * t258 + t12 * (t122 + t58) + t22 * (-t110 + t45)) * t158 - t31 * (t136 * t91 - t233 * t281) - t32 * (-t135 * t91 - t234 * t281) - t22 * ((-pkin(3) * t247 + t146) * t234 + (-pkin(3) * t246 + t147) * t233 + t135 * t76 - t136 * t77) - ((-t31 * t58 + t32 * t59) * t155 + (t31 * (t158 * t90 + t76) + t32 * (-t159 * t90 - t77) + t182) * t156) * qJD(4)) * m(5) + (g(1) * t121 + g(2) * t120 - g(3) * t138 - (t47 * (-t120 * t158 - t121 * t159) + t207 * t138) * qJD(3) + t33 * t194 - t47 * t193 + (-t158 * t61 - t159 * t60) * t137 + t207 * t132) * m(4); t156 * t6 * t220 + t2 * t246 / 0.2e1 + (t155 * t210 - t156 * t30) * t277 + (-t14 * t156 + (t10 * t159 + t158 * t9) * t155 + (t156 * t210 + t260) * qJD(3)) * t273 + t5 * t223 / 0.2e1 + t1 * t247 / 0.2e1 + (t155 * t211 - t156 * t29) * t276 + (-t13 * t156 + (t158 * t7 + t159 * t8) * t155 + (t156 * t211 + t261) * qJD(3)) * t271 + t11 * t236 / 0.2e1 - t156 * (t34 * t134 + t4 * t135 + t3 * t136 - t15 * t231 + t23 * t83 + t24 * t82) / 0.2e1 + (t155 * t209 - t34 * t156) * t275 + (-t15 * t156 + (t158 * t3 + t159 * t4) * t155 + (t34 * t155 + t156 * t209) * qJD(3)) * t219 + (t130 * t164 + t131 * t165 - t159 * t168) * t274 + (t128 * t164 + t129 * t165 - t158 * t168) * t272 + (t172 * t156 + (-t160 * t164 + t165 * t161) * t155) * t218 + ((t16 * t58 - t17 * t59 + t31 * t45 - t32 * t46 + (t182 + (t158 * t31 - t159 * t32) * t90) * qJD(3)) * t156 + (t31 * (-qJD(3) * t58 + t158 * t51) + t32 * (qJD(3) * t59 - t159 * t51) + t12 * t208 + t22 * (-t158 * t46 + t159 * t45) + t212 * t90) * t155 - t31 * (t127 * t136 + t231 * t68) - t32 * (-t127 * t135 - t231 * t69) - t22 * (t135 * t68 - t136 * t69) - g(1) * t69 - g(2) * t68 - g(3) * t127) * m(5);];
tau = t25;
