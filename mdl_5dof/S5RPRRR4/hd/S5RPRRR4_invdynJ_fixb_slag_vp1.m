% Calculate vector of inverse dynamics joint torques for
% S5RPRRR4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-05 18:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRR4_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR4_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR4_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR4_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR4_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR4_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:14:19
% EndTime: 2019-12-05 18:14:27
% DurationCPUTime: 4.35s
% Computational Cost: add. (10808->364), mult. (6100->453), div. (0->0), fcn. (4528->10), ass. (0->217)
t178 = qJD(1) ^ 2;
t173 = qJ(1) + pkin(9);
t167 = qJ(3) + t173;
t159 = cos(t167);
t172 = qJD(1) + qJD(3);
t248 = t159 * t172;
t237 = pkin(3) * t248;
t165 = cos(t173);
t177 = cos(qJ(1));
t287 = pkin(1) * t177;
t221 = pkin(2) * t165 + t287;
t301 = t221 * qJD(1);
t187 = t301 + t237;
t160 = qJ(4) + t167;
t154 = sin(t160);
t155 = cos(t160);
t106 = rSges(5,1) * t155 - t154 * rSges(5,2);
t166 = qJD(4) + t172;
t96 = t166 * t106;
t60 = -t187 - t96;
t174 = sin(qJ(5));
t176 = cos(qJ(5));
t145 = rSges(6,1) * t174 + rSges(6,2) * t176;
t158 = sin(t167);
t217 = rSges(4,1) * t158 + rSges(4,2) * t159;
t101 = t217 * t172;
t164 = sin(t173);
t286 = pkin(2) * t164;
t175 = sin(qJ(1));
t288 = pkin(1) * t175;
t222 = t286 + t288;
t204 = t222 * qJD(1);
t71 = t204 + t101;
t255 = t154 * t166;
t120 = pkin(4) * t255;
t272 = rSges(6,2) * t174;
t273 = rSges(6,1) * t176;
t147 = -t272 + t273;
t128 = t147 * qJD(5);
t171 = qJDD(1) + qJDD(3);
t163 = qJDD(4) + t171;
t170 = t172 ^ 2;
t162 = t178 * t288;
t191 = -qJDD(1) * t221 + t178 * t286 + t162;
t284 = pkin(3) * t159;
t285 = pkin(3) * t158;
t186 = t170 * t285 - t171 * t284 + t191;
t242 = qJD(5) * t154;
t252 = t155 * t166;
t108 = pkin(4) * t155 + pkin(8) * t154;
t251 = t155 * t174;
t136 = rSges(6,2) * t251;
t250 = t155 * t176;
t205 = -rSges(6,1) * t250 - rSges(6,3) * t154;
t79 = -t136 - t205;
t265 = -t108 - t79;
t253 = t154 * t176;
t236 = rSges(6,1) * t253;
t309 = t145 * qJD(5);
t232 = -t309 * t155 - t166 * t236;
t254 = t154 * t174;
t245 = rSges(6,2) * t254 + t155 * rSges(6,3);
t47 = t166 * t245 + t232;
t240 = qJD(5) * t166;
t87 = qJDD(5) * t154 + t155 * t240;
t13 = t128 * t242 + t87 * t145 + (-pkin(8) * t252 + t120 - t47) * t166 + t265 * t163 + t186;
t308 = -g(2) + t13;
t81 = rSges(5,1) * t255 + rSges(5,2) * t252;
t307 = -t106 * t163 + t166 * t81 - g(2) + t186;
t275 = rSges(4,1) * t159;
t114 = -t158 * rSges(4,2) + t275;
t312 = t101 * t172 - t114 * t171 - g(2) + t191;
t100 = t166 * t108;
t246 = t177 * t178;
t183 = (-qJDD(1) * t164 - t165 * t178) * pkin(2) + (-qJDD(1) * t175 - t246) * pkin(1);
t180 = (-t158 * t171 - t159 * t170) * pkin(3) + t183;
t241 = qJD(5) * t155;
t231 = t166 * t136 + t309 * t154;
t48 = t166 * t205 + t231;
t230 = -pkin(4) - t273;
t283 = t155 * pkin(8);
t57 = t154 * t230 + t245 + t283;
t88 = qJDD(5) * t155 - t154 * t240;
t12 = -t128 * t241 - t88 * t145 + t57 * t163 + (t48 - t100) * t166 + t180;
t306 = -g(3) + t12;
t216 = -rSges(5,1) * t154 - rSges(5,2) * t155;
t82 = -rSges(5,1) * t252 + rSges(5,2) * t255;
t305 = t163 * t216 + t166 * t82 - g(3) + t180;
t249 = t158 * t172;
t137 = rSges(4,2) * t249;
t102 = -rSges(4,1) * t248 + t137;
t311 = t102 * t172 - t171 * t217 - g(3) + t183;
t168 = Icges(6,4) * t176;
t210 = -Icges(6,2) * t174 + t168;
t303 = Icges(6,1) * t174 + t168;
t243 = t303 + t210;
t263 = Icges(6,4) * t174;
t141 = Icges(6,2) * t176 + t263;
t144 = Icges(6,1) * t176 - t263;
t244 = t141 - t144;
t310 = (t174 * t243 + t176 * t244) * t166;
t76 = Icges(6,4) * t250 - Icges(6,2) * t251 + Icges(6,6) * t154;
t134 = Icges(6,4) * t251;
t78 = Icges(6,1) * t250 + Icges(6,5) * t154 - t134;
t212 = t174 * t76 - t176 * t78;
t304 = t212 * t155;
t157 = t164 * rSges(3,2);
t276 = rSges(3,1) * t165;
t220 = -t276 - t287;
t110 = t157 + t220;
t218 = rSges(3,1) * t164 + rSges(3,2) * t165;
t198 = t218 + t288;
t238 = pkin(3) * t249;
t189 = t204 + t238;
t207 = t236 - t245;
t67 = t166 * t207;
t206 = t145 * t241 + t67 - t166 * (-pkin(4) * t154 + t283);
t30 = t189 + t206;
t107 = t145 * t242;
t31 = t166 * t265 + t107 - t187;
t302 = t154 * t31 + t155 * t30;
t200 = t210 * t166;
t297 = -Icges(6,6) * t166 + qJD(5) * t141;
t42 = t154 * t297 - t155 * t200;
t201 = t144 * t166;
t295 = -Icges(6,5) * t166 + qJD(5) * t303;
t44 = t154 * t295 - t155 * t201;
t75 = Icges(6,6) * t155 - t154 * t210;
t133 = Icges(6,4) * t254;
t77 = -Icges(6,1) * t253 + Icges(6,5) * t155 + t133;
t45 = t174 * t77 + t176 * t75;
t140 = Icges(6,5) * t176 - Icges(6,6) * t174;
t73 = Icges(6,3) * t155 - t140 * t154;
t300 = qJD(5) * t45 - t166 * t73 + t174 * t42 - t176 * t44;
t41 = -t154 * t200 - t155 * t297;
t43 = -t154 * t201 - t155 * t295;
t46 = t174 * t78 + t176 * t76;
t74 = Icges(6,5) * t250 - Icges(6,6) * t251 + Icges(6,3) * t154;
t299 = qJD(5) * t46 - t166 * t74 + t174 * t41 - t176 * t43;
t139 = Icges(6,5) * t174 + Icges(6,6) * t176;
t298 = -Icges(6,3) * t166 + qJD(5) * t139;
t126 = t210 * qJD(5);
t127 = t144 * qJD(5);
t209 = t176 * t141 + t174 * t303;
t296 = qJD(5) * t209 + t126 * t174 - t127 * t176 - t139 * t166;
t277 = -Icges(6,2) * t250 - t134 + t78;
t279 = t155 * t303 + t76;
t293 = t174 * t277 + t176 * t279;
t278 = Icges(6,2) * t253 + t133 + t77;
t280 = -t154 * t303 + t75;
t292 = -t174 * t278 - t176 * t280;
t291 = t87 / 0.2e1;
t290 = t88 / 0.2e1;
t289 = -rSges(6,3) - pkin(8);
t282 = t155 * t73 + t75 * t254;
t281 = t154 * t73 + t77 * t250;
t68 = t166 * t79;
t268 = t174 * t75;
t267 = t176 * t77;
t208 = t174 * t141 - t176 * t303;
t89 = t139 * t154;
t54 = -t155 * t208 + t89;
t266 = t54 * t166;
t258 = t139 * t155;
t257 = t140 * t166;
t97 = t145 * t154;
t256 = t145 * t155;
t247 = t172 * t114;
t239 = m(3) + m(4) + m(5);
t235 = -t100 - t68 + t107;
t229 = -t242 / 0.2e1;
t228 = t242 / 0.2e1;
t227 = -t241 / 0.2e1;
t226 = t241 / 0.2e1;
t224 = -t74 - t267;
t223 = t120 - t232;
t148 = rSges(2,1) * t177 - t175 * rSges(2,2);
t219 = rSges(2,1) * t175 + rSges(2,2) * t177;
t24 = -t253 * t77 + t282;
t25 = t155 * t74 - t253 * t78 + t76 * t254;
t215 = t25 * t154 + t24 * t155;
t26 = -t251 * t75 + t281;
t27 = t154 * t74 - t304;
t214 = t27 * t154 + t26 * t155;
t213 = -t267 + t268;
t86 = -t106 - t284;
t85 = t216 - t285;
t194 = -t257 * t154 - t155 * t298 + t166 * t212;
t193 = t154 * t298 - t257 * t155 + t166 * t213;
t192 = t140 * qJD(5) + t166 * t208;
t190 = t154 * t207 + t155 * t79;
t58 = t154 * t289 + t155 * t230 + t136;
t55 = t57 - t285;
t56 = t58 - t284;
t10 = qJD(5) * t214 + t266;
t16 = -qJD(5) * t213 + t174 * t44 + t176 * t42;
t17 = -qJD(5) * t212 + t174 * t43 + t176 * t41;
t20 = t192 * t154 - t155 * t296;
t21 = t154 * t296 + t192 * t155;
t53 = t154 * t208 + t258;
t52 = t53 * t166;
t9 = qJD(5) * t215 + t52;
t185 = (t52 + ((t27 + t282 + t304) * t155 + (-t26 + (t224 - t268) * t155 + t25 + t281) * t154) * qJD(5)) * t229 + (-qJD(5) * t208 + t126 * t176 + t127 * t174) * t166 + (t46 + t54) * t291 + (t45 + t53) * t290 + (-t266 + ((t25 + (-t74 + t268) * t155 - t281) * t155 + (t154 * t224 - t24 + t282) * t154) * qJD(5) + t10) * t227 + (t16 + t21) * t226 + (Icges(5,3) + t209) * t163 + (t17 + t20 + t9) * t228;
t184 = Icges(4,3) * t171 + t185;
t95 = t166 * t216;
t59 = -t95 + t189;
t182 = t60 * (t238 + t81) - t59 * (t82 - t237);
t181 = ((-t272 * t31 - t289 * t30) * t154 + (-t230 * t30 + t289 * t31) * t155) * t166;
t179 = t31 * (t223 + t238) + t181 - t30 * (t231 - t237);
t72 = -t301 - t247;
t34 = qJD(5) * t190 + qJD(2);
t11 = qJDD(2) + t88 * t79 + t87 * t207 + (-t154 * t48 + t155 * t47) * qJD(5);
t6 = t154 * t299 + t194 * t155;
t5 = t154 * t300 + t193 * t155;
t4 = t194 * t154 - t155 * t299;
t3 = t193 * t154 - t155 * t300;
t1 = [t184 + (Icges(2,3) + Icges(3,3)) * qJDD(1) + ((qJDD(1) * t110 + t218 * t178 - g(2) + t162) * t110 + (t246 * pkin(1) + t198 * qJDD(1) + g(3) + (-0.2e1 * t157 + t276 - t220 + t110) * t178) * t198) * m(3) + ((qJDD(1) * t219 + g(3)) * t219 + (qJDD(1) * t148 + g(2)) * t148) * m(2) + (-(t31 + t187 - t235) * t30 + (t221 * t30 + t222 * t31) * qJD(1) + t179 + t308 * (-t221 + t56) + t306 * (-t222 + t55)) * m(6) + ((t221 * t59 + t222 * t60) * qJD(1) + t182 + t307 * (-t221 + t86) + t305 * (-t222 + t85)) * m(5) + (t312 * (-t114 - t221) + t311 * (-t217 - t222) + (t172 * t275 - t137 - t247) * t71) * m(4); m(6) * t11 + t239 * qJDD(2) + (-m(6) - t239) * g(1); t184 + (-t31 * (t206 + t238) + t30 * (t235 - t237) + t179 + t308 * t56 + t306 * t55) * m(6) + (-t60 * (-t95 + t238) + t59 * (-t96 - t237) + t182 + t307 * t86 + t305 * t85) * m(5) + (t101 * t72 - t102 * t71 + (-t172 * t71 - t312) * t114 - (t172 * t72 + t311) * t217) * m(4); t185 + (t181 + t308 * t58 + t306 * t57 + (-t206 + t223) * t31 + (-t231 + t235) * t30) * m(6) + (-t59 * t82 + t60 * t81 + (t166 * t60 + t305) * t216 + (-t166 * t59 - t307) * t106) * m(5); t163 * (t46 * t154 + t45 * t155) / 0.2e1 + t166 * ((t166 * t46 + t16) * t155 + (-t166 * t45 + t17) * t154) / 0.2e1 - t9 * t255 / 0.2e1 + t155 * (t53 * t163 + t21 * t166 + t24 * t88 + t25 * t87 + (t154 * t6 + t155 * t5) * qJD(5)) / 0.2e1 + t215 * t290 + ((t166 * t25 + t5) * t155 + (-t166 * t24 + t6) * t154) * t226 + t10 * t252 / 0.2e1 + t154 * (t54 * t163 + t20 * t166 + t26 * t88 + t27 * t87 + (t154 * t4 + t3 * t155) * qJD(5)) / 0.2e1 + t214 * t291 + ((t166 * t27 + t3) * t155 + (-t166 * t26 + t4) * t154) * t228 - t166 * ((-t244 * t174 + t243 * t176) * t166 + ((t154 * t277 + t155 * t278) * t176 + (-t154 * t279 - t155 * t280) * t174) * qJD(5)) / 0.2e1 + ((t89 * t241 + t257) * t155 + (t310 + (t293 * t154 + (-t292 - t258) * t155) * qJD(5)) * t154) * t227 + ((-t242 * t258 + t257) * t154 + (-t310 + (t292 * t155 + (-t293 + t89) * t154) * qJD(5)) * t155) * t229 + (t11 * t190 + t34 * ((-t48 - t68) * t154 + (t47 + t67) * t155) + t13 * t97 - t12 * t256 + (t31 * t252 - t30 * t255) * t145 + t302 * t128 - (t256 * t31 - t30 * t97) * t166 - (t34 * (-t154 * t97 - t155 * t256) + t302 * t147) * qJD(5) - g(1) * t147 - g(2) * t97 + g(3) * t256) * m(6);];
tau = t1;
