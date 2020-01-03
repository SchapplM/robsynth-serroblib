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
% Datum: 2020-01-03 11:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:52:02
% EndTime: 2020-01-03 11:52:10
% DurationCPUTime: 4.62s
% Computational Cost: add. (10808->367), mult. (6100->460), div. (0->0), fcn. (4528->10), ass. (0->217)
t189 = qJ(1) + pkin(9);
t178 = sin(t189);
t179 = cos(t189);
t112 = rSges(3,1) * t179 - rSges(3,2) * t178;
t193 = cos(qJ(1));
t185 = t193 * pkin(1);
t332 = t112 + t185;
t181 = qJ(3) + t189;
t169 = sin(t181);
t188 = qJD(1) + qJD(3);
t273 = t169 * t188;
t254 = pkin(3) * t273;
t191 = sin(qJ(1));
t184 = t191 * pkin(1);
t316 = pkin(2) * t178 + t184;
t312 = t316 * qJD(1);
t202 = t312 + t254;
t180 = qJD(4) + t188;
t171 = qJ(4) + t181;
t161 = sin(t171);
t162 = cos(t171);
t318 = rSges(5,1) * t161 + rSges(5,2) * t162;
t329 = t318 * t180;
t59 = t202 + t329;
t170 = cos(t181);
t317 = rSges(4,1) * t169 + rSges(4,2) * t170;
t96 = t317 * t188;
t331 = -t312 - t96;
t190 = sin(qJ(5));
t192 = cos(qJ(5));
t182 = Icges(6,4) * t192;
t226 = -Icges(6,2) * t190 + t182;
t315 = Icges(6,1) * t190 + t182;
t267 = t315 + t226;
t288 = Icges(6,4) * t190;
t137 = Icges(6,2) * t192 + t288;
t140 = Icges(6,1) * t192 - t288;
t268 = t137 - t140;
t330 = (t190 * t267 + t192 * t268) * t180;
t141 = rSges(6,1) * t190 + rSges(6,2) * t192;
t259 = qJD(5) * t162;
t149 = t161 * pkin(4);
t278 = t161 * t190;
t252 = rSges(6,2) * t278;
t277 = t161 * t192;
t240 = rSges(6,1) * t277 - t252;
t75 = -rSges(6,3) * t162 + t240;
t290 = -pkin(8) * t162 + t149 + t75;
t210 = -t141 * t259 - t180 * t290;
t30 = t202 - t210;
t260 = qJD(5) * t161;
t248 = t141 * t260;
t272 = t170 * t188;
t134 = pkin(3) * t272;
t261 = qJD(1) * t179;
t295 = pkin(1) * qJD(1);
t266 = pkin(2) * t261 + t193 * t295;
t249 = t134 + t266;
t215 = -t248 + t249;
t104 = pkin(4) * t162 + pkin(8) * t161;
t274 = t162 * t192;
t253 = rSges(6,1) * t274;
t275 = t162 * t190;
t76 = -rSges(6,2) * t275 + rSges(6,3) * t161 + t253;
t58 = t104 + t76;
t31 = t180 * t58 + t215;
t93 = t141 * t161;
t94 = t141 * t162;
t327 = -t30 * t93 - t31 * t94;
t143 = rSges(6,1) * t192 - rSges(6,2) * t190;
t122 = t143 * qJD(5);
t187 = qJDD(1) + qJDD(3);
t177 = qJDD(4) + t187;
t283 = pkin(1) * qJDD(1);
t173 = t193 * t283;
t194 = qJD(1) ^ 2;
t282 = pkin(2) * qJDD(1);
t204 = t179 * t282 - t194 * t316 + t173;
t301 = pkin(3) * t187;
t302 = pkin(3) * t188 ^ 2;
t199 = -t169 * t302 + t170 * t301 + t204;
t276 = t162 * t180;
t279 = t161 * t180;
t241 = -pkin(4) * t279 + pkin(8) * t276;
t251 = t180 * t277;
t256 = qJD(5) * t192;
t257 = qJD(5) * t190;
t271 = -rSges(6,3) * t276 - t180 * t252;
t47 = rSges(6,2) * t162 * t256 + (t162 * t257 + t251) * rSges(6,1) + t271;
t258 = qJD(5) * t180;
t84 = -qJDD(5) * t161 - t162 * t258;
t13 = -t122 * t260 + t84 * t141 + (t241 - t47) * t180 + t58 * t177 + t199;
t326 = -g(2) + t13;
t102 = rSges(5,1) * t162 - rSges(5,2) * t161;
t325 = t102 * t177 - t180 * t329 - g(2) + t199;
t110 = rSges(4,1) * t170 - rSges(4,2) * t169;
t324 = t110 * t187 - t188 * t96 - g(2) + t204;
t168 = pkin(2) * t179;
t263 = t194 * t185 + t191 * t283;
t242 = t168 * t194 + t178 * t282 + t263;
t221 = t169 * t301 + t170 * t302 + t242;
t269 = pkin(4) * t276 + pkin(8) * t279;
t250 = t180 * t275;
t270 = rSges(6,3) * t279 + t180 * t253;
t48 = -rSges(6,1) * t161 * t257 + (-t161 * t256 - t250) * rSges(6,2) + t270;
t85 = -qJDD(5) * t162 + t161 * t258;
t12 = t122 * t259 - t85 * t141 + (t48 + t269) * t180 + t290 * t177 + t221;
t323 = -g(3) + t12;
t79 = rSges(5,1) * t276 - rSges(5,2) * t279;
t322 = t177 * t318 + t180 * t79 - g(3) + t221;
t97 = rSges(4,1) * t272 - rSges(4,2) * t273;
t321 = t187 * t317 + t188 * t97 - g(3) + t242;
t64 = t180 * t76;
t319 = t104 * t180 + t64;
t92 = t180 * t102;
t60 = t249 + t92;
t111 = rSges(3,1) * t178 + rSges(3,2) * t179;
t217 = t226 * t180;
t309 = -Icges(6,6) * t180 + qJD(5) * t137;
t42 = -t161 * t309 + t162 * t217;
t218 = t140 * t180;
t306 = -Icges(6,5) * t180 + qJD(5) * t315;
t44 = -t161 * t306 + t162 * t218;
t71 = -Icges(6,6) * t162 + t161 * t226;
t73 = -Icges(6,5) * t162 + t140 * t161;
t45 = t190 * t73 + t192 * t71;
t136 = Icges(6,5) * t192 - Icges(6,6) * t190;
t69 = -Icges(6,3) * t162 + t136 * t161;
t311 = qJD(5) * t45 - t180 * t69 + t190 * t42 - t192 * t44;
t135 = Icges(6,5) * t190 + Icges(6,6) * t192;
t310 = -Icges(6,3) * t180 + qJD(5) * t135;
t120 = t226 * qJD(5);
t121 = t140 * qJD(5);
t225 = t137 * t192 + t190 * t315;
t308 = qJD(5) * t225 + t120 * t190 - t121 * t192 - t135 * t180;
t72 = Icges(6,4) * t274 - Icges(6,2) * t275 + Icges(6,6) * t161;
t127 = Icges(6,4) * t275;
t74 = Icges(6,1) * t274 + Icges(6,5) * t161 - t127;
t229 = t190 * t74 + t192 * t72;
t41 = t161 * t217 + t162 * t309;
t43 = t161 * t218 + t162 * t306;
t70 = Icges(6,5) * t274 - Icges(6,6) * t275 + Icges(6,3) * t161;
t307 = qJD(5) * t229 - t180 * t70 - t190 * t41 + t192 * t43;
t304 = t84 / 0.2e1;
t303 = t85 / 0.2e1;
t299 = t161 * t315 + t71;
t298 = t162 * t315 + t72;
t297 = -t137 * t161 + t73;
t296 = -Icges(6,2) * t274 - t127 + t74;
t294 = t190 * t71;
t293 = t190 * t72;
t292 = t192 * t73;
t280 = t135 * t161;
t54 = -t137 * t275 + t274 * t315 + t280;
t291 = t54 * t180;
t87 = t135 * t162;
t216 = t136 * t180;
t264 = t168 + t185;
t262 = qJD(1) * t178;
t255 = m(3) + m(4) + m(5);
t159 = pkin(3) * t169;
t82 = t159 + t318;
t246 = -t260 / 0.2e1;
t245 = t260 / 0.2e1;
t244 = -t259 / 0.2e1;
t243 = t259 / 0.2e1;
t160 = pkin(3) * t170;
t83 = t102 + t160;
t68 = t110 * t188 + t266;
t144 = rSges(2,1) * t193 - rSges(2,2) * t191;
t142 = rSges(2,1) * t191 + rSges(2,2) * t193;
t61 = t73 * t277;
t24 = -t162 * t69 - t278 * t71 + t61;
t62 = t74 * t277;
t25 = t162 * t70 + t278 * t72 - t62;
t234 = -t161 * t25 - t162 * t24;
t63 = t71 * t275;
t26 = -t161 * t69 - t274 * t73 + t63;
t228 = -t192 * t74 + t293;
t27 = t161 * t70 - t228 * t162;
t233 = -t27 * t161 - t162 * t26;
t232 = -t161 * t31 + t162 * t30;
t231 = t161 * t75 + t162 * t76;
t230 = t292 - t294;
t224 = -t137 * t190 + t192 * t315;
t223 = t134 + t79;
t222 = -t248 + t319;
t212 = t190 * t297 + t192 * t299;
t211 = t190 * t296 + t192 * t298;
t208 = -t161 * t216 - t162 * t310 + t180 * t228;
t207 = t161 * t310 - t162 * t216 + t180 * t230;
t206 = -qJD(5) * t136 + t224 * t180;
t205 = -rSges(6,2) * t250 + t269 + t270;
t57 = t149 + (-rSges(6,3) - pkin(8)) * t162 + t240;
t203 = t134 + t205;
t56 = t160 + t58;
t55 = t159 + t57;
t200 = -rSges(6,1) * t251 + t241 - t271;
t10 = qJD(5) * t233 - t291;
t16 = qJD(5) * t230 + t190 * t44 + t192 * t42;
t17 = qJD(5) * t228 + t190 * t43 + t192 * t41;
t20 = t161 * t206 + t162 * t308;
t21 = -t161 * t308 + t162 * t206;
t53 = t161 * t224 - t87;
t52 = t53 * t180;
t9 = qJD(5) * t234 + t52;
t198 = (t52 + ((t26 + t62 - t63 + (t69 - t293) * t161) * t161 + (-t61 - t27 + (-t228 + t69) * t162 + (t292 + t294) * t161) * t162) * qJD(5)) * t245 - t229 * t304 - t84 * t54 / 0.2e1 + (t224 * qJD(5) + t120 * t192 + t121 * t190) * t180 + (t45 + t53) * t303 + (t291 + ((-t25 + t63 + (t70 - t292) * t162) * t162 + (t24 - t61 + (t70 + t294) * t161) * t161) * qJD(5) + t10) * t243 + (Icges(5,3) + t225) * t177 + (t17 + t20 + t9) * t246 + (t16 + t21) * t244;
t197 = t200 - t254;
t196 = Icges(4,3) * t187 + t198;
t195 = t327 * qJD(5);
t34 = qJD(5) * t231 + qJD(2);
t11 = -t75 * t84 - t76 * t85 + qJDD(2) + (t161 * t48 - t162 * t47) * qJD(5);
t6 = t161 * t307 + t162 * t208;
t5 = -t161 * t311 + t162 * t207;
t4 = t161 * t208 - t162 * t307;
t3 = t161 * t207 + t162 * t311;
t1 = [t196 + (Icges(2,3) + Icges(3,3)) * qJDD(1) + (t324 * (t110 + t264) + t321 * (t316 + t317) + (t68 - t97 - t266) * t331) * m(4) + ((qJDD(1) * t112 - g(2) + t173) * t332 + (-t194 * t332 + qJDD(1) * t111 - g(3) + t263 + (0.2e1 * rSges(3,1) * t261 - 0.2e1 * rSges(3,2) * t262 - qJD(1) * t112) * qJD(1)) * (t184 + t111)) * m(3) + (-g(2) * t144 - g(3) * t142 + (t142 ^ 2 + t144 ^ 2) * qJDD(1)) * m(2) + (t31 * (-pkin(2) * t262 - t191 * t295 + t197) + t195 + t326 * (t56 + t264) + t323 * (t55 + t316) + (t203 + t266 - t215 + t31 - t319) * t30) * m(6) + (t325 * (t83 + t264) + t322 * (t82 + t316) + (-t60 + t223 + t266) * t59) * m(5); m(6) * t11 + t255 * qJDD(2) + (-m(6) - t255) * g(1); t196 + (t195 + t326 * t56 + t323 * t55 + (t197 - t210 + t254) * t31 + (t203 - t134 - t222) * t30) * m(6) + (t325 * t83 + t322 * t82 + (-t134 - t92 + t223) * t59) * m(5) + (-t331 * t97 - t68 * t96 + (t188 * t331 + t324) * t110 + (t188 * t68 + t321) * t317) * m(4); t198 + (t195 + t326 * t58 + t323 * t57 + (t200 - t210) * t31 + (t205 - t222) * t30) * m(6) + (t59 * t79 - t60 * t329 + (-t180 * t59 + t325) * t102 + (t180 * t60 + t322) * t318) * m(5); t177 * (t161 * t229 - t162 * t45) / 0.2e1 + t180 * ((t180 * t229 - t16) * t162 + (t180 * t45 - t17) * t161) / 0.2e1 + t9 * t279 / 0.2e1 - t162 * (t53 * t177 + t21 * t180 + t24 * t85 + t25 * t84 + (-t161 * t6 - t162 * t5) * qJD(5)) / 0.2e1 + t234 * t303 + ((-t180 * t25 - t5) * t162 + (t180 * t24 - t6) * t161) * t244 - t10 * t276 / 0.2e1 - t161 * (-t54 * t177 + t20 * t180 + t26 * t85 + t27 * t84 + (-t161 * t4 - t162 * t3) * qJD(5)) / 0.2e1 + t233 * t304 + ((-t180 * t27 - t3) * t162 + (t180 * t26 - t4) * t161) * t246 - t180 * ((-t268 * t190 + t267 * t192) * t180 + ((t161 * t296 - t162 * t297) * t192 + (-t161 * t298 + t162 * t299) * t190) * qJD(5)) / 0.2e1 + ((-t259 * t280 - t216) * t162 + (-t330 + (-t211 * t161 + (t87 + t212) * t162) * qJD(5)) * t161) * t243 + ((t87 * t260 - t216) * t161 + (t330 + (-t212 * t162 + (-t280 + t211) * t161) * qJD(5)) * t162) * t245 + (t11 * t231 + t34 * ((t180 * t75 - t47) * t162 + (t48 - t64) * t161) + t232 * t122 + ((-t180 * t31 + t12) * t162 + (-t180 * t30 - t13) * t161) * t141 - t327 * t180 - (t34 * (-t161 * t93 - t162 * t94) + t232 * t143) * qJD(5) - g(1) * t143 + g(2) * t93 - g(3) * t94) * m(6);];
tau = t1;
