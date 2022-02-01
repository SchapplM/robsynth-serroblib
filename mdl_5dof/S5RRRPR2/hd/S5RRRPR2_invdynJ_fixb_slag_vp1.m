% Calculate vector of inverse dynamics joint torques for
% S5RRRPR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% m [6x1]
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
% Datum: 2022-01-20 11:31
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPR2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR2_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR2_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR2_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR2_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:30:33
% EndTime: 2022-01-20 11:30:40
% DurationCPUTime: 3.74s
% Computational Cost: add. (11392->389), mult. (6350->470), div. (0->0), fcn. (4656->10), ass. (0->235)
t185 = qJ(1) + qJ(2);
t177 = sin(t185);
t178 = cos(t185);
t184 = qJD(1) + qJD(2);
t182 = t184 ^ 2;
t183 = qJDD(1) + qJDD(2);
t187 = sin(qJ(1));
t189 = cos(qJ(1));
t190 = qJD(1) ^ 2;
t216 = (-qJDD(1) * t187 - t189 * t190) * pkin(1);
t199 = t216 + (-t177 * t183 - t178 * t182) * pkin(2);
t320 = t199 - g(1);
t180 = qJ(3) + t185;
t169 = pkin(9) + t180;
t163 = cos(t169);
t159 = t163 * pkin(8);
t162 = sin(t169);
t112 = pkin(4) * t162 - t159;
t176 = qJD(3) + t184;
t266 = t163 * t176;
t122 = pkin(8) * t266;
t186 = sin(qJ(5));
t268 = t162 * t186;
t138 = rSges(6,2) * t268;
t188 = cos(qJ(5));
t284 = rSges(6,2) * t188;
t245 = qJD(5) * t284;
t220 = rSges(6,3) * t266 + t176 * t138 - t163 * t245;
t253 = qJD(5) * t186;
t242 = t163 * t253;
t259 = t163 * rSges(6,3) + t138;
t267 = t162 * t188;
t79 = rSges(6,1) * t267 - t259;
t68 = t176 * t79;
t319 = -rSges(6,1) * t242 + t176 * t112 + t122 + t220 + t68;
t283 = pkin(1) * qJD(1);
t247 = t187 * t283;
t123 = rSges(3,1) * t177 + rSges(3,2) * t178;
t273 = t123 * t184;
t102 = -t247 - t273;
t285 = rSges(6,1) * t188;
t150 = -rSges(6,2) * t186 + t285;
t132 = t150 * qJD(5);
t148 = rSges(6,1) * t186 + t284;
t175 = qJDD(3) + t183;
t170 = sin(t180);
t292 = pkin(3) * t170;
t233 = -t112 - t79 - t292;
t255 = qJD(5) * t163;
t171 = cos(t180);
t174 = t176 ^ 2;
t262 = t171 * t174;
t160 = t163 * pkin(4);
t305 = -t162 * pkin(8) - t160;
t265 = t163 * t186;
t248 = rSges(6,2) * t265;
t244 = t176 * t248 + (rSges(6,1) * t253 + t245) * t162;
t264 = t163 * t188;
t139 = rSges(6,1) * t264;
t306 = t162 * rSges(6,3) + t139;
t48 = t176 * t306 - t244;
t254 = qJD(5) * t176;
t91 = -qJDD(5) * t163 + t162 * t254;
t11 = -pkin(3) * t262 - t132 * t255 + t91 * t148 + (t176 * t305 - t48) * t176 + t233 * t175 + t199;
t318 = -g(1) + t11;
t164 = pkin(3) * t171;
t168 = pkin(2) * t178;
t181 = t189 * pkin(1);
t294 = pkin(1) * t187;
t231 = qJDD(1) * t181 - t190 * t294;
t293 = pkin(2) * t177;
t209 = t183 * t168 - t182 * t293 + t231;
t200 = t175 * t164 - t174 * t292 + t209;
t256 = qJD(5) * t162;
t269 = t162 * t176;
t80 = -t248 + t306;
t278 = t80 - t305;
t47 = (-t176 * t267 - t242) * rSges(6,1) + t220;
t90 = qJDD(5) * t162 + t163 * t254;
t12 = -t132 * t256 - t90 * t148 + (-pkin(4) * t269 + t122 + t47) * t176 + t278 * t175 + t200;
t317 = -g(2) + t12;
t109 = rSges(5,1) * t162 + rSges(5,2) * t163;
t120 = rSges(5,2) * t269;
t316 = -t175 * t109 - t176 * (rSges(5,1) * t266 - t120) + (-t170 * t175 - t262) * pkin(3) + t320;
t157 = t163 * rSges(5,1);
t110 = -rSges(5,2) * t162 + t157;
t315 = -t109 * t174 + t175 * t110 - g(2) + t200;
t118 = rSges(4,1) * t170 + rSges(4,2) * t171;
t161 = t171 * rSges(4,1);
t263 = t170 * t176;
t89 = -rSges(4,2) * t263 + t161 * t176;
t314 = -t118 * t175 - t176 * t89 + t320;
t106 = t176 * t118;
t119 = -rSges(4,2) * t170 + t161;
t313 = -t106 * t176 + t119 * t175 - g(2) + t209;
t260 = t178 * t184;
t261 = t177 * t184;
t108 = rSges(3,1) * t260 - rSges(3,2) * t261;
t312 = -t108 * t184 - t123 * t183 - g(1) + t216;
t124 = t178 * rSges(3,1) - rSges(3,2) * t177;
t311 = t124 * t183 - t184 * t273 - g(2) + t231;
t307 = t110 + t164;
t309 = t176 * t307;
t304 = t164 + t278;
t179 = Icges(6,4) * t188;
t223 = -Icges(6,2) * t186 + t179;
t145 = Icges(6,1) * t186 + t179;
t251 = pkin(2) * t261;
t215 = -t247 - t251;
t303 = t215 + t247;
t302 = -t251 + t319;
t142 = Icges(6,5) * t188 - Icges(6,6) * t186;
t141 = Icges(6,5) * t186 + Icges(6,6) * t188;
t205 = Icges(6,3) * t176 - qJD(5) * t141;
t218 = t223 * t163;
t76 = Icges(6,6) * t162 + t218;
t280 = t186 * t76;
t276 = Icges(6,4) * t186;
t146 = Icges(6,1) * t188 - t276;
t219 = t146 * t163;
t78 = Icges(6,5) * t162 + t219;
t225 = -t188 * t78 + t280;
t301 = -t142 * t269 + t163 * t205 + t176 * t225;
t217 = t142 * t163;
t75 = Icges(6,4) * t267 - Icges(6,2) * t268 - Icges(6,6) * t163;
t281 = t186 * t75;
t137 = Icges(6,4) * t268;
t77 = Icges(6,1) * t267 - Icges(6,5) * t163 - t137;
t226 = -t188 * t77 + t281;
t300 = t162 * t205 + (t217 + t226) * t176;
t143 = Icges(6,2) * t188 + t276;
t221 = t186 * t143 - t188 * t145;
t299 = t142 * qJD(5) + t176 * t221;
t73 = Icges(6,5) * t267 - Icges(6,6) * t268 - Icges(6,3) * t163;
t24 = -t162 * t226 - t163 * t73;
t287 = -Icges(6,2) * t267 - t137 + t77;
t289 = t145 * t162 + t75;
t298 = -t186 * t287 - t188 * t289;
t111 = t148 * t256;
t241 = -pkin(4) - t285;
t246 = t189 * t283;
t250 = pkin(2) * t260;
t214 = t246 + t250;
t31 = t176 * t304 - t111 + t214;
t252 = t31 * t292;
t243 = t148 * t255;
t213 = -t243 - t251;
t203 = t213 - t247;
t30 = t176 * t233 + t203;
t191 = (t30 * (-t139 - t160 - t164) - t252 + (t30 * (-rSges(6,3) - pkin(8)) + t31 * t241) * t162) * t176;
t297 = t191 + (-t111 + t244) * t30 - (-t30 * t304 - t252) * t176;
t296 = t90 / 0.2e1;
t295 = t91 / 0.2e1;
t291 = -t162 * t73 - t77 * t264;
t74 = Icges(6,3) * t162 + t217;
t290 = t162 * t74 + t78 * t264;
t288 = -t145 * t163 - t76;
t286 = -t143 * t163 + t78;
t271 = t141 * t163;
t55 = -t162 * t221 - t271;
t279 = t55 * t176;
t274 = t119 * t176;
t272 = t141 * t162;
t270 = t142 * t176;
t258 = -t143 + t146;
t257 = t145 + t223;
t249 = pkin(3) * t263;
t240 = -t256 / 0.2e1;
t239 = t256 / 0.2e1;
t238 = -t255 / 0.2e1;
t237 = t255 / 0.2e1;
t84 = -t109 - t292;
t63 = t78 * t267;
t235 = t163 * t74 - t63;
t234 = -t73 + t280;
t105 = t119 + t168;
t151 = rSges(2,1) * t189 - rSges(2,2) * t187;
t149 = rSges(2,1) * t187 + rSges(2,2) * t189;
t25 = -t268 * t76 - t235;
t230 = t162 * t25 - t24 * t163;
t26 = -t265 * t75 - t291;
t27 = -t265 * t76 + t290;
t229 = t27 * t162 - t26 * t163;
t228 = -t31 * t162 - t30 * t163;
t227 = t162 * t79 + t163 * t80;
t45 = t186 * t77 + t188 * t75;
t46 = t186 * t78 + t188 * t76;
t222 = t143 * t188 + t145 * t186;
t83 = t168 + t307;
t104 = -t118 - t293;
t98 = t176 * t109;
t212 = -t249 - t98 - t251;
t211 = -t186 * t286 + t188 * t288;
t210 = -t89 - t250;
t82 = t84 - t293;
t54 = t168 + t304;
t208 = (-t186 * t257 + t188 * t258) * t176;
t207 = Icges(6,5) * t176 - qJD(5) * t145;
t206 = Icges(6,6) * t176 - qJD(5) * t143;
t57 = t162 * t241 + t159 + t259 - t292;
t53 = t57 - t293;
t56 = -t163 * t221 + t272;
t52 = t56 * t176;
t10 = qJD(5) * t229 + t52;
t129 = t223 * qJD(5);
t130 = t146 * qJD(5);
t42 = t162 * t206 + t176 * t218;
t44 = t162 * t207 + t176 * t219;
t16 = -qJD(5) * t226 + t186 * t44 + t188 * t42;
t41 = t163 * t206 - t223 * t269;
t43 = -t146 * t269 + t163 * t207;
t17 = -qJD(5) * t225 + t186 * t43 + t188 * t41;
t195 = -qJD(5) * t222 - t129 * t186 + t130 * t188 + t141 * t176;
t20 = t162 * t299 + t195 * t163;
t21 = t195 * t162 - t163 * t299;
t9 = qJD(5) * t230 + t279;
t198 = (t52 + ((t25 - t63 + (t74 + t281) * t163 + t291) * t163 + t290 * t162) * qJD(5)) * t237 + (-qJD(5) * t221 + t129 * t188 + t130 * t186) * t176 + (t46 + t56) * t296 + (t45 + t55) * t295 + (-t279 + ((t163 * t234 + t27 - t290) * t163 + (t162 * t234 + t235 + t26) * t162) * qJD(5) + t9) * t240 + (t17 + t20) * t239 + (t16 + t21 + t10) * t238 + (Icges(5,3) + Icges(4,3) + t222) * t175;
t197 = -qJD(5) * t45 + t176 * t73 - t186 * t42 + t188 * t44;
t196 = -qJD(5) * t46 + t176 * t74 - t186 * t41 + t188 * t43;
t194 = Icges(3,3) * t183 + t198;
t61 = t176 * t84 + t215;
t62 = t214 + t309;
t192 = (t61 * (-t157 - t164) + t62 * t84) * t176;
t103 = t124 * t184 + t246;
t100 = t148 * t163;
t99 = t148 * t162;
t70 = t214 + t274;
t69 = t215 - t106;
t34 = qJD(5) * t227 + qJD(4);
t13 = t79 * t90 - t80 * t91 + qJDD(4) + (t162 * t48 + t163 * t47) * qJD(5);
t6 = t196 * t162 - t163 * t301;
t5 = t197 * t162 - t163 * t300;
t4 = t162 * t301 + t196 * t163;
t3 = t162 * t300 + t197 * t163;
t1 = [Icges(2,3) * qJDD(1) + t194 + (t311 * (t124 + t181) + t312 * (-t123 - t294) + (-t108 - t246 + t103) * t102) * m(3) + ((t149 ^ 2 + t151 ^ 2) * qJDD(1) + g(1) * t149 - g(2) * t151) * m(2) + (t30 * (-t214 + t244) + t191 + t317 * (t181 + t54) + t318 * (t53 - t294) + (-t203 + t30 + t249 - t247 + t302) * t31) * m(6) + (t61 * (t120 - t214) + t192 + t315 * (t181 + t83) + t316 * (t82 - t294) + (-t212 + t61 + t303) * t62) * m(5) + (t69 * (t210 - t246) + t313 * (t105 + t181) + t314 * (t104 - t294) + (t69 - t303 - t251) * t70) * m(4); t194 + (t317 * t54 + t318 * t53 + (-t213 + t302) * t31 + t297) * m(6) + (t192 + t315 * t83 + t316 * t82 + (-t212 - t251) * t62 + (t120 + t309) * t61) * m(5) + ((t250 + t274 + t210) * t69 + t313 * t105 + t314 * t104) * m(4) + (-t102 * t108 - t103 * t273 + (t102 * t184 + t311) * t124 + (t103 * t184 - t312) * t123) * m(3); t198 + (t317 * t304 + t318 * t57 + (t243 + t319) * t31 + t297) * m(6) + (t62 * t98 - (-t292 * t62 - t307 * t61) * t176 + t61 * t120 + t192 + t315 * t307 + t316 * t84) * m(5) + (-t69 * t89 - t70 * t106 + (t176 * t69 + t313) * t119 + (t176 * t70 - t314) * t118) * m(4); (t13 - g(3)) * m(6) + (qJDD(4) - g(3)) * m(5); t10 * t266 / 0.2e1 + t162 * (t56 * t175 + t20 * t176 + t26 * t91 + t27 * t90 + (t162 * t4 - t163 * t3) * qJD(5)) / 0.2e1 + t229 * t296 + ((t176 * t27 - t3) * t163 + (t176 * t26 + t4) * t162) * t239 + t9 * t269 / 0.2e1 - t163 * (t55 * t175 + t21 * t176 + t24 * t91 + t25 * t90 + (t162 * t6 - t163 * t5) * qJD(5)) / 0.2e1 + t230 * t295 + ((t176 * t25 - t5) * t163 + (t176 * t24 + t6) * t162) * t238 + t175 * (t162 * t46 - t45 * t163) / 0.2e1 + t176 * ((t176 * t46 - t16) * t163 + (t176 * t45 + t17) * t162) / 0.2e1 + ((-t256 * t271 + t270) * t162 + (t208 + (-t298 * t163 + (t272 + t211) * t162) * qJD(5)) * t163) * t240 + ((-t255 * t272 - t270) * t163 + (t208 + (t211 * t162 + (-t298 + t271) * t163) * qJD(5)) * t162) * t237 - t176 * ((t258 * t186 + t257 * t188) * t176 + ((t162 * t286 - t163 * t287) * t188 + (t162 * t288 + t163 * t289) * t186) * qJD(5)) / 0.2e1 + (t13 * t227 + t34 * ((t47 + t68) * t163 + (-t176 * t80 + t48) * t162) + t228 * t132 + ((-t176 * t31 - t11) * t163 + (t176 * t30 - t12) * t162) * t148 - (-t100 * t31 + t30 * t99) * t176 - (t34 * (-t100 * t163 - t162 * t99) + t228 * t150) * qJD(5) + g(1) * t100 + g(2) * t99 - g(3) * t150) * m(6);];
tau = t1;
