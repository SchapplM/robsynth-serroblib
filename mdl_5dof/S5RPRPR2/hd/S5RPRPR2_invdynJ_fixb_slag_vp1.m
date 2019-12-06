% Calculate vector of inverse dynamics joint torques for
% S5RPRPR2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-05 17:50
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR2_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR2_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR2_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:49:20
% EndTime: 2019-12-05 17:49:29
% DurationCPUTime: 4.89s
% Computational Cost: add. (9387->394), mult. (5723->476), div. (0->0), fcn. (4240->10), ass. (0->229)
t189 = qJD(1) ^ 2;
t181 = pkin(9) + qJ(5);
t174 = sin(t181);
t176 = cos(t181);
t139 = rSges(6,1) * t174 + rSges(6,2) * t176;
t182 = qJD(1) + qJD(3);
t183 = qJ(1) + pkin(8);
t178 = qJ(3) + t183;
t169 = sin(t178);
t170 = cos(t178);
t229 = rSges(4,1) * t169 + rSges(4,2) * t170;
t102 = t229 * t182;
t175 = sin(t183);
t304 = pkin(2) * t175;
t187 = sin(qJ(1));
t306 = pkin(1) * t187;
t235 = t304 + t306;
t214 = t235 * qJD(1);
t78 = t214 + t102;
t254 = qJD(5) * t182;
t108 = qJDD(5) * t169 + t170 * t254;
t292 = rSges(6,1) * t176;
t140 = -rSges(6,2) * t174 + t292;
t118 = t140 * qJD(5);
t180 = qJDD(1) + qJDD(3);
t173 = t189 * t306;
t177 = cos(t183);
t188 = cos(qJ(1));
t305 = pkin(1) * t188;
t234 = pkin(2) * t177 + t305;
t198 = -qJDD(1) * t234 + t189 * t304 + t173;
t195 = qJDD(4) * t170 + t198;
t270 = t169 * t182;
t155 = pkin(3) * t270;
t267 = t170 * t182;
t233 = qJ(4) * t267 - t155;
t161 = qJD(4) * t169;
t236 = -0.2e1 * t161 - t233;
t303 = pkin(3) * t170;
t120 = qJ(4) * t169 + t303;
t185 = cos(pkin(9));
t171 = pkin(4) * t185 + pkin(3);
t148 = t170 * t171;
t186 = -pkin(7) - qJ(4);
t263 = qJ(4) + t186;
t74 = t169 * t263 - t148 + t303;
t269 = t170 * t174;
t150 = rSges(6,2) * t269;
t268 = t170 * t176;
t250 = rSges(6,1) * t268;
t216 = rSges(6,3) * t169 + t250;
t88 = -t150 + t216;
t252 = -t120 + t74 - t88;
t256 = qJD(5) * t169;
t265 = t182 * t186;
t261 = t170 * t265 + t171 * t270;
t271 = t169 * t176;
t251 = rSges(6,1) * t271;
t326 = t139 * qJD(5);
t247 = -t170 * t326 - t182 * t251;
t272 = t169 * t174;
t323 = rSges(6,2) * t272 + t170 * rSges(6,3);
t55 = t182 * t323 + t247;
t8 = t118 * t256 + t108 * t139 + t252 * t180 + (t233 + t236 - t55 + t261) * t182 + t195;
t333 = -g(2) + t8;
t109 = qJDD(5) * t170 - t169 * t254;
t142 = t169 * t265;
t162 = qJD(4) * t170;
t264 = t188 * t189;
t193 = (-qJDD(1) * t175 - t177 * t189) * pkin(2) + (-qJDD(1) * t187 - t264) * pkin(1);
t227 = -pkin(3) * t169 + qJ(4) * t170;
t257 = -t182 * t120 + t162;
t192 = qJDD(4) * t169 + t180 * t227 + t193 + (t162 + t257) * t182;
t237 = -t171 - t292;
t255 = qJD(5) * t170;
t324 = t263 * t170;
t246 = t182 * t150 + t169 * t326;
t56 = -t182 * t216 + t246;
t7 = -t109 * t139 - t118 * t255 + (t142 + t56 + (t120 - t148) * t182) * t182 + t192 + (t323 + (pkin(3) + t237) * t169 - t324) * t180;
t332 = -g(3) + t7;
t294 = rSges(5,1) * t185;
t133 = t270 * t294;
t184 = sin(pkin(9));
t291 = rSges(5,2) * t184;
t156 = t170 * t291;
t217 = -rSges(5,3) * t169 - t170 * t294;
t91 = -t156 - t217;
t282 = -t120 - t91;
t289 = rSges(5,3) * t170;
t24 = t282 * t180 + (t133 + (-t169 * t291 - t289) * t182 + t236) * t182 + t195;
t331 = -g(2) + t24;
t295 = rSges(4,1) * t170;
t121 = -t169 * rSges(4,2) + t295;
t330 = t102 * t182 - t121 * t180 - g(2) + t198;
t134 = t182 * t156;
t228 = t291 - t294;
t197 = t169 * t228 + t289;
t23 = t180 * t197 + (t182 * t217 + t134) * t182 + t192;
t329 = -g(3) + t23;
t154 = rSges(4,2) * t270;
t103 = -rSges(4,1) * t267 + t154;
t328 = t103 * t182 - t180 * t229 - g(3) + t193;
t167 = Icges(6,4) * t176;
t221 = -Icges(6,2) * t174 + t167;
t322 = Icges(6,1) * t174 + t167;
t259 = t322 + t221;
t280 = Icges(6,4) * t174;
t129 = Icges(6,2) * t176 + t280;
t132 = Icges(6,1) * t176 - t280;
t260 = t129 - t132;
t327 = (t174 * t259 + t176 * t260) * t182;
t320 = t234 * qJD(1);
t83 = Icges(6,4) * t268 - Icges(6,2) * t269 + Icges(6,6) * t169;
t147 = Icges(6,4) * t269;
t85 = Icges(6,1) * t268 + Icges(6,5) * t169 - t147;
t223 = t174 * t83 - t176 * t85;
t325 = t223 * t170;
t168 = t175 * rSges(3,2);
t296 = rSges(3,1) * t177;
t232 = -t296 - t305;
t114 = t168 + t232;
t230 = rSges(3,1) * t175 + rSges(3,2) * t177;
t208 = t230 + t306;
t262 = -t182 * t227 - t161;
t199 = t214 + t262;
t218 = t251 - t323;
t75 = t182 * t218;
t319 = t139 * t255 + t75 - t182 * (-t324 + (pkin(3) - t171) * t169);
t27 = t199 + t319;
t106 = t139 * t256;
t204 = t162 - t320;
t28 = t182 * t252 + t106 + t204;
t321 = t169 * t28 + t170 * t27;
t76 = t182 * t88;
t318 = t182 * t74 + t106 - t76;
t82 = Icges(6,6) * t170 - t169 * t221;
t146 = Icges(6,4) * t272;
t84 = -Icges(6,1) * t271 + Icges(6,5) * t170 + t146;
t40 = t174 * t84 + t176 * t82;
t210 = t221 * t182;
t314 = -Icges(6,6) * t182 + qJD(5) * t129;
t52 = t169 * t314 - t170 * t210;
t211 = t132 * t182;
t312 = -Icges(6,5) * t182 + qJD(5) * t322;
t54 = t169 * t312 - t170 * t211;
t128 = Icges(6,5) * t176 - Icges(6,6) * t174;
t80 = Icges(6,3) * t170 - t128 * t169;
t317 = qJD(5) * t40 + t174 * t52 - t176 * t54 - t182 * t80;
t41 = t174 * t85 + t176 * t83;
t51 = -t169 * t210 - t170 * t314;
t53 = -t169 * t211 - t170 * t312;
t81 = Icges(6,5) * t268 - Icges(6,6) * t269 + Icges(6,3) * t169;
t316 = qJD(5) * t41 + t174 * t51 - t176 * t53 - t182 * t81;
t127 = Icges(6,5) * t174 + Icges(6,6) * t176;
t315 = -Icges(6,3) * t182 + qJD(5) * t127;
t116 = t221 * qJD(5);
t117 = t132 * qJD(5);
t220 = t129 * t176 + t174 * t322;
t313 = qJD(5) * t220 + t116 * t174 - t117 * t176 - t127 * t182;
t297 = -Icges(6,2) * t268 - t147 + t85;
t299 = t170 * t322 + t83;
t310 = t174 * t297 + t176 * t299;
t298 = Icges(6,2) * t271 + t146 + t84;
t300 = -t169 * t322 + t82;
t309 = -t174 * t298 - t176 * t300;
t308 = t108 / 0.2e1;
t307 = t109 / 0.2e1;
t302 = t170 * t80 + t82 * t272;
t301 = t169 * t80 + t84 * t268;
t286 = t174 * t82;
t285 = t176 * t84;
t219 = t174 * t129 - t176 * t322;
t94 = t127 * t169;
t48 = -t170 * t219 + t94;
t284 = t48 * t182;
t283 = rSges(5,3) + qJ(4);
t275 = t127 * t170;
t274 = t128 * t182;
t100 = t139 * t169;
t273 = t139 * t170;
t266 = t182 * t121;
t253 = m(3) + m(4) + m(5);
t244 = -pkin(3) - t294;
t243 = -t256 / 0.2e1;
t242 = t256 / 0.2e1;
t241 = -t255 / 0.2e1;
t240 = t255 / 0.2e1;
t238 = -t81 - t285;
t158 = rSges(2,1) * t188 - t187 * rSges(2,2);
t231 = rSges(2,1) * t187 + rSges(2,2) * t188;
t29 = -t271 * t84 + t302;
t30 = t170 * t81 - t271 * t85 + t83 * t272;
t226 = t30 * t169 + t29 * t170;
t31 = -t269 * t82 + t301;
t32 = t169 * t81 - t325;
t225 = t32 * t169 + t31 * t170;
t224 = -t285 + t286;
t203 = -t274 * t169 - t170 * t315 + t182 * t223;
t202 = t169 * t315 - t274 * t170 + t182 * t224;
t201 = t128 * qJD(5) + t182 * t219;
t200 = t320 - t257;
t62 = -t250 - t148 + t150 + (-rSges(6,3) + t186) * t169;
t196 = t169 * t218 + t170 * t88;
t61 = t169 * t237 - t170 * t186 + t323;
t67 = -t169 * t283 + t170 * t244 + t156;
t66 = t283 * t170 + (-pkin(3) + t228) * t169;
t47 = t169 * t219 + t275;
t42 = t47 * t182;
t11 = qJD(5) * t226 + t42;
t12 = qJD(5) * t225 + t284;
t16 = -qJD(5) * t224 + t174 * t54 + t176 * t52;
t17 = -qJD(5) * t223 + t174 * t53 + t176 * t51;
t20 = t201 * t169 - t170 * t313;
t21 = t169 * t313 + t201 * t170;
t194 = (t42 + ((t302 + t32 + t325) * t170 + (-t31 + (t238 - t286) * t170 + t30 + t301) * t169) * qJD(5)) * t243 + (-qJD(5) * t219 + t116 * t176 + t117 * t174) * t182 + (t41 + t48) * t308 + (t40 + t47) * t307 + (t12 - t284 + ((t30 + (-t81 + t286) * t170 - t301) * t170 + (t169 * t238 - t29 + t302) * t169) * qJD(5)) * t241 + (t16 + t21) * t240 + (t17 + t20 + t11) * t242 + (Icges(4,3) + t220 + Icges(5,2) * t185 ^ 2 + (Icges(5,1) * t184 + 0.2e1 * Icges(5,4) * t185) * t184) * t180;
t86 = t182 * t197;
t43 = -t86 + t199;
t44 = t182 * t282 + t204;
t191 = t44 * (t133 + t155 - t161) + ((t283 * t43 - t291 * t44) * t169 + (-t244 * t43 - t283 * t44) * t170) * t182 - t43 * (t134 + t162);
t190 = t28 * (-t161 - t247 + t261) + (-t28 * t323 - t27 * (-t216 - t148)) * t182 - t27 * (t142 + t162 + t246);
t87 = t182 * t91;
t79 = -t320 - t266;
t39 = qJD(5) * t196 + qJD(2);
t13 = qJDD(2) + t109 * t88 + t108 * t218 + (-t169 * t56 + t170 * t55) * qJD(5);
t6 = t169 * t316 + t203 * t170;
t5 = t169 * t317 + t202 * t170;
t4 = t203 * t169 - t170 * t316;
t3 = t202 * t169 - t170 * t317;
t1 = [t194 + (Icges(2,3) + Icges(3,3)) * qJDD(1) + ((qJDD(1) * t114 + t189 * t230 - g(2) + t173) * t114 + (t264 * pkin(1) + t208 * qJDD(1) + g(3) + (-0.2e1 * t168 + t296 - t232 + t114) * t189) * t208) * m(3) + ((qJDD(1) * t231 + g(3)) * t231 + (qJDD(1) * t158 + g(2)) * t158) * m(2) + (-(t28 + t200 - t318) * t27 + (t234 * t27 + t235 * t28) * qJD(1) + t190 + t333 * (-t234 + t62) + t332 * (-t235 + t61)) * m(6) + (-(t44 + t87 + t200) * t43 + (t234 * t43 + t235 * t44) * qJD(1) + t191 + t331 * (-t234 + t67) + t329 * (-t235 + t66)) * m(5) + (t330 * (-t121 - t234) + t328 * (-t229 - t235) + (t295 * t182 - t154 - t266) * t78) * m(4); m(6) * t13 + t253 * qJDD(2) + (-m(6) - t253) * g(1); t194 + (t190 - t28 * (t262 + t319) + t27 * (t257 + t318) + t333 * t62 + t332 * t61) * m(6) + (t191 - t44 * (-t86 + t262) + t43 * (-t87 + t257) + t331 * t67 + t329 * t66) * m(5) + (t102 * t79 - t103 * t78 + (-t182 * t78 - t330) * t121 - (t182 * t79 + t328) * t229) * m(4); (-m(5) - m(6)) * (g(2) * t170 + g(3) * t169) + m(5) * (t169 * t23 + t170 * t24) + m(6) * (t169 * t7 + t170 * t8); t180 * (t41 * t169 + t40 * t170) / 0.2e1 + t182 * ((t182 * t41 + t16) * t170 + (-t182 * t40 + t17) * t169) / 0.2e1 - t11 * t270 / 0.2e1 + t170 * (t30 * t108 + t29 * t109 + t47 * t180 + t21 * t182 + (t169 * t6 + t170 * t5) * qJD(5)) / 0.2e1 + t226 * t307 + ((t182 * t30 + t5) * t170 + (-t182 * t29 + t6) * t169) * t240 + t12 * t267 / 0.2e1 + t169 * (t108 * t32 + t109 * t31 + t48 * t180 + t20 * t182 + (t169 * t4 + t170 * t3) * qJD(5)) / 0.2e1 + t225 * t308 + ((t182 * t32 + t3) * t170 + (-t182 * t31 + t4) * t169) * t242 - t182 * ((-t260 * t174 + t259 * t176) * t182 + ((t169 * t297 + t170 * t298) * t176 + (-t169 * t299 - t170 * t300) * t174) * qJD(5)) / 0.2e1 + ((t94 * t255 + t274) * t170 + (t327 + (t310 * t169 + (-t309 - t275) * t170) * qJD(5)) * t169) * t241 + ((-t256 * t275 + t274) * t169 + (-t327 + (t309 * t170 + (-t310 + t94) * t169) * qJD(5)) * t170) * t243 + (t13 * t196 + t39 * ((-t56 - t76) * t169 + (t55 + t75) * t170) + t8 * t100 - t7 * t273 + (t28 * t267 - t27 * t270) * t139 + t321 * t118 - (-t100 * t27 + t273 * t28) * t182 - (t39 * (-t100 * t169 - t170 * t273) + t321 * t140) * qJD(5) - g(1) * t140 - g(2) * t100 + g(3) * t273) * m(6);];
tau = t1;
