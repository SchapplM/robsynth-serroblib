% Calculate vector of inverse dynamics joint torques for
% S5RPPRR5
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 17:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR5_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR5_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR5_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR5_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR5_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR5_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR5_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR5_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR5_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:56:29
% EndTime: 2019-12-31 17:56:37
% DurationCPUTime: 6.34s
% Computational Cost: add. (9858->381), mult. (10822->478), div. (0->0), fcn. (10992->8), ass. (0->201)
t272 = qJ(1) + pkin(8);
t254 = sin(t272);
t255 = cos(t272);
t326 = sin(qJ(4));
t327 = cos(qJ(4));
t129 = -t254 * t327 + t255 * t326;
t181 = qJDD(1) - qJDD(4);
t182 = qJD(1) - qJD(4);
t177 = t255 * pkin(3);
t187 = qJD(1) ^ 2;
t169 = qJD(3) * t254;
t186 = cos(qJ(1));
t296 = pkin(1) * qJDD(1);
t179 = t186 * t296;
t233 = qJD(1) * t254;
t280 = t255 * pkin(2) + t254 * qJ(3);
t234 = qJD(1) * t255;
t283 = qJ(3) * t234 + t169;
t206 = qJDD(1) * t280 - qJDD(3) * t255 + t179 + (-pkin(2) * t233 + t169 + t283) * qJD(1);
t184 = sin(qJ(1));
t321 = t184 * pkin(1);
t212 = -pkin(3) * t254 - t321;
t190 = qJDD(1) * t177 + t187 * t212 + t206;
t183 = sin(qJ(5));
t185 = cos(qJ(5));
t244 = rSges(6,1) * t183 + rSges(6,2) * t185;
t158 = -rSges(6,1) * t185 + rSges(6,2) * t183;
t143 = t158 * qJD(5);
t277 = qJD(5) * t143;
t128 = -t254 * t326 - t255 * t327;
t293 = t128 * t185;
t294 = t128 * t183;
t75 = -rSges(6,1) * t293 + rSges(6,2) * t294 + t129 * rSges(6,3);
t305 = -t128 * pkin(4) + pkin(7) * t129 + t75;
t275 = qJD(5) * t183;
t99 = t182 * t129;
t219 = t128 * t275 + t185 * t99;
t274 = qJD(5) * t185;
t220 = t128 * t274 - t183 * t99;
t98 = t182 * t128;
t32 = rSges(6,1) * t219 + t220 * rSges(6,2) + t98 * rSges(6,3);
t367 = t99 * pkin(4) + t98 * pkin(7) + t32;
t78 = qJD(5) * t98 + qJDD(5) * t129;
t9 = -t129 * t277 + t305 * t181 + t182 * t367 + t78 * t244 + t190;
t369 = g(2) - t9;
t172 = t255 * qJ(3);
t248 = t254 * pkin(2);
t133 = t248 - t172;
t210 = -t133 + t212;
t180 = t186 * pkin(1);
t213 = -t177 - t180;
t170 = qJD(3) * t255;
t247 = qJDD(3) * t254 + (-qJD(1) * t280 + 0.2e1 * t170) * qJD(1);
t189 = qJDD(1) * t210 + t187 * t213 + t247;
t217 = t129 * t275 - t185 * t98;
t218 = t129 * t274 + t183 * t98;
t31 = rSges(6,1) * t217 + t218 * rSges(6,2) + t99 * rSges(6,3);
t363 = t98 * pkin(4) - t99 * pkin(7) - t31;
t290 = t129 * t185;
t291 = t129 * t183;
t74 = -rSges(6,1) * t290 + rSges(6,2) * t291 - t128 * rSges(6,3);
t41 = -t129 * pkin(4) - t128 * pkin(7) + t74;
t79 = qJD(5) * t99 - qJDD(5) * t128;
t8 = -t128 * t277 - t41 * t181 + t182 * t363 - t79 * t244 + t189;
t368 = t8 - g(1);
t286 = t129 * rSges(5,1) - t128 * rSges(5,2);
t52 = -t98 * rSges(5,1) - t99 * rSges(5,2);
t19 = t181 * t286 - t182 * t52 + t189;
t366 = t19 - g(1);
t365 = t182 * t41;
t361 = -qJD(1) * t133 + t169 - t283;
t303 = Icges(6,4) * t185;
t230 = -Icges(6,2) * t183 + t303;
t64 = -Icges(6,6) * t129 + t128 * t230;
t308 = t183 * t64;
t304 = Icges(6,4) * t183;
t232 = Icges(6,1) * t185 - t304;
t68 = -Icges(6,5) * t129 + t128 * t232;
t235 = t185 * t68 - t308;
t65 = Icges(6,6) * t128 + t129 * t230;
t309 = t183 * t65;
t69 = Icges(6,5) * t128 + t129 * t232;
t360 = -t185 * t69 + t309;
t228 = Icges(6,5) * t185 - Icges(6,6) * t183;
t61 = Icges(6,3) * t128 + t129 * t228;
t319 = -t128 * t61 - t290 * t69;
t60 = -Icges(6,3) * t129 + t128 * t228;
t359 = -(t60 - t309) * t129 + t319;
t276 = qJD(5) * t244;
t358 = t129 * t276 + t182 * t305;
t236 = -t183 * t68 - t185 * t64;
t238 = -t183 * t69 - t185 * t65;
t278 = qJD(5) * t129;
t279 = qJD(5) * t128;
t328 = -t182 / 0.2e1;
t357 = ((-t236 * t128 + t238 * t129) * qJD(5) + t236 * t279 - t238 * t278) * t328;
t100 = rSges(5,1) * t128 + rSges(5,2) * t129;
t53 = t99 * rSges(5,1) - t98 * rSges(5,2);
t20 = -t100 * t181 + t182 * t53 + t190;
t356 = t20 - g(2);
t355 = t129 * t60;
t352 = t61 * t129;
t350 = t61 + t308;
t349 = t182 * t286;
t194 = -t280 - t180;
t253 = t177 - t194;
t281 = t255 * rSges(4,1) + t254 * rSges(4,3);
t337 = t194 - t281;
t345 = -qJD(1) * t253 + t170;
t196 = -t248 + t212;
t343 = -t361 + (t196 - t212) * qJD(1);
t227 = Icges(6,5) * t183 + Icges(6,6) * t185;
t85 = t227 * t128;
t84 = t227 * t129;
t339 = t255 * rSges(3,1) - t254 * rSges(3,2);
t127 = t180 + t339;
t229 = Icges(6,2) * t185 + t304;
t231 = Icges(6,1) * t183 + t303;
t226 = -t183 * t229 + t185 * t231;
t44 = t129 * t226 + t85;
t42 = t44 * t182;
t45 = t128 * t226 - t84;
t43 = t45 * t182;
t333 = t128 * t276 - t365;
t313 = t229 * t129 - t69;
t315 = -t231 * t129 - t65;
t332 = t183 * t313 + t185 * t315;
t329 = m(4) + m(5);
t318 = t128 * t60 + t290 * t68;
t317 = t293 * t69 - t352;
t316 = t293 * t68 - t355;
t314 = -t231 * t128 - t64;
t312 = t229 * t128 - t68;
t289 = t228 * t182;
t285 = -t229 + t232;
t284 = -t230 - t231;
t273 = -m(6) - t329;
t271 = t187 * t321;
t259 = t279 / 0.2e1;
t258 = -t278 / 0.2e1;
t257 = t278 / 0.2e1;
t135 = rSges(3,1) * t254 + rSges(3,2) * t255;
t126 = -t135 - t321;
t252 = -t255 / 0.2e1;
t251 = t254 / 0.2e1;
t26 = Icges(6,4) * t219 + Icges(6,2) * t220 + Icges(6,6) * t98;
t28 = Icges(6,1) * t219 + Icges(6,4) * t220 + Icges(6,5) * t98;
t198 = qJD(5) * t236 + t183 * t26 - t185 * t28;
t25 = Icges(6,4) * t217 + Icges(6,2) * t218 + Icges(6,6) * t99;
t27 = Icges(6,1) * t217 + Icges(6,4) * t218 + Icges(6,5) * t99;
t199 = qJD(5) * t238 + t183 * t25 - t185 * t27;
t23 = Icges(6,5) * t217 + Icges(6,6) * t218 + Icges(6,3) * t99;
t24 = Icges(6,5) * t219 + Icges(6,6) * t220 + Icges(6,3) * t98;
t250 = -(-t128 * t23 + t129 * t199 - t360 * t98 - t61 * t99) * t128 + t129 * (-t128 * t24 + t129 * t198 + t235 * t98 - t60 * t99);
t249 = -t128 * (t128 * t199 + t129 * t23 + t360 * t99 - t61 * t98) + t129 * (t128 * t198 + t129 * t24 - t235 * t99 - t60 * t98);
t246 = t254 * rSges(4,1);
t159 = rSges(2,1) * t186 - rSges(2,2) * t184;
t157 = rSges(2,1) * t184 + rSges(2,2) * t186;
t15 = -t291 * t65 - t319;
t16 = -t291 * t64 + t318;
t243 = -t128 * t15 + t129 * t16;
t17 = -t294 * t65 + t317;
t18 = -t294 * t64 + t316;
t242 = -t128 * t17 + t129 * t18;
t202 = qJD(1) * t210 + t169;
t29 = t202 + t333;
t201 = (t280 - t213) * qJD(1) - t170;
t30 = t201 + t358;
t241 = -t128 * t29 - t129 * t30;
t240 = t128 * t32 + t129 * t31;
t239 = t128 * t75 + t129 * t74;
t107 = -t183 * t231 - t185 * t229;
t214 = -t180 * t187 - t184 * t296;
t211 = t183 * t312 + t185 * t314;
t207 = (t183 * t284 + t185 * t285) * t182;
t141 = t230 * qJD(5);
t142 = t232 * qJD(5);
t197 = qJD(5) * t107 - t141 * t183 + t142 * t185;
t195 = -t246 - t248 - t321;
t193 = t172 + t196;
t11 = -t360 * qJD(5) - t183 * t27 - t185 * t25;
t12 = qJD(5) * t235 - t183 * t28 - t185 * t26;
t140 = t228 * qJD(5);
t13 = t128 * t140 + t129 * t197 + t226 * t98 - t227 * t99;
t14 = t128 * t197 - t129 * t140 - t226 * t99 - t227 * t98;
t6 = qJD(5) * t243 + t42;
t7 = qJD(5) * t242 + t43;
t188 = t182 * (-t142 * t183 + t229 * t275 + (-qJD(5) * t231 - t141) * t185) + t257 * t6 - (-t236 + t45) * t78 / 0.2e1 - (-t238 + t44) * t79 / 0.2e1 + (t12 + t14) * t258 + (-Icges(5,3) + t107) * t181 + (t11 + t13 + t7) * t259;
t174 = t255 * rSges(4,3);
t167 = rSges(4,3) * t234;
t134 = t246 - t174;
t91 = t244 * t128;
t90 = t244 * t129;
t59 = -t100 * t182 + t201;
t58 = t202 + t349;
t40 = -t271 + qJDD(1) * t281 + qJD(1) * (-rSges(4,1) * t233 + t167) + t206;
t39 = -t187 * t281 + t214 + t247 + (-t133 - t134) * qJDD(1);
t33 = qJD(5) * t239 + qJD(2);
t10 = (-t128 * t65 - t129 * t64) * t183 + t317 + t318;
t5 = qJD(5) * t240 + t74 * t78 - t75 * t79 + qJDD(2);
t1 = [-m(2) * (-g(1) * t157 + g(2) * t159) + (-t42 + ((t17 + (-t235 + t61) * t129) * t129 + (t350 * t128 + t18 - t316 - t355) * t128) * qJD(5)) * t258 + (t43 + ((t129 * t360 + t15 + t316) * t129 + (-t350 * t129 - t10 + t16) * t128) * qJD(5)) * t259 - t188 + t357 + ((-t135 * t187 - g(2) + t179 - t271) * t127 + (-t187 * t339 - g(1) + t214) * t126) * m(3) + (m(3) * (-t135 * t126 + t339 * t127) + m(2) * (t157 ^ 2 + t159 ^ 2) + Icges(2,3) + Icges(3,3) + Icges(4,2)) * qJDD(1) + (-t369 * (t253 + t305) + (t345 + t363) * t29 + (-t244 * t279 + t29 + t343 + t365 + t367) * t30 + t368 * (t193 - t41)) * m(6) + (t356 * (-t100 + t253) + (t345 - t52) * t58 + (t343 + t53 + t58 - t349) * t59 + t366 * (t193 + t286)) * m(5) + (-(-g(2) + t40) * t337 + (-g(1) + t39) * (t172 + t174 + t195) + (-t167 - (t195 + t134 + t321) * qJD(1) + t361) * (t337 * qJD(1) + t170)) * m(4); m(6) * t5 + (m(3) + t329) * qJDD(2) + (-m(3) + t273) * g(3); t273 * (g(1) * t254 - g(2) * t255) + 0.2e1 * (t251 * t8 + t252 * t9) * m(6) + 0.2e1 * (t19 * t251 + t20 * t252) * m(5) + 0.2e1 * (t251 * t39 + t252 * t40) * m(4); (-t43 + ((-t15 - t359) * t129 + (-t16 + (-t360 + t60) * t128 - t352) * t128) * qJD(5)) * t259 + (t42 + ((t10 - t17) * t129 + (t128 * t235 - t18 + t359) * t128) * qJD(5)) * t258 + t188 + t357 + (t368 * t41 + (t333 - t367) * t30 + (-t358 - t363) * t29 + t369 * t305) * m(6) + (t52 * t58 - t53 * t59 - (-t182 * t59 + t366) * t286 + (t182 * t58 + t356) * t100) * m(5); t98 * t7 / 0.2e1 + t129 * (qJD(5) * t249 + t14 * t182 + t17 * t79 + t18 * t78 + t181 * t45) / 0.2e1 + t78 * t242 / 0.2e1 + (t17 * t99 + t18 * t98 + t249) * t257 + t99 * t6 / 0.2e1 - t128 * (qJD(5) * t250 + t13 * t182 + t15 * t79 + t16 * t78 + t181 * t44) / 0.2e1 + t79 * t243 / 0.2e1 - (t15 * t99 + t16 * t98 + t250) * t279 / 0.2e1 + t181 * (t128 * t238 - t129 * t236) / 0.2e1 + t182 * (-t11 * t128 + t12 * t129 - t236 * t98 - t238 * t99) / 0.2e1 + ((t85 * t278 - t289) * t129 + (t207 + (-t332 * t128 + (-t84 + t211) * t129) * qJD(5)) * t128) * t258 + ((t84 * t279 + t289) * t128 + (t207 + (t211 * t129 + (-t332 - t85) * t128) * qJD(5)) * t129) * t259 + ((t183 * t285 - t185 * t284) * t182 + ((t128 * t313 - t129 * t312) * t185 + (-t128 * t315 + t129 * t314) * t183) * qJD(5)) * t328 + (t5 * t239 + t33 * (t74 * t98 - t75 * t99 + t240) + t241 * t143 - (-t128 * t8 - t129 * t9 + t29 * t99 - t30 * t98) * t244 - (-t29 * t90 + t30 * t91) * t182 - (t33 * (t128 * t91 + t129 * t90) + t241 * t158) * qJD(5) - g(1) * t91 - g(2) * t90 - g(3) * t158) * m(6);];
tau = t1;
