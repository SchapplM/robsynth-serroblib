% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPR13_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR13_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR13_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR13_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR13_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:43:19
% EndTime: 2019-12-31 21:43:52
% DurationCPUTime: 15.15s
% Computational Cost: add. (5838->596), mult. (15689->819), div. (0->0), fcn. (11289->8), ass. (0->265)
t199 = cos(qJ(2));
t196 = sin(qJ(2));
t192 = sin(pkin(5));
t264 = qJD(1) * t192;
t245 = t196 * t264;
t193 = cos(pkin(5));
t263 = qJD(1) * t193;
t256 = pkin(1) * t263;
t155 = -pkin(7) * t245 + t199 * t256;
t236 = qJD(2) + t263;
t121 = -pkin(2) * t236 - t155;
t244 = t199 * t264;
t367 = qJD(3) - t244;
t297 = t367 / 0.2e1;
t195 = sin(qJ(3));
t198 = cos(qJ(3));
t140 = t195 * t236 + t198 * t245;
t299 = t140 / 0.2e1;
t134 = qJD(5) + t140;
t303 = t134 / 0.2e1;
t213 = t198 * t236;
t139 = t195 * t245 - t213;
t194 = sin(qJ(5));
t197 = cos(qJ(5));
t97 = t139 * t194 + t197 * t367;
t312 = t97 / 0.2e1;
t96 = t139 * t197 - t194 * t367;
t314 = t96 / 0.2e1;
t361 = -t367 / 0.2e1;
t302 = -t139 / 0.2e1;
t366 = Ifges(4,4) + Ifges(5,6);
t363 = t366 * t302;
t205 = -t140 * qJ(4) + t121;
t59 = t139 * pkin(3) + t205;
t311 = pkin(3) + pkin(9);
t158 = pkin(7) * t244 + t196 * t256;
t122 = pkin(8) * t236 + t158;
t152 = (-pkin(2) * t199 - pkin(8) * t196 - pkin(1)) * t192;
t131 = qJD(1) * t152;
t75 = t122 * t195 - t198 * t131;
t212 = pkin(4) * t140 + t75;
t350 = qJD(4) + t212;
t34 = -t311 * t367 + t350;
t37 = t139 * t311 + t205;
t8 = -t194 * t37 + t197 * t34;
t9 = t194 * t34 + t197 * t37;
t368 = t8 * mrSges(6,1) + t121 * mrSges(4,2) - t9 * mrSges(6,2) - t59 * mrSges(5,3) + Ifges(5,4) * t361 + Ifges(4,5) * t297 + 0.2e1 * Ifges(6,5) * t312 + 0.2e1 * Ifges(6,6) * t314 + 0.2e1 * Ifges(6,3) * t303 + t363 + (Ifges(4,1) + Ifges(5,2)) * t299;
t259 = qJD(3) * t195;
t237 = pkin(3) * t259 - qJD(4) * t195;
t247 = t195 * pkin(3) * t244 + t158;
t273 = qJ(4) * t198;
t365 = t237 - t247 + t367 * (pkin(9) * t195 - t273);
t310 = pkin(4) + pkin(8);
t181 = t310 * t198;
t211 = (pkin(2) * t196 - pkin(8) * t199) * t192;
t156 = qJD(1) * t211;
t92 = -t195 * t155 + t156 * t198;
t364 = qJD(3) * t181 - (pkin(4) * t198 * t199 - t196 * t311) * t264 + t92;
t360 = -qJD(2) / 0.2e1;
t359 = Ifges(5,4) - Ifges(4,5);
t337 = Ifges(5,5) - Ifges(4,6);
t358 = Ifges(4,3) + Ifges(5,1);
t238 = -qJ(4) * t195 - pkin(2);
t166 = -t198 * t311 + t238;
t180 = t310 * t195;
t111 = t166 * t197 + t180 * t194;
t357 = -qJD(5) * t111 - t365 * t194 + t364 * t197;
t110 = -t166 * t194 + t180 * t197;
t356 = qJD(5) * t110 + t364 * t194 + t365 * t197;
t268 = t195 * t199;
t93 = t198 * t155 + t195 * t156;
t355 = -t310 * t259 - (-pkin(4) * t268 + qJ(4) * t196) * t264 - t93;
t353 = Ifges(3,6) * t360;
t334 = -qJD(4) - t75;
t63 = -pkin(3) * t367 - t334;
t352 = t63 * mrSges(5,1) + t75 * mrSges(4,3);
t76 = t122 * t198 + t131 * t195;
t65 = -qJ(4) * t367 - t76;
t351 = t65 * mrSges(5,1) - t76 * mrSges(4,3);
t300 = -t140 / 0.2e1;
t301 = t139 / 0.2e1;
t349 = Ifges(4,1) * t299 + Ifges(4,4) * t302 - Ifges(5,2) * t300 - Ifges(5,6) * t301 - t297 * t359 + t368;
t348 = -t121 * mrSges(4,1) + t59 * mrSges(5,2) + Ifges(5,5) * t361 + Ifges(4,6) * t297 + (Ifges(4,2) + Ifges(5,3)) * t302 + t366 * t299;
t347 = -Ifges(5,4) / 0.2e1;
t262 = qJD(2) * t192;
t240 = qJD(1) * t262;
t233 = t199 * t240;
t272 = t192 * t196;
t248 = t195 * t272;
t235 = qJD(3) * t248;
t103 = qJD(1) * t235 - qJD(3) * t213 - t198 * t233;
t261 = qJD(2) * t196;
t243 = t192 * t261;
t214 = t311 * t243;
t157 = qJD(2) * t211;
t147 = qJD(1) * t157;
t186 = pkin(7) * t272;
t290 = pkin(1) * t199;
t164 = t193 * t290 - t186;
t159 = t164 * qJD(2);
t148 = qJD(1) * t159;
t258 = qJD(3) * t198;
t29 = -t122 * t258 - t131 * t259 + t147 * t198 - t195 * t148;
t12 = -pkin(4) * t103 - qJD(1) * t214 - t29;
t104 = qJD(3) * t140 + t195 * t233;
t271 = t192 * t199;
t165 = t193 * t196 * pkin(1) + pkin(7) * t271;
t160 = t165 * qJD(2);
t149 = qJD(1) * t160;
t203 = qJ(4) * t103 - qJD(4) * t140 + t149;
t17 = t104 * t311 + t203;
t1 = qJD(5) * t8 + t12 * t194 + t17 * t197;
t346 = t1 * mrSges(6,2);
t2 = -qJD(5) * t9 + t12 * t197 - t17 * t194;
t345 = t2 * mrSges(6,1);
t306 = -t103 / 0.2e1;
t344 = -t104 / 0.2e1;
t234 = t196 * t240;
t27 = -pkin(3) * t234 - t29;
t84 = -t103 * mrSges(5,1) + mrSges(5,2) * t234;
t335 = m(5) * t27 + t84;
t333 = t103 * t359 + t104 * t337 + t234 * t358;
t332 = t1 * t194 + t197 * t2;
t254 = Ifges(4,5) / 0.2e1 + t347;
t331 = t254 * t367 + t352 + t368;
t151 = pkin(8) * t193 + t165;
t208 = t151 * t259 - t152 * t258 - t195 * t157 - t198 * t159;
t30 = -t192 * (qJ(4) * t261 - qJD(4) * t199) + t208;
t280 = Ifges(3,6) * t193;
t285 = Ifges(3,4) * t196;
t328 = t63 * mrSges(5,2) + t353 - (t280 + (Ifges(3,2) * t199 + t285) * t192) * qJD(1) / 0.2e1 - t65 * mrSges(5,3) - t75 * mrSges(4,1) - t76 * mrSges(4,2);
t284 = Ifges(6,4) * t194;
t222 = Ifges(6,2) * t197 + t284;
t283 = Ifges(6,4) * t197;
t224 = Ifges(6,1) * t194 + t283;
t227 = mrSges(6,1) * t197 - mrSges(6,2) * t194;
t229 = t194 * t8 - t197 * t9;
t279 = Ifges(6,6) * t197;
t281 = Ifges(6,5) * t194;
t295 = -t197 / 0.2e1;
t296 = -t194 / 0.2e1;
t304 = -t134 / 0.2e1;
t313 = -t97 / 0.2e1;
t315 = -t96 / 0.2e1;
t293 = Ifges(6,4) * t97;
t32 = Ifges(6,2) * t96 + Ifges(6,6) * t134 + t293;
t95 = Ifges(6,4) * t96;
t33 = Ifges(6,1) * t97 + Ifges(6,5) * t134 + t95;
t289 = pkin(4) * t139;
t42 = -t65 - t289;
t327 = mrSges(6,3) * t229 + (t279 + t281) * t304 + t222 * t315 + t224 * t313 + t42 * t227 + t295 * t32 + t296 * t33;
t326 = -Ifges(4,4) * t299 - Ifges(4,2) * t302 + Ifges(5,6) * t300 + Ifges(5,3) * t301 + t297 * t337 - t348;
t40 = qJD(5) * t96 + t104 * t194 + t197 * t234;
t41 = -qJD(5) * t97 + t104 * t197 - t194 * t234;
t7 = t40 * Ifges(6,1) + t41 * Ifges(6,4) - t103 * Ifges(6,5);
t324 = t7 / 0.2e1;
t323 = t32 / 0.2e1;
t322 = t33 / 0.2e1;
t321 = t40 / 0.2e1;
t320 = t41 / 0.2e1;
t319 = Ifges(5,2) * t306 + Ifges(5,6) * t344 + t234 * t347;
t308 = pkin(1) * mrSges(3,1);
t307 = pkin(1) * mrSges(3,2);
t39 = Ifges(6,5) * t40;
t38 = Ifges(6,6) * t41;
t282 = Ifges(3,5) * t193;
t107 = mrSges(5,1) * t139 - mrSges(5,3) * t367;
t49 = -mrSges(6,1) * t96 + mrSges(6,2) * t97;
t277 = -t107 + t49;
t276 = -mrSges(3,1) * t236 + mrSges(4,1) * t139 + mrSges(4,2) * t140 + mrSges(3,3) * t245;
t274 = qJ(4) * t139;
t270 = t194 * t195;
t269 = t195 * t197;
t267 = t197 * t199;
t105 = -mrSges(4,2) * t367 - mrSges(4,3) * t139;
t266 = t105 - t107;
t106 = mrSges(4,1) * t367 - mrSges(4,3) * t140;
t108 = mrSges(5,1) * t140 + mrSges(5,2) * t367;
t265 = t106 - t108;
t91 = t198 * t151 + t195 * t152;
t260 = qJD(2) * t199;
t5 = -Ifges(6,3) * t103 + t38 + t39;
t255 = Ifges(4,4) / 0.2e1 + Ifges(5,6) / 0.2e1;
t253 = -Ifges(5,2) / 0.2e1 - Ifges(4,1) / 0.2e1;
t252 = Ifges(4,6) / 0.2e1 - Ifges(5,5) / 0.2e1;
t251 = Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1;
t250 = -Ifges(5,3) / 0.2e1 - Ifges(4,2) / 0.2e1;
t242 = t192 * t260;
t90 = -t195 * t151 + t152 * t198;
t231 = t345 - t346;
t79 = pkin(3) * t271 - t90;
t228 = -t149 * mrSges(3,1) - t148 * mrSges(3,2);
t226 = mrSges(6,1) * t194 + mrSges(6,2) * t197;
t225 = Ifges(6,1) * t197 - t284;
t223 = -Ifges(6,2) * t194 + t283;
t221 = Ifges(6,5) * t197 - Ifges(6,6) * t194;
t20 = -mrSges(6,1) * t103 - mrSges(6,3) * t40;
t21 = mrSges(6,2) * t103 + mrSges(6,3) * t41;
t218 = t194 * t21 + t197 * t20;
t163 = t193 * t195 + t198 * t272;
t52 = pkin(4) * t163 + pkin(9) * t271 + t79;
t162 = -t193 * t198 + t248;
t150 = t186 + (-pkin(2) - t290) * t193;
t207 = -qJ(4) * t163 + t150;
t60 = t162 * t311 + t207;
t15 = -t194 * t60 + t197 * t52;
t16 = t194 * t52 + t197 * t60;
t66 = -mrSges(6,2) * t134 + mrSges(6,3) * t96;
t67 = mrSges(6,1) * t134 - mrSges(6,3) * t97;
t217 = -t194 * t67 + t197 * t66;
t216 = t195 * t65 + t198 * t63;
t215 = -t195 * t76 + t198 * t75;
t78 = qJ(4) * t271 - t91;
t44 = -t151 * t258 - t152 * t259 + t157 * t198 - t195 * t159;
t210 = -t162 * t194 + t192 * t267;
t114 = t162 * t197 + t194 * t271;
t28 = -t122 * t259 + t131 * t258 + t195 * t147 + t198 * t148;
t113 = -t235 + (qJD(3) * t193 + t242) * t198;
t204 = -qJ(4) * t113 - qJD(4) * t163 + t160;
t22 = -qJ(4) * t234 - qJD(4) * t367 - t28;
t202 = t255 * t140 + t252 * t367 + t348 - t351;
t182 = Ifges(3,4) * t244;
t177 = Ifges(3,5) * t233;
t173 = -pkin(3) * t198 + t238;
t161 = -qJ(4) * t258 + t237;
t154 = -mrSges(3,2) * t236 + mrSges(3,3) * t244;
t127 = (t194 * t268 + t196 * t197) * t264;
t126 = (-t194 * t196 + t195 * t267) * t264;
t119 = Ifges(3,1) * t245 + Ifges(3,5) * t236 + t182;
t112 = qJD(3) * t163 + t195 * t242;
t94 = -t244 * t273 + t247;
t89 = -mrSges(5,2) * t139 - mrSges(5,3) * t140;
t87 = pkin(3) * t140 + t274;
t86 = -mrSges(4,2) * t234 - mrSges(4,3) * t104;
t85 = mrSges(4,1) * t234 + mrSges(4,3) * t103;
t83 = mrSges(5,1) * t104 - mrSges(5,3) * t234;
t82 = -pkin(3) * t245 - t92;
t80 = -qJ(4) * t245 - t93;
t77 = pkin(3) * t162 + t207;
t73 = Ifges(5,1) * t367 - t140 * Ifges(5,4) + t139 * Ifges(5,5);
t70 = t140 * Ifges(4,5) - t139 * Ifges(4,6) + Ifges(4,3) * t367;
t64 = t140 * t311 + t274;
t61 = -pkin(4) * t162 - t78;
t58 = qJD(5) * t114 + t112 * t194 + t197 * t243;
t57 = qJD(5) * t210 + t112 * t197 - t194 * t243;
t54 = mrSges(4,1) * t104 - mrSges(4,2) * t103;
t53 = -mrSges(5,2) * t104 + mrSges(5,3) * t103;
t51 = t76 - t289;
t48 = -t103 * Ifges(4,1) - t104 * Ifges(4,4) + Ifges(4,5) * t234;
t47 = -t103 * Ifges(4,4) - t104 * Ifges(4,2) + Ifges(4,6) * t234;
t45 = Ifges(5,5) * t234 + t103 * Ifges(5,6) + t104 * Ifges(5,3);
t36 = pkin(3) * t112 + t204;
t35 = -pkin(3) * t243 - t44;
t26 = pkin(3) * t104 + t203;
t25 = t112 * t311 + t204;
t19 = -pkin(4) * t112 - t30;
t18 = pkin(4) * t113 - t214 - t44;
t14 = t194 * t51 + t197 * t64;
t13 = -t194 * t64 + t197 * t51;
t11 = -pkin(4) * t104 - t22;
t10 = -mrSges(6,1) * t41 + mrSges(6,2) * t40;
t6 = Ifges(6,4) * t40 + Ifges(6,2) * t41 - Ifges(6,6) * t103;
t4 = -qJD(5) * t16 + t18 * t197 - t194 * t25;
t3 = qJD(5) * t15 + t18 * t194 + t197 * t25;
t23 = [t276 * t160 + t104 * (-Ifges(5,5) * t271 - Ifges(5,6) * t163 + Ifges(5,3) * t162) / 0.2e1 + t103 * (-Ifges(5,4) * t271 - Ifges(5,2) * t163 + Ifges(5,6) * t162) / 0.2e1 + t27 * (mrSges(5,1) * t163 - mrSges(5,2) * t271) + t29 * (-mrSges(4,1) * t271 - mrSges(4,3) * t163) + t22 * (mrSges(5,1) * t162 + mrSges(5,3) * t271) + t28 * (mrSges(4,2) * t271 - mrSges(4,3) * t162) + (t148 * t199 + t149 * t196 + (-t155 * t199 - t158 * t196) * qJD(2)) * t192 * mrSges(3,3) + ((-t164 * mrSges(3,3) + t282 + (-0.2e1 * t307 + 0.3e1 / 0.2e1 * Ifges(3,4) * t199) * t192) * t260 + (-t165 * mrSges(3,3) - 0.3e1 / 0.2e1 * t280 + (-0.2e1 * t308 - 0.3e1 / 0.2e1 * t285) * t192 + t254 * t163 - t252 * t162 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) - t251) * t271) * t261) * t264 + (t326 + t351) * t112 + (t349 + t352) * t113 + (Ifges(6,1) * t58 + Ifges(6,4) * t57) * t312 + (Ifges(6,5) * t58 + Ifges(6,6) * t57) * t303 - t163 * t346 + (Ifges(4,4) * t163 - Ifges(4,2) * t162 - Ifges(4,6) * t271) * t344 + t163 * t345 + (t48 + t5) * t163 / 0.2e1 + ((t73 + t70) * t196 + t199 * t119) * t262 / 0.2e1 + m(4) * (t121 * t160 + t149 * t150 - t208 * t76 + t28 * t91 + t29 * t90 - t44 * t75) - t208 * t105 + (-Ifges(4,4) * t162 - Ifges(4,5) * t271 - Ifges(6,5) * t210 + Ifges(6,6) * t114 + (Ifges(4,1) + Ifges(6,3)) * t163) * t306 + (t1 * t114 + t2 * t210 + t57 * t9 - t58 * t8) * mrSges(6,3) + t11 * (-mrSges(6,1) * t114 - mrSges(6,2) * t210) + (-Ifges(6,4) * t210 + Ifges(6,2) * t114 + Ifges(6,6) * t163) * t320 + (-Ifges(6,1) * t210 + Ifges(6,4) * t114 + Ifges(6,5) * t163) * t321 - t210 * t324 + qJD(2) ^ 2 * (Ifges(3,5) * t199 - Ifges(3,6) * t196) * t192 / 0.2e1 + (Ifges(6,4) * t58 + Ifges(6,2) * t57) * t314 + (Ifges(5,4) * t300 + Ifges(4,5) * t299 + Ifges(5,5) * t301 + Ifges(4,6) * t302 + t297 * t358 + t328) * t243 + m(6) * (t1 * t16 + t11 * t61 + t15 * t2 + t19 * t42 + t3 * t9 + t4 * t8) + m(5) * (t22 * t78 + t26 * t77 + t27 * t79 + t30 * t65 + t35 * t63 + t36 * t59) + m(3) * (t148 * t165 - t149 * t164 - t155 * t160 + t158 * t159) - t333 * t271 / 0.2e1 + (t177 / 0.2e1 + t228) * t193 + t26 * (-mrSges(5,2) * t162 - mrSges(5,3) * t163) + t162 * t45 / 0.2e1 - t162 * t47 / 0.2e1 + t149 * (mrSges(4,1) * t162 + mrSges(4,2) * t163) + t159 * t154 + t150 * t54 + t114 * t6 / 0.2e1 + t44 * t106 + t30 * t107 + t35 * t108 + t36 * t89 + t90 * t85 + t91 * t86 + t77 * t53 + t78 * t83 + t79 * t84 + t163 * t319 + t58 * t322 + t57 * t323 + t15 * t20 + t16 * t21 + t19 * t49 + t42 * (-mrSges(6,1) * t57 + mrSges(6,2) * t58) + t61 * t10 + t3 * t66 + t4 * t67; -t276 * t158 + t355 * t49 + (-m(4) * t149 - t54) * pkin(2) + (t42 * (-mrSges(6,1) * t269 + mrSges(6,2) * t270) + (Ifges(6,5) * t270 + Ifges(6,6) * t269) * t303 + (Ifges(6,1) * t270 + Ifges(6,4) * t269) * t312 + (Ifges(6,4) * t270 + Ifges(6,2) * t269) * t314 + t269 * t323 + t270 * t322 + t215 * mrSges(4,3) + t216 * mrSges(5,1) + (m(4) * t215 + m(5) * t216) * pkin(8) + (t269 * t9 - t270 * t8) * mrSges(6,3) + (-t266 * pkin(8) + t326) * t195 + (-t265 * pkin(8) + t349) * t198) * qJD(3) + t177 + (-t126 * t9 + t127 * t8) * mrSges(6,3) + (Ifges(6,5) * t127 + Ifges(6,6) * t126) * t304 + (Ifges(6,1) * t127 + Ifges(6,4) * t126) * t313 + (t27 * mrSges(5,1) - t29 * mrSges(4,3) - t26 * mrSges(5,3) + t319 + t149 * mrSges(4,2) + t39 / 0.2e1 + t38 / 0.2e1 + t5 / 0.2e1 + t48 / 0.2e1 - t255 * t104 + (-Ifges(6,3) / 0.2e1 + t253) * t103 + (-m(4) * t29 + t335 - t85) * pkin(8) + t231) * t195 + (t161 - t94) * t89 + t228 - m(5) * (t59 * t94 + t63 * t82 + t65 * t80) - m(4) * (t121 * t158 - t75 * t92 + t76 * t93) + t356 * t66 + t357 * t67 + (t1 * t111 + t11 * t181 + t110 * t2 + t355 * t42 + t356 * t9 + t357 * t8) * m(6) + (Ifges(6,4) * t127 + Ifges(6,2) * t126) * t315 + (-t22 * mrSges(5,1) + t28 * mrSges(4,3) + t26 * mrSges(5,2) - t149 * mrSges(4,1) + t47 / 0.2e1 - t45 / 0.2e1 - t40 * t224 / 0.2e1 - t41 * t222 / 0.2e1 + t11 * t227 + t7 * t296 + t6 * t295 + t250 * t104 + (t281 / 0.2e1 + t279 / 0.2e1 - t255) * t103 + (-t1 * t197 + t194 * t2) * mrSges(6,3) + (m(4) * t28 - m(5) * t22 - t83 + t86) * pkin(8) + (t221 * t304 + t225 * t313 + t223 * t315 - t42 * t226 + t33 * t295 + t194 * t323 + (t194 * t9 + t197 * t8) * mrSges(6,3)) * qJD(5)) * t198 + m(5) * (t161 * t59 + t173 * t26) + t181 * t10 + t173 * t53 - t155 * t154 - t42 * (-mrSges(6,1) * t126 + mrSges(6,2) * t127) - t127 * t33 / 0.2e1 - t126 * t32 / 0.2e1 - t93 * t105 - t92 * t106 - t80 * t107 - t82 * t108 + t110 * t20 + t111 * t21 + ((t353 - t70 / 0.2e1 - t73 / 0.2e1 + t158 * mrSges(3,3) - t251 * t367 - t254 * t140 + t252 * t139 + (t280 / 0.2e1 + (t308 + t285 / 0.2e1) * t192) * qJD(1) + (-t195 * t359 - t337 * t198) * qJD(2) / 0.2e1 - t328) * t196 + (-t182 / 0.2e1 - t119 / 0.2e1 + Ifges(3,5) * t360 + (-t282 / 0.2e1 + (t307 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t196) * t192) * qJD(1) + t155 * mrSges(3,3) + (t139 * t250 + t202) * t195 + (t139 * t255 + t140 * t253 - t331) * t198) * t199) * t264; t277 * qJD(4) + (t202 + (t250 - t253) * t139 + t327) * t140 + t327 * qJD(5) + (t331 + t363) * t139 + (t11 * qJ(4) - t13 * t8 - t14 * t9 + t350 * t42) * m(6) + (-t83 + t10) * qJ(4) + t221 * t306 + t6 * t296 + t333 + t212 * t49 + t266 * t75 + (-pkin(3) * t27 - qJ(4) * t22 + t334 * t65 - t59 * t87 - t63 * t76) * m(5) - ((-m(6) * t229 + t217) * qJD(5) + m(6) * t332 + t218) * t311 + t265 * t76 + t11 * t226 - t87 * t89 - pkin(3) * t84 + t223 * t320 + t225 * t321 + t197 * t324 - t22 * mrSges(5,3) + t27 * mrSges(5,2) - t28 * mrSges(4,2) + t29 * mrSges(4,1) - t332 * mrSges(6,3) - t14 * t66 - t13 * t67; -t277 * t367 + t217 * qJD(5) + (t217 + t89) * t140 - m(5) * (-t140 * t59 - t367 * t65) + t218 + (-t134 * t229 - t367 * t42 + t332) * m(6) + t335; -t42 * (mrSges(6,1) * t97 + mrSges(6,2) * t96) + (Ifges(6,1) * t96 - t293) * t313 + t32 * t312 + (Ifges(6,5) * t96 - Ifges(6,6) * t97) * t304 - t8 * t66 + t9 * t67 + (t8 * t96 + t9 * t97) * mrSges(6,3) + t231 + t5 + (-Ifges(6,2) * t97 + t33 + t95) * t315;];
tauc = t23(:);
