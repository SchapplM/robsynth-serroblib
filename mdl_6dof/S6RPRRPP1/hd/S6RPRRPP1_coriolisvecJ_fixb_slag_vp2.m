% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPP1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:27:55
% EndTime: 2019-03-09 04:28:16
% DurationCPUTime: 10.46s
% Computational Cost: add. (5453->557), mult. (13423->731), div. (0->0), fcn. (8477->8), ass. (0->247)
t201 = cos(qJ(4));
t246 = qJD(3) * qJD(4);
t199 = sin(qJ(4));
t200 = sin(qJ(3));
t249 = qJD(4) * t200;
t202 = cos(qJ(3));
t251 = qJD(3) * t202;
t320 = t199 * t249 - t201 * t251;
t122 = -qJD(1) * t320 + t201 * t246;
t248 = qJD(4) * t201;
t328 = t199 * t251 + t200 * t248;
t123 = -qJD(1) * t328 - t199 * t246;
t196 = sin(pkin(10));
t269 = cos(pkin(10));
t75 = t122 * t196 - t123 * t269;
t312 = -t75 / 0.2e1;
t76 = t122 * t269 + t196 * t123;
t310 = t76 / 0.2e1;
t337 = Ifges(6,1) + Ifges(7,1);
t336 = Ifges(6,4) - Ifges(7,5);
t323 = Ifges(7,4) + Ifges(6,5);
t227 = t269 * t199;
t162 = t196 * t201 + t227;
t255 = qJD(1) * t202;
t130 = t162 * t255;
t149 = t162 * qJD(4);
t260 = t130 - t149;
t226 = t269 * t201;
t265 = t196 * t199;
t209 = t226 - t265;
t207 = t202 * t209;
t131 = qJD(1) * t207;
t150 = t209 * qJD(4);
t259 = t131 - t150;
t342 = -Ifges(4,1) / 0.2e1;
t194 = Ifges(4,4) * t255;
t341 = -t194 / 0.2e1;
t253 = qJD(3) * t200;
t229 = qJD(1) * t253;
t340 = t337 * t310 + t336 * t312 + t323 * t229 / 0.2e1;
t191 = sin(pkin(9)) * pkin(1) + pkin(7);
t175 = t191 * qJD(1);
t139 = qJD(2) * t202 - t200 * t175;
t224 = pkin(3) * t200 - pkin(8) * t202;
t169 = t224 * qJD(1);
t100 = -t199 * t139 + t201 * t169;
t262 = t201 * t202;
t211 = pkin(4) * t200 - qJ(5) * t262;
t79 = qJD(1) * t211 + t100;
t101 = t201 * t139 + t199 * t169;
t239 = t199 * t255;
t84 = -qJ(5) * t239 + t101;
t30 = -t196 * t84 + t269 * t79;
t289 = -qJ(5) - pkin(8);
t228 = qJD(4) * t289;
t247 = qJD(5) * t201;
t145 = t199 * t228 + t247;
t208 = -qJD(5) * t199 + t201 * t228;
t93 = t145 * t196 - t208 * t269;
t339 = -t30 - t93;
t31 = t196 * t79 + t269 * t84;
t94 = t145 * t269 + t196 * t208;
t338 = -t31 + t94;
t252 = qJD(3) * t201;
t256 = qJD(1) * t200;
t166 = -t199 * t256 + t252;
t167 = qJD(3) * t199 + t201 * t256;
t107 = -t269 * t166 + t167 * t196;
t307 = t107 / 0.2e1;
t234 = Ifges(4,5) * qJD(3) / 0.2e1;
t335 = Ifges(6,6) - Ifges(7,6);
t333 = -qJ(6) * t256 + t338;
t332 = pkin(5) * t256 - t339;
t187 = qJD(4) - t255;
t210 = t196 * t166 + t167 * t269;
t331 = -t107 * t336 + t187 * t323 + t210 * t337;
t140 = t200 * qJD(2) + t202 * t175;
t117 = pkin(4) * t239 + t140;
t250 = qJD(4) * t199;
t330 = pkin(4) * t250 - pkin(5) * t260 + qJ(6) * t259 - qJD(6) * t162 - t117;
t128 = -qJD(3) * pkin(3) - t139;
t129 = qJD(3) * pkin(8) + t140;
t241 = -cos(pkin(9)) * pkin(1) - pkin(2);
t158 = -pkin(3) * t202 - t200 * pkin(8) + t241;
t134 = t158 * qJD(1);
t81 = -t129 * t199 + t201 * t134;
t82 = t129 * t201 + t134 * t199;
t215 = t82 * t199 + t81 * t201;
t280 = Ifges(5,4) * t201;
t219 = -Ifges(5,2) * t199 + t280;
t281 = Ifges(5,4) * t199;
t221 = Ifges(5,1) * t201 - t281;
t222 = mrSges(5,1) * t199 + mrSges(5,2) * t201;
t278 = Ifges(5,6) * t199;
t279 = Ifges(5,5) * t201;
t293 = t201 / 0.2e1;
t294 = -t199 / 0.2e1;
t295 = t187 / 0.2e1;
t297 = t167 / 0.2e1;
t275 = t167 * Ifges(5,4);
t97 = t166 * Ifges(5,2) + t187 * Ifges(5,6) + t275;
t159 = Ifges(5,4) * t166;
t98 = t167 * Ifges(5,1) + t187 * Ifges(5,5) + t159;
t203 = -t215 * mrSges(5,3) + t97 * t294 + t98 * t293 + t128 * t222 + t166 * t219 / 0.2e1 + t221 * t297 + (-t278 + t279) * t295;
t329 = t139 * mrSges(4,3) + t256 * t342 - t203 - t234 + t341;
t61 = qJ(5) * t166 + t82;
t272 = t196 * t61;
t60 = -qJ(5) * t167 + t81;
t51 = pkin(4) * t187 + t60;
t13 = t269 * t51 - t272;
t10 = -t187 * pkin(5) + qJD(6) - t13;
t102 = -t166 * pkin(4) + qJD(5) + t128;
t29 = t107 * pkin(5) - qJ(6) * t210 + t102;
t327 = -t102 * mrSges(6,2) - t10 * mrSges(7,2) + t13 * mrSges(6,3) + t29 * mrSges(7,3);
t54 = t269 * t61;
t14 = t196 * t51 + t54;
t11 = qJ(6) * t187 + t14;
t43 = Ifges(7,5) * t210 + Ifges(7,6) * t187 + Ifges(7,3) * t107;
t46 = Ifges(6,4) * t210 - Ifges(6,2) * t107 + Ifges(6,6) * t187;
t326 = -mrSges(6,1) * t102 - mrSges(7,1) * t29 + mrSges(7,2) * t11 + mrSges(6,3) * t14 + t46 / 0.2e1 - t43 / 0.2e1;
t325 = -Ifges(6,6) / 0.2e1;
t324 = Ifges(7,6) / 0.2e1;
t311 = t75 / 0.2e1;
t308 = -t107 / 0.2e1;
t304 = t210 / 0.2e1;
t233 = -Ifges(4,6) * qJD(3) / 0.2e1;
t25 = t75 * mrSges(7,1) - t76 * mrSges(7,3);
t26 = t75 * mrSges(6,1) + t76 * mrSges(6,2);
t321 = -t25 - t26;
t86 = -mrSges(6,2) * t187 - mrSges(6,3) * t107;
t89 = -mrSges(7,2) * t107 + mrSges(7,3) * t187;
t286 = t86 + t89;
t87 = mrSges(6,1) * t187 - mrSges(6,3) * t210;
t88 = -mrSges(7,1) * t187 + mrSges(7,2) * t210;
t285 = t87 - t88;
t168 = t191 * t262;
t111 = t199 * t158 + t168;
t132 = t139 * qJD(3);
t170 = t224 * qJD(3);
t157 = qJD(1) * t170;
t33 = -t129 * t250 + t201 * t132 + t134 * t248 + t199 * t157;
t34 = -qJD(4) * t82 - t132 * t199 + t201 * t157;
t216 = -t199 * t34 + t201 * t33;
t15 = pkin(4) * t229 - qJ(5) * t122 - qJD(5) * t167 + t34;
t20 = qJ(5) * t123 + qJD(5) * t166 + t33;
t4 = t196 * t15 + t269 * t20;
t1 = qJ(6) * t229 + qJD(6) * t187 + t4;
t3 = t15 * t269 - t196 * t20;
t2 = -pkin(5) * t229 - t3;
t318 = -t34 * mrSges(5,1) - t3 * mrSges(6,1) + t2 * mrSges(7,1) + t33 * mrSges(5,2) + t4 * mrSges(6,2) - t1 * mrSges(7,3) - Ifges(5,5) * t122 - Ifges(5,6) * t123;
t177 = t241 * qJD(1);
t225 = Ifges(6,3) / 0.2e1 + Ifges(7,2) / 0.2e1 + Ifges(5,3) / 0.2e1;
t244 = t324 + t325;
t245 = Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1;
t282 = Ifges(4,4) * t200;
t317 = t244 * t107 + t225 * t187 + t245 * t210 + t11 * mrSges(7,3) + t13 * mrSges(6,1) + t177 * mrSges(4,1) + t81 * mrSges(5,1) + t233 - (Ifges(4,2) * t202 + t282) * qJD(1) / 0.2e1 + Ifges(6,6) * t308 + Ifges(7,6) * t307 + t167 * Ifges(5,5) + t166 * Ifges(5,6) - t10 * mrSges(7,1) - t14 * mrSges(6,2) - t82 * mrSges(5,2) + t323 * t304 + (Ifges(6,3) + Ifges(7,2) + Ifges(5,3)) * t295;
t316 = Ifges(7,5) * t310 + Ifges(7,3) * t311 + t229 * t324;
t315 = -t76 * Ifges(6,4) / 0.2e1 + Ifges(6,2) * t311 + t229 * t325;
t305 = -t210 / 0.2e1;
t303 = t122 / 0.2e1;
t302 = t123 / 0.2e1;
t299 = -t166 / 0.2e1;
t298 = -t167 / 0.2e1;
t296 = -t187 / 0.2e1;
t292 = pkin(4) * t167;
t291 = pkin(4) * t196;
t266 = t191 * t199;
t257 = t201 * t170 + t253 * t266;
t36 = -t200 * t247 + t211 * qJD(3) + (-t168 + (qJ(5) * t200 - t158) * t199) * qJD(4) + t257;
t258 = t158 * t248 + t199 * t170;
t263 = t200 * t201;
t41 = (-qJ(5) * qJD(4) - qJD(3) * t191) * t263 + (-qJD(5) * t200 + (-qJ(5) * qJD(3) - qJD(4) * t191) * t202) * t199 + t258;
t9 = t196 * t36 + t269 * t41;
t56 = -mrSges(7,2) * t75 + mrSges(7,3) * t229;
t57 = -mrSges(6,2) * t229 - mrSges(6,3) * t75;
t288 = t56 + t57;
t58 = mrSges(6,1) * t229 - mrSges(6,3) * t76;
t59 = -mrSges(7,1) * t229 + t76 * mrSges(7,2);
t287 = t59 - t58;
t144 = t201 * t158;
t92 = -qJ(5) * t263 + t144 + (-pkin(4) - t266) * t202;
t264 = t199 * t200;
t99 = -qJ(5) * t264 + t111;
t40 = t196 * t92 + t269 * t99;
t284 = mrSges(5,3) * t166;
t283 = mrSges(5,3) * t167;
t273 = t177 * mrSges(4,2);
t243 = mrSges(4,3) * t256;
t261 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t166 + mrSges(5,2) * t167 + t243;
t146 = pkin(4) * t264 + t200 * t191;
t242 = mrSges(4,3) * t255;
t113 = pkin(4) * t328 + t191 * t251;
t193 = -pkin(4) * t201 - pkin(3);
t240 = t269 * pkin(4);
t230 = m(4) * t191 + mrSges(4,3);
t223 = mrSges(5,1) * t201 - mrSges(5,2) * t199;
t220 = Ifges(5,1) * t199 + t280;
t218 = Ifges(5,2) * t201 + t281;
t217 = Ifges(5,5) * t199 + Ifges(5,6) * t201;
t214 = t199 * t81 - t201 * t82;
t104 = mrSges(5,1) * t229 - mrSges(5,3) * t122;
t105 = -mrSges(5,2) * t229 + mrSges(5,3) * t123;
t213 = -t199 * t104 + t201 * t105;
t125 = -mrSges(5,2) * t187 + t284;
t126 = mrSges(5,1) * t187 - t283;
t212 = -t199 * t125 - t201 * t126;
t8 = -t196 * t41 + t269 * t36;
t39 = -t196 * t99 + t269 * t92;
t133 = t140 * qJD(3);
t83 = -t123 * pkin(4) + t133;
t192 = -t240 - pkin(5);
t189 = qJ(6) + t291;
t186 = Ifges(7,2) * t229;
t185 = Ifges(5,3) * t229;
t184 = Ifges(6,3) * t229;
t179 = t289 * t201;
t178 = -qJD(3) * mrSges(4,2) + t242;
t138 = t209 * t200;
t137 = t162 * t200;
t116 = -t179 * t269 + t265 * t289;
t115 = -t179 * t196 - t227 * t289;
t110 = -t202 * t266 + t144;
t103 = -pkin(5) * t209 - qJ(6) * t162 + t193;
t91 = qJD(3) * t207 - t162 * t249;
t90 = t196 * t320 - t226 * t249 - t227 * t251;
t78 = -mrSges(5,1) * t123 + mrSges(5,2) * t122;
t73 = pkin(5) * t137 - qJ(6) * t138 + t146;
t72 = Ifges(7,4) * t76;
t71 = Ifges(6,5) * t76;
t70 = Ifges(6,6) * t75;
t69 = Ifges(7,6) * t75;
t65 = t122 * Ifges(5,1) + t123 * Ifges(5,4) + Ifges(5,5) * t229;
t64 = t122 * Ifges(5,4) + t123 * Ifges(5,2) + Ifges(5,6) * t229;
t63 = -qJD(4) * t111 + t257;
t62 = (-t200 * t252 - t202 * t250) * t191 + t258;
t53 = mrSges(6,1) * t107 + mrSges(6,2) * t210;
t52 = mrSges(7,1) * t107 - mrSges(7,3) * t210;
t42 = pkin(5) * t210 + qJ(6) * t107 + t292;
t37 = t202 * pkin(5) - t39;
t35 = -qJ(6) * t202 + t40;
t19 = -pkin(5) * t90 - qJ(6) * t91 - qJD(6) * t138 + t113;
t18 = t269 * t60 - t272;
t17 = t196 * t60 + t54;
t7 = -pkin(5) * t253 - t8;
t6 = qJ(6) * t253 - qJD(6) * t202 + t9;
t5 = t75 * pkin(5) - t76 * qJ(6) - qJD(6) * t210 + t83;
t12 = [((0.3e1 / 0.2e1 * t194 + t234 + 0.2e1 * t273 + (-m(4) * t139 + m(5) * t128 + t261) * t191 - t329) * qJD(3) - t245 * t76 - t244 * t75 - t186 / 0.2e1 - t184 / 0.2e1 - t185 / 0.2e1 - t72 / 0.2e1 - t71 / 0.2e1 + t70 / 0.2e1 - t69 / 0.2e1 + t230 * t132 + t318) * t202 + (-t137 * t336 + t138 * t337) * t310 + (Ifges(6,2) * t308 - Ifges(7,3) * t307 + t295 * t335 + t304 * t336 + t326) * t90 + (-t137 * t4 - t138 * t3) * mrSges(6,3) + (-t1 * t137 + t138 * t2) * mrSges(7,2) + m(7) * (t1 * t35 + t10 * t7 + t11 * t6 + t19 * t29 + t2 * t37 + t5 * t73) + m(6) * (t102 * t113 + t13 * t8 + t14 * t9 + t146 * t83 + t3 * t39 + t4 * t40) + (Ifges(7,5) * t138 + Ifges(7,3) * t137) * t311 + (Ifges(6,4) * t138 - Ifges(6,2) * t137) * t312 + t137 * t315 + t137 * t316 + t138 * t340 + (t221 * t303 + t219 * t302 + t64 * t294 + t65 * t293 + (-t199 * t33 - t201 * t34) * mrSges(5,3) + (mrSges(4,3) + t222) * t133 + (-t201 * t97 / 0.2e1 + t98 * t294 + t128 * t223 + t218 * t299 + t220 * t298 + t217 * t296 + t214 * mrSges(5,3)) * qJD(4) + ((t241 * mrSges(4,1) + (t279 / 0.2e1 - t278 / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(4,4)) * t200 + t245 * t138 + t244 * t137 + (0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2) - t225) * t202) * qJD(1) - t230 * t140 + t233 + t317) * qJD(3) + (t78 + (m(4) + m(5)) * t133 - t178 * qJD(3)) * t191) * t200 + (t323 * t295 + t337 * t304 + t331 / 0.2e1 + Ifges(6,4) * t308 + Ifges(7,5) * t307 - t327) * t91 + m(5) * (t34 * t110 + t33 * t111 + t82 * t62 + t81 * t63) + t19 * t52 + t35 * t56 + t40 * t57 + t39 * t58 + t37 * t59 + t73 * t25 + t9 * t86 + t8 * t87 + t7 * t88 + t6 * t89 + t110 * t104 + t111 * t105 + t113 * t53 + t62 * t125 + t63 * t126 + t5 * (mrSges(7,1) * t137 - mrSges(7,3) * t138) + t83 * (mrSges(6,1) * t137 + mrSges(6,2) * t138) + t146 * t26; t286 * t91 + t285 * t90 + t288 * t138 + t287 * t137 + (-t78 + t321) * t202 + (t212 * qJD(4) + t213) * t200 + ((t125 * t201 - t126 * t199 + t178 - t242) * t202 + (t52 + t53 - t243 + t261) * t200) * qJD(3) + m(7) * (t1 * t138 - t10 * t90 + t11 * t91 + t2 * t137 - t202 * t5 + t253 * t29) + m(6) * (t102 * t253 + t13 * t90 - t3 * t137 + t4 * t138 + t14 * t91 - t202 * t83) + m(4) * (t132 * t200 - t133 * t202 + (-t139 * t200 + t140 * t202) * qJD(3)) + m(5) * ((-qJD(3) * t214 - t133) * t202 + (qJD(3) * t128 - t215 * qJD(4) + t216) * t200); t213 * pkin(8) + ((m(6) * t102 + t53) * t199 * pkin(4) + (-m(5) * t215 + t212) * pkin(8) + t203) * qJD(4) + t216 * mrSges(5,3) + (-t149 * t335 + t150 * t323) * t295 + (t162 * t337 + t209 * t336) * t310 + (-t130 * t336 + t131 * t337) * t305 + (-t149 * t336 + t150 * t337) * t304 + t330 * t52 + t331 * (-t131 / 0.2e1 + t150 / 0.2e1) + t332 * t88 + (t1 * t116 + t10 * t332 + t103 * t5 + t11 * t333 + t115 * t2 + t29 * t330) * m(7) + t333 * t89 + (-t130 * t335 + t131 * t323) * t296 + ((t234 - t273 + t341 + t329) * t202 + ((t282 / 0.2e1 + (t342 + Ifges(4,2) / 0.2e1) * t202) * qJD(1) + t140 * mrSges(4,3) + t233 - t317) * t200 + (t162 * t323 + t209 * t335 + t217) * t253 / 0.2e1) * qJD(1) + (-mrSges(7,1) * t260 + mrSges(7,3) * t259) * t29 + (-mrSges(6,1) * t260 - mrSges(6,2) * t259) * t102 - t261 * t140 + (Ifges(6,4) * t150 + Ifges(7,5) * t131 - Ifges(6,2) * t149 + Ifges(7,3) * t130) * t308 + t64 * t293 + m(6) * (-t115 * t3 + t116 * t4 - t13 * t93 + t14 * t94 + t193 * t83) - m(6) * (t102 * t117 + t13 * t30 + t14 * t31) - m(5) * (t100 * t81 + t101 * t82 + t128 * t140) + t287 * t115 + t288 * t116 + (t46 - t43) * (t130 / 0.2e1 - t149 / 0.2e1) + t162 * t340 + (Ifges(6,4) * t131 + Ifges(7,5) * t150 - Ifges(6,2) * t130 + Ifges(7,3) * t149) * t307 + t218 * t302 + t220 * t303 + (t13 * t259 + t14 * t260 - t162 * t3 + t209 * t4) * mrSges(6,3) + (t1 * t209 - t10 * t259 + t11 * t260 + t162 * t2) * mrSges(7,2) + (Ifges(7,5) * t162 - Ifges(7,3) * t209) * t311 + (Ifges(6,4) * t162 + Ifges(6,2) * t209) * t312 + t5 * (-mrSges(7,1) * t209 - mrSges(7,3) * t162) + t83 * (-mrSges(6,1) * t209 + mrSges(6,2) * t162) - t209 * t315 - t209 * t316 + m(5) * (-pkin(3) * t133 + pkin(8) * t216) + t338 * t86 + t339 * t87 - pkin(3) * t78 + t103 * t25 - t117 * t53 - t101 * t125 - t100 * t126 - t132 * mrSges(4,2) - t139 * t178 + t193 * t26 + t199 * t65 / 0.2e1 + (-mrSges(4,1) - t223) * t133; -(Ifges(6,4) * t307 + Ifges(7,5) * t308 + t323 * t296 + t337 * t305 + t327) * t107 + (-Ifges(6,2) * t307 + Ifges(7,3) * t308 - t335 * t296 - t336 * t305 + t326) * t210 + ((t196 * t4 + t269 * t3) * pkin(4) - t102 * t292 + t13 * t17 - t14 * t18) * m(6) + t57 * t291 + t186 - t53 * t292 + t184 + t185 + t58 * t240 + (t1 * t189 - t10 * t17 + t192 * t2 - t29 * t42 + (qJD(6) - t18) * t11) * m(7) + t72 + t71 - t70 + t69 + (Ifges(5,5) * t166 - Ifges(5,6) * t167) * t296 + t97 * t297 + (Ifges(5,1) * t166 - t275) * t298 + (-Ifges(5,2) * t167 + t159 + t98) * t299 + (t284 - t125) * t81 - t318 + (t283 + t126) * t82 + t285 * t17 - t286 * t18 + t331 * t307 - t42 * t52 + qJD(6) * t89 - t128 * (mrSges(5,1) * t167 + mrSges(5,2) * t166) + t189 * t56 + t192 * t59; t285 * t210 + t286 * t107 + (-t10 * t210 + t107 * t11 + t5) * m(7) + (t107 * t14 + t13 * t210 + t83) * m(6) - t321; t210 * t52 - t187 * t89 + 0.2e1 * (t2 / 0.2e1 + t29 * t304 + t11 * t296) * m(7) + t59;];
tauc  = t12(:);
