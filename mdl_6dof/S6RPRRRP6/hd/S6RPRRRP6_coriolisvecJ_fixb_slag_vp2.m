% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRP6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP6_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP6_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:14:16
% EndTime: 2019-03-09 06:14:42
% DurationCPUTime: 13.53s
% Computational Cost: add. (12244->555), mult. (31862->742), div. (0->0), fcn. (24180->8), ass. (0->234)
t231 = sin(pkin(10));
t309 = pkin(7) + qJ(2);
t221 = t309 * t231;
t211 = qJD(1) * t221;
t232 = cos(pkin(10));
t222 = t309 * t232;
t212 = qJD(1) * t222;
t235 = sin(qJ(3));
t238 = cos(qJ(3));
t168 = -t211 * t238 - t235 * t212;
t163 = -qJD(3) * pkin(3) - t168;
t210 = t231 * t238 + t232 * t235;
t201 = t210 * qJD(1);
t234 = sin(qJ(4));
t237 = cos(qJ(4));
t179 = qJD(3) * t237 - t201 * t234;
t116 = -pkin(4) * t179 + t163;
t180 = qJD(3) * t234 + t201 * t237;
t233 = sin(qJ(5));
t236 = cos(qJ(5));
t125 = t179 * t233 + t180 * t236;
t385 = -t231 * t235 + t238 * t232;
t200 = t385 * qJD(1);
t195 = qJD(4) - t200;
t268 = -pkin(2) * t232 - pkin(1);
t220 = qJD(1) * t268 + qJD(2);
t139 = -pkin(3) * t200 - pkin(8) * t201 + t220;
t169 = -t235 * t211 + t238 * t212;
t164 = qJD(3) * pkin(8) + t169;
t97 = t237 * t139 - t164 * t234;
t82 = -pkin(9) * t180 + t97;
t69 = pkin(4) * t195 + t82;
t98 = t139 * t234 + t164 * t237;
t83 = pkin(9) * t179 + t98;
t75 = t233 * t83;
t28 = t236 * t69 - t75;
t388 = qJ(6) * t125;
t18 = t28 - t388;
t190 = qJD(5) + t195;
t15 = pkin(5) * t190 + t18;
t77 = t236 * t83;
t29 = t233 * t69 + t77;
t261 = t236 * t179 - t180 * t233;
t357 = qJ(6) * t261;
t19 = t29 + t357;
t322 = -t190 / 0.2e1;
t330 = t125 / 0.2e1;
t331 = -t125 / 0.2e1;
t203 = t210 * qJD(3);
t192 = qJD(1) * t203;
t202 = t385 * qJD(3);
t191 = qJD(1) * t202;
t133 = qJD(4) * t179 + t191 * t237;
t134 = -qJD(4) * t180 - t191 * t234;
t57 = qJD(5) * t261 + t133 * t236 + t134 * t233;
t244 = t385 * qJD(2);
t126 = qJD(1) * t244 + qJD(3) * t168;
t147 = pkin(3) * t192 - pkin(8) * t191;
t47 = -qJD(4) * t98 - t126 * t234 + t237 * t147;
t26 = pkin(4) * t192 - pkin(9) * t133 + t47;
t277 = qJD(4) * t237;
t278 = qJD(4) * t234;
t46 = t237 * t126 + t139 * t277 + t234 * t147 - t164 * t278;
t32 = pkin(9) * t134 + t46;
t8 = -qJD(5) * t29 - t233 * t32 + t236 * t26;
t2 = pkin(5) * t192 - qJ(6) * t57 - qJD(6) * t125 + t8;
t58 = -qJD(5) * t125 - t133 * t233 + t134 * t236;
t275 = qJD(5) * t236;
t276 = qJD(5) * t233;
t7 = t233 * t26 + t236 * t32 + t69 * t275 - t276 * t83;
t3 = qJ(6) * t58 + qJD(6) * t261 + t7;
t345 = t8 * mrSges(6,1) + t2 * mrSges(7,1) - t7 * mrSges(6,2) - t3 * mrSges(7,2);
t376 = Ifges(6,6) + Ifges(7,6);
t378 = Ifges(6,5) + Ifges(7,5);
t389 = Ifges(6,3) + Ifges(7,3);
t350 = t192 * t389 + t376 * t58 + t378 * t57;
t380 = Ifges(6,1) + Ifges(7,1);
t379 = Ifges(6,4) + Ifges(7,4);
t395 = t379 * t261;
t372 = t380 * t125 + t378 * t190 + t395;
t377 = Ifges(6,2) + Ifges(7,2);
t393 = t125 * t379;
t373 = t376 * t190 + t377 * t261 + t393;
t381 = -t261 / 0.2e1;
t74 = -pkin(5) * t261 + qJD(6) + t116;
t398 = t345 + t350 + (-t125 * t376 + t378 * t261) * t322 + (t125 * t19 + t15 * t261) * mrSges(7,3) + (t125 * t29 + t261 * t28) * mrSges(6,3) - t116 * (mrSges(6,1) * t125 + mrSges(6,2) * t261) - t74 * (mrSges(7,1) * t125 + mrSges(7,2) * t261) + t373 * t330 + (-t125 * t377 + t372 + t395) * t381 + (t380 * t261 - t393) * t331;
t165 = pkin(3) * t201 - pkin(8) * t200;
t108 = t234 * t165 + t237 * t168;
t338 = -pkin(9) - pkin(8);
t267 = qJD(4) * t338;
t289 = t200 * t234;
t397 = pkin(9) * t289 + t234 * t267 - t108;
t107 = t237 * t165 - t168 * t234;
t312 = pkin(9) * t237;
t396 = -pkin(4) * t201 + t200 * t312 + t237 * t267 - t107;
t214 = t233 * t237 + t234 * t236;
t145 = t214 * t200;
t347 = qJD(4) + qJD(5);
t172 = t347 * t214;
t353 = t145 - t172;
t247 = t233 * t234 - t236 * t237;
t146 = t247 * t200;
t171 = t347 * t247;
t352 = t146 - t171;
t307 = Ifges(5,4) * t180;
t110 = t179 * Ifges(5,2) + t195 * Ifges(5,6) + t307;
t176 = Ifges(5,4) * t179;
t111 = t180 * Ifges(5,1) + t195 * Ifges(5,5) + t176;
t250 = t234 * t98 + t237 * t97;
t253 = Ifges(5,5) * t237 - Ifges(5,6) * t234;
t305 = Ifges(5,4) * t237;
t255 = -Ifges(5,2) * t234 + t305;
t306 = Ifges(5,4) * t234;
t257 = Ifges(5,1) * t237 - t306;
t258 = mrSges(5,1) * t234 + mrSges(5,2) * t237;
t314 = t237 / 0.2e1;
t315 = -t234 / 0.2e1;
t323 = t180 / 0.2e1;
t394 = t179 * t255 / 0.2e1 + t257 * t323 + t163 * t258 + t195 * t253 / 0.2e1 + t110 * t315 + t111 * t314 - t250 * mrSges(5,3);
t391 = -t200 * Ifges(4,2) / 0.2e1;
t390 = -t377 * t58 / 0.2e1 - t379 * t57 / 0.2e1 - t376 * t192 / 0.2e1;
t223 = t338 * t234;
t224 = t338 * t237;
t361 = t223 * t275 + t224 * t276 + t396 * t233 + t397 * t236;
t178 = t233 * t223 - t236 * t224;
t360 = -qJD(5) * t178 - t397 * t233 + t396 * t236;
t194 = Ifges(4,4) * t200;
t346 = t194 / 0.2e1 + t201 * Ifges(4,1) / 0.2e1;
t387 = t220 * mrSges(4,2) + Ifges(4,5) * qJD(3) + t346 + t394;
t386 = t202 * t234 + t210 * t277;
t374 = t378 * t192 + t379 * t58 + t380 * t57;
t371 = qJ(6) * t353 - qJD(6) * t247 + t361;
t370 = -pkin(5) * t201 - qJ(6) * t352 - qJD(6) * t214 + t360;
t129 = pkin(4) * t289 + t169;
t359 = pkin(4) * t278 - pkin(5) * t353 - t129;
t79 = -mrSges(6,1) * t261 + mrSges(6,2) * t125;
t358 = m(6) * t116 + t79;
t158 = t247 * t210;
t103 = -mrSges(7,2) * t190 + mrSges(7,3) * t261;
t104 = -mrSges(6,2) * t190 + mrSges(6,3) * t261;
t356 = -t104 - t103;
t105 = mrSges(7,1) * t190 - mrSges(7,3) * t125;
t106 = mrSges(6,1) * t190 - mrSges(6,3) * t125;
t355 = -t105 - t106;
t354 = Ifges(5,5) * t133 + Ifges(5,6) * t134;
t167 = -pkin(3) * t385 - pkin(8) * t210 + t268;
t175 = -t221 * t235 + t222 * t238;
t170 = t237 * t175;
t115 = t234 * t167 + t170;
t351 = -t238 * t221 - t222 * t235;
t349 = -t234 * t47 + t237 * t46;
t348 = t47 * mrSges(5,1) - t46 * mrSges(5,2);
t344 = (m(3) * qJ(2) + mrSges(3,3)) * (t231 ^ 2 + t232 ^ 2);
t343 = (Ifges(6,5) / 0.2e1 + Ifges(7,5) / 0.2e1) * t125 + (Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t190 + (Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1) * t261 - t19 * mrSges(7,2) - t29 * mrSges(6,2) - t98 * mrSges(5,2) + t195 * Ifges(5,3) + t180 * Ifges(5,5) + t179 * Ifges(5,6) - Ifges(4,6) * qJD(3) - t201 * Ifges(4,4) + t391 + t15 * mrSges(7,1) + t220 * mrSges(4,1) + t28 * mrSges(6,1) + t97 * mrSges(5,1) - t376 * t381 - t378 * t331 - t389 * t322;
t342 = t57 / 0.2e1;
t341 = t58 / 0.2e1;
t333 = t261 / 0.2e1;
t329 = t133 / 0.2e1;
t328 = t134 / 0.2e1;
t325 = -t179 / 0.2e1;
t324 = -t180 / 0.2e1;
t321 = t190 / 0.2e1;
t320 = t192 / 0.2e1;
t319 = -t195 / 0.2e1;
t50 = -mrSges(7,2) * t192 + mrSges(7,3) * t58;
t51 = -mrSges(6,2) * t192 + mrSges(6,3) * t58;
t308 = t51 + t50;
t35 = t236 * t82 - t75;
t114 = t237 * t167 - t175 * t234;
t92 = -pkin(4) * t385 - t210 * t312 + t114;
t287 = t210 * t234;
t99 = -pkin(9) * t287 + t115;
t45 = t233 * t92 + t236 * t99;
t245 = t210 * qJD(2);
t127 = qJD(1) * t245 + qJD(3) * t169;
t292 = t127 * t351;
t280 = qJD(3) * mrSges(4,1) + mrSges(5,1) * t179 - mrSges(5,2) * t180 - mrSges(4,3) * t201;
t269 = Ifges(5,3) * t192 + t354;
t228 = -pkin(4) * t237 - pkin(3);
t16 = -t58 * mrSges(7,1) + t57 * mrSges(7,2);
t34 = -t233 * t82 - t77;
t44 = -t233 * t99 + t236 * t92;
t263 = t192 * mrSges(4,1) + t191 * mrSges(4,2);
t140 = qJD(3) * t351 + t244;
t166 = pkin(3) * t203 - pkin(8) * t202;
t262 = -t140 * t234 + t237 * t166;
t177 = t236 * t223 + t224 * t233;
t142 = pkin(4) * t287 - t351;
t259 = mrSges(5,1) * t237 - mrSges(5,2) * t234;
t256 = Ifges(5,1) * t234 + t305;
t254 = Ifges(5,2) * t237 + t306;
t252 = Ifges(5,5) * t234 + Ifges(5,6) * t237;
t251 = -t234 * t46 - t237 * t47;
t249 = t234 * t97 - t237 * t98;
t137 = -mrSges(5,2) * t195 + mrSges(5,3) * t179;
t138 = mrSges(5,1) * t195 - mrSges(5,3) * t180;
t248 = t137 * t237 - t138 * t234;
t38 = -t202 * t312 + pkin(4) * t203 + (-t170 + (pkin(9) * t210 - t167) * t234) * qJD(4) + t262;
t59 = t237 * t140 + t234 * t166 + t167 * t277 - t175 * t278;
t43 = -pkin(9) * t386 + t59;
t9 = t233 * t38 + t236 * t43 + t92 * t275 - t276 * t99;
t10 = -qJD(5) * t45 - t233 * t43 + t236 * t38;
t141 = qJD(3) * t175 + t245;
t102 = pkin(4) * t386 + t141;
t88 = -pkin(4) * t134 + t127;
t227 = pkin(4) * t236 + pkin(5);
t189 = pkin(5) * t247 + t228;
t181 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t200;
t157 = t214 * t210;
t150 = -qJ(6) * t247 + t178;
t149 = -qJ(6) * t214 + t177;
t113 = -mrSges(5,2) * t192 + mrSges(5,3) * t134;
t112 = mrSges(5,1) * t192 - mrSges(5,3) * t133;
t101 = pkin(4) * t180 + pkin(5) * t125;
t100 = pkin(5) * t157 + t142;
t91 = -mrSges(5,1) * t134 + mrSges(5,2) * t133;
t78 = -mrSges(7,1) * t261 + mrSges(7,2) * t125;
t73 = t133 * Ifges(5,1) + t134 * Ifges(5,4) + t192 * Ifges(5,5);
t72 = t133 * Ifges(5,4) + t134 * Ifges(5,2) + t192 * Ifges(5,6);
t71 = t158 * t347 - t214 * t202;
t70 = -t172 * t210 - t202 * t247;
t60 = -qJD(4) * t115 + t262;
t49 = mrSges(6,1) * t192 - mrSges(6,3) * t57;
t48 = mrSges(7,1) * t192 - mrSges(7,3) * t57;
t41 = -pkin(5) * t71 + t102;
t33 = -qJ(6) * t157 + t45;
t30 = -pkin(5) * t385 + qJ(6) * t158 + t44;
t23 = -pkin(5) * t58 + t88;
t21 = t35 - t388;
t20 = t34 - t357;
t17 = -mrSges(6,1) * t58 + mrSges(6,2) * t57;
t5 = qJ(6) * t71 - qJD(6) * t157 + t9;
t4 = pkin(5) * t203 - qJ(6) * t70 + qJD(6) * t158 + t10;
t1 = [t372 * t70 / 0.2e1 + t373 * t71 / 0.2e1 + (Ifges(4,1) * t191 - Ifges(4,4) * t192 + t255 * t328 + t257 * t329 + t72 * t315 + t73 * t314 + (mrSges(4,3) + t258) * t127 + t251 * mrSges(5,3) + (t163 * t259 + t254 * t325 + t256 * t324 + t252 * t319 + t111 * t315 - t237 * t110 / 0.2e1 + t249 * mrSges(5,3)) * qJD(4)) * t210 + t268 * t263 + (-t157 * t377 - t158 * t379) * t341 + (t379 * t71 + t380 * t70) * t330 + (-t157 * t379 - t158 * t380) * t342 - t374 * t158 / 0.2e1 + (t376 * t71 + t378 * t70) * t321 + (-t157 * t376 - t158 * t378 + t210 * t253) * t320 + (t377 * t71 + t379 * t70) * t333 + 0.2e1 * t344 * qJD(2) * qJD(1) + m(5) * (t114 * t47 + t115 * t46 + t141 * t163 + t59 * t98 + t60 * t97 - t292) + m(4) * (t126 * t175 + t140 * t169 - t141 * t168 - t292) + t140 * t181 + (t391 + t343) * t203 + (-t168 * t202 - t169 * t203 - t175 * t192) * mrSges(4,3) - (mrSges(4,3) * t191 + t91) * t351 - (-mrSges(4,3) * t126 - Ifges(4,4) * t191 + t320 * t389 + t376 * t341 + t378 * t342 + t345 + t348) * t385 + (t346 + t387) * t202 - t280 * t141 - (t269 + t350 + t354) * t385 / 0.2e1 + m(7) * (t100 * t23 + t15 * t4 + t19 * t5 + t2 * t30 + t3 * t33 + t41 * t74) + m(6) * (t10 * t28 + t102 * t116 + t142 * t88 + t29 * t9 + t44 * t8 + t45 * t7) - (Ifges(4,2) + Ifges(5,3) / 0.2e1) * t192 * t385 + t142 * t17 + t59 * t137 + t60 * t138 + t114 * t112 + t115 * t113 + t116 * (-mrSges(6,1) * t71 + mrSges(6,2) * t70) + t4 * t105 + t10 * t106 + t100 * t16 + t102 * t79 + t5 * t103 + t9 * t104 + t74 * (-mrSges(7,1) * t71 + mrSges(7,2) * t70) + t41 * t78 + t45 * t51 + t30 * t48 + t44 * t49 + t33 * t50 + (-t157 * t7 + t158 * t8 - t28 * t70 + t29 * t71) * mrSges(6,3) + (-t15 * t70 - t157 * t3 + t158 * t2 + t19 * t71) * mrSges(7,3) + t23 * (mrSges(7,1) * t157 - mrSges(7,2) * t158) + t88 * (mrSges(6,1) * t157 - mrSges(6,2) * t158) + t157 * t390; t263 + (-t181 - t248) * t200 + t248 * qJD(4) + t308 * t214 - (t48 + t49) * t247 + (-t78 - t79 + t280) * t201 - m(4) * (-t168 * t201 + t169 * t200) + t234 * t113 + t237 * t112 - t355 * t353 - t356 * t352 - t344 * qJD(1) ^ 2 + (t15 * t353 + t19 * t352 - t2 * t247 - t201 * t74 + t214 * t3) * m(7) + (-t116 * t201 + t214 * t7 - t247 * t8 + t28 * t353 + t29 * t352) * m(6) + (-t163 * t201 - t195 * t249 - t251) * m(5); t370 * t105 + (t149 * t2 + t15 * t370 + t150 * t3 + t189 * t23 + t19 * t371 + t359 * t74) * m(7) + t371 * t103 + (-t145 * t379 - t146 * t380) * t331 + (-mrSges(4,1) - t259) * t127 + t374 * t214 / 0.2e1 + (-t145 * t376 - t146 * t378) * t322 + (pkin(4) * t234 * t358 + t394) * qJD(4) + (t169 * mrSges(4,3) - t343) * t201 + (-t112 * t234 + t113 * t237 + m(5) * t349 + (-m(5) * t250 - t234 * t137 - t237 * t138) * qJD(4)) * pkin(8) + t349 * mrSges(5,3) + t23 * (mrSges(7,1) * t247 + mrSges(7,2) * t214) + (t214 * t380 - t247 * t379) * t342 + (t214 * t378 - t247 * t376 + t252) * t320 + (t214 * t379 - t247 * t377) * t341 + t88 * (mrSges(6,1) * t247 + mrSges(6,2) * t214) + (-t145 * t377 - t146 * t379) * t381 + t372 * (t146 / 0.2e1 - t171 / 0.2e1) + t373 * (t145 / 0.2e1 - t172 / 0.2e1) + (-t171 * t380 - t172 * t379) * t330 + (-t171 * t378 - t172 * t376) * t321 + (-t171 * t379 - t172 * t377) * t333 + t359 * t78 + t360 * t106 + t361 * t104 + (-t116 * t129 + t177 * t8 + t178 * t7 + t228 * t88 + t28 * t360 + t29 * t361) * m(6) + t189 * t16 - t168 * t181 + t177 * t49 + t178 * t51 + (-t194 / 0.2e1 + t168 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1) * t201 - t387) * t200 + (-t214 * t8 - t247 * t7 - t28 * t352 + t29 * t353) * mrSges(6,3) + (-t15 * t352 + t19 * t353 - t2 * t214 - t247 * t3) * mrSges(7,3) + (-mrSges(6,1) * t353 + mrSges(6,2) * t352) * t116 + (-mrSges(7,1) * t353 + mrSges(7,2) * t352) * t74 + t280 * t169 + t254 * t328 + t256 * t329 + t72 * t314 + (-pkin(3) * t127 - t107 * t97 - t108 * t98 - t163 * t169) * m(5) + t149 * t48 + t150 * t50 - t108 * t137 - t107 * t138 - t129 * t79 - t126 * mrSges(4,2) - pkin(3) * t91 + Ifges(4,5) * t191 - Ifges(4,6) * t192 + t228 * t17 + t234 * t73 / 0.2e1 + t247 * t390; (-Ifges(5,2) * t180 + t111 + t176) * t325 + t269 + (t236 * t49 + t308 * t233 + (t233 * t355 - t236 * t356) * qJD(5) + m(6) * (t233 * t7 + t236 * t8 + t275 * t29 - t276 * t28) - t358 * t180) * pkin(4) + (t2 * t227 - t101 * t74 - t15 * t20 - t19 * t21 + (-t15 * t276 + t19 * t275 + t233 * t3) * pkin(4)) * m(7) - t163 * (mrSges(5,1) * t180 + mrSges(5,2) * t179) + t348 + t110 * t323 + (Ifges(5,1) * t179 - t307) * t324 + (Ifges(5,5) * t179 - Ifges(5,6) * t180) * t319 + (t179 * t97 + t180 * t98) * mrSges(5,3) - m(6) * (t28 * t34 + t29 * t35) - t97 * t137 + t98 * t138 - t20 * t105 - t34 * t106 - t101 * t78 - t21 * t103 - t35 * t104 + t227 * t48 + t398; (-(-t15 + t18) * t19 + (-t125 * t74 + t2) * pkin(5)) * m(7) + (-t125 * t78 + t48) * pkin(5) + t19 * t105 + t29 * t106 - t18 * t103 - t28 * t104 + t398; -t261 * t103 + t125 * t105 + 0.2e1 * (t23 / 0.2e1 + t19 * t381 + t15 * t330) * m(7) + t16;];
tauc  = t1(:);
