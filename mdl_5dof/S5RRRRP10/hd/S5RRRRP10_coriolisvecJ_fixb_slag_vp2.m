% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRP10_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP10_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP10_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP10_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP10_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:08:10
% EndTime: 2019-12-31 22:08:45
% DurationCPUTime: 15.28s
% Computational Cost: add. (7420->614), mult. (20057->855), div. (0->0), fcn. (14818->8), ass. (0->253)
t360 = Ifges(6,4) + Ifges(5,4);
t361 = Ifges(5,1) + Ifges(6,1);
t359 = Ifges(6,5) + Ifges(5,5);
t358 = Ifges(6,2) + Ifges(5,2);
t357 = Ifges(6,6) + Ifges(5,6);
t218 = cos(pkin(5));
t209 = qJD(1) * t218 + qJD(2);
t220 = sin(qJ(3));
t223 = cos(qJ(3));
t221 = sin(qJ(2));
t217 = sin(pkin(5));
t276 = qJD(1) * t217;
t263 = t221 * t276;
t164 = t209 * t220 + t223 * t263;
t224 = cos(qJ(2));
t262 = t224 * t276;
t203 = qJD(3) - t262;
t219 = sin(qJ(4));
t222 = cos(qJ(4));
t124 = -t164 * t219 + t203 * t222;
t163 = t209 * t223 - t220 * t263;
t158 = qJD(4) - t163;
t125 = t164 * t222 + t203 * t219;
t375 = t360 * t125;
t353 = t358 * t124 + t357 * t158 + t375;
t377 = -t353 / 0.2e1;
t376 = t360 * t124;
t352 = t361 * t125 + t359 * t158 + t376;
t281 = t223 * t224;
t152 = (t219 * t221 + t222 * t281) * t276;
t202 = -pkin(3) * t223 - pkin(9) * t220 - pkin(2);
t282 = t222 * t223;
t214 = pkin(8) * t282;
t269 = qJD(5) * t222;
t255 = pkin(3) * t220 - pkin(9) * t223;
t198 = t255 * qJD(3);
t273 = qJD(3) * t220;
t314 = pkin(8) * t219;
t278 = t222 * t198 + t273 * t314;
t287 = qJ(5) * t220;
t315 = pkin(4) * t220;
t317 = pkin(1) * t224;
t268 = t218 * t317;
t181 = -pkin(7) * t263 + qJD(1) * t268;
t236 = (pkin(2) * t221 - pkin(8) * t224) * t217;
t182 = qJD(1) * t236;
t120 = t223 * t181 + t220 * t182;
t106 = pkin(9) * t263 + t120;
t213 = t218 * t221 * pkin(1);
t284 = t217 * t224;
t123 = (t213 + (pkin(7) + t255) * t284) * qJD(1);
t53 = -t219 * t106 + t222 * t123;
t374 = t152 * qJ(5) - t262 * t315 - t53 - t220 * t269 + (-qJ(5) * t282 + t315) * qJD(3) + (-t214 + (-t202 + t287) * t219) * qJD(4) + t278;
t151 = (-t219 * t281 + t221 * t222) * t276;
t270 = qJD(4) * t222;
t279 = t219 * t198 + t202 * t270;
t283 = t220 * t222;
t54 = t222 * t106 + t219 * t123;
t373 = -qJ(5) * t151 - t54 + (-pkin(8) * qJD(3) - qJ(5) * qJD(4)) * t283 + (-qJD(5) * t220 + (-pkin(8) * qJD(4) - qJ(5) * qJD(3)) * t223) * t219 + t279;
t119 = -t220 * t181 + t182 * t223;
t105 = -pkin(3) * t263 - t119;
t272 = qJD(3) * t223;
t372 = pkin(8) * t272 - t105 + (t219 * t272 + t220 * t270 + t151) * pkin(4);
t371 = t360 * t222;
t370 = t360 * t219;
t293 = t163 * Ifges(4,2);
t146 = -t209 * pkin(2) - t181;
t83 = -t163 * pkin(3) - t164 * pkin(9) + t146;
t277 = pkin(7) * t284 + t213;
t184 = t277 * qJD(1);
t147 = t209 * pkin(8) + t184;
t178 = (-pkin(2) * t224 - pkin(8) * t221 - pkin(1)) * t217;
t156 = qJD(1) * t178;
t99 = t147 * t223 + t156 * t220;
t86 = pkin(9) * t203 + t99;
t24 = -t219 * t86 + t222 * t83;
t17 = -qJ(5) * t125 + t24;
t15 = pkin(4) * t158 + t17;
t25 = t219 * t83 + t222 * t86;
t18 = qJ(5) * t124 + t25;
t265 = Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1;
t266 = Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1;
t267 = Ifges(6,5) / 0.2e1 + Ifges(5,5) / 0.2e1;
t368 = -Ifges(4,6) * t203 / 0.2e1 - Ifges(4,4) * t164 / 0.2e1;
t367 = -t293 / 0.2e1 + t368;
t44 = Ifges(6,5) * t125 + Ifges(6,6) * t124 + Ifges(6,3) * t158;
t45 = Ifges(5,5) * t125 + Ifges(5,6) * t124 + Ifges(5,3) * t158;
t345 = t266 * t124 + t267 * t125 + t265 * t158 - t18 * mrSges(6,2) - t25 * mrSges(5,2) - t99 * mrSges(4,3) + t44 / 0.2e1 + t45 / 0.2e1 + t367 + t146 * mrSges(4,1) + t15 * mrSges(6,1) + t24 * mrSges(5,1) + t368;
t369 = -t345 + t293 / 0.2e1;
t274 = qJD(2) * t224;
t261 = t220 * t274;
t131 = t209 * t273 + (t221 * t272 + t261) * t276;
t260 = t223 * t274;
t130 = t209 * t272 + (-t221 * t273 + t260) * t276;
t275 = qJD(2) * t217;
t258 = qJD(1) * t275;
t256 = t221 * t258;
t63 = qJD(4) * t124 + t130 * t222 + t219 * t256;
t64 = -qJD(4) * t125 - t130 * t219 + t222 * t256;
t355 = t357 * t131 + t358 * t64 + t360 * t63;
t354 = t359 * t131 + t360 * t64 + t361 * t63;
t366 = -t357 * t219 + t359 * t222;
t365 = -t358 * t219 + t371;
t364 = t361 * t222 - t370;
t363 = -Ifges(3,6) * t209 / 0.2e1;
t10 = Ifges(5,5) * t63 + Ifges(5,6) * t64 + Ifges(5,3) * t131;
t9 = Ifges(6,5) * t63 + Ifges(6,6) * t64 + Ifges(6,3) * t131;
t362 = t10 + t9;
t356 = Ifges(5,3) + Ifges(6,3);
t285 = t217 * t221;
t210 = pkin(7) * t285;
t193 = -t210 + t268;
t172 = t219 * t202 + t214;
t157 = Ifges(4,4) * t163;
t301 = Ifges(4,5) * t203;
t308 = Ifges(4,1) * t164;
t95 = t157 + t301 + t308;
t98 = -t220 * t147 + t156 * t223;
t351 = t98 * mrSges(4,3) - t95 / 0.2e1 - t301 / 0.2e1 - t146 * mrSges(4,2) - t157 / 0.2e1;
t350 = t359 * t219 + t357 * t222;
t349 = t358 * t222 + t370;
t348 = t361 * t219 + t371;
t271 = qJD(4) * t219;
t183 = qJD(2) * t236;
t173 = qJD(1) * t183;
t185 = t193 * qJD(2);
t174 = qJD(1) * t185;
t42 = -t147 * t273 + t156 * t272 + t220 * t173 + t223 * t174;
t35 = pkin(9) * t256 + t42;
t186 = t277 * qJD(2);
t175 = qJD(1) * t186;
t57 = t131 * pkin(3) - t130 * pkin(9) + t175;
t5 = t219 * t57 + t222 * t35 + t83 * t270 - t271 * t86;
t6 = -qJD(4) * t25 - t219 * t35 + t222 * t57;
t347 = -t219 * t6 + t222 * t5;
t237 = t25 * t219 + t24 * t222;
t251 = mrSges(6,1) * t219 + mrSges(6,2) * t222;
t253 = mrSges(5,1) * t219 + mrSges(5,2) * t222;
t318 = t222 / 0.2e1;
t321 = -t219 / 0.2e1;
t326 = t158 / 0.2e1;
t335 = t125 / 0.2e1;
t337 = t124 / 0.2e1;
t85 = -pkin(3) * t203 - t98;
t41 = -pkin(4) * t124 + qJD(5) + t85;
t346 = t237 * mrSges(5,3) + (t15 * t222 + t18 * t219) * mrSges(6,3) - t251 * t41 - t253 * t85 - t365 * t337 - t364 * t335 - t366 * t326 - t353 * t321 - t352 * t318;
t344 = t63 / 0.2e1;
t343 = t64 / 0.2e1;
t340 = pkin(1) * mrSges(3,1);
t339 = pkin(1) * mrSges(3,2);
t338 = -t124 / 0.2e1;
t336 = -t125 / 0.2e1;
t334 = t131 / 0.2e1;
t327 = -t158 / 0.2e1;
t189 = -t218 * t223 + t220 * t285;
t325 = -t189 / 0.2e1;
t190 = t218 * t220 + t223 * t285;
t323 = t190 / 0.2e1;
t322 = t218 / 0.2e1;
t316 = pkin(4) * t219;
t311 = t42 * mrSges(4,2);
t43 = -t147 * t272 - t156 * t273 + t173 * t223 - t220 * t174;
t310 = t43 * mrSges(4,1);
t309 = -qJ(5) - pkin(9);
t307 = Ifges(3,4) * t221;
t296 = t130 * Ifges(4,1);
t295 = t130 * Ifges(4,4);
t294 = t131 * Ifges(4,4);
t176 = t210 + (-pkin(2) - t317) * t218;
t102 = t189 * pkin(3) - t190 * pkin(9) + t176;
t177 = pkin(8) * t218 + t277;
t116 = t223 * t177 + t220 * t178;
t104 = -pkin(9) * t284 + t116;
t38 = t219 * t102 + t222 * t104;
t133 = mrSges(4,1) * t203 - mrSges(4,3) * t164;
t74 = -mrSges(5,1) * t124 + mrSges(5,2) * t125;
t288 = t133 - t74;
t114 = pkin(3) * t164 - pkin(9) * t163;
t40 = t219 * t114 + t222 * t98;
t286 = t163 * t219;
t280 = -mrSges(3,1) * t209 - mrSges(4,1) * t163 + mrSges(4,2) * t164 + mrSges(3,3) * t263;
t264 = Ifges(4,5) * t130 - Ifges(4,6) * t131 + Ifges(4,3) * t256;
t259 = t221 * t275;
t19 = -t64 * mrSges(6,1) + t63 * mrSges(6,2);
t257 = qJD(4) * t309;
t37 = t222 * t102 - t104 * t219;
t39 = t222 * t114 - t219 * t98;
t115 = -t220 * t177 + t178 * t223;
t103 = pkin(3) * t284 - t115;
t254 = mrSges(5,1) * t222 - mrSges(5,2) * t219;
t252 = mrSges(6,1) * t222 - mrSges(6,2) * t219;
t66 = -t177 * t272 - t178 * t273 + t183 * t223 - t220 * t185;
t139 = -t190 * t219 - t222 * t284;
t235 = -t190 * t222 + t219 * t284;
t65 = -t177 * t273 + t178 * t272 + t220 * t183 + t223 * t185;
t51 = pkin(9) * t259 + t65;
t137 = qJD(3) * t190 + t217 * t261;
t138 = -qJD(3) * t189 + t217 * t260;
t75 = t137 * pkin(3) - t138 * pkin(9) + t186;
t7 = t102 * t270 - t104 * t271 + t219 * t75 + t222 * t51;
t206 = Ifges(3,4) * t262;
t233 = -t181 * mrSges(3,3) + Ifges(3,1) * t263 / 0.2e1 + t206 / 0.2e1 + t209 * Ifges(3,5);
t1 = pkin(4) * t131 - qJ(5) * t63 - qJD(5) * t125 + t6;
t2 = qJ(5) * t64 + qJD(5) * t124 + t5;
t232 = -t6 * mrSges(5,1) - t1 * mrSges(6,1) + t5 * mrSges(5,2) + t2 * mrSges(6,2);
t52 = -pkin(3) * t259 - t66;
t8 = -qJD(4) * t38 - t219 * t51 + t222 * t75;
t36 = -pkin(3) * t256 - t43;
t230 = -t308 / 0.2e1 + t351;
t229 = t98 * mrSges(4,1) + t363 - (Ifges(3,2) * t224 + t307) * t276 / 0.2e1 + t203 * Ifges(4,3) + t164 * Ifges(4,5) + t163 * Ifges(4,6) - t184 * mrSges(3,3) - t99 * mrSges(4,2);
t216 = -pkin(4) * t222 - pkin(3);
t205 = t309 * t222;
t204 = t309 * t219;
t201 = Ifges(3,5) * t224 * t258;
t199 = (pkin(8) + t316) * t220;
t197 = t222 * t202;
t188 = -qJD(5) * t219 + t222 * t257;
t187 = t219 * t257 + t269;
t180 = -t209 * mrSges(3,2) + mrSges(3,3) * t262;
t171 = -t223 * t314 + t197;
t143 = -t219 * t287 + t172;
t134 = -qJ(5) * t283 + t197 + (-pkin(4) - t314) * t223;
t132 = -mrSges(4,2) * t203 + mrSges(4,3) * t163;
t113 = -qJD(4) * t172 + t278;
t112 = (-t222 * t273 - t223 * t271) * pkin(8) + t279;
t110 = -mrSges(4,2) * t256 - mrSges(4,3) * t131;
t109 = mrSges(4,1) * t256 - mrSges(4,3) * t130;
t90 = mrSges(5,1) * t158 - mrSges(5,3) * t125;
t89 = mrSges(6,1) * t158 - mrSges(6,3) * t125;
t88 = -mrSges(5,2) * t158 + mrSges(5,3) * t124;
t87 = -mrSges(6,2) * t158 + mrSges(6,3) * t124;
t82 = qJD(4) * t139 + t138 * t222 + t219 * t259;
t81 = qJD(4) * t235 - t138 * t219 + t222 * t259;
t79 = mrSges(4,1) * t131 + mrSges(4,2) * t130;
t73 = -mrSges(6,1) * t124 + mrSges(6,2) * t125;
t70 = pkin(4) * t286 + t99;
t69 = -pkin(4) * t139 + t103;
t68 = Ifges(4,5) * t256 - t294 + t296;
t67 = -t131 * Ifges(4,2) + Ifges(4,6) * t256 + t295;
t32 = -mrSges(5,2) * t131 + mrSges(5,3) * t64;
t31 = -mrSges(6,2) * t131 + mrSges(6,3) * t64;
t30 = mrSges(5,1) * t131 - mrSges(5,3) * t63;
t29 = mrSges(6,1) * t131 - mrSges(6,3) * t63;
t27 = -qJ(5) * t286 + t40;
t26 = qJ(5) * t139 + t38;
t23 = -qJ(5) * t163 * t222 + pkin(4) * t164 + t39;
t22 = pkin(4) * t189 + qJ(5) * t235 + t37;
t21 = -pkin(4) * t81 + t52;
t20 = -mrSges(5,1) * t64 + mrSges(5,2) * t63;
t16 = -pkin(4) * t64 + t36;
t4 = qJ(5) * t81 + qJD(5) * t139 + t7;
t3 = pkin(4) * t137 - qJ(5) * t82 + qJD(5) * t235 + t8;
t11 = [t2 * (-mrSges(6,2) * t189 + mrSges(6,3) * t139) + t5 * (-mrSges(5,2) * t189 + mrSges(5,3) * t139) + t176 * t79 + t185 * t180 + t146 * (mrSges(4,1) * t137 + mrSges(4,2) * t138) + t24 * (mrSges(5,1) * t137 - mrSges(5,3) * t82) + t15 * (mrSges(6,1) * t137 - mrSges(6,3) * t82) + t25 * (-mrSges(5,2) * t137 + mrSges(5,3) * t81) + t18 * (-mrSges(6,2) * t137 + mrSges(6,3) * t81) + t138 * t95 / 0.2e1 + t65 * t132 + t66 * t133 + t163 * (Ifges(4,4) * t138 - Ifges(4,2) * t137) / 0.2e1 + t137 * t367 + m(6) * (t1 * t22 + t15 * t3 + t16 * t69 + t18 * t4 + t2 * t26 + t21 * t41) + m(5) * (t103 * t36 + t24 * t8 + t25 * t7 + t37 * t6 + t38 * t5 + t52 * t85) + m(4) * (t115 * t43 + t116 * t42 + t146 * t186 + t175 * t176 + t65 * t99 + t66 * t98) + (t137 * t359 + t360 * t81 + t361 * t82) * t335 + t352 * t82 / 0.2e1 + t353 * t81 / 0.2e1 + t164 * (Ifges(4,1) * t138 - Ifges(4,4) * t137) / 0.2e1 + (t45 + t44) * t137 / 0.2e1 + m(3) * (t174 * t277 - t175 * t193 - t181 * t186 + t184 * t185) + ((Ifges(3,5) * t322 - t193 * mrSges(3,3) + (-0.2e1 * t339 + 0.3e1 / 0.2e1 * Ifges(3,4) * t224) * t217) * t224 + (Ifges(4,5) * t323 + Ifges(4,6) * t325 - Ifges(3,6) * t218 - t277 * mrSges(3,3) + (-0.2e1 * t340 - 0.3e1 / 0.2e1 * t307 + (-Ifges(4,3) / 0.2e1 + 0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t224) * t217) * t221) * t258 + (-t137 * t99 - t138 * t98 - t189 * t42 - t190 * t43) * mrSges(4,3) + (t139 * t360 + t189 * t359 - t235 * t361) * t344 + t1 * (mrSges(6,1) * t189 + mrSges(6,3) * t235) + t6 * (mrSges(5,1) * t189 + mrSges(5,3) * t235) + t36 * (-mrSges(5,1) * t139 - mrSges(5,2) * t235) + t16 * (-mrSges(6,1) * t139 - mrSges(6,2) * t235) - t354 * t235 / 0.2e1 + (t139 * t357 + t189 * t356 - t235 * t359) * t334 + (t139 * t358 + t189 * t357 - t235 * t360) * t343 + t115 * t109 + t116 * t110 + t103 * t20 + t3 * t89 + t8 * t90 + t41 * (-mrSges(6,1) * t81 + mrSges(6,2) * t82) + t85 * (-mrSges(5,1) * t81 + mrSges(5,2) * t82) + t4 * t87 + t7 * t88 + t69 * t19 + t21 * t73 + t52 * t74 + t37 * t30 + t38 * t32 + t22 * t29 + t26 * t31 + t355 * t139 / 0.2e1 + (t137 * t356 + t357 * t81 + t359 * t82) * t326 + (t137 * t357 + t358 * t81 + t360 * t82) * t337 - t284 * t310 - t264 * t284 / 0.2e1 - t131 * (Ifges(4,4) * t190 - Ifges(4,2) * t189 - Ifges(4,6) * t284) / 0.2e1 + t130 * (Ifges(4,1) * t190 - Ifges(4,4) * t189 - Ifges(4,5) * t284) / 0.2e1 + t174 * (-t218 * mrSges(3,2) + mrSges(3,3) * t284) + (-mrSges(3,1) * t218 + mrSges(4,1) * t189 + mrSges(4,2) * t190 + mrSges(3,3) * t285) * t175 + t280 * t186 + t203 * (Ifges(4,5) * t138 - Ifges(4,6) * t137) / 0.2e1 + t362 * t189 / 0.2e1 + (t233 * t224 + (t363 + t229) * t221) * t275 + t284 * t311 + t201 * t322 + t68 * t323 + t67 * t325; t199 * t19 - t174 * mrSges(3,2) - t175 * mrSges(3,1) - t181 * t180 + t171 * t30 + t172 * t32 - t85 * (-mrSges(5,1) * t151 + mrSges(5,2) * t152) - t41 * (-mrSges(6,1) * t151 + mrSges(6,2) * t152) + t143 * t31 - t120 * t132 - t119 * t133 + t134 * t29 - m(5) * (t105 * t85 + t24 * t53 + t25 * t54) - m(4) * (t119 * t98 + t120 * t99 + t146 * t184) + (t15 * t152 - t151 * t18) * mrSges(6,3) + (t113 - t53) * t90 + (-t25 * t151 + t24 * t152) * mrSges(5,3) - t105 * t74 + t373 * t87 + (t1 * t134 + t143 * t2 + t374 * t15 + t16 * t199 + t373 * t18 + t372 * t41) * m(6) + t374 * t89 - t352 * t152 / 0.2e1 + t372 * t73 + (t112 - t54) * t88 + (((-m(4) * t99 - t132) * pkin(8) - t369) * t220 + ((-m(4) * t98 + m(5) * t85 - t288) * pkin(8) - t346 - t230) * t223) * qJD(3) + ((qJD(2) * (Ifges(4,5) * t220 + Ifges(4,6) * t223) / 0.2e1 + (t340 + t307 / 0.2e1) * t276 + (-qJD(2) + t209 / 0.2e1) * Ifges(3,6) - t229) * t221 + ((t339 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t221) * t276 - t206 / 0.2e1 + t230 * t223 + t369 * t220 - t233) * t224) * t276 + (-m(4) * t175 - t79) * pkin(2) + (t42 * mrSges(4,3) + t67 / 0.2e1 - t9 / 0.2e1 - t10 / 0.2e1 - t175 * mrSges(4,1) + t295 / 0.2e1 - t266 * t64 - t267 * t63 + (m(4) * t42 + t110) * pkin(8) + (-Ifges(4,2) / 0.2e1 - t265) * t131 + t232) * t223 - t280 * t184 + t201 + (t68 / 0.2e1 + t36 * t253 + t16 * t251 + t175 * mrSges(4,2) - t43 * mrSges(4,3) - t294 / 0.2e1 + t296 / 0.2e1 + (-t1 * t222 - t2 * t219) * mrSges(6,3) + (-t219 * t5 - t222 * t6) * mrSges(5,3) + (-m(4) * t43 + m(5) * t36 - t109 + t20) * pkin(8) + (t85 * t254 + t41 * t252 + (t15 * t219 - t18 * t222) * mrSges(6,3) + (t219 * t24 - t222 * t25) * mrSges(5,3) + t349 * t338 + t348 * t336 + t350 * t327 + t222 * t377) * qJD(4) + t364 * t344 + t365 * t343 + t366 * t334 + (qJD(4) * t352 + t355) * t321 + t354 * t318) * t220 + t151 * t377 + m(5) * (t112 * t25 + t113 * t24 + t171 * t6 + t172 * t5) + (Ifges(6,5) * t152 + Ifges(6,6) * t151) * t327 + (Ifges(5,5) * t152 + Ifges(5,6) * t151) * t327 + (Ifges(6,1) * t152 + Ifges(6,4) * t151) * t336 + (Ifges(5,1) * t152 + Ifges(5,4) * t151) * t336 + (Ifges(6,4) * t152 + Ifges(6,2) * t151) * t338 + (Ifges(5,4) * t152 + Ifges(5,2) * t151) * t338; t216 * t19 + t204 * t29 - t205 * t31 - t98 * t132 + (t188 - t23) * t89 + t264 + (-t219 * t30 + t222 * t32 + (-m(5) * t237 - t219 * t88 - t222 * t90) * qJD(4) + m(5) * t347) * pkin(9) + t347 * mrSges(5,3) + t348 * t344 + t349 * t343 + t350 * t334 + ((-Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t164 + t346 + t351) * t163 + (t187 - t27) * t87 - m(6) * (t15 * t23 + t18 * t27 + t41 * t70) + (-t1 * t219 + t2 * t222) * mrSges(6,3) + (-t346 + (m(6) * t41 + t73) * t316) * qJD(4) - t345 * t164 + t310 - t311 + m(6) * (t1 * t204 + t15 * t188 + t16 * t216 + t18 * t187 - t2 * t205) - t39 * t90 - t40 * t88 - t70 * t73 - pkin(3) * t20 + t354 * t219 / 0.2e1 + t355 * t318 + (-pkin(3) * t36 - t24 * t39 - t25 * t40 - t85 * t99) * m(5) + t288 * t99 - t16 * t252 - t36 * t254; -t85 * (mrSges(5,1) * t125 + mrSges(5,2) * t124) - t41 * (mrSges(6,1) * t125 + mrSges(6,2) * t124) + t362 - t232 + (-(-t15 + t17) * t18 + (-t125 * t41 + t1) * pkin(4)) * m(6) + t25 * t90 - t17 * t87 - t24 * t88 + t18 * t89 + (-t125 * t73 + t29) * pkin(4) + (t124 * t24 + t125 * t25) * mrSges(5,3) + (t124 * t15 + t125 * t18) * mrSges(6,3) + (t124 * t361 - t375) * t336 + t353 * t335 + (t124 * t359 - t125 * t357) * t327 + (-t125 * t358 + t352 + t376) * t338; -t124 * t87 + t125 * t89 + 0.2e1 * (t16 / 0.2e1 + t18 * t338 + t15 * t335) * m(6) + t19;];
tauc = t11(:);
