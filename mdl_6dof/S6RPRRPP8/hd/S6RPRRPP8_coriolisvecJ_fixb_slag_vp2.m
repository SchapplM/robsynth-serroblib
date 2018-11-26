% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
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
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 16:15
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPP8_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP8_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP8_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:15:09
% EndTime: 2018-11-23 16:15:15
% DurationCPUTime: 6.58s
% Computational Cost: add. (3364->514), mult. (7422->620), div. (0->0), fcn. (3861->4), ass. (0->228)
t325 = Ifges(6,6) - Ifges(7,6);
t323 = Ifges(5,1) + Ifges(7,3);
t316 = Ifges(7,4) + Ifges(6,5);
t315 = Ifges(5,5) + Ifges(7,5);
t322 = Ifges(7,2) + Ifges(6,3);
t154 = -pkin(1) - pkin(7);
t137 = qJD(1) * t154 + qJD(2);
t153 = cos(qJ(3));
t102 = -qJD(3) * pkin(3) - t137 * t153;
t150 = sin(qJ(4));
t152 = cos(qJ(4));
t151 = sin(qJ(3));
t120 = t151 * t137;
t101 = qJD(3) * pkin(8) + t120;
t121 = pkin(3) * t151 - pkin(8) * t153 + qJ(2);
t163 = t121 * qJD(1);
t47 = t101 * t150 - t152 * t163;
t48 = t152 * t101 + t150 * t163;
t172 = t150 * t48 - t152 * t47;
t235 = qJD(1) * t151;
t144 = qJD(4) + t235;
t298 = -qJD(5) - t47;
t32 = -pkin(4) * t144 - t298;
t33 = -t144 * qJ(5) - t48;
t174 = t150 * t33 + t152 * t32;
t268 = pkin(4) + qJ(6);
t234 = qJD(1) * t153;
t114 = qJD(3) * t150 + t152 * t234;
t169 = pkin(5) * t114 + t47;
t303 = qJD(5) + t169;
t12 = -t144 * t268 + t303;
t226 = t152 * qJD(3);
t113 = t150 * t234 - t226;
t302 = -t113 * pkin(5) + qJD(6);
t14 = -t33 + t302;
t176 = t12 * t152 - t14 * t150;
t251 = Ifges(7,6) * t152;
t181 = Ifges(7,2) * t150 + t251;
t196 = mrSges(7,2) * t152 - mrSges(7,3) * t150;
t197 = mrSges(6,2) * t150 + mrSges(6,3) * t152;
t199 = mrSges(5,1) * t150 + mrSges(5,2) * t152;
t162 = -qJ(5) * t114 + t102;
t21 = t113 * t268 + t162;
t276 = t152 / 0.2e1;
t277 = -t152 / 0.2e1;
t278 = t150 / 0.2e1;
t279 = -t150 / 0.2e1;
t280 = t144 / 0.2e1;
t281 = -t144 / 0.2e1;
t282 = t114 / 0.2e1;
t284 = t113 / 0.2e1;
t285 = -t113 / 0.2e1;
t318 = (-Ifges(5,4) - t325) * t150;
t294 = (Ifges(6,2) + t323) * t152 + t318;
t319 = (Ifges(5,4) + Ifges(6,6)) * t152;
t297 = (-Ifges(5,2) - Ifges(6,3)) * t150 + t319;
t106 = Ifges(7,6) * t114;
t246 = t114 * Ifges(6,6);
t299 = t322 * t113 + t316 * t144 + t106 - t246;
t108 = Ifges(5,4) * t113;
t253 = Ifges(7,6) * t113;
t300 = t323 * t114 + t315 * t144 - t108 + t253;
t36 = pkin(4) * t113 + t162;
t107 = Ifges(6,6) * t113;
t42 = t144 * Ifges(6,4) - t114 * Ifges(6,2) + t107;
t247 = t114 * Ifges(5,4);
t43 = -t113 * Ifges(5,2) + t144 * Ifges(5,6) + t247;
t324 = t300 * t276 + t277 * t42 + t299 * t278 + t279 * t43 + t102 * t199 + t181 * t284 + (Ifges(6,4) * t152 - Ifges(6,5) * t150) * t281 + (t315 * t152 + (Ifges(7,4) - Ifges(5,6)) * t150) * t280 + t297 * t285 + t294 * t282 - t21 * t196 - t36 * t197 + t174 * mrSges(6,1) + t176 * mrSges(7,1) - t172 * mrSges(5,3);
t244 = Ifges(4,5) * qJD(3);
t259 = Ifges(4,4) * t151;
t321 = t244 / 0.2e1 + (t153 * Ifges(4,1) - t259) * qJD(1) / 0.2e1 + t324;
t230 = qJD(4) * t150;
t241 = qJ(5) * t152;
t320 = pkin(4) * t230 - qJD(5) * t150 - t235 * t241 - t120;
t317 = -qJD(3) / 0.2e1;
t232 = qJD(3) * t153;
t212 = qJD(1) * t232;
t228 = qJD(4) * t153;
t217 = t150 * t228;
t305 = t151 * t226 + t217;
t77 = t305 * qJD(1) - qJD(4) * t226;
t233 = qJD(3) * t151;
t220 = t150 * t233;
t78 = -qJD(1) * t220 + qJD(4) * t114;
t314 = t316 * t212 + t322 * t78 + t325 * t77;
t313 = (-Ifges(5,4) + Ifges(7,6)) * t78 - t323 * t77 + t315 * t212;
t170 = qJ(6) * t150 - t241;
t213 = t268 * t151;
t312 = qJD(1) * t150 * t213 + qJD(4) * t170 - qJD(6) * t152 + t320;
t229 = qJD(4) * t152;
t311 = pkin(4) * t150 * t235 - qJ(5) * t229 + t320;
t274 = pkin(5) * t151;
t286 = pkin(5) + pkin(8);
t204 = pkin(3) * t153 + pkin(8) * t151;
t117 = t204 * qJD(1);
t237 = t152 * t153;
t60 = t150 * t117 + t137 * t237;
t310 = -t286 * t230 - (qJ(5) * t153 + t150 * t274) * qJD(1) - t60;
t129 = t286 * t152;
t164 = -t152 * t274 - t153 * t268;
t239 = t150 * t153;
t59 = t117 * t152 - t137 * t239;
t309 = -qJD(1) * t164 + qJD(4) * t129 + t59;
t221 = t137 * t232;
t308 = qJD(4) * t163 + t221;
t307 = t319 + (Ifges(5,1) + Ifges(6,2)) * t150;
t306 = t316 * t152 + (Ifges(6,4) - Ifges(7,5)) * t150;
t304 = qJD(1) * (qJ(2) * (m(3) + m(4)) + mrSges(3,3));
t301 = -t318 + (Ifges(5,2) + t322) * t152;
t283 = -t114 / 0.2e1;
t214 = Ifges(4,6) * t317;
t112 = qJD(3) * t204 + qJD(2);
t93 = t112 * qJD(1);
t10 = -t101 * t229 - t308 * t150 + t152 * t93;
t9 = -t101 * t230 + t150 * t93 + t308 * t152;
t201 = -t10 * t150 + t152 * t9;
t4 = -qJ(5) * t212 - qJD(5) * t144 - t9;
t6 = -pkin(4) * t212 - t10;
t202 = t150 * t6 - t152 * t4;
t296 = -t48 - t302;
t292 = -m(5) * t172 + m(6) * t174;
t205 = Ifges(7,5) / 0.2e1 + Ifges(5,5) / 0.2e1 - Ifges(6,4) / 0.2e1;
t206 = Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1 - Ifges(5,6) / 0.2e1;
t207 = Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1 + Ifges(5,3) / 0.2e1;
t258 = Ifges(4,4) * t153;
t291 = t206 * t113 + t205 * t114 + t207 * t144 + t14 * mrSges(7,2) + t32 * mrSges(6,2) + Ifges(5,6) * t285 + Ifges(6,4) * t283 + t214 - (-t151 * Ifges(4,2) + t258) * qJD(1) / 0.2e1 - t12 * mrSges(7,3) - t33 * mrSges(6,3) - t47 * mrSges(5,1) - t48 * mrSges(5,2) + t316 * t284 + t315 * t282 + (Ifges(5,3) + Ifges(7,1) + Ifges(6,1)) * t280;
t289 = -t77 / 0.2e1;
t288 = -t78 / 0.2e1;
t275 = pkin(4) * t152;
t54 = mrSges(6,1) * t78 - mrSges(6,3) * t212;
t58 = -mrSges(5,2) * t212 - mrSges(5,3) * t78;
t267 = -t54 + t58;
t56 = -t77 * mrSges(6,1) + mrSges(6,2) * t212;
t57 = mrSges(5,1) * t212 + mrSges(5,3) * t77;
t266 = t56 - t57;
t61 = -mrSges(7,2) * t114 + mrSges(7,3) * t113;
t64 = -mrSges(6,2) * t113 - mrSges(6,3) * t114;
t265 = t61 + t64;
t261 = mrSges(5,3) * t113;
t79 = -mrSges(5,2) * t144 - t261;
t82 = mrSges(6,1) * t113 - mrSges(6,3) * t144;
t264 = t79 - t82;
t260 = mrSges(5,3) * t114;
t80 = mrSges(5,1) * t144 - t260;
t84 = mrSges(6,1) * t114 + mrSges(6,2) * t144;
t263 = -t80 + t84;
t83 = -mrSges(7,1) * t113 + mrSges(7,2) * t144;
t262 = t82 - t83;
t250 = qJ(2) * mrSges(4,1);
t249 = qJ(2) * mrSges(4,2);
t245 = -qJD(3) * mrSges(4,1) + mrSges(5,1) * t113 + mrSges(5,2) * t114 + mrSges(4,3) * t234;
t242 = qJ(5) * t113;
t240 = qJD(3) * mrSges(4,2);
t238 = t151 * t154;
t86 = t150 * t121 + t152 * t238;
t231 = qJD(3) * t154;
t227 = qJD(4) * t154;
t225 = qJD(1) * qJD(2);
t224 = t79 - t262;
t81 = mrSges(7,1) * t114 - mrSges(7,3) * t144;
t223 = t81 + t263;
t218 = t153 * t231;
t222 = t150 * t112 + t121 * t229 + t152 * t218;
t115 = t137 * t233;
t216 = t151 * t227;
t215 = -t244 / 0.2e1;
t55 = -t78 * mrSges(7,1) + mrSges(7,2) * t212;
t211 = -qJ(5) * t150 - pkin(3);
t210 = qJD(3) * t268;
t134 = t150 * t238;
t85 = t121 * t152 - t134;
t208 = t305 * qJ(5) + t151 * t231 + t228 * t275;
t1 = -pkin(5) * t77 - qJD(6) * t144 - t210 * t234 - t10;
t3 = -pkin(5) * t78 - t4;
t203 = t1 * t150 + t152 * t3;
t68 = -qJ(5) * t151 - t86;
t200 = mrSges(5,1) * t152 - mrSges(5,2) * t150;
t198 = mrSges(6,2) * t152 - mrSges(6,3) * t150;
t195 = mrSges(7,2) * t150 + mrSges(7,3) * t152;
t185 = Ifges(5,5) * t150 + Ifges(5,6) * t152;
t178 = -Ifges(7,3) * t150 + t251;
t175 = t12 * t150 + t14 * t152;
t173 = t150 * t32 - t152 * t33;
t171 = t150 * t47 + t152 * t48;
t53 = -t77 * mrSges(7,1) - mrSges(7,3) * t212;
t27 = -t121 * t230 - t150 * t218 + (t112 - t216) * t152;
t168 = qJ(5) * t77 - qJD(5) * t114 + t115;
t160 = -t150 * t224 + t152 * t223;
t159 = t10 * mrSges(5,1) - t9 * mrSges(5,2) + t6 * mrSges(6,2) + t3 * mrSges(7,2) - t4 * mrSges(6,3) - t1 * mrSges(7,3);
t145 = pkin(4) * t239;
t142 = Ifges(6,1) * t212;
t141 = Ifges(7,1) * t212;
t140 = Ifges(5,3) * t212;
t128 = t286 * t150;
t126 = -mrSges(4,3) * t235 - t240;
t122 = t211 - t275;
t116 = (mrSges(4,1) * t151 + mrSges(4,2) * t153) * qJD(1);
t105 = -t152 * t268 + t211;
t88 = t145 + (-t154 - t241) * t153;
t76 = Ifges(6,4) * t77;
t75 = Ifges(7,4) * t78;
t74 = Ifges(5,5) * t77;
t73 = Ifges(6,5) * t78;
t72 = Ifges(7,5) * t77;
t71 = Ifges(5,6) * t78;
t69 = -pkin(4) * t151 - t85;
t67 = t145 + (-t154 + t170) * t153;
t62 = pkin(4) * t114 + t242;
t52 = -pkin(5) * t239 - t68;
t51 = -pkin(4) * t234 - t59;
t50 = -qJ(5) * t234 - t60;
t37 = t134 + (pkin(5) * t153 - t121) * t152 - t213;
t34 = t114 * t268 + t242;
t31 = -pkin(4) * t220 - qJD(5) * t237 + t208;
t26 = -t150 * t216 + t222;
t25 = mrSges(7,2) * t77 + mrSges(7,3) * t78;
t24 = mrSges(5,1) * t78 - mrSges(5,2) * t77;
t23 = -mrSges(6,2) * t78 + mrSges(6,3) * t77;
t22 = -pkin(4) * t232 - t27;
t19 = -t77 * Ifges(5,4) - t78 * Ifges(5,2) + Ifges(5,6) * t212;
t18 = Ifges(6,4) * t212 + t77 * Ifges(6,2) + t78 * Ifges(6,6);
t13 = -qJ(5) * t232 + (t150 * t227 - qJD(5)) * t151 - t222;
t11 = (qJ(6) * qJD(4) - qJD(5)) * t237 + (qJD(6) * t153 - t151 * t210) * t150 + t208;
t8 = (-pkin(5) * t229 + qJ(5) * qJD(3)) * t153 + (qJD(5) + (pkin(5) * qJD(3) - t227) * t150) * t151 + t222;
t7 = pkin(4) * t78 + t168;
t5 = -pkin(5) * t217 + qJD(3) * t164 - qJD(6) * t151 - t27;
t2 = qJD(6) * t113 + t268 * t78 + t168;
t15 = [t22 * t84 + t85 * t57 + t86 * t58 + t88 * t23 + t26 * t79 + t27 * t80 + t5 * t81 + t13 * t82 + t8 * t83 + t68 * t54 + t69 * t56 + t11 * t61 + t31 * t64 + t67 * t25 + t37 * t53 + t52 * t55 + m(6) * (t13 * t33 + t22 * t32 + t31 * t36 + t4 * t68 + t6 * t69 + t7 * t88) + m(5) * (t10 * t85 + t26 * t48 - t27 * t47 + t86 * t9) + m(7) * (t1 * t37 + t11 * t21 + t12 * t5 + t14 * t8 + t2 * t67 + t3 * t52) + (t116 + 0.2e1 * t304) * qJD(2) + (t140 / 0.2e1 + t141 / 0.2e1 + t142 / 0.2e1 + t76 / 0.2e1 - t74 / 0.2e1 + t75 / 0.2e1 - t71 / 0.2e1 - t72 / 0.2e1 + t73 / 0.2e1 - t205 * t77 + mrSges(4,1) * t225 + t159 + ((m(5) * t102 + t245) * t154 + (0.3e1 / 0.2e1 * t259 - 0.2e1 * t249) * qJD(1) + t215 - t321) * qJD(3) + t206 * t78) * t151 + (mrSges(4,2) * t225 - t154 * t24 + t18 * t277 - t7 * t197 - t2 * t196 + t78 * t181 / 0.2e1 + (-t10 * t152 - t150 * t9) * mrSges(5,3) + (t1 * t152 - t150 * t3) * mrSges(7,1) + (t150 * t4 + t152 * t6) * mrSges(6,1) + ((0.2e1 * t250 - 0.3e1 / 0.2e1 * t258 + t205 * t237 + t206 * t239 + (0.3e1 / 0.2e1 * Ifges(4,2) - 0.3e1 / 0.2e1 * Ifges(4,1) + t207) * t151) * qJD(1) + t154 * t126 + (-m(5) * t154 + t199) * t120 + t214 + t291) * qJD(3) + t297 * t288 + t19 * t279 + t314 * t278 + t294 * t289 + t313 * t276 + (-t173 * mrSges(6,1) - t175 * mrSges(7,1) - t171 * mrSges(5,3) + t102 * t200 + t178 * t282 + t185 * t281 + t21 * t195 - t36 * t198 + t299 * t276 + t43 * t277 + t42 * t278 + t300 * t279 + t306 * t280 + t307 * t283 + t301 * t284) * qJD(4)) * t153; (-t23 - t24 - t25 - m(6) * t7 - m(7) * t2 + (m(5) * t171 + m(6) * t173 + m(7) * t175 + t150 * t223 + t152 * t224 + t126) * qJD(3)) * t153 + ((t55 + t267) * t152 + (t53 + t266) * t150 + (t245 + t265) * qJD(3) + t160 * qJD(4) + m(5) * (qJD(3) * t102 + t229 * t47 - t230 * t48 + t201 - t221) + m(6) * (qJD(3) * t36 + t229 * t32 + t230 * t33 + t202) + m(7) * (qJD(3) * t21 + t12 * t229 - t14 * t230 + t203)) * t151 + (m(7) * t176 - t116 + t160 + t292 - t304) * qJD(1); (-pkin(3) * t115 + pkin(8) * t201 - t102 * t120 + t47 * t59 - t48 * t60) * m(5) + (t202 * pkin(8) + t122 * t7 + t311 * t36 - t32 * t51 - t33 * t50) * m(6) + t311 * t64 + (t1 * t128 + t105 * t2 + t309 * t12 + t129 * t3 + t310 * t14 + t312 * t21) * m(7) + t312 * t61 + t313 * t278 + t314 * t277 + t309 * t81 + t310 * t83 + t307 * t289 + t301 * t288 + ((-t150 * t264 + t152 * t263 + t292) * pkin(8) + t324) * qJD(4) + ((-t126 - t240) * t153 + ((-mrSges(4,1) - t200) * qJD(3) - t245) * t151) * t137 + t201 * mrSges(5,3) + t202 * mrSges(6,1) + t203 * mrSges(7,1) - t2 * t195 + t7 * t198 + t77 * t178 / 0.2e1 + t122 * t23 + t128 * t53 + t129 * t55 + t105 * t25 - t51 * t84 - t60 * t79 - t59 * t80 - t50 * t82 - pkin(3) * t24 + (t150 * t266 + t152 * t267) * pkin(8) + ((qJD(3) * t185 / 0.2e1 + t214 + (t258 / 0.2e1 - t250) * qJD(1) + t306 * t317 - t291) * t153 + ((t249 - t259 / 0.2e1 + (Ifges(4,1) / 0.2e1 - Ifges(4,2) / 0.2e1) * t153) * qJD(1) + t215 + t321) * t151) * qJD(1) + t19 * t276 + t18 * t279; (Ifges(6,2) * t113 + t246 + t43) * t282 + ((-Ifges(5,6) + t316) * t114 + (Ifges(6,4) - t315) * t113) * t281 + (qJ(5) * t3 - t1 * t268 + t296 * t12 + t303 * t14 - t21 * t34) * m(7) + t169 * t83 - t268 * t53 + t140 + t141 + t142 + t76 - t74 + t75 - t71 - t72 + t73 - t36 * (-mrSges(6,2) * t114 + mrSges(6,3) * t113) - t21 * (mrSges(7,2) * t113 + mrSges(7,3) * t114) - t102 * (mrSges(5,1) * t114 - mrSges(5,2) * t113) - pkin(4) * t56 - t34 * t61 - t62 * t64 + (-t54 + t55) * qJ(5) + (-t323 * t113 + t106 - t247 + t299) * t283 + (t322 * t114 + t107 - t253 + t42) * t285 + (t261 + t264) * t47 - t262 * qJD(5) + (t260 - t263) * t48 + t159 + (t113 * t12 + t114 * t14) * mrSges(7,1) + (t113 * t32 - t114 * t33) * mrSges(6,1) + (-Ifges(5,2) * t114 - t108 + t300) * t284 + t296 * t81 + (-pkin(4) * t6 - qJ(5) * t4 + t298 * t33 - t32 * t48 - t36 * t62) * m(6); t265 * t114 + t262 * t144 + t53 + t56 + (t114 * t21 - t14 * t144 + t1) * m(7) + (t114 * t36 + t144 * t33 + t6) * m(6); -t113 * t61 + t144 * t81 + 0.2e1 * (t3 / 0.2e1 + t21 * t285 + t12 * t280) * m(7) + t55;];
tauc  = t15(:);
