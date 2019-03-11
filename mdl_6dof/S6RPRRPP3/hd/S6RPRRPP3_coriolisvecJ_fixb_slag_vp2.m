% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2019-03-09 04:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPP3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP3_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP3_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:35:09
% EndTime: 2019-03-09 04:35:23
% DurationCPUTime: 7.23s
% Computational Cost: add. (3420->524), mult. (8193->627), div. (0->0), fcn. (4508->6), ass. (0->233)
t318 = -Ifges(5,4) + Ifges(7,6);
t316 = Ifges(5,1) + Ifges(7,3);
t312 = Ifges(7,4) + Ifges(6,5);
t311 = Ifges(5,5) + Ifges(7,5);
t315 = Ifges(7,2) + Ifges(6,3);
t162 = cos(qJ(3));
t235 = qJD(1) * t162;
t151 = Ifges(4,4) * t235;
t317 = -t151 / 0.2e1;
t159 = sin(qJ(4));
t314 = t318 * t159;
t160 = sin(qJ(3));
t236 = qJD(1) * t160;
t313 = t236 / 0.2e1;
t214 = Ifges(4,5) * qJD(3) / 0.2e1;
t232 = qJD(3) * t160;
t211 = qJD(1) * t232;
t161 = cos(qJ(4));
t230 = qJD(4) * t159;
t216 = t160 * t230;
t218 = t161 * t235;
t227 = qJD(3) * qJD(4);
t84 = qJD(1) * t216 - qJD(3) * t218 - t161 * t227;
t229 = qJD(4) * t161;
t231 = qJD(3) * t162;
t302 = t159 * t231 + t160 * t229;
t85 = qJD(1) * t302 + t159 * t227;
t310 = t311 * t211 - t316 * t84 + t318 * t85;
t309 = t315 * t85 + (Ifges(6,6) - Ifges(7,6)) * t84 + t312 * t211;
t173 = -qJ(5) * t161 + qJ(6) * t159;
t171 = t173 * t162;
t207 = pkin(4) * t230 - qJD(5) * t159;
t149 = sin(pkin(9)) * pkin(1) + pkin(7);
t130 = t149 * qJD(1);
t103 = t160 * qJD(2) + t162 * t130;
t220 = t159 * pkin(4) * t235 + t103;
t308 = -qJD(1) * t171 + qJD(4) * t173 - qJD(6) * t161 + t207 - t220;
t241 = t159 * t162;
t224 = pkin(5) * t241;
t284 = pkin(5) + pkin(8);
t102 = qJD(2) * t162 - t160 * t130;
t202 = pkin(3) * t160 - pkin(8) * t162;
t124 = t202 * qJD(1);
t54 = t161 * t102 + t159 * t124;
t307 = -t284 * t230 - (qJ(5) * t160 - t224) * qJD(1) - t54;
t136 = t284 * t161;
t239 = t161 * t162;
t223 = pkin(5) * t239;
t268 = pkin(4) + qJ(6);
t53 = -t159 * t102 + t161 * t124;
t306 = qJD(4) * t136 - (-t160 * t268 + t223) * qJD(1) + t53;
t219 = -cos(pkin(9)) * pkin(1) - pkin(2);
t113 = -pkin(3) * t162 - t160 * pkin(8) + t219;
t169 = t113 * qJD(1);
t98 = qJD(3) * t102;
t305 = qJD(4) * t169 + t98;
t147 = -qJD(4) + t235;
t233 = qJD(3) * t159;
t121 = t161 * t236 + t233;
t97 = qJD(3) * pkin(8) + t103;
t36 = t159 * t97 - t161 * t169;
t172 = pkin(5) * t121 + t36;
t301 = qJD(5) + t172;
t11 = t147 * t268 + t301;
t37 = t159 * t169 + t161 * t97;
t28 = t147 * qJ(5) - t37;
t228 = t161 * qJD(3);
t120 = t159 * t236 - t228;
t300 = -t120 * pkin(5) + qJD(6);
t13 = -t28 + t300;
t96 = -qJD(3) * pkin(3) - t102;
t170 = -t121 * qJ(5) + t96;
t16 = t120 * t268 + t170;
t174 = t159 * t37 - t161 * t36;
t296 = -qJD(5) - t36;
t27 = pkin(4) * t147 - t296;
t175 = t159 * t28 + t161 * t27;
t256 = Ifges(6,6) * t159;
t184 = -Ifges(6,2) * t161 + t256;
t257 = Ifges(5,4) * t161;
t189 = -Ifges(5,2) * t159 + t257;
t193 = -mrSges(7,2) * t161 + mrSges(7,3) * t159;
t194 = -mrSges(6,2) * t159 - mrSges(6,3) * t161;
t196 = mrSges(5,1) * t159 + mrSges(5,2) * t161;
t274 = t161 / 0.2e1;
t275 = -t161 / 0.2e1;
t276 = t159 / 0.2e1;
t277 = -t159 / 0.2e1;
t279 = -t147 / 0.2e1;
t280 = t121 / 0.2e1;
t281 = -t121 / 0.2e1;
t282 = t120 / 0.2e1;
t283 = -t120 / 0.2e1;
t290 = -Ifges(5,6) + t312;
t291 = Ifges(6,4) - t311;
t252 = Ifges(7,6) * t161;
t255 = Ifges(6,6) * t161;
t293 = t315 * t159 + t252 - t255;
t295 = t316 * t161 + t314;
t116 = Ifges(5,4) * t120;
t254 = Ifges(7,6) * t120;
t297 = t316 * t121 - t311 * t147 - t116 + t254;
t114 = Ifges(7,6) * t121;
t249 = t121 * Ifges(6,6);
t298 = t315 * t120 - t312 * t147 + t114 - t249;
t35 = t120 * pkin(4) + t170;
t115 = Ifges(6,6) * t120;
t48 = -t147 * Ifges(6,4) - t121 * Ifges(6,2) + t115;
t250 = t121 * Ifges(5,4);
t49 = -t120 * Ifges(5,2) - t147 * Ifges(5,6) + t250;
t163 = t175 * mrSges(6,1) + (t11 * t161 - t13 * t159) * mrSges(7,1) - t174 * mrSges(5,3) + t16 * t193 + t184 * t281 + t189 * t283 + t35 * t194 + t96 * t196 + t48 * t275 + t49 * t277 + t293 * t282 + t295 * t280 + t298 * t276 + t297 * t274 + (t159 * t290 - t161 * t291) * t279;
t304 = t102 * mrSges(4,3) - Ifges(4,1) * t313 - t163 - t214 + t317;
t303 = t255 + t257 + (Ifges(5,1) + Ifges(6,2)) * t159;
t99 = qJD(3) * t103;
t299 = t256 - t314 + (Ifges(5,2) + t315) * t161;
t213 = -Ifges(4,6) * qJD(3) / 0.2e1;
t294 = -t37 - t300;
t127 = t202 * qJD(3);
t112 = qJD(1) * t127;
t7 = t159 * t112 + t161 * t305 - t230 * t97;
t8 = t112 * t161 - t159 * t305 - t97 * t229;
t199 = -t159 * t8 + t161 * t7;
t4 = -qJ(5) * t211 + qJD(5) * t147 - t7;
t5 = -pkin(4) * t211 - t8;
t200 = t159 * t5 - t161 * t4;
t132 = t219 * qJD(1);
t203 = -Ifges(5,6) / 0.2e1 + Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1;
t204 = Ifges(6,4) / 0.2e1 - Ifges(7,5) / 0.2e1 - Ifges(5,5) / 0.2e1;
t205 = Ifges(6,1) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(7,1) / 0.2e1;
t259 = Ifges(4,4) * t160;
t289 = -t203 * t120 + t204 * t121 + t205 * t147 - t13 * mrSges(7,2) - t132 * mrSges(4,1) - t27 * mrSges(6,2) - t213 + (Ifges(4,2) * t162 + t259) * qJD(1) / 0.2e1 - Ifges(5,6) * t283 - Ifges(6,4) * t281 + t11 * mrSges(7,3) + t28 * mrSges(6,3) + t36 * mrSges(5,1) + t37 * mrSges(5,2) - t312 * t282 - t311 * t280 - (Ifges(5,3) + Ifges(7,1) + Ifges(6,1)) * t279;
t288 = -t84 / 0.2e1;
t287 = t84 / 0.2e1;
t286 = -t85 / 0.2e1;
t278 = t147 / 0.2e1;
t59 = mrSges(6,1) * t85 - mrSges(6,3) * t211;
t63 = -mrSges(5,2) * t211 - mrSges(5,3) * t85;
t267 = -t59 + t63;
t61 = -t84 * mrSges(6,1) + mrSges(6,2) * t211;
t62 = mrSges(5,1) * t211 + mrSges(5,3) * t84;
t266 = t61 - t62;
t67 = -mrSges(7,2) * t121 + mrSges(7,3) * t120;
t70 = -mrSges(6,2) * t120 - mrSges(6,3) * t121;
t265 = t67 + t70;
t261 = mrSges(5,3) * t120;
t88 = mrSges(5,2) * t147 - t261;
t91 = mrSges(6,1) * t120 + mrSges(6,3) * t147;
t264 = t88 - t91;
t260 = mrSges(5,3) * t121;
t89 = -mrSges(5,1) * t147 - t260;
t93 = mrSges(6,1) * t121 - mrSges(6,2) * t147;
t263 = -t89 + t93;
t92 = -mrSges(7,1) * t120 - mrSges(7,2) * t147;
t262 = -t91 + t92;
t248 = t132 * mrSges(4,2);
t222 = mrSges(4,3) * t236;
t247 = -qJD(3) * mrSges(4,1) + mrSges(5,1) * t120 + mrSges(5,2) * t121 + t222;
t244 = qJ(5) * t120;
t243 = t149 * t159;
t242 = t159 * t160;
t240 = t160 * t161;
t238 = t113 * t229 + t159 * t127;
t123 = t149 * t239;
t66 = t159 * t113 + t123;
t237 = pkin(4) * t242 + t160 * t149;
t226 = t88 + t262;
t90 = mrSges(7,1) * t121 + mrSges(7,3) * t147;
t225 = t90 + t263;
t221 = mrSges(4,3) * t235;
t122 = t149 * t241;
t212 = m(4) * t149 + mrSges(4,3);
t60 = -t85 * mrSges(7,1) + mrSges(7,2) * t211;
t210 = -qJ(5) * t159 - pkin(3);
t209 = -pkin(4) - t243;
t208 = t149 * t161 - qJ(5);
t65 = t161 * t113 - t122;
t206 = pkin(4) * t302 + qJ(5) * t216 + t149 * t231;
t1 = -pkin(5) * t84 + qJD(6) * t147 - t211 * t268 - t8;
t2 = -pkin(5) * t85 - t4;
t201 = t1 * t159 + t161 * t2;
t56 = qJ(5) * t162 - t66;
t198 = qJD(4) * t123 + t113 * t230 - t161 * t127;
t197 = mrSges(5,1) * t161 - mrSges(5,2) * t159;
t195 = mrSges(6,2) * t161 - mrSges(6,3) * t159;
t192 = mrSges(7,2) * t159 + mrSges(7,3) * t161;
t187 = Ifges(6,4) * t159 + Ifges(6,5) * t161;
t186 = Ifges(7,4) * t161 - Ifges(7,5) * t159;
t185 = Ifges(5,5) * t159 + Ifges(5,6) * t161;
t178 = -Ifges(7,3) * t159 + t252;
t58 = -t84 * mrSges(7,1) - mrSges(7,3) * t211;
t167 = t84 * qJ(5) - t121 * qJD(5) + t99;
t166 = -t8 * mrSges(5,1) + t7 * mrSges(5,2) - t5 * mrSges(6,2) - t2 * mrSges(7,2) + t4 * mrSges(6,3) + t1 * mrSges(7,3);
t155 = t162 * pkin(4);
t146 = Ifges(6,1) * t211;
t145 = Ifges(7,1) * t211;
t144 = Ifges(5,3) * t211;
t135 = t284 * t159;
t133 = -qJD(3) * mrSges(4,2) + t221;
t129 = -pkin(4) * t161 + t210;
t111 = -t161 * t268 + t210;
t106 = -qJ(5) * t229 + t207;
t86 = -qJ(5) * t240 + t237;
t79 = Ifges(6,4) * t84;
t78 = Ifges(7,4) * t85;
t77 = Ifges(5,5) * t84;
t76 = Ifges(6,5) * t85;
t75 = Ifges(7,5) * t84;
t74 = Ifges(5,6) * t85;
t68 = pkin(4) * t121 + t244;
t64 = t160 * t173 + t237;
t57 = t155 - t65;
t55 = -qJ(5) * t218 + t220;
t43 = -pkin(5) * t242 - t56;
t42 = -pkin(4) * t236 - t53;
t41 = -qJ(5) * t236 - t54;
t39 = t121 * t268 + t244;
t38 = qJ(6) * t162 + t122 + t155 + (pkin(5) * t160 - t113) * t161;
t33 = (-qJ(5) * t231 - qJD(5) * t160) * t161 + t206;
t32 = mrSges(7,2) * t84 + mrSges(7,3) * t85;
t31 = mrSges(5,1) * t85 - mrSges(5,2) * t84;
t30 = -mrSges(6,2) * t85 + mrSges(6,3) * t84;
t25 = -t84 * Ifges(5,4) - t85 * Ifges(5,2) + Ifges(5,6) * t211;
t24 = Ifges(6,4) * t211 + t84 * Ifges(6,2) + t85 * Ifges(6,6);
t20 = t232 * t243 - t198;
t19 = (-t160 * t228 - t162 * t230) * t149 + t238;
t15 = t209 * t232 + t198;
t14 = (t149 * t230 + qJD(5)) * t162 + t208 * t232 - t238;
t12 = qJD(3) * t171 + (qJD(6) * t159 + (qJ(6) * qJD(4) - qJD(5)) * t161) * t160 + t206;
t10 = t85 * pkin(4) + t167;
t9 = -qJD(5) * t162 + (-pkin(5) * t240 - t122) * qJD(4) + (-t160 * t208 - t224) * qJD(3) + t238;
t6 = -pkin(5) * t216 + qJD(6) * t162 + (t223 + (-qJ(6) + t209) * t160) * qJD(3) + t198;
t3 = t120 * qJD(6) + t268 * t85 + t167;
t17 = [t86 * t30 + t19 * t88 + t20 * t89 + t6 * t90 + t14 * t91 + t9 * t92 + t15 * t93 + t65 * t62 + t66 * t63 + t12 * t67 + t33 * t70 + t38 * t58 + t56 * t59 + t43 * t60 + t57 * t61 + t64 * t32 + m(6) * (t10 * t86 + t14 * t28 + t15 * t27 + t33 * t35 + t4 * t56 + t5 * t57) + m(7) * (t1 * t38 + t11 * t6 + t12 * t16 + t13 * t9 + t2 * t43 + t3 * t64) + m(5) * (t37 * t19 - t36 * t20 + t8 * t65 + t7 * t66) + (-t144 / 0.2e1 - t145 / 0.2e1 - t146 / 0.2e1 - t78 / 0.2e1 - t79 / 0.2e1 + t75 / 0.2e1 - t76 / 0.2e1 + t77 / 0.2e1 + t74 / 0.2e1 - t203 * t85 + t166 - t204 * t84 + t212 * t98 + (0.3e1 / 0.2e1 * t151 + 0.2e1 * t248 + t214 + (-m(4) * t102 + m(5) * t96 + t247) * t149 - t304) * qJD(3)) * t162 + (t189 * t286 + t184 * t287 + t3 * t193 + t10 * t194 + t24 * t275 + (-t159 * t7 - t161 * t8) * mrSges(5,3) + (t1 * t161 - t159 * t2) * mrSges(7,1) + (t159 * t4 + t161 * t5) * mrSges(6,1) + (mrSges(4,3) + t196) * t99 + ((-0.3e1 / 0.2e1 * t259 + t219 * mrSges(4,1) - t204 * t240 + t203 * t242 + (0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2) - t205) * t162) * qJD(1) + t213 - t212 * t103 - t289) * qJD(3) + t295 * t288 + t293 * t85 / 0.2e1 + t25 * t277 + t309 * t276 + t310 * t274 + (t31 + (m(5) + m(4)) * t99 - qJD(3) * t133) * t149 + (t185 * t278 + t178 * t280 + t96 * t197 - t35 * t195 + t16 * t192 + t49 * t275 + (-t159 * t36 - t161 * t37) * mrSges(5,3) + (-t11 * t159 - t13 * t161) * mrSges(7,1) + (-t159 * t27 + t161 * t28) * mrSges(6,1) + t303 * t281 + (t187 + t186) * t279 + t299 * t282 + t297 * t277 + t48 * t276 + t298 * t274) * qJD(4)) * t160; (-t30 - t31 - t32 + (t159 * t225 + t161 * t226 + t133 - t221) * qJD(3) + m(5) * (t228 * t37 + t233 * t36 - t99) + m(6) * (-t228 * t28 + t233 * t27 - t10) + m(7) * (t11 * t233 + t13 * t228 - t3)) * t162 + ((t60 + t267) * t161 + (t58 + t266) * t159 + (-t222 + t247 + t265) * qJD(3) + (-t159 * t226 + t161 * t225) * qJD(4) + m(5) * (qJD(3) * t96 + t229 * t36 - t230 * t37 + t199) + m(6) * (qJD(3) * t35 + t229 * t27 + t230 * t28 + t200) + m(7) * (qJD(3) * t16 + t11 * t229 - t13 * t230 + t201)) * t160; t135 * t58 + t136 * t60 + t129 * t30 - t102 * t133 + t111 * t32 - t98 * mrSges(4,2) - t54 * t88 - t53 * t89 - t41 * t91 - t42 * t93 + (t1 * t135 + t11 * t306 + t111 * t3 + t13 * t307 + t136 * t2 + t16 * t308) * m(7) + t308 * t67 + t309 * t275 + t310 * t276 - pkin(3) * t31 + t306 * t90 + t307 * t92 + t303 * t288 + t299 * t286 + t199 * mrSges(5,3) + t200 * mrSges(6,1) + t201 * mrSges(7,1) + (-mrSges(4,1) - t197) * t99 - t3 * t192 + t10 * t195 + t178 * t287 + t24 * t277 + t25 * t274 - m(5) * (t103 * t96 - t36 * t53 + t37 * t54) - m(6) * (t27 * t42 + t28 * t41 + t35 * t55) + (t106 - t55) * t70 + ((t103 * mrSges(4,3) + Ifges(4,4) * t313 + t213 + t289) * t160 + (t317 + (Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1) * t236 - t248 + t214 + t304) * t162 + (t185 / 0.2e1 - t186 / 0.2e1 - t187 / 0.2e1) * t232) * qJD(1) - t247 * t103 + (t163 + (-m(5) * t174 + m(6) * t175 - t159 * t264 + t161 * t263) * pkin(8)) * qJD(4) + (t159 * t266 + t161 * t267) * pkin(8) + m(6) * (t200 * pkin(8) + t10 * t129 + t106 * t35) + m(5) * (-pkin(3) * t99 + t199 * pkin(8)); t144 + t145 + t146 - t96 * (mrSges(5,1) * t121 - mrSges(5,2) * t120) - t16 * (mrSges(7,2) * t120 + mrSges(7,3) * t121) - t35 * (-mrSges(6,2) * t121 + mrSges(6,3) * t120) + t78 + t79 - t75 + t76 - t77 - t74 - t166 - t39 * t67 - t68 * t70 - pkin(4) * t61 + (qJ(5) * t2 - t1 * t268 + t294 * t11 + t13 * t301 - t16 * t39) * m(7) + (-t59 + t60) * qJ(5) + (Ifges(6,2) * t120 + t249 + t49) * t280 + t172 * t92 - t268 * t58 + (t11 * t120 + t121 * t13) * mrSges(7,1) + (t120 * t27 - t121 * t28) * mrSges(6,1) + (t315 * t121 + t115 - t254 + t48) * t283 + (-t316 * t120 + t114 - t250 + t298) * t281 + t262 * qJD(5) + (t260 - t263) * t37 + (t261 + t264) * t36 + (t120 * t291 + t121 * t290) * t278 + t294 * t90 + (-pkin(4) * t5 - qJ(5) * t4 - t27 * t37 + t28 * t296 - t35 * t68) * m(6) + (-Ifges(5,2) * t121 - t116 + t297) * t282; t265 * t121 + t262 * t147 + t58 + t61 + (t121 * t16 + t13 * t147 + t1) * m(7) + (t121 * t35 - t147 * t28 + t5) * m(6); -t120 * t67 - t147 * t90 + 0.2e1 * (t2 / 0.2e1 + t11 * t279 + t16 * t283) * m(7) + t60;];
tauc  = t17(:);
