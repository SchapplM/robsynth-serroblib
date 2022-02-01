% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% m [6x1]
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
% Datum: 2022-01-20 12:09
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRR6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR6_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:07:40
% EndTime: 2022-01-20 12:07:51
% DurationCPUTime: 3.60s
% Computational Cost: add. (7050->365), mult. (11986->520), div. (0->0), fcn. (7919->8), ass. (0->198)
t221 = sin(qJ(4));
t222 = sin(qJ(3));
t225 = cos(qJ(4));
t226 = cos(qJ(3));
t189 = -t221 * t222 + t225 * t226;
t296 = -pkin(8) - pkin(7);
t252 = qJD(3) * t296;
t194 = t222 * t252;
t195 = t226 * t252;
t203 = t296 * t222;
t215 = t226 * pkin(8);
t204 = pkin(7) * t226 + t215;
t227 = cos(qJ(2));
t279 = pkin(1) * qJD(1);
t255 = t227 * t279;
t261 = qJD(4) * t225;
t262 = qJD(4) * t221;
t298 = -t189 * t255 + t225 * t194 + t221 * t195 + t203 * t261 - t204 * t262;
t147 = t221 * t203 + t225 * t204;
t190 = t221 * t226 + t222 * t225;
t297 = -qJD(4) * t147 + t190 * t255 - t194 * t221 + t225 * t195;
t216 = qJD(3) + qJD(4);
t129 = t216 * t189;
t288 = pkin(9) * t129;
t314 = -t288 + t297;
t130 = t216 * t190;
t126 = t130 * pkin(9);
t313 = t126 - t298;
t217 = qJD(1) + qJD(2);
t166 = t189 * t217;
t167 = t190 * t217;
t220 = sin(qJ(5));
t224 = cos(qJ(5));
t103 = t166 * t220 + t167 * t224;
t211 = -t226 * pkin(3) - pkin(2);
t169 = t211 * t217 - t255;
t111 = -t166 * pkin(4) + t169;
t214 = qJD(5) + t216;
t245 = t224 * t166 - t167 * t220;
t280 = Ifges(6,4) * t103;
t107 = t130 * t217;
t278 = pkin(1) * qJD(2);
t253 = qJD(1) * t278;
t244 = t227 * t253;
t202 = t226 * t244;
t223 = sin(qJ(2));
t256 = t223 * t279;
t199 = pkin(7) * t217 + t256;
t248 = pkin(8) * t217 + t199;
t238 = qJD(3) * t248;
t120 = -t222 * t238 + t202;
t237 = t222 * t244;
t121 = -t226 * t238 - t237;
t154 = t248 * t222;
t143 = qJD(3) * pkin(3) - t154;
t155 = t248 * t226;
t40 = t225 * t120 + t221 * t121 + t143 * t261 - t155 * t262;
t22 = -pkin(9) * t107 + t40;
t106 = t129 * t217;
t142 = t225 * t155;
t86 = t143 * t221 + t142;
t41 = -t86 * qJD(4) - t120 * t221 + t225 * t121;
t23 = -pkin(9) * t106 + t41;
t287 = pkin(9) * t166;
t72 = t86 + t287;
t275 = t220 * t72;
t160 = t167 * pkin(9);
t140 = t221 * t155;
t85 = t225 * t143 - t140;
t71 = -t160 + t85;
t67 = pkin(4) * t216 + t71;
t24 = t224 * t67 - t275;
t3 = t24 * qJD(5) + t22 * t224 + t220 * t23;
t37 = qJD(5) * t245 + t106 * t224 - t107 * t220;
t38 = -t103 * qJD(5) - t106 * t220 - t107 * t224;
t274 = t224 * t72;
t25 = t220 * t67 + t274;
t4 = -t25 * qJD(5) - t22 * t220 + t224 * t23;
t98 = Ifges(6,4) * t245;
t51 = Ifges(6,1) * t103 + Ifges(6,5) * t214 + t98;
t312 = t4 * mrSges(6,1) - t3 * mrSges(6,2) + Ifges(6,5) * t37 + Ifges(6,6) * t38 - (Ifges(6,5) * t245 - Ifges(6,6) * t103) * t214 / 0.2e1 - (-Ifges(6,2) * t103 + t51 + t98) * t245 / 0.2e1 - t111 * (mrSges(6,1) * t103 + mrSges(6,2) * t245) - (Ifges(6,1) * t245 - t280) * t103 / 0.2e1;
t311 = t129 / 0.2e1;
t310 = -t166 / 0.2e1;
t309 = t245 / 0.2e1;
t277 = t103 * t25;
t264 = qJD(3) * t222;
t138 = -t199 * t264 + t202;
t263 = qJD(3) * t226;
t139 = -t199 * t263 - t237;
t236 = t138 * t226 - t139 * t222;
t50 = Ifges(6,2) * t245 + Ifges(6,6) * t214 + t280;
t307 = t50 / 0.2e1;
t146 = t225 * t203 - t204 * t221;
t286 = pkin(9) * t190;
t112 = t146 - t286;
t185 = t189 * pkin(9);
t113 = t185 + t147;
t62 = t112 * t224 - t113 * t220;
t306 = t62 * qJD(5) + t314 * t220 - t313 * t224;
t63 = t112 * t220 + t113 * t224;
t305 = -t63 * qJD(5) + t313 * t220 + t314 * t224;
t301 = t24 * t245;
t209 = pkin(3) * t225 + pkin(4);
t259 = qJD(5) * t224;
t260 = qJD(5) * t220;
t265 = t221 * t224;
t89 = t154 * t221 - t142;
t75 = t89 - t287;
t90 = -t225 * t154 - t140;
t76 = -t160 + t90;
t300 = t220 * t76 - t224 * t75 - t209 * t260 + (-t221 * t259 + (-t220 * t225 - t265) * qJD(4)) * pkin(3);
t266 = t220 * t221;
t299 = -t220 * t75 - t224 * t76 + t209 * t259 + (-t221 * t260 + (t224 * t225 - t266) * qJD(4)) * pkin(3);
t208 = pkin(1) * t223 + pkin(7);
t285 = -pkin(8) - t208;
t186 = t285 * t222;
t187 = t208 * t226 + t215;
t123 = t221 * t186 + t225 * t187;
t293 = t103 / 0.2e1;
t291 = t167 / 0.2e1;
t289 = pkin(1) * t227;
t284 = mrSges(5,3) * t166;
t283 = Ifges(4,4) * t222;
t281 = Ifges(5,4) * t167;
t276 = t167 * mrSges(5,3);
t273 = Ifges(4,5) * qJD(3);
t272 = Ifges(4,6) * qJD(3);
t268 = t217 * t222;
t267 = t217 * t226;
t207 = t223 * t253;
t251 = t217 * t264;
t176 = pkin(3) * t251 + t207;
t258 = -qJD(1) - t217;
t257 = -qJD(2) + t217;
t254 = t227 * t278;
t212 = pkin(3) * t264;
t117 = pkin(4) * t130 + t212;
t247 = qJD(3) * t285;
t246 = t199 * (t222 ^ 2 + t226 ^ 2);
t122 = t225 * t186 - t187 * t221;
t125 = t189 * t220 + t190 * t224;
t124 = t189 * t224 - t190 * t220;
t52 = t124 * qJD(5) + t129 * t224 - t130 * t220;
t242 = -t125 * t4 - t24 * t52;
t240 = -mrSges(4,1) * t226 + mrSges(4,2) * t222;
t239 = -t129 * t85 - t190 * t41;
t93 = t122 - t286;
t94 = t185 + t123;
t47 = -t220 * t94 + t224 * t93;
t48 = t220 * t93 + t224 * t94;
t197 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t268;
t198 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t267;
t235 = t197 * t226 + t198 * t222;
t234 = t222 * t197 - t226 * t198;
t163 = -t189 * pkin(4) + t211;
t233 = (Ifges(4,2) * t226 + t283) * t217;
t232 = (mrSges(4,1) * t222 + mrSges(4,2) * t226) * qJD(3);
t148 = t222 * t247 + t226 * t254;
t149 = -t222 * t254 + t226 * t247;
t56 = t225 * t148 + t221 * t149 + t186 * t261 - t187 * t262;
t57 = -qJD(4) * t123 - t148 * t221 + t225 * t149;
t159 = Ifges(5,4) * t166;
t96 = Ifges(5,2) * t166 + Ifges(5,6) * t216 + t281;
t97 = Ifges(5,1) * t167 + Ifges(5,5) * t216 + t159;
t229 = mrSges(6,3) * t301 - t40 * mrSges(5,2) + t103 * t307 - t169 * (mrSges(5,1) * t167 + mrSges(5,2) * t166) - Ifges(5,6) * t107 + Ifges(5,5) * t106 + t41 * mrSges(5,1) + t85 * t284 + t96 * t291 - t167 * (Ifges(5,1) * t166 - t281) / 0.2e1 - t216 * (Ifges(5,5) * t166 - Ifges(5,6) * t167) / 0.2e1 + (-Ifges(5,2) * t167 + t159 + t97) * t310 + t312;
t164 = t233 + t272;
t206 = Ifges(4,4) * t267;
t165 = Ifges(4,1) * t268 + t206 + t273;
t200 = -t217 * pkin(2) - t255;
t53 = -t125 * qJD(5) - t129 * t220 - t130 * t224;
t82 = pkin(4) * t107 + t176;
t228 = (-t189 * t107 + t130 * t310) * Ifges(5,2) + (t189 * t106 - t107 * t190 - t130 * t291 + t166 * t311) * Ifges(5,4) + ((0.3e1 * Ifges(4,4) * t226 + (Ifges(4,1) - 0.2e1 * Ifges(4,2)) * t222) * t217 + t165) * t263 / 0.2e1 + (t124 * t38 + t53 * t309) * Ifges(6,2) + (t124 * t37 + t38 * t125 + t53 * t293 + t52 * t309) * Ifges(6,4) + t236 * mrSges(4,3) - (t233 + t164) * t264 / 0.2e1 + (-t86 * t130 + t40 * t189) * mrSges(5,3) + t169 * (mrSges(5,1) * t130 + mrSges(5,2) * t129) + t216 * (Ifges(5,5) * t129 - Ifges(5,6) * t130) / 0.2e1 + t200 * t232 + t97 * t311 + t240 * t207 + qJD(3) ^ 2 * (Ifges(4,5) * t226 - Ifges(4,6) * t222) / 0.2e1 + (t37 * t125 + t52 * t293) * Ifges(6,1) + (t106 * t190 + t129 * t291) * Ifges(5,1) + t53 * t307 + (t3 * t124 + t25 * t53) * mrSges(6,3) + (Ifges(4,1) * t226 - t283) * t251 + t52 * t51 / 0.2e1 + t111 * (-mrSges(6,1) * t53 + mrSges(6,2) * t52) + t82 * (-mrSges(6,1) * t124 + mrSges(6,2) * t125) - t130 * t96 / 0.2e1 + t176 * (-mrSges(5,1) * t189 + mrSges(5,2) * t190) + t214 * (Ifges(6,5) * t52 + Ifges(6,6) * t53) / 0.2e1;
t213 = t223 * t278;
t210 = -pkin(2) - t289;
t201 = t211 - t289;
t196 = t213 + t212;
t179 = t240 * t217;
t178 = pkin(3) * t265 + t209 * t220;
t177 = -pkin(3) * t266 + t209 * t224;
t168 = t217 * t232;
t151 = t163 - t289;
t137 = mrSges(5,1) * t216 - t276;
t136 = -mrSges(5,2) * t216 + t284;
t131 = pkin(3) * t268 + pkin(4) * t167;
t110 = -mrSges(5,1) * t166 + mrSges(5,2) * t167;
t108 = t117 + t213;
t88 = mrSges(6,1) * t214 - mrSges(6,3) * t103;
t87 = -mrSges(6,2) * t214 + mrSges(6,3) * t245;
t58 = mrSges(5,1) * t107 + mrSges(5,2) * t106;
t55 = -mrSges(6,1) * t245 + mrSges(6,2) * t103;
t44 = t57 - t288;
t43 = -t126 + t56;
t27 = t224 * t71 - t275;
t26 = -t220 * t71 - t274;
t9 = -mrSges(6,1) * t38 + mrSges(6,2) * t37;
t6 = -t48 * qJD(5) - t220 * t43 + t224 * t44;
t5 = t47 * qJD(5) + t220 * t44 + t224 * t43;
t1 = [t228 + (-t37 * t47 + t38 * t48 + t242) * mrSges(6,3) + (-t106 * t122 - t107 * t123 + t239) * mrSges(5,3) + m(5) * (t122 * t41 + t123 * t40 + t169 * t196 + t176 * t201 + t56 * t86 + t57 * t85) + m(6) * (t108 * t111 + t151 * t82 + t24 * t6 + t25 * t5 + t3 * t48 + t4 * t47) + ((t179 + m(4) * (qJD(1) * t210 + t200) + t258 * mrSges(3,1)) * t223 + (m(4) * t246 + t258 * mrSges(3,2) - t234) * t227) * t278 + t5 * t87 + t6 * t88 + t108 * t55 + t56 * t136 + t57 * t137 + t151 * t9 + t196 * t110 + (m(4) * t236 - t235 * qJD(3)) * t208 + t201 * t58 + t210 * t168; t228 + ((t257 * mrSges(3,2) + t234) * t227 + (t257 * mrSges(3,1) - t110 - t179 - t55) * t223) * t279 + t305 * t88 + t306 * t87 + (-t37 * t62 + t38 * t63 + t242) * mrSges(6,3) + (-t106 * t146 - t107 * t147 + t239) * mrSges(5,3) + (t222 * pkin(3) * t110 - t235 * pkin(7)) * qJD(3) + t297 * t137 + t298 * t136 + t117 * t55 + t163 * t9 - pkin(2) * t168 + t211 * t58 + (t163 * t82 + t3 * t63 + t4 * t62 + t306 * t25 + t305 * t24 + (t117 - t256) * t111) * m(6) + (t146 * t41 + t147 * t40 + t176 * t211 + t298 * t86 + t297 * t85 + (t212 - t256) * t169) * m(5) + (-(t200 * t223 + t227 * t246) * t279 - pkin(2) * t207 + t236 * pkin(7)) * m(4); t86 * t276 + t235 * t199 - m(5) * (t85 * t89 + t86 * t90) + (-t177 * t37 + t178 * t38 + t277) * mrSges(6,3) + t300 * t88 + t299 * t87 + t229 - t131 * t55 - t90 * t136 - t89 * t137 - t138 * mrSges(4,2) + t139 * mrSges(4,1) + ((t273 / 0.2e1 - t165 / 0.2e1 - t206 / 0.2e1 - t200 * mrSges(4,2)) * t226 + (-t272 / 0.2e1 + t164 / 0.2e1 - t200 * mrSges(4,1) + (t283 / 0.2e1 + (Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1) * t226) * t217 + (-m(5) * t169 - t110) * pkin(3)) * t222) * t217 + (t136 * t261 - t137 * t262 + m(5) * (t221 * t40 + t225 * t41 + t86 * t261 - t85 * t262) + (-t106 * t225 - t107 * t221) * mrSges(5,3)) * pkin(3) + (-t111 * t131 + t4 * t177 + t3 * t178 + t24 * t300 + t25 * t299) * m(6); mrSges(6,3) * t277 + (t137 + t276) * t86 + t229 - t27 * t87 - t26 * t88 - t85 * t136 + (-t167 * t55 + (-t220 * t88 + t224 * t87) * qJD(5) + (t220 * t38 - t224 * t37) * mrSges(6,3) + (-t111 * t167 + t220 * t3 + t224 * t4 - t24 * t260 + t25 * t259) * m(6)) * pkin(4) - m(6) * (t24 * t26 + t25 * t27); t50 * t293 - t24 * t87 + t25 * t88 + (t301 + t277) * mrSges(6,3) + t312;];
tauc = t1(:);
