% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRRR7
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
% Datum: 2019-12-31 22:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRR7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR7_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR7_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR7_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR7_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:21:23
% EndTime: 2019-12-31 22:21:35
% DurationCPUTime: 5.45s
% Computational Cost: add. (9852->465), mult. (25467->655), div. (0->0), fcn. (18354->8), ass. (0->227)
t208 = sin(qJ(3));
t209 = sin(qJ(2));
t212 = cos(qJ(3));
t213 = cos(qJ(2));
t181 = t208 * t213 + t212 * t209;
t171 = t181 * qJD(1);
t205 = qJD(2) + qJD(3);
t204 = qJD(4) + t205;
t206 = sin(qJ(5));
t210 = cos(qJ(5));
t257 = t212 * t213;
t180 = -t208 * t209 + t257;
t170 = t180 * qJD(1);
t207 = sin(qJ(4));
t211 = cos(qJ(4));
t229 = t170 * t207 + t211 * t171;
t111 = t204 * t210 - t206 * t229;
t146 = t205 * t180;
t135 = t146 * qJD(1);
t147 = t205 * t181;
t136 = t147 * qJD(1);
t242 = t211 * t170 - t171 * t207;
t60 = t242 * qJD(4) + t135 * t211 - t136 * t207;
t40 = qJD(5) * t111 + t210 * t60;
t112 = t204 * t206 + t210 * t229;
t41 = -qJD(5) * t112 - t206 * t60;
t10 = -mrSges(6,1) * t41 + mrSges(6,2) * t40;
t296 = t135 * pkin(8);
t313 = -pkin(7) - pkin(6);
t193 = t313 * t213;
t186 = qJD(1) * t193;
t172 = t208 * t186;
t192 = t313 * t209;
t185 = qJD(1) * t192;
t178 = qJD(2) * pkin(2) + t185;
t139 = t212 * t178 + t172;
t164 = t171 * pkin(8);
t109 = t139 - t164;
t101 = pkin(3) * t205 + t109;
t175 = t212 * t186;
t140 = t178 * t208 - t175;
t297 = pkin(8) * t170;
t110 = t140 + t297;
t259 = t211 * t110;
t63 = t101 * t207 + t259;
t245 = qJD(2) * t313;
t252 = qJD(3) * t212;
t253 = qJD(3) * t208;
t93 = t171 * t245 + t178 * t252 + t186 * t253;
t65 = -pkin(8) * t136 + t93;
t182 = t208 * t192;
t220 = (t257 * t313 - t182) * qJD(2) * qJD(1);
t94 = -t140 * qJD(3) + t220;
t14 = t207 * t65 - t211 * (t94 - t296) + t63 * qJD(4);
t325 = -m(6) * t14 - t10;
t124 = qJD(5) - t242;
t238 = mrSges(6,1) * t206 + mrSges(6,2) * t210;
t266 = t110 * t207;
t62 = t101 * t211 - t266;
t53 = -pkin(4) * t204 - t62;
t226 = t53 * t238;
t233 = Ifges(6,5) * t210 - Ifges(6,6) * t206;
t284 = Ifges(6,4) * t210;
t235 = -Ifges(6,2) * t206 + t284;
t285 = Ifges(6,4) * t206;
t237 = Ifges(6,1) * t210 - t285;
t300 = t210 / 0.2e1;
t301 = -t206 / 0.2e1;
t307 = t112 / 0.2e1;
t286 = Ifges(6,4) * t112;
t47 = Ifges(6,2) * t111 + Ifges(6,6) * t124 + t286;
t106 = Ifges(6,4) * t111;
t48 = Ifges(6,1) * t112 + Ifges(6,5) * t124 + t106;
t324 = t124 * t233 / 0.2e1 + t47 * t301 + t48 * t300 + t226 + t111 * t235 / 0.2e1 + t237 * t307;
t287 = Ifges(5,4) * t229;
t121 = Ifges(5,4) * t242;
t273 = mrSges(5,1) * t204 + mrSges(6,1) * t111 - mrSges(6,2) * t112 - mrSges(5,3) * t229;
t199 = pkin(2) * t212 + pkin(3);
t250 = qJD(4) * t211;
t251 = qJD(4) * t207;
t262 = t207 * t208;
t133 = t199 * t250 + (-t208 * t251 + (t211 * t212 - t262) * qJD(3)) * pkin(2);
t145 = t212 * t185 + t172;
t113 = -t164 + t145;
t144 = -t185 * t208 + t175;
t227 = t144 - t297;
t71 = t211 * t113 + t207 * t227;
t321 = t133 - t71;
t261 = t208 * t211;
t320 = -t113 * t207 + t211 * t227 + t199 * t251 + (t208 * t250 + (t207 * t212 + t261) * qJD(3)) * pkin(2);
t150 = -t212 * t193 + t182;
t54 = pkin(9) * t204 + t63;
t200 = -pkin(2) * t213 - pkin(1);
t191 = qJD(1) * t200;
t148 = -t170 * pkin(3) + t191;
t68 = -pkin(4) * t242 - pkin(9) * t229 + t148;
t21 = -t206 * t54 + t210 * t68;
t22 = t206 * t68 + t210 * t54;
t319 = -t206 * t21 + t210 * t22;
t13 = t62 * qJD(4) + t211 * t65 + (-t178 * t253 + t186 * t252 + t220 - t296) * t207;
t256 = qJD(1) * t209;
t202 = pkin(2) * t256;
t116 = pkin(3) * t136 + qJD(2) * t202;
t61 = qJD(4) * t229 + t135 * t207 + t211 * t136;
t17 = pkin(4) * t61 - pkin(9) * t60 + t116;
t2 = t21 * qJD(5) + t13 * t210 + t17 * t206;
t267 = qJD(5) * t22;
t3 = -t13 * t206 + t17 * t210 - t267;
t318 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t40 + Ifges(6,6) * t41;
t91 = pkin(4) * t229 - pkin(9) * t242;
t231 = t206 * t22 + t21 * t210;
t283 = Ifges(5,5) * t204;
t290 = Ifges(5,1) * t229;
t83 = t121 + t283 + t290;
t317 = t148 * mrSges(5,2) + t83 / 0.2e1 + t283 / 0.2e1 + t121 / 0.2e1 - t231 * mrSges(6,3) + t324;
t316 = t40 / 0.2e1;
t315 = t41 / 0.2e1;
t314 = t61 / 0.2e1;
t311 = pkin(1) * mrSges(3,1);
t310 = pkin(1) * mrSges(3,2);
t309 = -t111 / 0.2e1;
t308 = -t112 / 0.2e1;
t306 = -t124 / 0.2e1;
t304 = t170 / 0.2e1;
t303 = -t171 / 0.2e1;
t302 = t171 / 0.2e1;
t299 = m(4) * t191;
t298 = pkin(3) * t171;
t149 = t212 * t192 + t193 * t208;
t122 = -pkin(8) * t181 + t149;
t123 = pkin(8) * t180 + t150;
t230 = t211 * t122 - t123 * t207;
t295 = t14 * t230;
t294 = t2 * t210;
t293 = t3 * t206;
t292 = t62 * mrSges(5,3);
t291 = mrSges(4,3) * t170;
t289 = Ifges(3,4) * t209;
t288 = Ifges(4,4) * t171;
t282 = Ifges(6,5) * t112;
t281 = Ifges(5,2) * t242;
t280 = Ifges(5,6) * t204;
t279 = Ifges(6,6) * t111;
t278 = Ifges(6,3) * t124;
t277 = t229 * t63;
t276 = t171 * mrSges(4,3);
t271 = Ifges(3,5) * qJD(2);
t270 = Ifges(3,6) * qJD(2);
t269 = qJD(2) * mrSges(3,1);
t268 = qJD(2) * mrSges(3,2);
t265 = t242 * t206;
t264 = t242 * t210;
t166 = pkin(2) * t261 + t207 * t199;
t255 = qJD(1) * t213;
t254 = qJD(2) * t209;
t249 = qJD(5) * t206;
t248 = qJD(5) * t210;
t244 = t271 / 0.2e1;
t243 = -t270 / 0.2e1;
t130 = pkin(2) * t254 + pkin(3) * t147;
t239 = mrSges(6,1) * t210 - mrSges(6,2) * t206;
t236 = Ifges(6,1) * t206 + t284;
t234 = Ifges(6,2) * t210 + t285;
t232 = Ifges(6,5) * t206 + Ifges(6,6) * t210;
t143 = t180 * t207 + t181 * t211;
t154 = -t180 * pkin(3) + t200;
t228 = t211 * t180 - t181 * t207;
t80 = -pkin(4) * t228 - t143 * pkin(9) + t154;
t85 = t122 * t207 + t123 * t211;
t34 = -t206 * t85 + t210 * t80;
t35 = t206 * t80 + t210 * t85;
t165 = -pkin(2) * t262 + t199 * t211;
t187 = t209 * t245;
t188 = t213 * t245;
t97 = t212 * t187 + t208 * t188 + t192 * t252 + t193 * t253;
t75 = t298 + t91;
t222 = -qJD(5) * t231 - t293;
t221 = t222 * mrSges(6,3);
t98 = -t150 * qJD(3) - t187 * t208 + t212 * t188;
t219 = -pkin(8) * t146 + t98;
t8 = t40 * Ifges(6,4) + t41 * Ifges(6,2) + t61 * Ifges(6,6);
t9 = t40 * Ifges(6,1) + t41 * Ifges(6,4) + t61 * Ifges(6,5);
t218 = -t13 * mrSges(5,2) + mrSges(6,3) * t294 + t206 * t9 / 0.2e1 + t236 * t316 + t234 * t315 + t232 * t314 - Ifges(5,6) * t61 + Ifges(5,5) * t60 + t8 * t300 + (-t239 - mrSges(5,1)) * t14 + t324 * qJD(5);
t46 = t278 + t279 + t282;
t82 = t280 + t281 + t287;
t217 = t148 * mrSges(5,1) + t21 * mrSges(6,1) + t46 / 0.2e1 - t82 / 0.2e1 - t287 / 0.2e1 + t282 / 0.2e1 - t280 / 0.2e1 + t279 / 0.2e1 + t278 / 0.2e1 - t22 * mrSges(6,2);
t15 = mrSges(6,1) * t61 - mrSges(6,3) * t40;
t16 = -mrSges(6,2) * t61 + mrSges(6,3) * t41;
t76 = -mrSges(6,2) * t124 + mrSges(6,3) * t111;
t77 = mrSges(6,1) * t124 - mrSges(6,3) * t112;
t216 = -t206 * t15 + m(6) * (-t21 * t248 - t22 * t249 - t293 + t294) + t210 * t16 - t76 * t249 - t77 * t248;
t119 = Ifges(4,2) * t170 + t205 * Ifges(4,6) + t288;
t163 = Ifges(4,4) * t170;
t120 = t171 * Ifges(4,1) + t205 * Ifges(4,5) + t163;
t214 = -(-Ifges(4,2) * t171 + t120 + t163) * t170 / 0.2e1 - (-Ifges(5,2) * t229 + t121 + t83) * t242 / 0.2e1 - (Ifges(5,1) * t242 - t287 + t46) * t229 / 0.2e1 - t191 * (mrSges(4,1) * t171 + mrSges(4,2) * t170) - t22 * (-mrSges(6,2) * t229 - mrSges(6,3) * t265) - t21 * (mrSges(6,1) * t229 - mrSges(6,3) * t264) + t229 * t82 / 0.2e1 - t148 * (mrSges(5,1) * t229 + mrSges(5,2) * t242) + (Ifges(6,6) * t229 + t235 * t242) * t309 + (Ifges(6,3) * t229 + t233 * t242) * t306 + (Ifges(6,5) * t229 + t237 * t242) * t308 - t204 * (Ifges(5,5) * t242 - Ifges(5,6) * t229) / 0.2e1 + t242 * t292 - t242 * t226 - t48 * t264 / 0.2e1 + t47 * t265 / 0.2e1 + t218 + t119 * t302 + (Ifges(4,1) * t170 - t288) * t303 + t139 * t291 - t205 * (Ifges(4,5) * t170 - Ifges(4,6) * t171) / 0.2e1 + Ifges(4,5) * t135 - Ifges(4,6) * t136 - t93 * mrSges(4,2) + t94 * mrSges(4,1);
t201 = Ifges(3,4) * t255;
t190 = mrSges(3,3) * t255 - t268;
t189 = -mrSges(3,3) * t256 + t269;
t169 = Ifges(3,1) * t256 + t201 + t271;
t168 = t270 + (t213 * Ifges(3,2) + t289) * qJD(1);
t162 = pkin(9) + t166;
t161 = -pkin(4) - t165;
t153 = mrSges(4,1) * t205 - t276;
t152 = -mrSges(4,2) * t205 + t291;
t151 = t202 + t298;
t138 = -mrSges(4,1) * t170 + mrSges(4,2) * t171;
t114 = -mrSges(5,2) * t204 + mrSges(5,3) * t242;
t90 = -mrSges(5,1) * t242 + mrSges(5,2) * t229;
t79 = -pkin(8) * t147 + t97;
t74 = qJD(4) * t143 + t146 * t207 + t211 * t147;
t73 = qJD(4) * t228 + t146 * t211 - t147 * t207;
t72 = t202 + t75;
t67 = t109 * t211 - t266;
t66 = t109 * t207 + t259;
t57 = Ifges(6,3) * t61;
t28 = t206 * t91 + t210 * t62;
t27 = -t206 * t62 + t210 * t91;
t26 = t206 * t72 + t210 * t71;
t25 = -t206 * t71 + t210 * t72;
t24 = t206 * t75 + t210 * t67;
t23 = -t206 * t67 + t210 * t75;
t20 = pkin(4) * t74 - pkin(9) * t73 + t130;
t19 = qJD(4) * t85 + t207 * t79 - t211 * t219;
t18 = qJD(4) * t230 + t207 * t219 + t211 * t79;
t5 = -qJD(5) * t35 - t18 * t206 + t20 * t210;
t4 = qJD(5) * t34 + t18 * t210 + t20 * t206;
t1 = [(-t230 * t60 - t61 * t85 - t62 * t73 - t63 * t74) * mrSges(5,3) - t230 * t10 - (-t13 * mrSges(5,3) + t57 / 0.2e1 - Ifges(5,4) * t60 + t116 * mrSges(5,1) + (Ifges(6,3) / 0.2e1 + Ifges(5,2)) * t61 + t318) * t228 + t205 * (Ifges(4,5) * t146 - Ifges(4,6) * t147) / 0.2e1 + t191 * (mrSges(4,1) * t147 + mrSges(4,2) * t146) + (t8 * t301 + t237 * t316 + t235 * t315 + t233 * t314 + t116 * mrSges(5,2) + Ifges(5,1) * t60 - Ifges(5,4) * t61 + t9 * t300 + (mrSges(5,3) + t238) * t14 + (-t2 * t206 - t210 * t3) * mrSges(6,3) + (t232 * t306 + t234 * t309 + t236 * t308 + t53 * t239 + t48 * t301 - t210 * t47 / 0.2e1 - t319 * mrSges(6,3)) * qJD(5)) * t143 + m(5) * (t116 * t154 + t13 * t85 + t130 * t148 + t18 * t63 - t19 * t62 - t295) + m(6) * (t19 * t53 + t2 * t35 + t21 * t5 + t22 * t4 + t3 * t34 - t295) + (-t281 / 0.2e1 + t217) * t74 + (-t135 * t149 - t136 * t150 - t139 * t146 - t140 * t147 + t180 * t93 - t181 * t94) * mrSges(4,3) + (-t180 * t136 - t147 * t304) * Ifges(4,2) + (t180 * t135 - t136 * t181 + t146 * t304 - t147 * t302) * Ifges(4,4) + t200 * (mrSges(4,1) * t136 + mrSges(4,2) * t135) + m(4) * (t139 * t98 + t140 * t97 + t149 * t94 + t150 * t93) + (t169 / 0.2e1 - pkin(6) * t189 + t244 + (-0.2e1 * t310 + 0.3e1 / 0.2e1 * Ifges(3,4) * t213) * qJD(1)) * t213 * qJD(2) - t273 * t19 + (t290 / 0.2e1 + t317) * t73 + t97 * t152 + t98 * t153 + t154 * (mrSges(5,1) * t61 + mrSges(5,2) * t60) + t146 * t120 / 0.2e1 - t147 * t119 / 0.2e1 + t130 * t90 + t18 * t114 + t4 * t76 + t5 * t77 + (t135 * t181 + t146 * t302) * Ifges(4,1) + (-t168 / 0.2e1 - pkin(6) * t190 + t243 + (-0.2e1 * t311 - 0.3e1 / 0.2e1 * t289 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t213) * qJD(1) + (qJD(1) * (-mrSges(4,1) * t180 + mrSges(4,2) * t181) + 0.2e1 * t299 + t138) * pkin(2)) * t254 + t34 * t15 + t35 * t16; t321 * t114 + t140 * t276 + (-t133 * t77 + (-qJD(5) * t76 - t15) * t162 + (-t3 - t267) * mrSges(6,3)) * t206 + ((t152 * t212 - t153 * t208) * qJD(3) + (-t135 * t212 - t136 * t208) * mrSges(4,3)) * pkin(2) + (t133 * t76 + t162 * t16 + (-mrSges(6,3) * t21 - t162 * t77) * qJD(5)) * t210 + t214 + (-t165 * t60 - t166 * t61 + t277) * mrSges(5,3) - t151 * t90 - t145 * t152 - t144 * t153 + t161 * t10 - t26 * t76 - t25 * t77 + ((t244 - t169 / 0.2e1 - t201 / 0.2e1 + qJD(1) * t310 + (t189 - t269) * pkin(6)) * t213 + (t243 + t168 / 0.2e1 + (t311 + t289 / 0.2e1 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t213) * qJD(1) + (t190 + t268) * pkin(6) + (-t138 - t299) * pkin(2)) * t209) * qJD(1) - t273 * t320 + (t13 * t166 - t14 * t165 - t148 * t151 - t320 * t62 + t321 * t63) * m(5) + ((t208 * t93 + t212 * t94 + (-t139 * t208 + t140 * t212) * qJD(3)) * pkin(2) - t139 * t144 - t140 * t145) * m(4) + (t14 * t161 + (t222 + t294) * t162 - t21 * t25 - t22 * t26 + t320 * t53 + t319 * t133) * m(6); -m(5) * (-t62 * t66 + t63 * t67) - m(6) * (t21 * t23 + t22 * t24 + t53 * t66) + t273 * t66 + mrSges(5,3) * t277 + t214 + t221 - t139 * t152 - t67 * t114 - t24 * t76 - t23 * t77 + (-t171 * t90 + (-t207 * t61 - t211 * t60) * mrSges(5,3) + (-t273 * t207 + (-t206 * t77 + t210 * t76 + t114) * t211 + m(6) * (t207 * t53 + t319 * t211)) * qJD(4) + (0.2e1 * t148 * t303 + t13 * t207 - t14 * t211 + (-t207 * t62 + t211 * t63) * qJD(4)) * m(5)) * pkin(3) + t216 * (pkin(3) * t207 + pkin(9)) + (t153 + t276) * t140 - t325 * (-pkin(3) * t211 - pkin(4)); (t292 + (Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t229 - t317) * t242 + t216 * pkin(9) + (t63 * mrSges(5,3) - t217) * t229 - m(6) * (t21 * t27 + t22 * t28 + t53 * t63) + t218 + t221 + t273 * t63 - t62 * t114 - t28 * t76 - t27 * t77 + t325 * pkin(4); t57 - t53 * (mrSges(6,1) * t112 + mrSges(6,2) * t111) + (Ifges(6,1) * t111 - t286) * t308 + t47 * t307 + (Ifges(6,5) * t111 - Ifges(6,6) * t112) * t306 - t21 * t76 + t22 * t77 + (t111 * t21 + t112 * t22) * mrSges(6,3) + (-Ifges(6,2) * t112 + t106 + t48) * t309 + t318;];
tauc = t1(:);
