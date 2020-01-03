% Calculate vector of inverse dynamics joint torques for
% S5RRPPR9
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR9_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR9_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR9_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR9_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR9_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR9_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR9_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:40:42
% EndTime: 2019-12-31 19:40:58
% DurationCPUTime: 10.18s
% Computational Cost: add. (2200->495), mult. (4517->635), div. (0->0), fcn. (2259->6), ass. (0->233)
t329 = Ifges(5,1) + Ifges(4,3);
t328 = Ifges(4,4) + Ifges(3,5);
t312 = Ifges(3,6) - Ifges(4,6);
t327 = Ifges(4,6) - Ifges(5,5);
t326 = -Ifges(5,5) - t312;
t283 = -pkin(3) - pkin(7);
t132 = -pkin(2) + t283;
t144 = -pkin(2) - pkin(3);
t325 = -m(5) * t144 - m(6) * t132 - mrSges(5,2) + mrSges(6,3);
t138 = sin(qJ(5));
t141 = cos(qJ(5));
t139 = sin(qJ(2));
t142 = cos(qJ(2));
t226 = qJD(1) * qJD(2);
t84 = qJDD(1) * t139 + t142 * t226;
t70 = t84 * pkin(6);
t201 = qJDD(3) + t70;
t230 = t139 * qJD(1);
t147 = -qJ(4) * t84 - qJD(4) * t230 + t201;
t16 = qJDD(2) * t132 + t147;
t119 = t139 * qJD(3);
t135 = qJDD(1) * pkin(1);
t210 = t139 * t226;
t224 = t142 * qJDD(1);
t83 = t210 - t224;
t27 = t83 * pkin(2) - t84 * qJ(3) - qJD(1) * t119 - t135;
t167 = qJDD(4) - t27;
t7 = pkin(4) * t84 + t283 * t83 + t167;
t193 = t139 * pkin(4) + t142 * pkin(7);
t229 = t142 * qJD(1);
t59 = -qJD(1) * pkin(1) - pkin(2) * t229 - qJ(3) * t230;
t39 = pkin(3) * t229 + qJD(4) - t59;
t31 = qJD(1) * t193 + t39;
t105 = qJ(4) * t230;
t114 = pkin(6) * t230;
t228 = qJD(3) + t114;
t212 = -t105 + t228;
t37 = qJD(2) * t132 + t212;
t9 = -t138 * t37 + t141 * t31;
t1 = qJD(5) * t9 + t138 * t7 + t141 * t16;
t324 = t1 * mrSges(6,2);
t10 = t138 * t31 + t141 * t37;
t2 = -qJD(5) * t10 - t138 * t16 + t141 * t7;
t323 = t2 * mrSges(6,1);
t314 = t142 / 0.2e1;
t140 = sin(qJ(1));
t322 = g(2) * t140;
t321 = Ifges(5,6) + t328;
t187 = t142 * mrSges(4,1) + t139 * mrSges(4,3);
t189 = mrSges(3,1) * t142 - mrSges(3,2) * t139;
t320 = t187 + t189;
t122 = t139 * mrSges(5,1);
t307 = -t142 * mrSges(5,2) + t122;
t319 = -t142 * mrSges(6,3) - t307;
t235 = qJD(2) * t141;
t76 = t138 * t229 - t235;
t28 = qJD(5) * t76 - qJDD(2) * t138 + t141 * t83;
t318 = -t28 / 0.2e1;
t77 = qJD(2) * t138 + t141 * t229;
t29 = qJD(5) * t77 - qJDD(2) * t141 - t138 * t83;
t317 = -t29 / 0.2e1;
t74 = qJDD(5) + t84;
t316 = -t74 / 0.2e1;
t315 = -m(6) - m(5);
t215 = mrSges(5,3) * t229;
t270 = qJD(2) * mrSges(5,1) - mrSges(6,1) * t76 - mrSges(6,2) * t77 - t215;
t218 = mrSges(4,2) * t230;
t311 = mrSges(3,3) * t230 + t218 + (-mrSges(3,1) - mrSges(4,1)) * qJD(2);
t217 = mrSges(4,2) * t229;
t93 = qJD(2) * mrSges(4,3) + t217;
t310 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t229 + t93;
t184 = mrSges(6,1) * t138 + mrSges(6,2) * t141;
t134 = qJD(2) * qJ(3);
t115 = pkin(6) * t229;
t81 = -qJ(4) * t229 + t115;
t52 = -t134 - t81;
t42 = qJD(2) * pkin(4) - t52;
t309 = t184 * t42;
t269 = Ifges(3,4) * t139;
t179 = t142 * Ifges(3,2) + t269;
t101 = qJD(5) + t230;
t280 = Ifges(6,4) * t77;
t23 = Ifges(6,2) * t76 + Ifges(6,6) * t101 - t280;
t308 = Ifges(3,6) * qJD(2) + qJD(1) * t179 + t138 * t23;
t306 = qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t305 = -t139 * t312 + t142 * t328;
t111 = pkin(6) * t224;
t69 = -pkin(6) * t210 + t111;
t304 = t139 * t70 + t142 * t69;
t40 = t69 + t306;
t43 = -qJDD(2) * pkin(2) + t201;
t303 = t139 * t43 + t142 * t40;
t302 = t1 * t141 - t138 * t2;
t143 = cos(qJ(1));
t301 = g(1) * t143 + t322;
t227 = m(4) - t315;
t113 = Ifges(3,4) * t229;
t262 = Ifges(4,5) * t142;
t183 = t139 * Ifges(4,1) - t262;
t300 = Ifges(3,1) * t230 - t77 * Ifges(6,5) + t76 * Ifges(6,6) + t101 * Ifges(6,3) + qJD(1) * t183 + qJD(2) * t328 + t113;
t112 = Ifges(4,5) * t230;
t267 = Ifges(5,4) * t139;
t182 = -t142 * Ifges(5,1) - t267;
t67 = Ifges(6,4) * t76;
t24 = -Ifges(6,1) * t77 + Ifges(6,5) * t101 + t67;
t299 = -Ifges(4,3) * t229 + qJD(1) * t182 + qJD(2) * t327 + t141 * t24 + t112;
t298 = mrSges(2,1) + t320;
t185 = mrSges(6,1) * t141 - mrSges(6,2) * t138;
t253 = t142 * mrSges(4,3);
t255 = t142 * mrSges(5,1);
t297 = -(m(6) * pkin(4) + t185) * t142 - t255 - t253 + (m(4) * pkin(2) + mrSges(4,1) + t325) * t139;
t266 = Ifges(5,4) * t142;
t296 = t139 * t267 + (t262 - t266 + (-Ifges(5,2) + t329) * t139) * t142;
t295 = mrSges(2,2) - mrSges(3,3) - mrSges(4,2) + mrSges(5,3);
t294 = qJ(3) * t227;
t236 = qJD(2) * t139;
t220 = pkin(6) * t236;
t163 = qJD(4) * t142 + t220;
t19 = -qJ(4) * t83 + qJD(1) * t163 - t111 - t306;
t11 = mrSges(6,1) * t74 - mrSges(6,3) * t28;
t12 = -mrSges(6,2) * t74 + mrSges(6,3) * t29;
t191 = t10 * t138 + t9 * t141;
t293 = m(6) * (-qJD(5) * t191 + t302) + t141 * t12 - t138 * t11;
t292 = m(5) * t52 - m(6) * t42 - t270;
t289 = Ifges(6,1) * t318 + Ifges(6,4) * t317 + Ifges(6,5) * t316;
t287 = -t77 / 0.2e1;
t286 = t77 / 0.2e1;
t281 = -t142 / 0.2e1;
t279 = pkin(6) * t139;
t278 = pkin(6) * t142;
t127 = t142 * pkin(2);
t273 = -qJD(1) / 0.2e1;
t272 = -qJD(5) / 0.2e1;
t271 = pkin(6) - qJ(4);
t137 = qJ(3) + pkin(4);
t268 = Ifges(3,4) * t142;
t265 = Ifges(6,4) * t138;
t264 = Ifges(6,4) * t141;
t263 = Ifges(4,5) * t139;
t47 = qJDD(2) * mrSges(5,2) - t84 * mrSges(5,3);
t248 = t138 * t139;
t247 = t138 * t140;
t246 = t138 * t142;
t120 = t139 * qJ(3);
t245 = t139 * t141;
t244 = t139 * t143;
t243 = t140 * t141;
t242 = t141 * t142;
t241 = t141 * t143;
t240 = t142 * t143;
t234 = qJD(2) * t142;
t239 = qJ(3) * t234 + t119;
t238 = t127 + t120;
t237 = t143 * pkin(1) + t140 * pkin(6);
t233 = qJD(5) * t138;
t232 = qJD(5) * t141;
t231 = qJD(5) * t142;
t221 = Ifges(6,5) * t28 + Ifges(6,6) * t29 + Ifges(6,3) * t74;
t216 = mrSges(5,3) * t230;
t214 = t142 * pkin(3) + t238;
t213 = t144 * qJD(2);
t211 = t84 * mrSges(5,1) + t83 * mrSges(5,2);
t200 = -pkin(1) - t120;
t199 = -t226 / 0.2e1;
t198 = t226 / 0.2e1;
t46 = -qJDD(2) * mrSges(4,1) + t84 * mrSges(4,2);
t195 = pkin(2) * t240 + qJ(3) * t244 + t237;
t192 = g(1) * t140 - g(2) * t143;
t188 = mrSges(3,1) * t139 + mrSges(3,2) * t142;
t181 = Ifges(6,1) * t141 - t265;
t180 = Ifges(6,1) * t138 + t264;
t177 = -t139 * Ifges(5,2) - t266;
t176 = -Ifges(6,2) * t138 + t264;
t175 = Ifges(6,2) * t141 + t265;
t173 = Ifges(5,5) * t139 - Ifges(5,6) * t142;
t172 = Ifges(6,5) * t141 - Ifges(6,6) * t138;
t171 = Ifges(6,5) * t138 + Ifges(6,6) * t141;
t34 = -mrSges(6,2) * t101 + mrSges(6,3) * t76;
t35 = mrSges(6,1) * t101 + mrSges(6,3) * t77;
t170 = -t138 * t35 + t141 * t34;
t169 = -t138 * t34 - t141 * t35;
t166 = t214 + t193;
t36 = pkin(1) + t166;
t94 = t271 * t139;
t20 = -t138 * t94 + t141 * t36;
t21 = t138 * t36 + t141 * t94;
t85 = -qJD(2) * pkin(2) + t228;
t90 = t115 + t134;
t168 = -t90 * t139 + t85 * t142;
t164 = pkin(1) * t188;
t162 = t39 * (t139 * mrSges(5,2) + t255);
t161 = t59 * (mrSges(4,1) * t139 - t253);
t160 = t139 * (Ifges(3,1) * t142 - t269);
t156 = pkin(4) * t142 + t132 * t139;
t154 = t138 * t231 + t139 * t235;
t153 = t138 * t236 - t141 * t231;
t150 = Ifges(6,5) * t142 + t139 * t181;
t149 = Ifges(6,6) * t142 + t139 * t176;
t148 = Ifges(6,3) * t142 + t139 * t172;
t128 = t143 * pkin(6);
t107 = qJ(3) * t229;
t87 = qJD(2) * mrSges(5,2) - t216;
t86 = -pkin(1) - t238;
t82 = t307 * qJD(1);
t80 = pkin(2) * t230 - t107;
t79 = t187 * qJD(1);
t78 = t114 - t105;
t66 = pkin(1) + t214;
t65 = t184 * t142;
t63 = t139 * t241 - t247;
t62 = -t138 * t244 - t243;
t61 = -t138 * t143 - t139 * t243;
t60 = t139 * t247 - t241;
t54 = -Ifges(5,6) * qJD(2) + qJD(1) * t177;
t51 = -qJD(4) * t139 + t234 * t271;
t50 = pkin(2) * t236 - t239;
t48 = -mrSges(4,2) * t83 + qJDD(2) * mrSges(4,3);
t45 = -qJDD(2) * mrSges(5,1) - mrSges(5,3) * t83;
t44 = t144 * t230 + t107;
t41 = t213 + t212;
t38 = t139 * t213 + t239;
t33 = qJD(1) * t156 + t107;
t30 = qJD(2) * t156 + t239;
t18 = qJDD(2) * t144 + t147;
t17 = qJDD(2) * pkin(4) - t19;
t15 = t138 * t33 + t141 * t81;
t14 = -t138 * t81 + t141 * t33;
t13 = -pkin(3) * t83 + t167;
t8 = -mrSges(6,1) * t29 + mrSges(6,2) * t28;
t5 = t28 * Ifges(6,4) + t29 * Ifges(6,2) + t74 * Ifges(6,6);
t4 = -qJD(5) * t21 - t138 * t51 + t141 * t30;
t3 = qJD(5) * t20 + t138 * t30 + t141 * t51;
t6 = [(-m(3) * t237 - m(4) * t195 - t63 * mrSges(6,1) - t62 * mrSges(6,2) + t315 * (pkin(3) * t240 - qJ(4) * t140 + t195) + t295 * t140 + (-m(6) * t193 - t298 + t319) * t143) * g(2) + (-t54 / 0.2e1 - t10 * mrSges(6,2) + t9 * mrSges(6,1) + t311 * pkin(6) - t41 * mrSges(5,3) + t85 * mrSges(4,2) + t300 / 0.2e1) * t234 + (t139 * (Ifges(4,1) * t142 + t263) + t142 * (-Ifges(3,2) * t139 + t268) + t160) * t198 + (-t308 / 0.2e1 - t52 * mrSges(5,3) - t90 * mrSges(4,2) + t299 / 0.2e1 - t292 * qJ(4)) * t236 + (-mrSges(3,1) * t279 - mrSges(3,2) * t278 + t321 * t139 + t327 * t281 + (Ifges(3,6) - t326) * t314) * qJDD(2) + ((-Ifges(5,4) + Ifges(4,5)) * t84 + t329 * t83) * t281 - t139 * (Ifges(5,4) * t83 - Ifges(5,2) * t84) / 0.2e1 + (t221 + (Ifges(3,1) + Ifges(4,1)) * t84 + (-Ifges(3,4) + Ifges(4,5)) * t83) * t139 / 0.2e1 - t310 * t220 + t5 * t246 / 0.2e1 + (-t173 / 0.2e1 + t305 / 0.2e1) * qJD(2) ^ 2 + t189 * t135 + t76 * (qJD(2) * t149 + t175 * t231) / 0.2e1 + t101 * (qJD(2) * t148 + t171 * t231) / 0.2e1 - t164 * t226 + t94 * t47 - pkin(1) * (mrSges(3,1) * t83 + mrSges(3,2) * t84) + (-t139 * t18 + t142 * t19) * mrSges(5,3) + t86 * (mrSges(4,1) * t83 - mrSges(4,3) * t84) + t303 * mrSges(4,2) + t51 * t87 - t50 * t79 + t38 * t82 - t17 * t65 - t139 * t324 + t13 * t307 + t66 * t211 + t3 * t34 + t4 * t35 + t21 * t12 + t20 * t11 + (qJD(2) * t150 + t180 * t231) * t287 + t242 * t289 + t48 * t278 + t46 * t279 + qJD(2) * t161 + qJD(2) * t162 + (t1 * t246 - t10 * t153 - t154 * t9 + t2 * t242) * mrSges(6,3) + (-t61 * mrSges(6,1) - t60 * mrSges(6,2) + t315 * (-qJ(4) * t143 + t128) + (-m(4) - m(3)) * t128 + t295 * t143 + (-m(4) * (t200 - t127) - m(6) * (-t137 * t139 - pkin(1)) - m(5) * t200 + t122 + m(3) * pkin(1) + t325 * t142 + t298) * t140) * g(1) + m(5) * (t13 * t66 + t18 * t94 + t38 * t39 + t41 * t51) + m(6) * (t1 * t21 + t10 * t3 + t2 * t20 + t4 * t9) + t42 * (mrSges(6,1) * t153 + mrSges(6,2) * t154) + m(4) * (t27 * t86 + t50 * t59 + (qJD(2) * t168 + t303) * pkin(6)) + (-t278 * t83 + t279 * t84 + t304) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(6) * t304) + t296 * t199 + t139 * t323 + t292 * t163 + Ifges(2,3) * qJDD(1) + (Ifges(3,4) * t84 - Ifges(3,2) * t83) * t314 + t74 * (Ifges(6,3) * t139 - t142 * t172) / 0.2e1 + t29 * (Ifges(6,6) * t139 - t142 * t176) / 0.2e1 - t84 * t177 / 0.2e1 - t83 * t179 / 0.2e1 + t28 * (Ifges(6,5) * t139 - t142 * t181) / 0.2e1 - t27 * t187 + (-m(5) * t19 + m(6) * t17 - t45 + t8) * (-t142 * qJ(4) + t278) + (t138 * t24 + t141 * t23) * t231 / 0.2e1 + (t139 * Ifges(3,1) + t183 + t268) * t84 / 0.2e1 + (-t142 * Ifges(4,3) + t182 + t263) * t83 / 0.2e1; (-m(4) * t238 - m(5) * t214 - m(6) * t166 - t139 * t185 + t319 - t320) * g(3) + (t149 * t273 + t176 * t272) * t76 + (-t45 + t48) * qJ(3) + t41 * t215 + t52 * t216 + t90 * t218 + t308 * t230 / 0.2e1 - t230 * t309 + (t181 * t286 - t309) * qJD(5) + t310 * t114 - t311 * t115 + t270 * t78 + t321 * t84 + t144 * t47 + t23 * t233 / 0.2e1 - t24 * t232 / 0.2e1 - t141 * t5 / 0.2e1 + t54 * t229 / 0.2e1 + t137 * t8 - t85 * t217 - t81 * t87 + t80 * t79 - t44 * t82 - t69 * mrSges(3,2) - t70 * mrSges(3,1) + t326 * t83 + (t148 * t273 + t172 * t272) * t101 + (t143 * t297 - t240 * t294) * g(1) - t43 * mrSges(4,1) - pkin(2) * t46 + t40 * mrSges(4,3) - t15 * t34 - t14 * t35 + t18 * mrSges(5,2) - t19 * mrSges(5,1) + (-pkin(2) * t43 + qJ(3) * t40 - t59 * t80) * m(4) + (-t19 * qJ(3) + t144 * t18 - t39 * t44 - t41 * t81 - t52 * t78) * m(5) + t138 * t289 + (-t10 * t15 + t137 * t17 - t14 * t9 + t42 * t78) * m(6) + t301 * t188 + (t10 * t233 + t232 * t9 - t302) * mrSges(6,3) - (Ifges(4,1) * t229 + t112 + t299) * t230 / 0.2e1 - (-Ifges(3,2) * t230 + t113 + t300) * t229 / 0.2e1 + (-t142 * t294 + t297) * t322 + (-t232 * t35 - t233 * t34 + t293) * t132 + (m(4) * t90 - t292 + t93) * qJD(3) + (Ifges(5,3) + Ifges(4,2) + Ifges(3,3)) * qJDD(2) + (t150 * t286 - t9 * (mrSges(6,1) * t142 - mrSges(6,3) * t245) - t10 * (-mrSges(6,2) * t142 - mrSges(6,3) * t248) - t161 - t162 - m(4) * pkin(6) * t168 + (t296 / 0.2e1 + t164 - t160 / 0.2e1) * qJD(1)) * qJD(1) + t173 * t198 + t171 * t316 + t175 * t317 + t180 * t318 + t17 * t185 + t305 * t199; t169 * qJD(5) + (-t93 - t270) * qJD(2) + t227 * t142 * g(3) + ((t169 - t79 - t82) * qJD(1) - t301 * t227) * t139 - m(6) * (qJD(2) * t42 + t191 * t230) + t46 + t47 + t293 + (qJD(2) * t52 - t230 * t39 + t18) * m(5) + (-qJD(2) * t90 + t230 * t59 + t43) * m(4); t141 * t11 + t138 * t12 + t170 * qJD(5) + (t1 * t138 + t2 * t141 + (t10 * t141 - t9 * t138) * qJD(5) + t192) * m(6) + (t13 + t192) * m(5) + (t270 * t142 + (t170 + t87) * t139 - m(5) * (-t139 * t41 + t142 * t52) - m(6) * (-t10 * t245 - t142 * t42 + t248 * t9)) * qJD(1) + t211; -t324 + t323 - t42 * (-mrSges(6,1) * t77 + mrSges(6,2) * t76) + (Ifges(6,1) * t76 + t280) * t286 + t23 * t287 - t101 * (Ifges(6,5) * t76 + Ifges(6,6) * t77) / 0.2e1 - t9 * t34 + t10 * t35 - g(1) * (mrSges(6,1) * t62 - mrSges(6,2) * t63) - g(2) * (-mrSges(6,1) * t60 + mrSges(6,2) * t61) - g(3) * t65 + (-t10 * t77 + t76 * t9) * mrSges(6,3) + t221 - (Ifges(6,2) * t77 + t24 + t67) * t76 / 0.2e1;];
tau = t6;
