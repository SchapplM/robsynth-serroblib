% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2018-11-23 16:01
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRP8_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP8_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP8_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:01:17
% EndTime: 2018-11-23 16:01:25
% DurationCPUTime: 8.19s
% Computational Cost: add. (5042->476), mult. (11370->627), div. (0->0), fcn. (7307->6), ass. (0->224)
t317 = Ifges(6,1) + Ifges(7,1);
t326 = Ifges(6,4) - Ifges(7,5);
t316 = Ifges(7,4) + Ifges(6,5);
t315 = Ifges(6,6) - Ifges(7,6);
t325 = Ifges(6,3) + Ifges(7,2);
t147 = sin(pkin(9));
t149 = sin(qJ(3));
t241 = cos(pkin(9));
t272 = cos(qJ(3));
t155 = -t147 * t272 - t149 * t241;
t113 = t155 * qJD(1);
t321 = qJD(5) - t113;
t109 = Ifges(5,4) * t113;
t104 = qJD(3) * t113;
t148 = sin(qJ(5));
t187 = t241 * t272;
t166 = qJD(1) * t187;
t228 = qJD(1) * t149;
t112 = t147 * t228 - t166;
t150 = cos(qJ(5));
t91 = qJD(3) * t148 - t112 * t150;
t59 = qJD(5) * t91 + t104 * t148;
t287 = -t59 / 0.2e1;
t324 = Ifges(6,2) * t287;
t225 = qJD(5) * t150;
t226 = qJD(5) * t148;
t201 = t272 * qJD(4);
t151 = -pkin(1) - pkin(7);
t131 = qJD(1) * t151 + qJD(2);
t227 = qJD(3) * t149;
t210 = t131 * t227;
t153 = -t210 + (qJ(4) * t227 - t201) * qJD(1);
t202 = qJD(3) * t272;
t191 = t131 * t202;
t224 = t149 * qJD(4);
t87 = t191 + (-qJ(4) * t202 - t224) * qJD(1);
t44 = t147 * t153 + t241 * t87;
t223 = qJD(1) * qJD(3);
t200 = t149 * t223;
t103 = -qJD(3) * t166 + t147 * t200;
t145 = qJD(1) * qJD(2);
t189 = qJD(1) * t202;
t124 = pkin(3) * t189 + t145;
t49 = -pkin(4) * t103 - pkin(8) * t104 + t124;
t107 = -qJ(4) * t228 + t131 * t149;
t195 = t241 * t107;
t123 = t272 * t131;
t203 = qJD(1) * t272;
t108 = -qJ(4) * t203 + t123;
t99 = qJD(3) * pkin(3) + t108;
t63 = t147 * t99 + t195;
t57 = qJD(3) * pkin(8) + t63;
t125 = pkin(3) * t228 + qJD(1) * qJ(2) + qJD(4);
t68 = -pkin(4) * t113 + pkin(8) * t112 + t125;
t3 = t148 * t49 + t150 * t44 + t68 * t225 - t226 * t57;
t22 = t148 * t68 + t150 * t57;
t4 = -qJD(5) * t22 - t148 * t44 + t150 * t49;
t185 = -t4 * t148 + t3 * t150;
t21 = -t148 * t57 + t150 * t68;
t323 = -t21 * t225 - t22 * t226 + t185;
t304 = qJD(6) - t21;
t15 = -pkin(5) * t321 + t304;
t16 = qJ(6) * t321 + t22;
t1 = -qJ(6) * t103 + qJD(6) * t321 + t3;
t2 = pkin(5) * t103 - t4;
t186 = t1 * t150 + t2 * t148;
t322 = t15 * t225 - t16 * t226 + t186;
t320 = qJ(2) * (m(3) + m(4)) + mrSges(3,3);
t39 = mrSges(6,2) * t103 - mrSges(6,3) * t59;
t40 = -mrSges(7,2) * t59 - mrSges(7,3) * t103;
t258 = t39 + t40;
t167 = t150 * qJD(3) + t112 * t148;
t58 = qJD(5) * t167 + t104 * t150;
t37 = -mrSges(6,1) * t103 - mrSges(6,3) * t58;
t38 = t103 * mrSges(7,1) + t58 * mrSges(7,2);
t259 = -t37 + t38;
t319 = t259 * t148 + t258 * t150;
t288 = t58 / 0.2e1;
t281 = -t103 / 0.2e1;
t313 = -t316 * t103 + t317 * t58 - t326 * t59;
t312 = t167 * t315 + t316 * t91 + t321 * t325;
t268 = Ifges(7,5) * t167;
t89 = Ifges(6,4) * t167;
t311 = t316 * t321 + t317 * t91 - t268 + t89;
t247 = t112 * mrSges(5,3);
t310 = -qJD(3) * mrSges(5,1) - mrSges(6,1) * t167 + mrSges(6,2) * t91 - t247;
t64 = mrSges(7,2) * t167 + mrSges(7,3) * t321;
t271 = mrSges(6,3) * t167;
t65 = -mrSges(6,2) * t321 + t271;
t257 = t64 + t65;
t270 = mrSges(6,3) * t91;
t66 = mrSges(6,1) * t321 - t270;
t67 = -mrSges(7,1) * t321 + mrSges(7,2) * t91;
t256 = t66 - t67;
t309 = t113 * Ifges(5,2);
t20 = mrSges(6,1) * t59 + mrSges(6,2) * t58;
t248 = t104 * mrSges(5,3);
t308 = t248 + t20;
t175 = pkin(5) * t148 - qJ(6) * t150;
t73 = t108 * t147 + t195;
t307 = -qJD(6) * t148 + t175 * t321 - t73;
t306 = Ifges(5,5) * qJD(3);
t305 = Ifges(5,6) * qJD(3);
t229 = qJ(4) - t151;
t183 = mrSges(7,1) * t148 - mrSges(7,3) * t150;
t184 = mrSges(6,1) * t148 + mrSges(6,2) * t150;
t96 = t147 * t107;
t62 = t241 * t99 - t96;
t56 = -qJD(3) * pkin(4) - t62;
t23 = -pkin(5) * t167 - t91 * qJ(6) + t56;
t303 = t23 * t183 + t56 * t184;
t302 = -t148 * t315 + t150 * t316;
t251 = Ifges(7,5) * t148;
t253 = Ifges(6,4) * t148;
t301 = t150 * t317 + t251 - t253;
t298 = -t103 * t325 - t315 * t59 + t316 * t58;
t297 = Ifges(6,4) * t288 + t324 + Ifges(6,6) * t281 - t58 * Ifges(7,5) / 0.2e1 + t103 * Ifges(7,6) / 0.2e1 + Ifges(7,3) * t287;
t114 = -qJD(3) * t187 + t147 * t227;
t115 = t155 * qJD(3);
t120 = t147 * t149 - t187;
t43 = t147 * t87 - t241 * t153;
t242 = t120 * t43;
t295 = -t114 * t63 + t115 * t62 + t242;
t139 = t149 * pkin(3) + qJ(2);
t82 = -pkin(4) * t155 + pkin(8) * t120 + t139;
t126 = t229 * t149;
t127 = t229 * t272;
t86 = -t126 * t241 - t147 * t127;
t255 = t148 * t82 + t150 * t86;
t105 = t227 * t229 - t201;
t106 = -qJD(3) * t127 - t224;
t71 = t147 * t105 + t106 * t241;
t132 = pkin(3) * t202 + qJD(2);
t75 = -pkin(4) * t114 - pkin(8) * t115 + t132;
t9 = -qJD(5) * t255 - t148 * t71 + t150 * t75;
t273 = t148 / 0.2e1;
t88 = Ifges(7,5) * t91;
t31 = Ifges(7,6) * t321 - Ifges(7,3) * t167 + t88;
t269 = Ifges(6,4) * t91;
t34 = Ifges(6,2) * t167 + Ifges(6,6) * t321 + t269;
t294 = t34 * t273 - t148 * t31 / 0.2e1 - t306 / 0.2e1 - t125 * mrSges(5,2);
t293 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t292 = t125 * mrSges(5,1) + t16 * mrSges(7,3) + t21 * mrSges(6,1) - t305 / 0.2e1 - t15 * mrSges(7,1) - t22 * mrSges(6,2);
t291 = qJD(1) ^ 2;
t286 = t59 / 0.2e1;
t285 = t167 / 0.2e1;
t284 = -t167 / 0.2e1;
t283 = -t91 / 0.2e1;
t282 = t91 / 0.2e1;
t280 = -t321 / 0.2e1;
t279 = t321 / 0.2e1;
t278 = -t112 / 0.2e1;
t277 = t112 / 0.2e1;
t276 = -t113 / 0.2e1;
t267 = pkin(3) * t147;
t6 = pkin(5) * t59 - qJ(6) * t58 - qJD(6) * t91 + t43;
t265 = t120 * t6;
t85 = -t126 * t147 + t241 * t127;
t261 = t43 * t85;
t74 = t108 * t241 - t96;
t192 = pkin(3) * t203;
t76 = -t112 * pkin(4) - t113 * pkin(8) + t192;
t25 = t148 * t76 + t150 * t74;
t254 = Ifges(4,4) * t149;
t252 = Ifges(6,4) * t150;
t250 = Ifges(7,5) * t150;
t249 = t103 * mrSges(5,3);
t246 = t112 * Ifges(5,4);
t245 = t113 * mrSges(5,3);
t239 = t114 * t148;
t238 = t114 * t150;
t237 = t115 * t150;
t236 = t120 * t148;
t235 = t120 * t150;
t232 = t148 * t113;
t129 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t203;
t231 = t149 * t129;
t230 = t150 * t113;
t47 = -mrSges(7,1) * t167 - mrSges(7,3) * t91;
t222 = t47 + t310;
t218 = Ifges(4,4) * t272;
t211 = t241 * pkin(3);
t197 = t226 / 0.2e1;
t196 = t225 / 0.2e1;
t194 = -t103 * mrSges(5,1) + t104 * mrSges(5,2);
t70 = -t241 * t105 + t106 * t147;
t138 = -t211 - pkin(4);
t180 = -Ifges(6,2) * t148 + t252;
t177 = Ifges(7,3) * t148 + t250;
t176 = -t150 * pkin(5) - t148 * qJ(6);
t174 = t16 * t148 - t15 * t150;
t172 = t22 * t148 + t21 * t150;
t24 = -t148 * t74 + t150 * t76;
t29 = -t148 * t86 + t150 * t82;
t165 = -Ifges(4,5) * t149 - Ifges(4,6) * t272;
t8 = t148 * t75 + t150 * t71 + t82 * t225 - t226 * t86;
t162 = -t148 * t115 + t120 * t225;
t161 = t120 * t226 + t237;
t160 = qJ(2) * (mrSges(4,1) * t272 - mrSges(4,2) * t149);
t159 = t149 * (-Ifges(4,2) * t272 - t254);
t158 = (Ifges(4,1) * t272 - t254) * qJD(1);
t157 = (-Ifges(4,2) * t149 + t218) * qJD(1);
t122 = (t149 * mrSges(4,1) + mrSges(4,2) * t272) * qJD(1);
t156 = -t148 * t257 - t150 * t256;
t154 = (-Ifges(4,1) * t149 - t218) * t272;
t128 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t228;
t118 = t176 + t138;
t117 = Ifges(4,5) * qJD(3) + t158;
t116 = Ifges(4,6) * qJD(3) + t157;
t94 = -qJD(3) * mrSges(5,2) + t245;
t81 = -mrSges(5,1) * t113 - mrSges(5,2) * t112;
t79 = -t112 * Ifges(5,1) + t109 + t306;
t78 = -t246 + t305 + t309;
t46 = pkin(5) * t91 - qJ(6) * t167;
t41 = -t120 * t175 + t85;
t27 = pkin(5) * t155 - t29;
t26 = -qJ(6) * t155 + t255;
t19 = mrSges(7,1) * t59 - mrSges(7,3) * t58;
t18 = pkin(5) * t112 - t24;
t17 = -qJ(6) * t112 + t25;
t10 = t175 * t115 + (qJD(5) * t176 + qJD(6) * t150) * t120 + t70;
t7 = pkin(5) * t114 - t9;
t5 = -qJ(6) * t114 - qJD(6) * t155 + t8;
t11 = [m(5) * (t124 * t139 + t125 * t132 + t44 * t86 + t63 * t71 + t261) + (t128 * t202 - t129 * t227) * t151 - t183 * t265 - t184 * t242 + t139 * t194 + (t109 / 0.2e1 + Ifges(5,1) * t278 + t79 / 0.2e1 - t294) * t115 + qJD(3) ^ 2 * t165 / 0.2e1 + t56 * (-mrSges(6,1) * t162 + mrSges(6,2) * t161) + t23 * (-mrSges(7,1) * t162 - mrSges(7,3) * t161) + (-t124 * mrSges(5,2) - Ifges(5,1) * t104 - Ifges(5,4) * t103 - t177 * t286 - t180 * t287 + t196 * t34 + t197 * t311 - t281 * t302 - t288 * t301) * t120 + (0.2e1 * t160 - t159 + t154) * t223 + (-m(5) * t62 + m(6) * t56 + t310) * t70 + t311 * t237 / 0.2e1 - (qJD(5) * t31 + t313) * t235 / 0.2e1 - (t157 + t116) * t202 / 0.2e1 - (t158 + t117) * t227 / 0.2e1 + m(6) * (t21 * t9 + t22 * t8 + t255 * t3 + t29 * t4 + t261) + t255 * t39 + (t317 * t161 + t162 * t326) * t282 + (-t161 * t21 + t162 * t22 + t235 * t4 + t236 * t3) * mrSges(6,3) + (t1 * t236 + t15 * t161 + t16 * t162 - t2 * t235) * mrSges(7,2) + (t278 * Ifges(5,4) - t312 / 0.2e1 - Ifges(6,6) * t285 - Ifges(7,6) * t284 + t309 / 0.2e1 + t78 / 0.2e1 - t316 * t282 - t325 * t279 - t292) * t114 + (Ifges(7,5) * t161 - Ifges(7,3) * t162) * t284 + (t161 * t316 + t162 * t315) * t279 + (Ifges(6,4) * t161 + Ifges(6,2) * t162) * t285 + m(7) * (t1 * t26 + t10 * t23 + t15 * t7 + t16 * t5 + t2 * t27 + t41 * t6) + (Ifges(5,4) * t104 - t124 * mrSges(5,1) + Ifges(5,2) * t103 - Ifges(6,6) * t287 - Ifges(7,6) * t286 - t325 * t281 - t316 * t288 - t293 + t44 * mrSges(5,3) - t298 / 0.2e1) * t155 + t86 * t249 + t297 * t236 + t308 * t85 + 0.2e1 * t320 * t145 - t295 * mrSges(5,3) + t29 * t37 + t27 * t38 + t26 * t40 + t41 * t19 + t10 * t47 + t5 * t64 + t8 * t65 + t9 * t66 + t7 * t67 + t71 * t94 + 0.2e1 * qJD(2) * t122 + t132 * t81; (t128 * t272 - t231) * qJD(3) - t320 * t291 + (t19 + t308) * t120 - t222 * t115 + (t148 * t256 - t150 * t257 - t94) * t114 + m(5) * t295 + m(6) * (-t115 * t56 + t21 * t239 - t22 * t238 + t242) + m(7) * (-t115 * t23 - t15 * t239 - t16 * t238 + t265) - (m(5) * t44 + m(6) * t323 + m(7) * t322 + t156 * qJD(5) + t249 + t319) * t155 + (-m(5) * t125 - m(6) * t172 - m(7) * t174 - t122 + t156 - t81) * qJD(1); -t211 * t248 - t63 * t247 + t117 * t228 / 0.2e1 - t165 * t223 / 0.2e1 - t34 * t226 / 0.2e1 - t128 * t123 - mrSges(4,1) * t210 + t116 * t203 / 0.2e1 - Ifges(4,5) * t200 - t81 * t192 + (-t43 * mrSges(6,1) - t6 * mrSges(7,1) - Ifges(7,3) * t286 + t281 * t315 + t297 + t324) * t150 - mrSges(4,2) * t191 - Ifges(4,6) * t189 + (t301 * t91 + t302 * t321) * qJD(5) / 0.2e1 + t78 * t278 + (-t250 + t252) * t288 - t310 * t73 + (t246 + t312) * t277 + t313 * t273 + t251 * t286 + t253 * t287 + (-(-t180 / 0.2e1 + t177 / 0.2e1) * t167 + t303) * qJD(5) + (-t230 / 0.2e1 + t196) * t311 + (-t256 * t225 - t257 * t226 + m(6) * (-qJD(5) * t172 + t185) + m(7) * (-qJD(5) * t174 + t186) + t319) * (pkin(8) + t267) + (Ifges(5,2) * t276 - Ifges(6,6) * t284 - Ifges(7,6) * t285 - t280 * t325 - t283 * t316 + t292) * t112 + (t43 * mrSges(6,2) - t6 * mrSges(7,3) + t281 * t316 + t288 * t317) * t148 + (t138 * t43 - t21 * t24 - t22 * t25 - t56 * t73) * m(6) + (t118 * t6 - t15 * t18 - t16 * t17 + t307 * t23) * m(7) + (-t125 * t192 + t62 * t73 - t63 * t74 + (t147 * t44 - t241 * t43) * pkin(3)) * m(5) + t31 * t197 + (t159 / 0.2e1 - t160 - t154 / 0.2e1) * t291 + t62 * t245 + t249 * t267 + t131 * t231 + t307 * t47 + (-t15 * t230 + t16 * t232 + t322) * mrSges(7,2) + (t21 * t230 + t22 * t232 + t323) * mrSges(6,3) + (Ifges(5,1) * t277 + t177 * t285 + t180 * t284 + t280 * t302 + t283 * t301 + t294 - t303) * t113 - t43 * mrSges(5,1) - t44 * mrSges(5,2) - t17 * t64 - t25 * t65 - t24 * t66 - t18 * t67 - t74 * t94 + Ifges(5,6) * t103 + Ifges(5,5) * t104 + t118 * t19 + t138 * t20 + (t79 + t109) * t276; -t113 * t94 + t222 * t112 + (t257 * t321 - t259) * t150 + (-t256 * t321 + t258) * t148 + t194 + (t1 * t148 + t112 * t23 - t150 * t2 + t321 * (t148 * t15 + t150 * t16)) * m(7) + (t112 * t56 + t148 * t3 + t150 * t4 + t321 * (-t148 * t21 + t150 * t22)) * m(6) + (-t112 * t62 - t113 * t63 + t124) * m(5); (t256 + t270) * t22 + (t167 * t316 - t315 * t91) * t280 + (t167 * t317 - t269 + t31 + t88) * t283 + (Ifges(7,3) * t91 + t268) * t285 + (-t257 + t271) * t21 + t34 * t282 + (-t15 * t167 + t16 * t91) * mrSges(7,2) - t23 * (mrSges(7,1) * t91 - mrSges(7,3) * t167) - t56 * (mrSges(6,1) * t91 + mrSges(6,2) * t167) + (-pkin(5) * t2 + qJ(6) * t1 - t15 * t22 + t16 * t304 - t23 * t46) * m(7) + (-Ifges(6,2) * t91 + t311 + t89) * t284 + t293 - pkin(5) * t38 + qJ(6) * t40 - t46 * t47 + qJD(6) * t64 + t298; -t321 * t64 + t91 * t47 + 0.2e1 * (t2 / 0.2e1 + t16 * t280 + t23 * t282) * m(7) + t38;];
tauc  = t11(:);
