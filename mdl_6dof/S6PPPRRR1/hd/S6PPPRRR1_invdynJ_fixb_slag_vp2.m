% Calculate vector of inverse dynamics joint torques for
% S6PPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PPPRRR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPPRRR1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPPRRR1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPPRRR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_invdynJ_fixb_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPPRRR1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPPRRR1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPPRRR1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:39:16
% EndTime: 2019-03-08 18:39:29
% DurationCPUTime: 7.03s
% Computational Cost: add. (8147->485), mult. (21328->734), div. (0->0), fcn. (21337->18), ass. (0->240)
t320 = m(6) + m(7);
t157 = sin(qJ(5));
t326 = pkin(10) * t157;
t160 = cos(qJ(5));
t218 = pkin(5) * t157 - pkin(11) * t160;
t267 = cos(pkin(6));
t145 = qJD(1) * t267 + qJD(2);
t154 = cos(pkin(13));
t153 = sin(pkin(6));
t263 = sin(pkin(7));
t227 = t153 * t263;
t219 = t154 * t227;
t266 = cos(pkin(7));
t115 = -qJD(1) * t219 + t145 * t266 + qJD(3);
t152 = sin(pkin(8));
t158 = sin(qJ(4));
t155 = cos(pkin(8));
t150 = sin(pkin(14));
t151 = sin(pkin(13));
t225 = t154 * t266;
t264 = cos(pkin(14));
t181 = (-t150 * t151 + t225 * t264) * t153;
t215 = t264 * t263;
t88 = qJD(1) * t181 + t145 * t215;
t270 = t155 * t88;
t284 = cos(qJ(4));
t223 = t264 * t151;
t182 = t153 * (t150 * t225 + t223);
t228 = t150 * t263;
t89 = qJD(1) * t182 + t145 * t228;
t47 = t284 * t89 + (t115 * t152 + t270) * t158;
t325 = -pkin(10) * qJD(6) * t160 + t218 * qJD(5) - t47;
t156 = sin(qJ(6));
t159 = cos(qJ(6));
t248 = qJD(5) * t159;
t252 = qJD(4) * t157;
t132 = -t156 * t252 + t248;
t324 = -t132 / 0.2e1;
t133 = qJD(5) * t156 + t159 * t252;
t323 = -t133 / 0.2e1;
t250 = qJD(4) * t160;
t146 = qJD(6) - t250;
t322 = -t146 / 0.2e1;
t238 = mrSges(6,3) * t252;
t312 = -qJD(5) * mrSges(6,1) - mrSges(7,1) * t132 + mrSges(7,2) * t133 + t238;
t213 = -t159 * mrSges(7,1) + t156 * mrSges(7,2);
t196 = m(7) * pkin(5) - t213;
t214 = -mrSges(6,1) * t160 + mrSges(6,2) * t157;
t239 = m(7) * pkin(11) + mrSges(7,3);
t321 = pkin(4) * t320 + t157 * t239 + t160 * t196 + mrSges(5,1) - t214;
t242 = qJD(4) * qJD(5);
t138 = qJDD(4) * t160 - t157 * t242;
t319 = t138 / 0.2e1;
t139 = qJDD(4) * t157 + t160 * t242;
t318 = t139 / 0.2e1;
t140 = -pkin(5) * t160 - pkin(11) * t157 - pkin(4);
t235 = t155 * t284;
t221 = t88 * t235;
t236 = t152 * t284;
t188 = -t115 * t236 + t158 * t89 - t221;
t246 = qJD(6) * t156;
t249 = qJD(5) * t157;
t257 = t156 * t160;
t317 = t156 * t249 * pkin(10) - t140 * t246 + t325 * t159 - t257 * t188;
t244 = qJD(6) * t159;
t255 = t159 * t160;
t316 = t140 * t244 + t325 * t156 + t255 * t188 - t248 * t326;
t315 = mrSges(6,1) + t196;
t314 = mrSges(6,2) - t239;
t100 = qJD(6) * t132 + qJDD(5) * t156 + t139 * t159;
t101 = -qJD(6) * t133 + qJDD(5) * t159 - t139 * t156;
t63 = -mrSges(7,1) * t101 + mrSges(7,2) * t100;
t313 = -qJDD(5) * mrSges(6,1) + mrSges(6,3) * t139 + t63;
t148 = Ifges(6,4) * t250;
t130 = Ifges(7,4) * t132;
t93 = Ifges(7,1) * t133 + Ifges(7,5) * t146 + t130;
t311 = Ifges(6,1) * t252 + Ifges(6,5) * qJD(5) + t159 * t93 + t148;
t265 = cos(pkin(12));
t217 = t267 * t265;
t262 = sin(pkin(12));
t184 = t151 * t217 + t154 * t262;
t185 = -t151 * t262 + t154 * t217;
t307 = t185 * t266 - t265 * t227;
t164 = t150 * t184 - t264 * t307;
t172 = t153 * t265 * t266 + t185 * t263;
t310 = t172 * t152 + t164 * t155;
t216 = t267 * t262;
t186 = -t151 * t216 + t154 * t265;
t187 = -t151 * t265 - t154 * t216;
t226 = t153 * t262;
t306 = t187 * t266 + t263 * t226;
t165 = t150 * t186 - t264 * t306;
t173 = t187 * t263 - t226 * t266;
t309 = t173 * t152 + t165 * t155;
t175 = t215 * t267 + t181;
t183 = -t266 * t267 + t219;
t308 = -t183 * t152 + t175 * t155;
t305 = t266 * t152 + t155 * t215;
t113 = -mrSges(6,1) * t138 + mrSges(6,2) * t139;
t304 = mrSges(5,1) * qJDD(4) - t113;
t134 = t214 * qJD(4);
t303 = qJD(4) * mrSges(5,1) - t134;
t144 = qJDD(1) * t267 + qJDD(2);
t114 = -qJDD(1) * t219 + t144 * t266 + qJDD(3);
t232 = qJD(4) * t284;
t220 = t152 * t232;
t251 = qJD(4) * t158;
t259 = t152 * t158;
t86 = qJDD(1) * t181 + t144 * t215;
t87 = qJDD(1) * t182 + t144 * t228;
t25 = t155 * t158 * t86 + qJD(4) * t221 + t114 * t259 + t115 * t220 - t251 * t89 + t284 * t87;
t23 = qJDD(4) * pkin(10) + t25;
t247 = qJD(5) * t160;
t45 = qJD(4) * pkin(10) + t47;
t66 = t114 * t155 - t152 * t86;
t67 = t115 * t155 - t152 * t88;
t7 = t157 * t66 + t160 * t23 + t67 * t247 - t249 * t45;
t34 = t157 * t67 + t160 * t45;
t8 = -qJD(5) * t34 - t157 * t23 + t160 * t66;
t301 = -t157 * t8 + t160 * t7;
t234 = t152 * t251;
t26 = t114 * t236 - t115 * t234 - t158 * t87 - t89 * t232 + t235 * t86 - t251 * t270;
t24 = -qJDD(4) * pkin(4) - t26;
t11 = -t138 * pkin(5) - t139 * pkin(11) + t24;
t5 = qJDD(5) * pkin(11) + t7;
t32 = qJD(5) * pkin(11) + t34;
t43 = qJD(4) * t140 + t188;
t9 = -t156 * t32 + t159 * t43;
t1 = qJD(6) * t9 + t11 * t156 + t159 * t5;
t10 = t156 * t43 + t159 * t32;
t2 = -qJD(6) * t10 + t11 * t159 - t156 * t5;
t300 = t1 * t159 - t156 * t2;
t275 = Ifges(6,4) * t157;
t209 = Ifges(6,2) * t160 + t275;
t297 = Ifges(6,6) * qJD(5) / 0.2e1 + qJD(4) * t209 / 0.2e1 + Ifges(7,5) * t323 + Ifges(7,6) * t324 + Ifges(7,3) * t322;
t296 = g(1) * t173 + g(2) * t172 + g(3) * t183;
t295 = m(3) + m(4) + m(5) + t320;
t212 = t156 * mrSges(7,1) + t159 * mrSges(7,2);
t294 = -t320 * pkin(10) + mrSges(5,2) - t212;
t33 = -t157 * t45 + t160 * t67;
t31 = -qJD(5) * pkin(5) - t33;
t293 = -m(7) * t31 - t312;
t292 = -t2 * mrSges(7,1) + t1 * mrSges(7,2);
t161 = qJD(4) ^ 2;
t273 = Ifges(7,4) * t133;
t92 = Ifges(7,2) * t132 + Ifges(7,6) * t146 + t273;
t291 = -t92 / 0.2e1;
t290 = t100 / 0.2e1;
t289 = t101 / 0.2e1;
t131 = qJDD(6) - t138;
t288 = t131 / 0.2e1;
t286 = t133 / 0.2e1;
t6 = -qJDD(5) * pkin(5) - t8;
t280 = t157 * t6;
t274 = Ifges(6,4) * t160;
t272 = Ifges(7,4) * t156;
t271 = Ifges(7,4) * t159;
t269 = t157 * t188;
t260 = qJD(4) * mrSges(5,2);
t258 = t156 * t157;
t256 = t157 * t159;
t253 = mrSges(5,2) * qJDD(4);
t245 = qJD(6) * t157;
t241 = Ifges(7,5) * t100 + Ifges(7,6) * t101 + Ifges(7,3) * t131;
t237 = mrSges(6,3) * t250;
t233 = t156 * t247;
t222 = t242 / 0.2e1;
t211 = Ifges(7,1) * t159 - t272;
t210 = Ifges(7,1) * t156 + t271;
t208 = -Ifges(7,2) * t156 + t271;
t207 = Ifges(7,2) * t159 + t272;
t206 = Ifges(6,5) * t160 - Ifges(6,6) * t157;
t205 = Ifges(7,5) * t159 - Ifges(7,6) * t156;
t204 = Ifges(7,5) * t156 + Ifges(7,6) * t159;
t166 = -t152 * t175 - t155 * t183;
t102 = t153 * t223 + (t153 * t225 + t263 * t267) * t150;
t52 = t102 * t284 + t158 * t308;
t38 = t157 * t166 + t52 * t160;
t51 = t102 * t158 - t284 * t308;
t15 = t156 * t51 + t159 * t38;
t14 = -t156 * t38 + t159 * t51;
t174 = -t158 * t228 + t284 * t305;
t106 = t158 * t305 + t284 * t228;
t124 = -t152 * t215 + t155 * t266;
t74 = t106 * t160 + t124 * t157;
t57 = -t156 * t74 - t159 * t174;
t58 = -t156 * t174 + t159 * t74;
t73 = t106 * t157 - t160 * t124;
t126 = t155 * t157 + t160 * t259;
t125 = -t160 * t155 + t157 * t259;
t44 = -qJD(4) * pkin(4) + t188;
t198 = t44 * (mrSges(6,1) * t157 + mrSges(6,2) * t160);
t197 = t157 * (Ifges(6,1) * t160 - t275);
t111 = -t156 * t126 - t159 * t236;
t194 = -t159 * t126 + t156 * t236;
t193 = -t156 * t245 + t159 * t247;
t192 = t157 * t244 + t233;
t191 = Ifges(7,5) * t157 + t160 * t211;
t190 = Ifges(7,6) * t157 + t160 * t208;
t189 = Ifges(7,3) * t157 + t160 * t205;
t37 = t52 * t157 - t160 * t166;
t142 = -qJD(5) * mrSges(6,2) + t237;
t136 = t218 * qJD(4);
t122 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t138;
t120 = pkin(10) * t255 + t140 * t156;
t119 = -pkin(10) * t257 + t140 * t159;
t117 = mrSges(7,1) * t146 - mrSges(7,3) * t133;
t116 = -mrSges(7,2) * t146 + mrSges(7,3) * t132;
t110 = qJD(5) * t126 + t157 * t220;
t109 = -qJD(5) * t125 + t160 * t220;
t99 = t106 * qJD(4);
t98 = t174 * qJD(4);
t76 = -mrSges(7,2) * t131 + mrSges(7,3) * t101;
t75 = mrSges(7,1) * t131 - mrSges(7,3) * t100;
t72 = t150 * t306 + t186 * t264;
t71 = t150 * t307 + t184 * t264;
t65 = qJD(6) * t194 - t156 * t109 + t159 * t234;
t64 = qJD(6) * t111 + t159 * t109 + t156 * t234;
t60 = t100 * Ifges(7,1) + t101 * Ifges(7,4) + t131 * Ifges(7,5);
t59 = t100 * Ifges(7,4) + t101 * Ifges(7,2) + t131 * Ifges(7,6);
t56 = t152 * t165 - t155 * t173;
t55 = t152 * t164 - t155 * t172;
t54 = qJD(5) * t74 + t157 * t98;
t53 = -qJD(5) * t73 + t160 * t98;
t49 = t52 * qJD(4);
t48 = t51 * qJD(4);
t42 = -t158 * t309 + t72 * t284;
t41 = t158 * t72 + t284 * t309;
t40 = -t158 * t310 + t71 * t284;
t39 = t158 * t71 + t284 * t310;
t30 = t136 * t156 + t159 * t33;
t29 = t136 * t159 - t156 * t33;
t28 = -qJD(6) * t58 - t156 * t53 + t159 * t99;
t27 = qJD(6) * t57 + t156 * t99 + t159 * t53;
t19 = t157 * t56 + t160 * t42;
t17 = t157 * t55 + t160 * t40;
t13 = -qJD(5) * t37 - t48 * t160;
t4 = qJD(6) * t14 + t13 * t159 + t156 * t49;
t3 = -qJD(6) * t15 - t13 * t156 + t159 * t49;
t12 = [-t52 * t253 + m(5) * (t166 * t66 + t25 * t52 - t47 * t48) + m(3) * (t144 * t267 + (t151 ^ 2 + t154 ^ 2) * t153 ^ 2 * qJDD(1)) + t48 * t260 + m(7) * (t1 * t15 + t10 * t4 + t14 * t2 + t3 * t9) + m(6) * (t13 * t34 + t38 * t7) + m(4) * (t87 * t102 - t114 * t183 + t175 * t86) + t14 * t75 + t15 * t76 + t4 * t116 + t3 * t117 + t38 * t122 + t13 * t142 + m(2) * qJDD(1) + (-m(5) * t26 + m(6) * t24 - t304) * t51 + (m(5) * t188 + m(6) * t44 - t303) * t49 + (-m(6) * t8 + m(7) * t6 + t313) * t37 + (-m(6) * t33 - t293) * (qJD(5) * t38 - t48 * t157) + (-m(2) - t295) * g(3); -t174 * t113 + t27 * t116 + t28 * t117 + t74 * t122 + t99 * t134 + t53 * t142 + t57 * t75 + t58 * t76 + t313 * t73 + t312 * t54 + (-qJD(4) * t98 - qJDD(4) * t106) * mrSges(5,2) + (-qJD(4) * t99 + qJDD(4) * t174) * mrSges(5,1) + m(4) * (t114 * t266 + t215 * t86 + t228 * t87) + m(3) * t144 + m(5) * (t106 * t25 + t124 * t66 + t174 * t26 + t188 * t99 + t47 * t98) + m(6) * (-t174 * t24 - t33 * t54 + t34 * t53 + t44 * t99 + t7 * t74 - t73 * t8) + m(7) * (t1 * t58 + t10 * t27 + t2 * t57 + t28 * t9 + t31 * t54 + t6 * t73) + (-t267 * g(3) + (-g(1) * t262 + g(2) * t265) * t153) * t295; t109 * t142 + t111 * t75 - t194 * t76 + t64 * t116 + t65 * t117 + t126 * t122 + t134 * t234 + (-mrSges(5,1) * t161 - t253) * t259 + t313 * t125 + t312 * t110 + (-mrSges(5,2) * t161 + t304) * t236 + (-t1 * t194 + t10 * t64 + t110 * t31 + t111 * t2 + t125 * t6 + t65 * t9 + t296) * m(7) + (t34 * t109 - t33 * t110 - t8 * t125 + t7 * t126 + (-t24 * t284 + t251 * t44) * t152 + t296) * m(6) + (t66 * t155 + (t284 * t26 + t158 * t25 + (t158 * t188 + t284 * t47) * qJD(4)) * t152 + t296) * m(5) + (t114 + t296) * m(4); -t59 * t258 / 0.2e1 + t60 * t256 / 0.2e1 + (-t1 * t258 - t10 * t192 - t193 * t9 - t2 * t256) * mrSges(7,3) + (-pkin(4) * t24 + ((-t157 * t34 - t160 * t33) * qJD(5) + t301) * pkin(10) - t44 * t47 + (-t157 * t33 + t160 * t34) * t188) * m(6) + (Ifges(6,4) * t318 + Ifges(6,2) * t319 - t241 / 0.2e1 + t188 * t142 + pkin(10) * t122 + (-Ifges(6,2) * t157 + t274) * t222 - Ifges(7,3) * t288 - Ifges(7,6) * t289 - Ifges(7,5) * t290 + t292) * t160 - t188 * t260 + (t9 * mrSges(7,1) - t10 * mrSges(7,2) - pkin(10) * t142 - t297) * t249 + (-g(1) * t42 - g(2) * t40 - g(3) * t52 - t247 * t33 - t249 * t34 + t301) * mrSges(6,3) + t303 * t47 + t132 * (qJD(5) * t190 - t207 * t245) / 0.2e1 + t146 * (qJD(5) * t189 - t204 * t245) / 0.2e1 + (Ifges(6,1) * t139 + Ifges(6,4) * t319 + t205 * t288 + t208 * t289 + t211 * t290) * t157 + (t294 * t42 + t321 * t41) * g(1) + (t294 * t52 + t321 * t51) * g(3) + (t294 * t40 + t321 * t39) * g(2) + t312 * (pkin(10) * t247 + t269) + t311 * t247 / 0.2e1 + t316 * t116 + (t1 * t120 + t119 * t2 + (t247 * t31 + t280) * pkin(10) + t269 * t31 + t317 * t9 + t316 * t10) * m(7) + t317 * t117 + Ifges(5,3) * qJDD(4) + t24 * t214 + qJD(5) ^ 2 * t206 / 0.2e1 - t25 * mrSges(5,2) + t26 * mrSges(5,1) + t197 * t222 + qJD(5) * t198 + t212 * t280 + (qJD(5) * t191 - t210 * t245) * t286 + t233 * t291 + t313 * t326 + t274 * t318 + t209 * t319 - pkin(4) * t113 + t119 * t75 + t120 * t76 + t31 * (mrSges(7,1) * t192 + mrSges(7,2) * t193) + qJDD(5) * (Ifges(6,5) * t157 + Ifges(6,6) * t160) - (t156 * t93 + t159 * t92) * t245 / 0.2e1; (t237 - t142) * t33 - t161 * t197 / 0.2e1 + t297 * t252 + t146 * t31 * t212 + (m(7) * ((-t10 * t156 - t159 * t9) * qJD(6) + t300) - t156 * t75 + t159 * t76 - t117 * t244 - t116 * t246) * pkin(11) + (-t10 * t246 - t244 * t9 + t300) * mrSges(7,3) + t93 * t244 / 0.2e1 + (-t9 * (mrSges(7,1) * t157 - mrSges(7,3) * t255) - t10 * (-mrSges(7,2) * t157 - mrSges(7,3) * t257) - t198) * qJD(4) - (-Ifges(6,2) * t252 + t148 + t311) * t250 / 0.2e1 + (t314 * t38 + t315 * t37) * g(3) + (t314 * t17 - t315 * (-t157 * t40 + t160 * t55)) * g(2) + (t314 * t19 - t315 * (-t157 * t42 + t160 * t56)) * g(1) + (t238 + t293) * t34 - t206 * t242 / 0.2e1 + Ifges(6,3) * qJDD(5) + t6 * t213 - t7 * mrSges(6,2) + t8 * mrSges(6,1) + t204 * t288 + t207 * t289 + t210 * t290 + t246 * t291 - pkin(5) * t63 - t30 * t116 - t29 * t117 + Ifges(6,6) * t138 + Ifges(6,5) * t139 + t159 * t59 / 0.2e1 + (t132 * t208 + t133 * t211 + t146 * t205) * qJD(6) / 0.2e1 - (t132 * t190 + t133 * t191 + t146 * t189) * qJD(4) / 0.2e1 + (t250 * t92 + t60) * t156 / 0.2e1 + (-pkin(5) * t6 - t10 * t30 - t29 * t9) * m(7); -t31 * (mrSges(7,1) * t133 + mrSges(7,2) * t132) + (Ifges(7,1) * t132 - t273) * t323 + t92 * t286 + (Ifges(7,5) * t132 - Ifges(7,6) * t133) * t322 - t9 * t116 + t10 * t117 - g(1) * ((-t156 * t19 + t159 * t41) * mrSges(7,1) + (-t156 * t41 - t159 * t19) * mrSges(7,2)) - g(2) * ((-t156 * t17 + t159 * t39) * mrSges(7,1) + (-t156 * t39 - t159 * t17) * mrSges(7,2)) - g(3) * (mrSges(7,1) * t14 - mrSges(7,2) * t15) + (t10 * t133 + t132 * t9) * mrSges(7,3) + t241 + (-Ifges(7,2) * t133 + t130 + t93) * t324 - t292;];
tau  = t12;
