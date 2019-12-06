% Calculate vector of inverse dynamics joint torques for
% S5RPPRR4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-12-05 17:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR4_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR4_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR4_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:44:05
% EndTime: 2019-12-05 17:44:34
% DurationCPUTime: 10.22s
% Computational Cost: add. (4551->427), mult. (11556->586), div. (0->0), fcn. (8496->14), ass. (0->209)
t179 = sin(pkin(8));
t181 = cos(pkin(8));
t232 = t179 ^ 2 + t181 ^ 2;
t223 = qJD(1) * qJD(2);
t158 = qJDD(1) * qJ(2) + t223;
t178 = sin(pkin(9));
t180 = cos(pkin(9));
t184 = sin(qJ(4));
t187 = cos(qJ(4));
t140 = t178 * t187 + t180 * t184;
t132 = t140 * qJD(4);
t196 = qJD(1) * t140;
t269 = t181 * t196 - t132;
t236 = t180 * t187;
t201 = t178 * t184 - t236;
t230 = qJD(1) * t181;
t267 = t201 * qJD(4);
t268 = t201 * t230 - t267;
t68 = (-qJD(1) * t132 - qJDD(1) * t201) * t179;
t286 = Ifges(5,5) * t68;
t183 = sin(qJ(5));
t186 = cos(qJ(5));
t107 = t179 * t196;
t231 = qJD(1) * t179;
t213 = t178 * t231;
t109 = -t184 * t213 + t231 * t236;
t206 = -t186 * t107 - t109 * t183;
t69 = (qJD(1) * t267 - qJDD(1) * t140) * t179;
t21 = qJD(5) * t206 + t183 * t69 + t186 * t68;
t285 = Ifges(6,5) * t21;
t284 = Ifges(5,6) * t69;
t60 = -t107 * t183 + t109 * t186;
t22 = -qJD(5) * t60 - t183 * t68 + t186 * t69;
t283 = Ifges(6,6) * t22;
t224 = m(4) + m(5) + m(6);
t282 = -m(3) - t224;
t220 = qJDD(1) * t181;
t159 = qJDD(4) - t220;
t281 = Ifges(5,3) * t159;
t152 = qJDD(5) + t159;
t280 = Ifges(6,3) * t152;
t203 = mrSges(4,1) * t178 + mrSges(4,2) * t180;
t197 = t203 * t179;
t279 = -mrSges(5,1) * t107 - mrSges(5,2) * t109 - qJD(1) * t197;
t278 = t203 - mrSges(2,2);
t207 = -qJ(3) * t179 - pkin(1);
t145 = -pkin(2) * t181 + t207;
t125 = qJD(1) * t145 + qJD(2);
t114 = t180 * t125;
t238 = t179 * t180;
t219 = pkin(6) * t238;
t193 = -t219 + (-qJ(2) * t178 - pkin(3)) * t181;
t70 = qJD(1) * t193 + t114;
t214 = qJ(2) * t230;
t88 = t178 * t125 + t180 * t214;
t77 = -pkin(6) * t213 + t88;
t38 = t184 * t70 + t187 * t77;
t239 = t178 * t181;
t227 = qJD(3) * t179;
t94 = -qJD(1) * t227 + qJDD(1) * t145 + qJDD(2);
t73 = -t158 * t239 + t180 * t94;
t48 = (-pkin(3) * t181 - t219) * qJDD(1) + t73;
t221 = qJDD(1) * t179;
t209 = t178 * t221;
t237 = t180 * t181;
t74 = t158 * t237 + t178 * t94;
t55 = -pkin(6) * t209 + t74;
t13 = -qJD(4) * t38 - t184 * t55 + t187 * t48;
t6 = pkin(4) * t159 - pkin(7) * t68 + t13;
t225 = qJD(4) * t187;
t226 = qJD(4) * t184;
t12 = t184 * t48 + t187 * t55 + t70 * t225 - t226 * t77;
t7 = pkin(7) * t69 + t12;
t29 = -pkin(7) * t107 + t38;
t243 = t183 * t29;
t161 = qJD(4) - t230;
t37 = -t184 * t77 + t187 * t70;
t28 = -pkin(7) * t109 + t37;
t25 = pkin(4) * t161 + t28;
t8 = t186 * t25 - t243;
t2 = qJD(5) * t8 + t183 * t6 + t186 * t7;
t242 = t186 * t29;
t9 = t183 * t25 + t242;
t3 = -qJD(5) * t9 - t183 * t7 + t186 * t6;
t276 = -t3 * mrSges(6,1) + t2 * mrSges(6,2);
t275 = -t13 * mrSges(5,1) + t12 * mrSges(5,2);
t166 = t180 * pkin(3) + pkin(2);
t177 = pkin(9) + qJ(4);
t169 = cos(t177);
t182 = -pkin(6) - qJ(3);
t204 = -mrSges(3,1) * t181 + mrSges(3,2) * t179;
t274 = -m(4) * t207 - (-m(4) * pkin(2) - mrSges(4,1) * t180 + mrSges(4,2) * t178) * t181 + m(3) * pkin(1) - t204 - m(5) * (-t166 * t181 - pkin(1)) - m(6) * (-(pkin(4) * t169 + t166) * t181 - pkin(1)) + mrSges(2,1) + (mrSges(4,3) - m(5) * t182 + mrSges(5,3) - m(6) * (-pkin(7) + t182) + mrSges(6,3)) * t179;
t83 = -t140 * t183 - t186 * t201;
t273 = qJD(5) * t83 + t183 * t269 + t268 * t186;
t84 = t140 * t186 - t183 * t201;
t272 = -qJD(5) * t84 - t268 * t183 + t186 * t269;
t266 = m(6) * pkin(4);
t271 = mrSges(5,1) + t266;
t270 = -m(3) * qJ(2) - mrSges(3,3);
t265 = -t206 / 0.2e1;
t264 = -t60 / 0.2e1;
t263 = t60 / 0.2e1;
t262 = t8 * mrSges(6,3);
t261 = t9 * mrSges(6,3);
t259 = t109 / 0.2e1;
t155 = qJD(5) + t161;
t258 = -t155 / 0.2e1;
t255 = Ifges(6,4) * t60;
t254 = pkin(3) * t178;
t253 = pkin(4) * t109;
t168 = sin(t177);
t252 = pkin(4) * t168;
t251 = g(1) * t179;
t136 = t180 * t145;
t80 = t136 + t193;
t101 = qJ(2) * t237 + t178 * t145;
t240 = t178 * t179;
t89 = -pkin(6) * t240 + t101;
t42 = t184 * t80 + t187 * t89;
t171 = qJ(5) + t177;
t164 = sin(t171);
t165 = cos(t171);
t188 = cos(qJ(1));
t185 = sin(qJ(1));
t235 = t181 * t185;
t103 = t164 * t235 + t165 * t188;
t104 = -t164 * t188 + t165 * t235;
t250 = t103 * mrSges(6,1) + t104 * mrSges(6,2);
t234 = t181 * t188;
t105 = t164 * t234 - t165 * t185;
t106 = -t164 * t185 - t165 * t234;
t249 = -t105 * mrSges(6,1) + t106 * mrSges(6,2);
t248 = mrSges(5,3) * t107;
t247 = mrSges(5,3) * t109;
t246 = Ifges(5,4) * t109;
t170 = t179 * qJ(2);
t208 = t180 * t221;
t233 = mrSges(4,1) * t209 + mrSges(4,2) * t208;
t141 = pkin(3) * t240 + t170;
t229 = qJD(2) * t179;
t228 = qJD(2) * t181;
t156 = qJ(2) * t231 + qJD(3);
t138 = t179 * t158 + qJDD(3);
t218 = t280 + t283 + t285;
t217 = t281 + t284 + t286;
t126 = pkin(3) * t213 + t156;
t211 = -t69 * mrSges(5,1) + t68 * mrSges(5,2);
t210 = -t22 * mrSges(6,1) + t21 * mrSges(6,2);
t102 = pkin(3) * t209 + t138;
t41 = -t184 * t89 + t187 * t80;
t205 = -mrSges(3,1) * t220 + mrSges(3,2) * t221;
t202 = -mrSges(6,1) * t164 - mrSges(6,2) * t165;
t123 = t201 * t179;
t35 = -pkin(4) * t181 + pkin(7) * t123 + t41;
t122 = t140 * t179;
t36 = -pkin(7) * t122 + t42;
t14 = -t183 * t36 + t186 * t35;
t15 = t183 * t35 + t186 * t36;
t71 = -t122 * t186 + t123 * t183;
t72 = -t122 * t183 - t123 * t186;
t200 = t218 - t276;
t199 = -mrSges(4,1) * t181 - mrSges(4,3) * t238;
t198 = mrSges(4,2) * t181 - mrSges(4,3) * t240;
t117 = t168 * t234 - t169 * t185;
t115 = t168 * t235 + t169 * t188;
t127 = -t178 * t228 - t180 * t227;
t128 = -t178 * t227 + t180 * t228;
t30 = t184 * t127 + t187 * t128 + t80 * t225 - t226 * t89;
t31 = -qJD(4) * t42 + t187 * t127 - t128 * t184;
t167 = -qJDD(1) * pkin(1) + qJDD(2);
t143 = t252 + t254;
t134 = t199 * qJD(1);
t133 = t198 * qJD(1);
t130 = t199 * qJDD(1);
t129 = t198 * qJDD(1);
t118 = -t168 * t185 - t169 * t234;
t116 = -t168 * t188 + t169 * t235;
t112 = t179 * t267;
t111 = t179 * t132;
t100 = -qJ(2) * t239 + t136;
t99 = Ifges(5,4) * t107;
t90 = -pkin(4) * t112 + t229;
t87 = -t178 * t214 + t114;
t86 = mrSges(5,1) * t161 - t247;
t85 = -mrSges(5,2) * t161 - t248;
t82 = pkin(4) * t122 + t141;
t75 = pkin(4) * t107 + t126;
t54 = Ifges(6,4) * t206;
t52 = Ifges(5,1) * t109 + Ifges(5,5) * t161 - t99;
t51 = -Ifges(5,2) * t107 + Ifges(5,6) * t161 + t246;
t50 = -mrSges(5,2) * t159 + mrSges(5,3) * t69;
t49 = mrSges(5,1) * t159 - mrSges(5,3) * t68;
t45 = mrSges(6,1) * t155 - mrSges(6,3) * t60;
t44 = -mrSges(6,2) * t155 + mrSges(6,3) * t206;
t43 = -pkin(4) * t69 + t102;
t34 = -qJD(5) * t72 + t111 * t183 + t112 * t186;
t33 = qJD(5) * t71 - t111 * t186 + t112 * t183;
t32 = -mrSges(6,1) * t206 + mrSges(6,2) * t60;
t27 = Ifges(6,1) * t60 + Ifges(6,5) * t155 + t54;
t26 = Ifges(6,2) * t206 + Ifges(6,6) * t155 + t255;
t24 = pkin(7) * t111 + t31;
t23 = pkin(7) * t112 + t30;
t17 = -mrSges(6,2) * t152 + mrSges(6,3) * t22;
t16 = mrSges(6,1) * t152 - mrSges(6,3) * t21;
t11 = t186 * t28 - t243;
t10 = -t183 * t28 - t242;
t5 = -qJD(5) * t15 - t183 * t23 + t186 * t24;
t4 = qJD(5) * t14 + t183 * t24 + t186 * t23;
t1 = [t34 * t261 + (Ifges(6,1) * t33 + Ifges(6,4) * t34) * t263 + t233 * t170 + t138 * t197 + t82 * t210 + t141 * t211 + 0.2e1 * t232 * t158 * mrSges(3,3) + m(6) * (t14 * t3 + t15 * t2 + t4 * t9 + t43 * t82 + t5 * t8 + t75 * t90) + (-mrSges(6,1) * t43 + mrSges(6,3) * t2 + Ifges(6,4) * t21 + Ifges(6,2) * t22 + Ifges(6,6) * t152) * t71 + m(4) * (t100 * t73 + t101 * t74 + t127 * t87 + t128 * t88) + t167 * t204 - pkin(1) * t205 + (-Ifges(5,1) * t111 + Ifges(5,4) * t112) * t259 + (t111 * t37 + t112 * t38) * mrSges(5,3) - t107 * (-Ifges(5,4) * t111 + Ifges(5,2) * t112) / 0.2e1 + t126 * (-mrSges(5,1) * t112 - mrSges(5,2) * t111) + t161 * (-Ifges(5,5) * t111 + Ifges(5,6) * t112) / 0.2e1 + m(5) * (t102 * t141 + t12 * t42 + t13 * t41 + t30 * t38 + t31 * t37) + t206 * (Ifges(6,4) * t33 + Ifges(6,2) * t34) / 0.2e1 - t33 * t262 + (mrSges(5,1) * t102 - mrSges(5,3) * t12 - Ifges(5,4) * t68 - Ifges(5,2) * t69 - Ifges(5,6) * t159) * t122 + t74 * t198 + t73 * t199 + (Ifges(3,1) * t221 - (Ifges(4,4) * t180 - Ifges(4,2) * t178) * t209 + (Ifges(4,1) * t180 - Ifges(4,4) * t178) * t208 + m(4) * (qJ(2) * t138 + qJD(2) * t156) + (-Ifges(4,5) * t180 + Ifges(4,6) * t178 + Ifges(3,4)) * t220) * t179 - (t218 + t217) * t181 / 0.2e1 + (mrSges(6,2) * t43 - mrSges(6,3) * t3 + Ifges(6,1) * t21 + Ifges(6,4) * t22 + Ifges(6,5) * t152) * t72 + (-mrSges(5,2) * t102 + mrSges(5,3) * t13 - Ifges(5,1) * t68 - Ifges(5,4) * t69 - Ifges(5,5) * t159) * t123 + m(3) * (-pkin(1) * t167 + (t158 + t223) * qJ(2) * t232) + (-t118 * mrSges(5,1) - t106 * mrSges(6,1) - t117 * mrSges(5,2) - t105 * mrSges(6,2) + t274 * t188 + (m(4) * qJ(2) - t270 - m(5) * (-qJ(2) - t254) - m(6) * (-qJ(2) - t143) + t278) * t185) * g(2) + (m(5) * t126 - t279) * t229 + (Ifges(3,4) * t221 + Ifges(4,6) * t209 - Ifges(4,5) * t208 - t285 / 0.2e1 - t283 / 0.2e1 - t286 / 0.2e1 - t284 / 0.2e1 - t280 / 0.2e1 - t281 / 0.2e1 + (Ifges(3,2) + Ifges(4,3)) * t220 + t275 + t276) * t181 + t14 * t16 + t15 * t17 + t33 * t27 / 0.2e1 + t34 * t26 / 0.2e1 + t4 * t44 + t5 * t45 + t41 * t49 + t42 * t50 + t75 * (-mrSges(6,1) * t34 + mrSges(6,2) * t33) + t30 * t85 + t31 * t86 + t90 * t32 + (t116 * mrSges(5,1) + t104 * mrSges(6,1) - t115 * mrSges(5,2) - t103 * mrSges(6,2) + t274 * t185 + (-m(5) * t254 - m(6) * t143 + qJ(2) * t282 - mrSges(3,3) - t278) * t188) * g(3) - t111 * t52 / 0.2e1 + t112 * t51 / 0.2e1 + t101 * t129 + t100 * t130 + t128 * t133 + t127 * t134 + t155 * (Ifges(6,5) * t33 + Ifges(6,6) * t34) / 0.2e1 + Ifges(2,3) * qJDD(1); m(3) * t167 + t269 * t86 + t268 * t85 + t205 + t272 * t45 + t273 * t44 + t83 * t16 + t84 * t17 - t201 * t49 + t140 * t50 + t178 * t129 + t180 * t130 + (g(2) * t188 + g(3) * t185) * t282 + (t2 * t84 - t75 * t231 + t272 * t8 + t273 * t9 + t3 * t83) * m(6) + (t12 * t140 - t126 * t231 - t13 * t201 + t268 * t38 + t269 * t37) * m(5) + m(4) * (t178 * t74 + t180 * t73) + ((-t133 * t180 + t134 * t178) * t181 + (-t32 + t279) * t179 + (-t156 * t179 - (-t178 * t87 + t180 * t88) * t181) * m(4) + t270 * t232 * qJD(1)) * qJD(1); t107 * t85 + t109 * t86 - t206 * t44 + t60 * t45 + t224 * t181 * g(1) + m(4) * t138 + ((t178 * t133 + t180 * t134 - m(4) * (-t178 * t88 - t180 * t87)) * qJD(1) + (g(2) * t185 - g(3) * t188) * t224) * t179 + t210 + t211 + t233 + (-t206 * t9 + t60 * t8 + t43) * m(6) + (t107 * t38 + t109 * t37 + t102) * m(5); (t247 + t86) * t38 + t51 * t259 + (t183 * t2 + t186 * t3 + (-t183 * t8 + t186 * t9) * qJD(5)) * t266 - t275 + ((-t183 * t45 + t186 * t44) * qJD(5) + t186 * t16 + t183 * t17) * pkin(4) + (-mrSges(5,2) * t116 - t115 * t271 - t250) * g(2) + t217 - t109 * (-Ifges(5,1) * t107 - t246) / 0.2e1 - t126 * (mrSges(5,1) * t109 - mrSges(5,2) * t107) - t161 * (-Ifges(5,5) * t107 - Ifges(5,6) * t109) / 0.2e1 + (-t248 - t85) * t37 + (-Ifges(5,2) * t109 + t52 - t99) * t107 / 0.2e1 + (Ifges(6,5) * t258 + t262 + Ifges(6,1) * t264 + Ifges(6,4) * t265 - t27 / 0.2e1 - t75 * mrSges(6,2)) * t206 + (-mrSges(5,2) * t118 + t117 * t271 - t249) * g(3) - t32 * t253 - m(6) * (t10 * t8 + t11 * t9 + t253 * t75) + t200 + (m(6) * t252 + mrSges(5,1) * t168 + mrSges(5,2) * t169 - t202) * t251 - (Ifges(6,6) * t258 + Ifges(6,4) * t264 + Ifges(6,2) * t265 - t261 - t26 / 0.2e1 + t75 * mrSges(6,1)) * t60 - t11 * t44 - t10 * t45; -t75 * (mrSges(6,1) * t60 + mrSges(6,2) * t206) + (Ifges(6,1) * t206 - t255) * t264 + t26 * t263 + (Ifges(6,5) * t206 - Ifges(6,6) * t60) * t258 - t8 * t44 + t9 * t45 - g(2) * t250 - g(3) * t249 - t202 * t251 + (t206 * t8 + t60 * t9) * mrSges(6,3) + t200 + (-Ifges(6,2) * t60 + t27 + t54) * t265;];
tau = t1;
