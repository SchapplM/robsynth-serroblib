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
% Datum: 2020-01-03 11:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:30:27
% EndTime: 2020-01-03 11:30:57
% DurationCPUTime: 10.02s
% Computational Cost: add. (4551->426), mult. (11556->587), div. (0->0), fcn. (8496->14), ass. (0->209)
t224 = qJD(1) * qJD(2);
t158 = qJDD(1) * qJ(2) + t224;
t180 = sin(pkin(9));
t182 = cos(pkin(9));
t186 = sin(qJ(4));
t189 = cos(qJ(4));
t140 = t180 * t189 + t182 * t186;
t132 = t140 * qJD(4);
t183 = cos(pkin(8));
t197 = qJD(1) * t140;
t270 = t183 * t197 - t132;
t238 = t182 * t189;
t202 = t180 * t186 - t238;
t231 = qJD(1) * t183;
t268 = t202 * qJD(4);
t269 = t202 * t231 - t268;
t181 = sin(pkin(8));
t233 = t181 ^ 2 + t183 ^ 2;
t288 = m(5) + m(6);
t68 = (-qJD(1) * t132 - qJDD(1) * t202) * t181;
t287 = Ifges(5,5) * t68;
t185 = sin(qJ(5));
t188 = cos(qJ(5));
t107 = t181 * t197;
t232 = qJD(1) * t181;
t214 = t180 * t232;
t109 = -t186 * t214 + t232 * t238;
t208 = -t188 * t107 - t109 * t185;
t69 = (qJD(1) * t268 - qJDD(1) * t140) * t181;
t21 = qJD(5) * t208 + t185 * t69 + t188 * t68;
t286 = Ifges(6,5) * t21;
t285 = Ifges(5,6) * t69;
t60 = -t107 * t185 + t109 * t188;
t22 = -qJD(5) * t60 - t185 * t68 + t188 * t69;
t284 = Ifges(6,6) * t22;
t221 = qJDD(1) * t183;
t159 = qJDD(4) - t221;
t283 = Ifges(5,3) * t159;
t152 = qJDD(5) + t159;
t282 = Ifges(6,3) * t152;
t204 = mrSges(4,1) * t180 + mrSges(4,2) * t182;
t198 = t204 * t181;
t281 = -mrSges(5,1) * t107 - mrSges(5,2) * t109 - qJD(1) * t198;
t280 = -m(3) - t288;
t145 = -pkin(2) * t183 - qJ(3) * t181 - pkin(1);
t125 = qJD(1) * t145 + qJD(2);
t114 = t182 * t125;
t240 = t181 * t182;
t220 = pkin(6) * t240;
t195 = -t220 + (-qJ(2) * t180 - pkin(3)) * t183;
t70 = qJD(1) * t195 + t114;
t215 = qJ(2) * t231;
t88 = t180 * t125 + t182 * t215;
t77 = -pkin(6) * t214 + t88;
t38 = t186 * t70 + t189 * t77;
t241 = t180 * t183;
t228 = qJD(3) * t181;
t94 = -qJD(1) * t228 + qJDD(1) * t145 + qJDD(2);
t73 = -t158 * t241 + t182 * t94;
t48 = (-pkin(3) * t183 - t220) * qJDD(1) + t73;
t222 = qJDD(1) * t181;
t210 = t180 * t222;
t239 = t182 * t183;
t74 = t158 * t239 + t180 * t94;
t55 = -pkin(6) * t210 + t74;
t13 = -qJD(4) * t38 - t186 * t55 + t189 * t48;
t6 = pkin(4) * t159 - pkin(7) * t68 + t13;
t226 = qJD(4) * t189;
t227 = qJD(4) * t186;
t12 = t186 * t48 + t189 * t55 + t70 * t226 - t227 * t77;
t7 = pkin(7) * t69 + t12;
t29 = -pkin(7) * t107 + t38;
t245 = t185 * t29;
t161 = qJD(4) - t231;
t37 = -t186 * t77 + t189 * t70;
t28 = -pkin(7) * t109 + t37;
t25 = pkin(4) * t161 + t28;
t8 = t188 * t25 - t245;
t2 = qJD(5) * t8 + t185 * t6 + t188 * t7;
t244 = t188 * t29;
t9 = t185 * t25 + t244;
t3 = -qJD(5) * t9 - t185 * t7 + t188 * t6;
t278 = -t3 * mrSges(6,1) + t2 * mrSges(6,2);
t277 = -t13 * mrSges(5,1) + t12 * mrSges(5,2);
t166 = t182 * pkin(3) + pkin(2);
t179 = pkin(9) + qJ(4);
t169 = cos(t179);
t184 = -pkin(6) - qJ(3);
t205 = mrSges(3,1) * t183 - mrSges(3,2) * t181;
t276 = -t205 - mrSges(2,1) + (-m(4) * pkin(2) - mrSges(4,1) * t182 + mrSges(4,2) * t180 - m(6) * (pkin(4) * t169 + t166) - m(5) * t166) * t183 + (-m(4) * qJ(3) - mrSges(4,3) - m(6) * (pkin(7) - t184) - mrSges(6,3) + m(5) * t184 - mrSges(5,3)) * t181;
t168 = sin(t179);
t254 = pkin(4) * t168;
t256 = pkin(3) * t180;
t275 = -m(5) * t256 - m(6) * (t254 + t256) + mrSges(2,2) - mrSges(3,3) - t204;
t83 = -t140 * t185 - t188 * t202;
t274 = qJD(5) * t83 + t185 * t270 + t269 * t188;
t84 = t140 * t188 - t185 * t202;
t273 = -qJD(5) * t84 - t269 * t185 + t188 * t270;
t267 = m(6) * pkin(4);
t272 = -mrSges(5,1) - t267;
t271 = qJD(1) ^ 2 * t233;
t266 = -t208 / 0.2e1;
t265 = -t60 / 0.2e1;
t264 = t60 / 0.2e1;
t263 = t8 * mrSges(6,3);
t262 = t9 * mrSges(6,3);
t260 = t109 / 0.2e1;
t155 = qJD(5) + t161;
t259 = -t155 / 0.2e1;
t257 = Ifges(6,4) * t60;
t255 = pkin(4) * t109;
t253 = g(1) * t181;
t136 = t182 * t145;
t80 = t136 + t195;
t101 = qJ(2) * t239 + t180 * t145;
t242 = t180 * t181;
t89 = -pkin(6) * t242 + t101;
t42 = t186 * t80 + t189 * t89;
t171 = qJ(5) + t179;
t164 = sin(t171);
t165 = cos(t171);
t190 = cos(qJ(1));
t187 = sin(qJ(1));
t237 = t183 * t187;
t103 = -t164 * t237 - t165 * t190;
t104 = -t164 * t190 + t165 * t237;
t252 = t103 * mrSges(6,1) - t104 * mrSges(6,2);
t236 = t183 * t190;
t105 = t164 * t236 - t165 * t187;
t106 = t164 * t187 + t165 * t236;
t251 = t105 * mrSges(6,1) + t106 * mrSges(6,2);
t250 = mrSges(5,3) * t107;
t249 = mrSges(5,3) * t109;
t248 = Ifges(5,4) * t109;
t170 = t181 * qJ(2);
t209 = t182 * t222;
t235 = mrSges(4,1) * t210 + mrSges(4,2) * t209;
t141 = pkin(3) * t242 + t170;
t230 = qJD(2) * t181;
t229 = qJD(2) * t183;
t156 = qJ(2) * t232 + qJD(3);
t225 = m(4) + t288;
t138 = t181 * t158 + qJDD(3);
t219 = t282 + t284 + t286;
t218 = t283 + t285 + t287;
t126 = pkin(3) * t214 + t156;
t212 = -t69 * mrSges(5,1) + t68 * mrSges(5,2);
t211 = -t22 * mrSges(6,1) + t21 * mrSges(6,2);
t102 = pkin(3) * t210 + t138;
t41 = -t186 * t89 + t189 * t80;
t206 = -mrSges(3,1) * t221 + mrSges(3,2) * t222;
t203 = -mrSges(6,1) * t164 - mrSges(6,2) * t165;
t123 = t202 * t181;
t35 = -pkin(4) * t183 + pkin(7) * t123 + t41;
t122 = t140 * t181;
t36 = -pkin(7) * t122 + t42;
t14 = -t185 * t36 + t188 * t35;
t15 = t185 * t35 + t188 * t36;
t71 = -t122 * t188 + t123 * t185;
t72 = -t122 * t185 - t123 * t188;
t201 = t219 - t278;
t200 = -mrSges(4,1) * t183 - mrSges(4,3) * t240;
t199 = mrSges(4,2) * t183 - mrSges(4,3) * t242;
t117 = t168 * t236 - t169 * t187;
t115 = -t168 * t237 - t169 * t190;
t127 = -t180 * t229 - t182 * t228;
t128 = -t180 * t228 + t182 * t229;
t30 = t186 * t127 + t189 * t128 + t80 * t226 - t227 * t89;
t31 = -qJD(4) * t42 + t189 * t127 - t128 * t186;
t174 = t187 * pkin(1);
t167 = -qJDD(1) * pkin(1) + qJDD(2);
t134 = t200 * qJD(1);
t133 = t199 * qJD(1);
t130 = t200 * qJDD(1);
t129 = t199 * qJDD(1);
t118 = t168 * t187 + t169 * t236;
t116 = -t168 * t190 + t169 * t237;
t112 = t181 * t268;
t111 = t181 * t132;
t100 = -qJ(2) * t241 + t136;
t99 = Ifges(5,4) * t107;
t90 = -pkin(4) * t112 + t230;
t87 = -t180 * t215 + t114;
t86 = mrSges(5,1) * t161 - t249;
t85 = -mrSges(5,2) * t161 - t250;
t82 = pkin(4) * t122 + t141;
t75 = pkin(4) * t107 + t126;
t54 = Ifges(6,4) * t208;
t52 = Ifges(5,1) * t109 + Ifges(5,5) * t161 - t99;
t51 = -Ifges(5,2) * t107 + Ifges(5,6) * t161 + t248;
t50 = -mrSges(5,2) * t159 + mrSges(5,3) * t69;
t49 = mrSges(5,1) * t159 - mrSges(5,3) * t68;
t45 = mrSges(6,1) * t155 - mrSges(6,3) * t60;
t44 = -mrSges(6,2) * t155 + mrSges(6,3) * t208;
t43 = -pkin(4) * t69 + t102;
t34 = -qJD(5) * t72 + t111 * t185 + t112 * t188;
t33 = qJD(5) * t71 - t111 * t188 + t112 * t185;
t32 = -mrSges(6,1) * t208 + mrSges(6,2) * t60;
t27 = Ifges(6,1) * t60 + Ifges(6,5) * t155 + t54;
t26 = Ifges(6,2) * t208 + Ifges(6,6) * t155 + t257;
t24 = pkin(7) * t111 + t31;
t23 = pkin(7) * t112 + t30;
t17 = -mrSges(6,2) * t152 + mrSges(6,3) * t22;
t16 = mrSges(6,1) * t152 - mrSges(6,3) * t21;
t11 = t188 * t28 - t245;
t10 = -t185 * t28 - t244;
t5 = -qJD(5) * t15 - t185 * t23 + t188 * t24;
t4 = qJD(5) * t14 + t185 * t24 + t188 * t23;
t1 = [-pkin(1) * t206 + (Ifges(4,6) * t210 - Ifges(4,5) * t209 + Ifges(3,4) * t222 - t286 / 0.2e1 - t284 / 0.2e1 - t287 / 0.2e1 - t285 / 0.2e1 - t282 / 0.2e1 - t283 / 0.2e1 + (Ifges(4,3) + Ifges(3,2)) * t221 + t277 + t278) * t183 - t167 * t205 + m(5) * (t102 * t141 + t12 * t42 + t13 * t41 + t30 * t38 + t31 * t37) + (-(Ifges(4,4) * t182 - Ifges(4,2) * t180) * t210 + (Ifges(4,1) * t182 - Ifges(4,4) * t180) * t209 + Ifges(3,1) * t222 + m(4) * (qJ(2) * t138 + qJD(2) * t156) + (-Ifges(4,5) * t182 + Ifges(4,6) * t180 + Ifges(3,4)) * t221) * t181 + t82 * t211 + t141 * t212 + m(4) * (t100 * t73 + t101 * t74 + t127 * t87 + t128 * t88) + m(6) * (t14 * t3 + t15 * t2 + t4 * t9 + t43 * t82 + t5 * t8 + t75 * t90) - (t219 + t218) * t183 / 0.2e1 + (t37 * t111 + t38 * t112) * mrSges(5,3) + (-Ifges(5,1) * t111 + Ifges(5,4) * t112) * t260 - t107 * (-Ifges(5,4) * t111 + Ifges(5,2) * t112) / 0.2e1 + t126 * (-mrSges(5,1) * t112 - mrSges(5,2) * t111) + t161 * (-Ifges(5,5) * t111 + Ifges(5,6) * t112) / 0.2e1 + (mrSges(6,2) * t43 - mrSges(6,3) * t3 + Ifges(6,1) * t21 + Ifges(6,4) * t22 + Ifges(6,5) * t152) * t72 + t74 * t199 + t73 * t200 + (m(5) * t126 - t281) * t230 + (-t118 * mrSges(5,1) - t106 * mrSges(6,1) + t117 * mrSges(5,2) + t105 * mrSges(6,2) + (-m(4) + t280) * (t190 * pkin(1) + t187 * qJ(2)) + t275 * t187 + t276 * t190) * g(2) + (-m(4) * t174 - t116 * mrSges(5,1) - t104 * mrSges(6,1) - t115 * mrSges(5,2) - t103 * mrSges(6,2) + t280 * (-qJ(2) * t190 + t174) + (m(4) * qJ(2) - t275) * t190 + t276 * t187) * g(3) + (mrSges(5,1) * t102 - mrSges(5,3) * t12 - Ifges(5,4) * t68 - Ifges(5,2) * t69 - Ifges(5,6) * t159) * t122 + (-mrSges(5,2) * t102 + mrSges(5,3) * t13 - Ifges(5,1) * t68 - Ifges(5,4) * t69 - Ifges(5,5) * t159) * t123 + (-mrSges(6,1) * t43 + mrSges(6,3) * t2 + Ifges(6,4) * t21 + Ifges(6,2) * t22 + Ifges(6,6) * t152) * t71 + 0.2e1 * t233 * t158 * mrSges(3,3) + t14 * t16 + t15 * t17 + t33 * t27 / 0.2e1 + t34 * t26 / 0.2e1 + t4 * t44 + t5 * t45 + t41 * t49 + t42 * t50 + t75 * (-mrSges(6,1) * t34 + mrSges(6,2) * t33) + t30 * t85 + t31 * t86 + t90 * t32 + t138 * t198 + t235 * t170 + t34 * t262 + (Ifges(6,1) * t33 + Ifges(6,4) * t34) * t264 - t111 * t52 / 0.2e1 - t33 * t263 + t112 * t51 / 0.2e1 + t101 * t129 + t100 * t130 + t128 * t133 + t127 * t134 + t155 * (Ifges(6,5) * t33 + Ifges(6,6) * t34) / 0.2e1 + m(3) * (-pkin(1) * t167 + (t158 + t224) * qJ(2) * t233) + t208 * (Ifges(6,4) * t33 + Ifges(6,2) * t34) / 0.2e1 + Ifges(2,3) * qJDD(1); t270 * t86 + t269 * t85 + t273 * t45 + t274 * t44 - mrSges(3,3) * t271 + ((-t133 * t182 + t134 * t180) * t183 + (-t32 + t281) * t181) * qJD(1) + t206 + t83 * t16 + t84 * t17 - t202 * t49 + t140 * t50 + t180 * t129 + t182 * t130 + (g(2) * t190 + g(3) * t187) * (m(3) + t225) + (t2 * t84 - t75 * t232 + t273 * t8 + t274 * t9 + t3 * t83) * m(6) + (t12 * t140 - t126 * t232 - t13 * t202 + t269 * t38 + t270 * t37) * m(5) + (-(t156 * t181 + (-t180 * t87 + t182 * t88) * t183) * qJD(1) + t180 * t74 + t182 * t73) * m(4) + (-qJ(2) * t271 + t167) * m(3); t107 * t85 + t109 * t86 - t208 * t44 + t60 * t45 + t225 * t183 * g(1) + m(4) * t138 + ((-m(4) * (-t180 * t88 - t182 * t87) + t180 * t133 + t182 * t134) * qJD(1) + (-g(2) * t187 + g(3) * t190) * t225) * t181 + t211 + t212 + t235 + (-t208 * t9 + t60 * t8 + t43) * m(6) + (t107 * t38 + t109 * t37 + t102) * m(5); -(-t26 / 0.2e1 + t75 * mrSges(6,1) + Ifges(6,6) * t259 + Ifges(6,4) * t265 + Ifges(6,2) * t266 - t262) * t60 - t109 * (-Ifges(5,1) * t107 - t248) / 0.2e1 - t126 * (mrSges(5,1) * t109 - mrSges(5,2) * t107) - t161 * (-Ifges(5,5) * t107 - Ifges(5,6) * t109) / 0.2e1 + (-t85 - t250) * t37 + (t86 + t249) * t38 + (m(6) * t254 + mrSges(5,1) * t168 + mrSges(5,2) * t169 - t203) * t253 - t277 + (-t118 * mrSges(5,2) + t117 * t272 - t251) * g(3) + (mrSges(5,2) * t116 + t115 * t272 - t252) * g(2) + ((-t185 * t45 + t188 * t44) * qJD(5) + t188 * t16 + t185 * t17) * pkin(4) + t218 + t201 - t11 * t44 - t10 * t45 - t32 * t255 - m(6) * (t10 * t8 + t11 * t9 + t255 * t75) + t51 * t260 + (t185 * t2 + t188 * t3 + (-t185 * t8 + t188 * t9) * qJD(5)) * t267 + (-Ifges(5,2) * t109 + t52 - t99) * t107 / 0.2e1 + (-t27 / 0.2e1 - t75 * mrSges(6,2) + Ifges(6,5) * t259 + t263 + Ifges(6,1) * t265 + Ifges(6,4) * t266) * t208; -t75 * (mrSges(6,1) * t60 + mrSges(6,2) * t208) + (Ifges(6,1) * t208 - t257) * t265 + t26 * t264 + (Ifges(6,5) * t208 - Ifges(6,6) * t60) * t259 - t8 * t44 + t9 * t45 - g(2) * t252 - g(3) * t251 - t203 * t253 + (t208 * t8 + t60 * t9) * mrSges(6,3) + t201 + (-Ifges(6,2) * t60 + t27 + t54) * t266;];
tau = t1;
