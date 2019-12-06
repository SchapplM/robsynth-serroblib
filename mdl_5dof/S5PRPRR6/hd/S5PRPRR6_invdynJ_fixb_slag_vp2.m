% Calculate vector of inverse dynamics joint torques for
% S5PRPRR6
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPRR6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR6_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR6_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR6_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR6_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR6_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:56:28
% EndTime: 2019-12-05 15:56:45
% DurationCPUTime: 8.01s
% Computational Cost: add. (3402->424), mult. (8294->604), div. (0->0), fcn. (6449->14), ass. (0->200)
t147 = sin(qJ(5));
t150 = cos(qJ(5));
t177 = -t150 * mrSges(6,1) + t147 * mrSges(6,2);
t267 = -m(6) * pkin(4) - mrSges(5,1) + t177;
t266 = -m(6) * pkin(8) + mrSges(5,2) - mrSges(6,3);
t141 = sin(pkin(10));
t144 = cos(pkin(10));
t148 = sin(qJ(4));
t151 = cos(qJ(4));
t115 = t141 * t148 - t151 * t144;
t143 = sin(pkin(5));
t152 = cos(qJ(2));
t210 = t143 * t152;
t156 = t115 * t210;
t235 = pkin(7) + qJ(3);
t122 = t235 * t141;
t123 = t235 * t144;
t165 = -t151 * t122 - t123 * t148;
t269 = qJD(1) * t156 - qJD(3) * t115 + qJD(4) * t165;
t110 = t115 * qJD(4);
t116 = t141 * t151 + t144 * t148;
t111 = t116 * qJD(4);
t149 = sin(qJ(2));
t207 = qJD(1) * t143;
t188 = t149 * t207;
t282 = pkin(4) * t111 + pkin(8) * t110 - t188;
t281 = qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t109 = t116 * qJD(2);
t226 = t109 * mrSges(5,3);
t82 = qJD(4) * t150 - t109 * t147;
t83 = qJD(4) * t147 + t109 * t150;
t232 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t82 + mrSges(6,2) * t83 + t226;
t120 = qJD(2) * qJ(3) + t188;
t145 = cos(pkin(5));
t206 = qJD(1) * t145;
t131 = t144 * t206;
t227 = pkin(7) * qJD(2);
t73 = t131 + (-t120 - t227) * t141;
t86 = t144 * t120 + t141 * t206;
t74 = t144 * t227 + t86;
t35 = -t148 * t74 + t151 * t73;
t31 = -qJD(4) * pkin(4) - t35;
t280 = -m(6) * t31 - t232;
t278 = t141 ^ 2 + t144 ^ 2;
t140 = pkin(10) + qJ(4);
t136 = sin(t140);
t137 = cos(t140);
t277 = t266 * t136 + t267 * t137;
t68 = -qJD(2) * t111 - qJDD(2) * t115;
t108 = t115 * qJD(2);
t101 = qJD(5) + t108;
t67 = -qJD(2) * t110 + qJDD(2) * t116;
t29 = qJD(5) * t82 + qJDD(4) * t147 + t150 * t67;
t251 = t29 / 0.2e1;
t30 = -qJD(5) * t83 + qJDD(4) * t150 - t147 * t67;
t250 = t30 / 0.2e1;
t62 = qJDD(5) - t68;
t249 = t62 / 0.2e1;
t276 = -m(6) - m(5);
t135 = pkin(3) * t144 + pkin(2);
t65 = pkin(4) * t115 - pkin(8) * t116 - t135;
t78 = -t122 * t148 + t123 * t151;
t27 = -t147 * t78 + t150 * t65;
t275 = qJD(5) * t27 + t282 * t147 + t150 * t269;
t36 = t148 * t73 + t151 * t74;
t32 = qJD(4) * pkin(8) + t36;
t172 = -t152 * t207 + qJD(3);
t97 = -qJD(2) * t135 + t172;
t37 = pkin(4) * t108 - pkin(8) * t109 + t97;
t14 = -t147 * t32 + t150 * t37;
t274 = t14 * mrSges(6,1);
t15 = t147 * t37 + t150 * t32;
t273 = t15 * mrSges(6,2);
t28 = t147 * t65 + t150 * t78;
t272 = -qJD(5) * t28 - t147 * t269 + t282 * t150;
t11 = -mrSges(6,1) * t30 + mrSges(6,2) * t29;
t271 = -qJDD(4) * mrSges(5,1) + mrSges(5,3) * t67 + t11;
t268 = mrSges(4,3) * t278;
t176 = t147 * mrSges(6,1) + t150 * mrSges(6,2);
t265 = -mrSges(5,3) - t176;
t196 = qJDD(2) * t144;
t197 = qJDD(2) * t141;
t114 = -mrSges(4,1) * t196 + mrSges(4,2) * t197;
t21 = -t68 * mrSges(5,1) + t67 * mrSges(5,2);
t264 = t114 + t21;
t179 = -mrSges(4,1) * t144 + mrSges(4,2) * t141;
t263 = mrSges(5,1) * t108 + mrSges(5,2) * t109 + t179 * qJD(2);
t170 = t144 * t86 - t141 * (-t120 * t141 + t131);
t262 = t152 * t170;
t202 = qJD(5) * t150;
t161 = -t147 * t110 + t116 * t202;
t199 = qJDD(1) * t145;
t129 = t144 * t199;
t184 = qJD(2) * t207;
t127 = t152 * t184;
t200 = qJDD(1) * t143;
t100 = t149 * t200 + t127;
t84 = t100 + t281;
t58 = -t141 * t84 + t129;
t59 = t141 * t199 + t144 * t84;
t261 = -t141 * t58 + t144 * t59;
t126 = t149 * t184;
t99 = t152 * t200 - t126;
t163 = qJDD(3) - t99;
t75 = -qJDD(2) * t135 + t163;
t20 = -pkin(4) * t68 - pkin(8) * t67 + t75;
t47 = t129 + (-pkin(7) * qJDD(2) - t84) * t141;
t48 = pkin(7) * t196 + t59;
t7 = qJD(4) * t35 + t148 * t47 + t151 * t48;
t5 = qJDD(4) * pkin(8) + t7;
t1 = qJD(5) * t14 + t147 * t20 + t150 * t5;
t2 = -qJD(5) * t15 - t147 * t5 + t150 * t20;
t260 = t1 * t150 - t147 * t2;
t8 = -qJD(4) * t36 - t148 * t48 + t151 * t47;
t185 = m(4) * qJ(3) + mrSges(4,3);
t256 = mrSges(3,2) - t185 + t265;
t159 = m(4) * pkin(2) - t179;
t255 = mrSges(3,1) + t159 - t277;
t254 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t253 = Ifges(6,1) * t251 + Ifges(6,4) * t250 + Ifges(6,5) * t249;
t248 = -t82 / 0.2e1;
t247 = -t83 / 0.2e1;
t246 = t83 / 0.2e1;
t245 = -t101 / 0.2e1;
t242 = t109 / 0.2e1;
t241 = t150 / 0.2e1;
t240 = Ifges(6,4) * t83;
t239 = g(3) * t143;
t231 = mrSges(5,3) * t108;
t230 = Ifges(5,4) * t109;
t229 = Ifges(6,4) * t147;
t228 = Ifges(6,4) * t150;
t219 = cos(pkin(9));
t218 = t108 * t147;
t217 = t108 * t150;
t216 = t116 * t147;
t215 = t116 * t150;
t142 = sin(pkin(9));
t214 = t142 * t143;
t213 = t142 * t149;
t212 = t142 * t152;
t211 = t143 * t149;
t205 = qJD(2) * t149;
t203 = qJD(5) * t147;
t195 = Ifges(6,5) * t29 + Ifges(6,6) * t30 + Ifges(6,3) * t62;
t193 = t147 * t210;
t192 = t150 * t210;
t81 = Ifges(6,4) * t82;
t24 = t83 * Ifges(6,1) + t101 * Ifges(6,5) + t81;
t191 = t24 * t241;
t187 = t143 * t205;
t183 = -t203 / 0.2e1;
t182 = t143 * t219;
t181 = t219 * t149;
t180 = t219 * t152;
t175 = Ifges(6,1) * t150 - t229;
t174 = -Ifges(6,2) * t147 + t228;
t173 = Ifges(6,5) * t150 - Ifges(6,6) * t147;
t43 = -mrSges(6,2) * t101 + mrSges(6,3) * t82;
t44 = mrSges(6,1) * t101 - mrSges(6,3) * t83;
t169 = -t147 * t44 + t150 * t43;
t102 = -t141 * t211 + t144 * t145;
t103 = t141 * t145 + t144 * t211;
t166 = t151 * t102 - t103 * t148;
t52 = t102 * t148 + t103 * t151;
t39 = -t147 * t52 - t192;
t164 = -t150 * t52 + t193;
t162 = t31 * t176;
t160 = t150 * t110 + t116 * t203;
t157 = t116 * t210;
t104 = -t145 * t180 + t213;
t106 = t145 * t212 + t181;
t155 = -g(1) * t106 - g(2) * t104 + g(3) * t210;
t153 = qJD(2) ^ 2;
t118 = -qJD(2) * pkin(2) + t172;
t107 = -t145 * t213 + t180;
t105 = t145 * t181 + t212;
t98 = Ifges(5,4) * t108;
t93 = t136 * t145 + t137 * t211;
t90 = -qJD(4) * mrSges(5,2) - t231;
t89 = -qJDD(2) * pkin(2) + t163;
t72 = t107 * t137 + t136 * t214;
t70 = t105 * t137 - t136 * t182;
t63 = pkin(4) * t109 + pkin(8) * t108;
t56 = t109 * Ifges(5,1) + Ifges(5,5) * qJD(4) - t98;
t55 = -t108 * Ifges(5,2) + Ifges(5,6) * qJD(4) + t230;
t54 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t68;
t34 = qJD(2) * t157 + qJD(4) * t52;
t33 = -qJD(2) * t156 + qJD(4) * t166;
t23 = t82 * Ifges(6,2) + t101 * Ifges(6,6) + t240;
t22 = t83 * Ifges(6,5) + t82 * Ifges(6,6) + t101 * Ifges(6,3);
t19 = -mrSges(6,2) * t62 + mrSges(6,3) * t30;
t18 = mrSges(6,1) * t62 - mrSges(6,3) * t29;
t17 = t147 * t63 + t150 * t35;
t16 = -t147 * t35 + t150 * t63;
t13 = qJD(5) * t164 - t147 * t33 + t150 * t187;
t12 = qJD(5) * t39 + t147 * t187 + t150 * t33;
t6 = -qJDD(4) * pkin(4) - t8;
t3 = t29 * Ifges(6,4) + t30 * Ifges(6,2) + t62 * Ifges(6,6);
t4 = [t12 * t43 + t13 * t44 + t39 * t18 - t164 * t19 + t33 * t90 + t52 * t54 - t271 * t166 + t232 * t34 + (-t102 * t141 + t103 * t144) * qJDD(2) * mrSges(4,3) + (-m(2) - m(3) - m(4) + t276) * g(3) + m(5) * (t166 * t8 + t33 * t36 - t34 * t35 + t52 * t7) + m(4) * (t102 * t58 + t103 * t59) + m(6) * (-t1 * t164 + t12 * t15 + t13 * t14 - t166 * t6 + t2 * t39 + t31 * t34) + ((mrSges(3,1) * qJDD(2) - t264) * t152 + (-qJDD(2) * mrSges(3,2) + t263 * qJD(2)) * t149 + m(5) * (-t152 * t75 + t205 * t97) + m(4) * (qJD(2) * t262 + t118 * t205 - t152 * t89) + m(3) * (t100 * t149 + t152 * t99) + (-t149 * mrSges(3,1) + (-mrSges(3,2) + t268) * t152) * t153) * t143 + (m(3) * t145 ^ 2 + m(2)) * qJDD(1); ((t277 * t152 + (t235 * t276 + t265) * t149) * t143 + t276 * t135 * t210) * g(3) - t111 * t273 + t275 * t43 + t269 * t90 + t272 * t44 - t263 * t188 - t161 * t23 / 0.2e1 + (-pkin(2) * t89 + t170 * qJD(3) + t261 * qJ(3) - (t118 * t149 + t262) * t207) * m(4) + (t195 / 0.2e1 - t7 * mrSges(5,3) - Ifges(5,4) * t67 - Ifges(5,2) * t68 - Ifges(5,6) * qJDD(4) + t75 * mrSges(5,1) + Ifges(6,3) * t249 + Ifges(6,6) * t250 + Ifges(6,5) * t251 + t254) * t115 - t135 * t21 - pkin(2) * t114 - t110 * t56 / 0.2e1 + t111 * t22 / 0.2e1 - t111 * t55 / 0.2e1 + t78 * t54 + t27 * t18 + t28 * t19 - t3 * t216 / 0.2e1 + (-t1 * t216 + t14 * t160 - t15 * t161 - t2 * t215) * mrSges(6,3) + t89 * t179 + t31 * (mrSges(6,1) * t161 - mrSges(6,2) * t160) + t82 * (-Ifges(6,4) * t160 - Ifges(6,2) * t161 + Ifges(6,6) * t111) / 0.2e1 + t101 * (-Ifges(6,5) * t160 - Ifges(6,6) * t161 + Ifges(6,3) * t111) / 0.2e1 + (t75 * mrSges(5,2) - t8 * mrSges(5,3) + Ifges(5,1) * t67 + Ifges(5,4) * t68 + Ifges(5,5) * qJDD(4) + t173 * t249 + t174 * t250 + t175 * t251 + t176 * t6 + t183 * t24) * t116 + (-t152 * t239 + t126 + t99) * mrSges(3,1) + (t110 * t35 - t111 * t36) * mrSges(5,3) - t108 * (-Ifges(5,4) * t110 - Ifges(5,2) * t111) / 0.2e1 + t97 * (mrSges(5,1) * t111 - mrSges(5,2) * t110) + qJD(4) * (-Ifges(5,5) * t110 - Ifges(5,6) * t111) / 0.2e1 + (-Ifges(5,1) * t110 - Ifges(5,4) * t111) * t242 - t110 * t191 - t271 * t165 + (t1 * t28 + t272 * t14 + t275 * t15 - t165 * t6 + t2 * t27) * m(6) + (-t135 * t75 + t165 * t8 - t188 * t97 + t269 * t36 + t7 * t78) * m(5) + (t276 * (-t106 * t135 + t107 * t235) + t256 * t107 + t255 * t106) * g(1) + (t276 * (-t104 * t135 + t105 * t235) + t256 * t105 + t255 * t104) * g(2) + t215 * t253 + (-Ifges(6,1) * t160 - Ifges(6,4) * t161 + Ifges(6,5) * t111) * t246 + (t149 * t239 - t100 + t127) * mrSges(3,2) + (t35 * m(5) + t280) * (qJD(1) * t157 - qJD(3) * t116 - qJD(4) * t78) + (t278 * (-t127 + t281) + t261) * mrSges(4,3) + Ifges(3,3) * qJDD(2) + t111 * t274 + (Ifges(4,4) * t141 + Ifges(4,2) * t144) * t196 + (Ifges(4,1) * t141 + Ifges(4,4) * t144) * t197 - (t149 * t185 + t152 * t159) * t239; t147 * t19 + t150 * t18 - t232 * t109 + t169 * qJD(5) - t153 * t268 - (-t169 - t90) * t108 + (t1 * t147 - t109 * t31 + t150 * t2 + t155 + t101 * (-t14 * t147 + t15 * t150)) * m(6) + (t108 * t36 + t109 * t35 + t155 + t75) * m(5) + (-qJD(2) * t170 + t155 + t89) * m(4) + t264; -t109 * t274 + (t266 * t70 + t267 * (-t105 * t136 - t137 * t182)) * g(2) + (t266 * t72 + t267 * (-t107 * t136 + t137 * t214)) * g(1) + (t266 * t93 + t267 * (-t136 * t211 + t137 * t145)) * g(3) + ((-t203 - t218) * t15 + (-t202 - t217) * t14 + t260) * mrSges(6,3) + (t150 * t19 - t147 * t18 - t44 * t202 - t43 * t203 + m(6) * ((-t14 * t150 - t15 * t147) * qJD(5) + t260)) * pkin(8) + (-pkin(4) * t6 - t14 * t16 - t15 * t17) * m(6) + Ifges(5,5) * t67 + Ifges(5,6) * t68 - t17 * t43 - t16 * t44 + t24 * t217 / 0.2e1 + t6 * t177 + t108 * t162 - pkin(4) * t11 - t7 * mrSges(5,2) + t8 * mrSges(5,1) + (-t218 / 0.2e1 + t183) * t23 + (-t90 - t231) * t35 + (-Ifges(5,2) * t109 + t56 - t98) * t108 / 0.2e1 + (t191 + t162) * qJD(5) + (Ifges(6,5) * t147 + Ifges(6,6) * t150) * t249 + (Ifges(6,2) * t150 + t229) * t250 + (Ifges(6,1) * t147 + t228) * t251 + t147 * t253 + t3 * t241 + t55 * t242 + (t226 + t280) * t36 + (t101 * t173 + t174 * t82 + t175 * t83) * qJD(5) / 0.2e1 - t97 * (mrSges(5,1) * t109 - mrSges(5,2) * t108) - qJD(4) * (-Ifges(5,5) * t108 - Ifges(5,6) * t109) / 0.2e1 + (Ifges(6,3) * t109 - t108 * t173) * t245 + (Ifges(6,5) * t109 - t108 * t175) * t247 + (Ifges(6,6) * t109 - t108 * t174) * t248 - (-Ifges(5,1) * t108 + t22 - t230) * t109 / 0.2e1 + Ifges(5,3) * qJDD(4) + t109 * t273; -t31 * (mrSges(6,1) * t83 + mrSges(6,2) * t82) + (Ifges(6,1) * t82 - t240) * t247 + t23 * t246 + (Ifges(6,5) * t82 - Ifges(6,6) * t83) * t245 - t14 * t43 + t15 * t44 - g(1) * ((t106 * t150 - t147 * t72) * mrSges(6,1) + (-t106 * t147 - t150 * t72) * mrSges(6,2)) - g(2) * ((t104 * t150 - t147 * t70) * mrSges(6,1) + (-t104 * t147 - t150 * t70) * mrSges(6,2)) - g(3) * ((-t147 * t93 - t192) * mrSges(6,1) + (-t150 * t93 + t193) * mrSges(6,2)) + (t14 * t82 + t15 * t83) * mrSges(6,3) + t195 + (-Ifges(6,2) * t83 + t24 + t81) * t248 + t254;];
tau = t4;
