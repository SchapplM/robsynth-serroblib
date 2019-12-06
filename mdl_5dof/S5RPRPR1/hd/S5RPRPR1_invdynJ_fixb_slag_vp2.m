% Calculate vector of inverse dynamics joint torques for
% S5RPRPR1
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR1_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR1_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR1_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR1_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR1_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:47:03
% EndTime: 2019-12-05 17:47:21
% DurationCPUTime: 7.06s
% Computational Cost: add. (3345->377), mult. (6549->507), div. (0->0), fcn. (4232->12), ass. (0->171)
t166 = sin(qJ(3));
t271 = -t166 / 0.2e1;
t169 = cos(qJ(3));
t239 = t169 / 0.2e1;
t162 = sin(pkin(8));
t163 = cos(pkin(8));
t205 = qJD(3) * t169;
t206 = qJD(3) * t166;
t106 = t162 * t206 - t163 * t205;
t113 = -t162 * t169 - t163 * t166;
t107 = t113 * qJD(3);
t165 = sin(qJ(5));
t168 = cos(qJ(5));
t114 = -t162 * t166 + t163 * t169;
t63 = t113 * t168 - t114 * t165;
t173 = qJD(5) * t63 + t106 * t165 + t168 * t107;
t256 = t113 * t165 + t168 * t114;
t26 = qJD(5) * t256 - t106 * t168 + t107 * t165;
t207 = qJD(1) * t169;
t208 = qJD(1) * t166;
t104 = t162 * t208 - t163 * t207;
t235 = pkin(7) * t104;
t171 = -pkin(1) - pkin(6);
t133 = qJD(1) * t171 + qJD(2);
t96 = -qJ(4) * t208 + t133 * t166;
t84 = t162 * t96;
t97 = -qJ(4) * t207 + t169 * t133;
t88 = qJD(3) * pkin(3) + t97;
t43 = t163 * t88 - t84;
t30 = qJD(3) * pkin(4) + t235 + t43;
t105 = t113 * qJD(1);
t234 = pkin(7) * t105;
t224 = t163 * t96;
t44 = t162 * t88 + t224;
t31 = t44 + t234;
t6 = -t165 * t31 + t168 * t30;
t7 = t165 * t30 + t168 * t31;
t273 = -t173 * t6 - t26 * t7;
t157 = qJDD(3) + qJDD(5);
t159 = qJD(3) + qJD(5);
t193 = t104 * t165 + t168 * t105;
t201 = qJD(1) * qJD(3);
t119 = qJDD(1) * t169 - t166 * t201;
t120 = -qJDD(1) * t166 - t169 * t201;
t67 = -t119 * t162 + t120 * t163;
t68 = t119 * t163 + t120 * t162;
t17 = qJD(5) * t193 + t165 * t67 + t168 * t68;
t55 = t104 * t168 - t105 * t165;
t18 = qJD(5) * t55 - t165 * t68 + t168 * t67;
t200 = qJD(1) * qJD(4);
t132 = qJDD(1) * t171 + qJDD(2);
t73 = t169 * t132 - t133 * t206;
t40 = qJDD(3) * pkin(3) - qJ(4) * t119 - t169 * t200 + t73;
t74 = t166 * t132 + t133 * t205;
t45 = qJ(4) * t120 - t166 * t200 + t74;
t19 = -t162 * t45 + t163 * t40;
t8 = qJDD(3) * pkin(4) - pkin(7) * t68 + t19;
t20 = t162 * t40 + t163 * t45;
t9 = pkin(7) * t67 + t20;
t2 = qJD(5) * t6 + t165 * t8 + t168 * t9;
t238 = Ifges(6,4) * t55;
t23 = Ifges(6,2) * t193 + Ifges(6,6) * t159 - t238;
t52 = Ifges(6,4) * t193;
t24 = -Ifges(6,1) * t55 + Ifges(6,5) * t159 + t52;
t244 = -t55 / 0.2e1;
t3 = -qJD(5) * t7 - t165 * t9 + t168 * t8;
t124 = pkin(3) * t208 + qJD(1) * qJ(2) + qJD(4);
t71 = -pkin(4) * t105 + t124;
t272 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t17 + Ifges(6,6) * t18 + Ifges(6,3) * t157 - (Ifges(6,5) * t193 + Ifges(6,6) * t55) * t159 / 0.2e1 + (t193 * t6 - t55 * t7) * mrSges(6,3) - (Ifges(6,2) * t55 + t24 + t52) * t193 / 0.2e1 - t71 * (-mrSges(6,1) * t55 + mrSges(6,2) * t193) + t23 * t244 + (Ifges(6,1) * t193 + t238) * t55 / 0.2e1;
t160 = qJ(3) + pkin(8);
t149 = sin(t160);
t150 = cos(t160);
t155 = t166 * pkin(3);
t151 = qJ(5) + t160;
t140 = sin(t151);
t141 = cos(t151);
t188 = -t140 * mrSges(6,1) - t141 * mrSges(6,2);
t189 = t166 * mrSges(4,1) + t169 * mrSges(4,2);
t270 = t189 + m(5) * t155 + t149 * mrSges(5,1) + t150 * mrSges(5,2) + m(6) * (pkin(4) * t149 + t155) - t188;
t243 = -m(3) - m(4);
t254 = -m(5) - m(6) + t243;
t181 = t169 * (-qJD(3) * mrSges(4,2) - mrSges(4,3) * t208) - t166 * (qJD(3) * mrSges(4,1) - mrSges(4,3) * t207);
t269 = t181 * qJD(3);
t142 = pkin(3) * t163 + pkin(4);
t237 = pkin(3) * t162;
t101 = t142 * t168 - t165 * t237;
t50 = -t162 * t97 - t224;
t34 = t50 - t234;
t51 = t163 * t97 - t84;
t35 = t51 + t235;
t264 = t101 * qJD(5) - t165 * t34 - t168 * t35;
t102 = t142 * t165 + t168 * t237;
t263 = -t102 * qJD(5) + t165 * t35 - t168 * t34;
t202 = qJD(1) * qJD(2);
t134 = qJDD(1) * qJ(2) + t202;
t184 = t74 * t166 + t73 * t169;
t229 = mrSges(4,2) * t166;
t236 = pkin(3) * t169;
t247 = m(5) * pkin(3);
t253 = -t169 * (mrSges(4,1) + t247) - m(6) * (pkin(4) * t150 + t236) - mrSges(5,1) * t150 + mrSges(5,2) * t149 + t229;
t252 = t44 * t106 - t43 * t107 + t20 * t113 - t19 * t114;
t226 = Ifges(4,4) * t169;
t227 = Ifges(4,4) * t166;
t251 = qJ(2) * (mrSges(4,1) * t169 - t229) + (-Ifges(4,2) * t169 - t227) * t271 + (-Ifges(4,1) * t166 - t226) * t239;
t250 = mrSges(2,1) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3) - mrSges(3,2);
t249 = m(4) * t184 + t166 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t120) + t169 * (qJDD(3) * mrSges(4,1) - mrSges(4,3) * t119);
t248 = mrSges(2,2) - mrSges(3,3) - t270;
t242 = -t104 / 0.2e1;
t167 = sin(qJ(1));
t233 = g(1) * t167;
t170 = cos(qJ(1));
t232 = g(2) * t170;
t231 = -qJD(3) / 0.2e1;
t164 = -qJ(4) - pkin(6);
t210 = qJ(4) - t171;
t94 = -qJD(4) * t169 + t206 * t210;
t126 = t210 * t169;
t95 = -qJD(3) * t126 - qJD(4) * t166;
t47 = t162 * t94 + t163 * t95;
t230 = mrSges(6,1) * t141;
t228 = mrSges(6,2) * t140;
t225 = Ifges(5,4) * t104;
t143 = qJ(2) + t155;
t125 = t210 * t166;
t70 = -t163 * t125 - t162 * t126;
t203 = qJDD(1) * mrSges(3,2);
t135 = pkin(3) * t205 + qJD(2);
t61 = -mrSges(5,1) * t105 - mrSges(5,2) * t104;
t196 = -m(5) * t124 - t61;
t195 = -t67 * mrSges(5,1) + t68 * mrSges(5,2);
t194 = -t18 * mrSges(6,1) + t17 * mrSges(6,2);
t46 = -t162 * t95 + t163 * t94;
t192 = (t134 + t202) * qJ(2);
t69 = t125 * t162 - t163 * t126;
t190 = -g(1) * t170 - g(2) * t167;
t187 = t169 * Ifges(4,1) - t227;
t186 = -t166 * Ifges(4,2) + t226;
t185 = -Ifges(4,5) * t166 - Ifges(4,6) * t169;
t48 = -pkin(7) * t114 + t69;
t49 = pkin(7) * t113 + t70;
t21 = -t165 * t49 + t168 * t48;
t22 = t165 * t48 + t168 * t49;
t78 = -pkin(3) * t120 + qJDD(4) + t134;
t158 = -pkin(7) + t164;
t148 = -qJDD(1) * pkin(1) + qJDD(2);
t128 = t170 * t228;
t127 = t167 * t230;
t116 = t189 * qJD(1);
t109 = Ifges(4,5) * qJD(3) + qJD(1) * t187;
t108 = Ifges(4,6) * qJD(3) + qJD(1) * t186;
t98 = Ifges(5,4) * t105;
t89 = -pkin(4) * t113 + t143;
t83 = qJD(3) * mrSges(5,1) + mrSges(5,3) * t104;
t82 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t105;
t77 = pkin(3) * t207 - pkin(4) * t104;
t72 = -pkin(4) * t106 + t135;
t60 = qJDD(3) * mrSges(5,1) - mrSges(5,3) * t68;
t59 = -qJDD(3) * mrSges(5,2) + mrSges(5,3) * t67;
t54 = -t104 * Ifges(5,1) + Ifges(5,5) * qJD(3) + t98;
t53 = t105 * Ifges(5,2) + Ifges(5,6) * qJD(3) - t225;
t42 = mrSges(6,1) * t159 + mrSges(6,3) * t55;
t41 = -mrSges(6,2) * t159 + mrSges(6,3) * t193;
t36 = -pkin(4) * t67 + t78;
t33 = pkin(7) * t106 + t47;
t32 = -pkin(7) * t107 + t46;
t25 = -mrSges(6,1) * t193 - mrSges(6,2) * t55;
t13 = -mrSges(6,2) * t157 + mrSges(6,3) * t18;
t12 = mrSges(6,1) * t157 - mrSges(6,3) * t17;
t5 = -qJD(5) * t22 - t165 * t33 + t168 * t32;
t4 = qJD(5) * t21 + t165 * t32 + t168 * t33;
t1 = [(t254 * (t170 * pkin(1) + t167 * qJ(2)) + (-m(4) * pkin(6) + m(5) * t164 + m(6) * t158 - t250) * t170 + t248 * t167) * g(2) - t184 * mrSges(4,3) + t251 * t201 + t252 * mrSges(5,3) + t273 * mrSges(6,3) + t148 * mrSges(3,2) + t135 * t61 + t78 * (-mrSges(5,1) * t113 + mrSges(5,2) * t114) + t67 * (Ifges(5,4) * t114 + Ifges(5,2) * t113) + qJD(2) * t116 + qJ(2) * (-mrSges(4,1) * t120 + mrSges(4,2) * t119) + t124 * (-mrSges(5,1) * t106 + mrSges(5,2) * t107) + t107 * t54 / 0.2e1 + t68 * (Ifges(5,1) * t114 + Ifges(5,4) * t113) + t106 * t53 / 0.2e1 + t105 * (Ifges(5,4) * t107 + Ifges(5,2) * t106) / 0.2e1 + qJD(3) * (Ifges(5,5) * t107 + Ifges(5,6) * t106) / 0.2e1 + t47 * t82 + t46 * t83 + t69 * t60 + t70 * t59 + t72 * t25 + (0.2e1 * mrSges(3,3) + t189) * t134 + t173 * t24 / 0.2e1 + (Ifges(4,4) * t119 + Ifges(4,2) * t120) * t271 + t4 * t41 + t5 * t42 + t22 * t13 + t21 * t12 + (0.2e1 * Ifges(4,5) * t239 + Ifges(5,5) * t114 + 0.2e1 * Ifges(4,6) * t271 + Ifges(5,6) * t113) * qJDD(3) + t249 * t171 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (mrSges(6,2) * t36 - mrSges(6,3) * t3 + Ifges(6,1) * t17 + Ifges(6,4) * t18 + Ifges(6,5) * t157) * t256 + (Ifges(4,1) * t119 + Ifges(4,4) * t120) * t239 + t171 * t269 + m(5) * (t124 * t135 + t143 * t78 + t19 * t69 + t20 * t70 + t43 * t46 + t44 * t47) + m(6) * (t2 * t22 + t21 * t3 + t36 * t89 + t4 * t7 + t5 * t6 + t71 * t72) + m(4) * t192 + ((m(3) * pkin(1) - m(5) * (-pkin(1) + t164) - m(6) * (-pkin(1) + t158) - m(4) * t171 + t250) * t167 + (t254 * qJ(2) + t248) * t170) * g(1) + (Ifges(5,1) * t107 + Ifges(5,4) * t106) * t242 + (-mrSges(6,1) * t36 + mrSges(6,3) * t2 + Ifges(6,4) * t17 + Ifges(6,2) * t18 + Ifges(6,6) * t157) * t63 - t26 * t23 / 0.2e1 + t159 * (Ifges(6,5) * t173 - Ifges(6,6) * t26) / 0.2e1 + t71 * (mrSges(6,1) * t26 + mrSges(6,2) * t173) + (Ifges(6,1) * t173 - Ifges(6,4) * t26) * t244 + t193 * (Ifges(6,4) * t173 - Ifges(6,2) * t26) / 0.2e1 + qJD(3) ^ 2 * t185 / 0.2e1 + t120 * t186 / 0.2e1 + t119 * t187 / 0.2e1 + m(3) * (-pkin(1) * t148 + t192) + t89 * t194 + t143 * t195 - pkin(1) * t203 - t108 * t205 / 0.2e1 - t109 * t206 / 0.2e1; t203 - t106 * t82 + t107 * t83 - t113 * t59 + t114 * t60 + t256 * t12 - t63 * t13 + t26 * t41 + t173 * t42 + t269 + m(3) * t148 - m(5) * t252 + m(6) * (-t2 * t63 + t256 * t3 - t273) - (t232 - t233) * t254 + (-m(6) * t71 - t116 + t196 - t25 + (qJ(2) * t243 - mrSges(3,3)) * qJD(1)) * qJD(1) + t249; -g(1) * t127 - g(2) * t128 + (t185 * t231 + t108 * t239 + t166 * t109 / 0.2e1 - t251 * qJD(1) + t196 * t236) * qJD(1) + (t230 - t253) * t232 + (t228 + t253) * t233 + Ifges(4,5) * t119 + Ifges(4,6) * t120 - t124 * (-mrSges(5,1) * t104 + mrSges(5,2) * t105) + t101 * t12 + t102 * t13 - t51 * t82 - t50 * t83 + Ifges(5,5) * t68 + t73 * mrSges(4,1) - t74 * mrSges(4,2) - t77 * t25 + Ifges(5,6) * t67 - (Ifges(5,2) * t104 + t54 + t98) * t105 / 0.2e1 + t19 * mrSges(5,1) - t20 * mrSges(5,2) + (t162 * t59 + t163 * t60) * pkin(3) - m(5) * (t43 * t50 + t44 * t51) + t272 + t270 * g(3) + (Ifges(5,3) + Ifges(4,3)) * qJDD(3) + (Ifges(5,5) * t105 + Ifges(5,6) * t104) * t231 + t263 * t42 + t264 * t41 + (t101 * t3 + t102 * t2 + t263 * t6 + t264 * t7 - t71 * t77) * m(6) + (-t104 * t44 + t105 * t43) * mrSges(5,3) + t53 * t242 + (t162 * t20 + t163 * t19) * t247 - t181 * t133 + t104 * (Ifges(5,1) * t105 + t225) / 0.2e1; -t104 * t83 - t105 * t82 - t193 * t41 - t55 * t42 + t194 + t195 + (-t193 * t7 - t55 * t6 + t190 + t36) * m(6) + (-t104 * t43 - t105 * t44 + t190 + t78) * m(5); -t6 * t41 + t7 * t42 - g(1) * (-t167 * t228 + t127) - g(2) * (-t170 * t230 + t128) - g(3) * t188 + t272;];
tau = t1;
