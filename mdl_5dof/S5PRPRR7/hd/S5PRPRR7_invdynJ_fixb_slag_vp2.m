% Calculate vector of inverse dynamics joint torques for
% S5PRPRR7
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
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
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPRR7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR7_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR7_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR7_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR7_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR7_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:59:45
% EndTime: 2019-12-05 15:59:55
% DurationCPUTime: 5.19s
% Computational Cost: add. (1752->311), mult. (3316->422), div. (0->0), fcn. (1949->10), ass. (0->150)
t189 = mrSges(3,1) - mrSges(4,2);
t225 = mrSges(5,3) + t189;
t108 = sin(qJ(4));
t111 = cos(qJ(4));
t142 = mrSges(5,1) * t108 + mrSges(5,2) * t111;
t198 = pkin(4) * t108;
t104 = qJ(4) + qJ(5);
t96 = sin(t104);
t97 = cos(t104);
t224 = -m(6) * t198 - mrSges(6,1) * t96 - mrSges(6,2) * t97 + mrSges(3,2) - t142;
t105 = sin(pkin(8));
t106 = cos(pkin(8));
t212 = g(1) * t106 + g(2) * t105;
t203 = -m(4) - m(5);
t153 = m(6) - t203;
t107 = sin(qJ(5));
t110 = cos(qJ(5));
t109 = sin(qJ(2));
t163 = t110 * t111;
t131 = t107 * t108 - t163;
t124 = t131 * t109;
t114 = -pkin(2) - pkin(6);
t187 = pkin(7) - t114;
t77 = t187 * t108;
t78 = t187 * t111;
t37 = -t107 * t78 - t110 * t77;
t159 = qJD(4) * t108;
t67 = t187 * t159;
t68 = qJD(4) * t78;
t219 = qJD(1) * t124 - qJD(5) * t37 + t107 * t68 + t110 * t67;
t132 = t107 * t111 + t108 * t110;
t125 = t132 * t109;
t36 = t107 * t77 - t110 * t78;
t218 = -qJD(1) * t125 + qJD(5) * t36 + t107 * t67 - t110 * t68;
t112 = cos(qJ(2));
t194 = g(3) * t112;
t64 = t132 * qJD(2);
t154 = t111 * qJD(2);
t160 = qJD(2) * t108;
t65 = -t107 * t160 + t110 * t154;
t32 = mrSges(6,1) * t64 + mrSges(6,2) * t65;
t72 = t142 * qJD(2);
t186 = t32 + t72;
t206 = m(6) * pkin(4);
t217 = -mrSges(5,1) - t206;
t55 = t132 * t112;
t151 = qJD(2) * qJD(4);
t73 = qJDD(2) * t111 - t108 * t151;
t74 = -qJDD(2) * t108 - t111 * t151;
t216 = t111 * (qJDD(4) * mrSges(5,1) - mrSges(5,3) * t73) + t108 * (-qJDD(4) * mrSges(5,2) + mrSges(5,3) * t74);
t80 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t160;
t81 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t154;
t215 = -t108 * t81 + t111 * t80;
t214 = -t108 * t80 - t111 * t81;
t152 = qJD(1) * qJD(2);
t91 = t109 * t152;
t75 = qJDD(1) * t112 - t91;
t133 = qJDD(3) - t75;
t52 = qJDD(2) * t114 + t133;
t161 = qJD(1) * t112;
t145 = qJD(3) - t161;
t69 = qJD(2) * t114 + t145;
t26 = t111 * t52 - t159 * t69;
t158 = qJD(4) * t111;
t27 = t108 * t52 + t158 * t69;
t213 = -t108 * t27 - t111 * t26;
t101 = qJD(4) + qJD(5);
t208 = qJD(4) * t215 + t216;
t44 = -pkin(7) * t160 + t108 * t69;
t180 = t107 * t44;
t45 = -pkin(7) * t154 + t111 * t69;
t39 = qJD(4) * pkin(4) + t45;
t15 = t110 * t39 - t180;
t174 = t110 * t44;
t16 = t107 * t39 + t174;
t12 = qJDD(4) * pkin(4) - pkin(7) * t73 + t26;
t17 = pkin(7) * t74 + t27;
t2 = qJD(5) * t15 + t107 * t12 + t110 * t17;
t3 = -qJD(5) * t16 - t107 * t17 + t110 * t12;
t156 = qJD(5) * t107;
t34 = t101 * t163 - t107 * t159 - t108 * t156;
t123 = t132 * qJD(5);
t35 = -qJD(4) * t132 - t123;
t207 = -t131 * t3 + t132 * t2 + t15 * t35 + t16 * t34;
t115 = qJD(2) ^ 2;
t204 = t65 / 0.2e1;
t200 = mrSges(6,1) * t97;
t199 = mrSges(6,3) * t64;
t191 = t65 * mrSges(6,3);
t190 = t65 * Ifges(6,4);
t188 = -mrSges(3,2) + mrSges(4,3);
t166 = t106 * t109;
t185 = (-t105 * t96 + t166 * t97) * mrSges(6,1) + (-t105 * t97 - t166 * t96) * mrSges(6,2);
t167 = t105 * t109;
t184 = (t106 * t96 + t167 * t97) * mrSges(6,1) + (t106 * t97 - t167 * t96) * mrSges(6,2);
t182 = Ifges(5,4) * t108;
t181 = Ifges(5,4) * t111;
t149 = qJDD(2) * qJ(3);
t150 = qJDD(1) * t109;
t53 = t149 + t150 + (qJD(3) + t161) * qJD(2);
t175 = t109 * t53;
t162 = qJD(1) * t109;
t83 = qJD(2) * qJ(3) + t162;
t169 = t112 * t83;
t165 = t108 * t109;
t164 = t109 * t111;
t155 = qJDD(2) * mrSges(4,2);
t148 = pkin(4) * t154;
t147 = t112 * t152;
t90 = qJ(3) + t198;
t146 = -t151 / 0.2e1;
t143 = mrSges(5,1) * t111 - mrSges(5,2) * t108;
t141 = t111 * Ifges(5,1) - t182;
t140 = -t108 * Ifges(5,2) + t181;
t139 = -Ifges(5,5) * t108 - Ifges(5,6) * t111;
t135 = t109 * (-qJD(2) * pkin(2) + t145) + t169;
t134 = qJ(3) * t53 + qJD(3) * t83;
t128 = t83 * t143;
t127 = t108 * (-Ifges(5,2) * t111 - t182);
t126 = t111 * (-Ifges(5,1) * t108 - t181);
t122 = t131 * qJD(5);
t120 = t169 + (t108 ^ 2 + t111 ^ 2) * t69 * t109;
t100 = qJDD(4) + qJDD(5);
t24 = -qJD(2) * t123 + t107 * t74 + t110 * t73;
t25 = qJD(2) * t122 - t107 * t73 + t110 * t74;
t29 = -t64 * Ifges(6,2) + t101 * Ifges(6,6) + t190;
t59 = Ifges(6,4) * t64;
t30 = t65 * Ifges(6,1) + t101 * Ifges(6,5) - t59;
t66 = qJD(2) * t90 + t162;
t118 = t3 * mrSges(6,1) - t2 * mrSges(6,2) - t15 * t199 + t29 * t204 - t66 * (mrSges(6,1) * t65 - mrSges(6,2) * t64) - t65 * (-Ifges(6,1) * t64 - t190) / 0.2e1 + Ifges(6,6) * t25 + Ifges(6,5) * t24 - t101 * (-Ifges(6,5) * t64 - Ifges(6,6) * t65) / 0.2e1 + Ifges(6,3) * t100 + (-Ifges(6,2) * t65 + t30 - t59) * t64 / 0.2e1;
t113 = -pkin(7) - pkin(6);
t84 = pkin(4) * t158 + qJD(3);
t82 = t112 * t96 * mrSges(6,2);
t76 = t147 + t150;
t63 = Ifges(5,5) * qJD(4) + qJD(2) * t141;
t62 = Ifges(5,6) * qJD(4) + qJD(2) * t140;
t60 = -qJDD(2) * pkin(2) + t133;
t54 = t131 * t112;
t41 = mrSges(6,1) * t101 - t191;
t40 = -mrSges(6,2) * t101 - t199;
t38 = -mrSges(5,1) * t74 + mrSges(5,2) * t73;
t33 = -pkin(4) * t74 + t53;
t21 = t110 * t45 - t180;
t20 = -t107 * t45 - t174;
t14 = -mrSges(6,2) * t100 + mrSges(6,3) * t25;
t13 = mrSges(6,1) * t100 - mrSges(6,3) * t24;
t9 = -qJD(2) * t124 + t101 * t55;
t8 = qJD(2) * t125 + (qJD(4) * t131 + t122) * t112;
t4 = -mrSges(6,1) * t25 + mrSges(6,2) * t24;
t1 = [m(2) * qJDD(1) + t54 * t13 - t55 * t14 + t8 * t40 + t9 * t41 + (-m(2) - m(3) - t153) * g(3) + (-qJD(2) * t214 + t188 * qJDD(2) - t189 * t115 + t38 + t4) * t109 + (qJD(2) * t186 + qJDD(2) * t189 + t115 * t188 - t208) * t112 + m(5) * (qJD(2) * t120 + t112 * t213 + t175) + m(3) * (t109 * t76 + t112 * t75) + m(4) * (qJD(2) * t135 - t112 * t60 + t175) + m(6) * (qJD(2) * t112 * t66 + t109 * t33 + t15 * t9 + t16 * t8 - t2 * t55 + t3 * t54); -t108 * (Ifges(5,4) * t73 + Ifges(5,2) * t74) / 0.2e1 + t53 * t142 + t74 * t140 / 0.2e1 + t73 * t141 / 0.2e1 + (-g(3) * t109 + qJD(2) * qJD(3) - t147 + t149 + t53) * mrSges(4,3) + t111 * (Ifges(5,1) * t73 + Ifges(5,4) * t74) / 0.2e1 + (Ifges(4,1) + Ifges(3,3)) * qJDD(2) + (-t91 + t60) * mrSges(4,2) + (t147 - t76) * mrSges(3,2) + (t91 + t75) * mrSges(3,1) + qJDD(4) * (Ifges(5,5) * t111 - Ifges(5,6) * t108) + t101 * (Ifges(6,5) * t35 - Ifges(6,6) * t34) / 0.2e1 - pkin(2) * t155 - t62 * t158 / 0.2e1 - t63 * t159 / 0.2e1 + t126 * t151 / 0.2e1 + (t215 * t114 + t128 + t139 * qJD(4) / 0.2e1) * qJD(4) + t84 * t32 + t90 * t4 + qJD(3) * t72 - t64 * (Ifges(6,4) * t35 - Ifges(6,2) * t34) / 0.2e1 + t66 * (mrSges(6,1) * t34 + mrSges(6,2) * t35) - t34 * t29 / 0.2e1 + t35 * t30 / 0.2e1 + t36 * t13 + t37 * t14 + qJ(3) * t38 + (-qJD(1) * t120 - t114 * t213 + t134) * m(5) + ((-qJ(3) * t153 - mrSges(4,3) + t224) * t112 + (-m(6) * (-pkin(2) + t113) + mrSges(6,3) - m(5) * t114 + m(4) * pkin(2) + t225) * t109) * t212 + (-t153 * (pkin(2) * t112 + qJ(3) * t109) + (-m(5) * pkin(6) + m(6) * t113 - t225) * t112 + t224 * t109) * g(3) + (-pkin(2) * t60 - qJD(1) * t135 + t134) * m(4) - (-mrSges(6,1) * t33 + Ifges(6,4) * t24 + Ifges(6,2) * t25 + Ifges(6,6) * t100) * t132 - (mrSges(6,2) * t33 + Ifges(6,1) * t24 + Ifges(6,4) * t25 + Ifges(6,5) * t100) * t131 + t127 * t146 + (Ifges(6,1) * t35 - Ifges(6,4) * t34) * t204 + t213 * mrSges(5,3) + t214 * t162 + t216 * t114 - t186 * t161 + (-t207 - t194) * mrSges(6,3) + t218 * t40 + t219 * t41 + (t2 * t37 + t3 * t36 + t33 * t90 + (-t161 + t84) * t66 + t218 * t16 + t219 * t15) * m(6); t155 - t115 * mrSges(4,3) - t131 * t13 + t132 * t14 + t34 * t40 + t35 * t41 - m(5) * t213 + m(6) * t207 + m(4) * t60 + (-m(6) * t66 + t203 * t83 - t186) * qJD(2) + (-t109 * t212 + t194) * t153 + t208; -t32 * t148 - m(6) * (t148 * t66 + t15 * t20 + t16 * t21) + t63 * t160 / 0.2e1 + t62 * t154 / 0.2e1 + Ifges(5,5) * t73 + Ifges(5,6) * t74 - t21 * t40 - t20 * t41 + t26 * mrSges(5,1) - t27 * mrSges(5,2) + t118 - g(3) * (t82 + (-t111 * t206 - t200) * t112) - qJD(2) * t128 + Ifges(5,3) * qJDD(4) + t139 * t146 + t16 * t191 + t143 * t194 + (t107 * t2 + t110 * t3 + (-t107 * t15 + t110 * t16) * qJD(5)) * t206 - t215 * t69 + (-t126 / 0.2e1 + t127 / 0.2e1) * t115 + (-(-t105 * t165 + t106 * t111) * mrSges(5,2) - t184 + t217 * (t105 * t164 + t106 * t108)) * g(2) + (-(-t105 * t111 - t106 * t165) * mrSges(5,2) - t185 + t217 * (-t105 * t108 + t106 * t164)) * g(1) + ((qJD(5) * t40 + t13) * t110 + t107 * t14 - t156 * t41) * pkin(4); (t41 + t191) * t16 - g(3) * (-t112 * t200 + t82) - t15 * t40 - g(1) * t185 - g(2) * t184 + t118;];
tau = t1;
