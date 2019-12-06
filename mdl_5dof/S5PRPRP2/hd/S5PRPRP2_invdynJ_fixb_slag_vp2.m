% Calculate vector of inverse dynamics joint torques for
% S5PRPRP2
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPRP2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP2_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP2_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP2_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP2_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP2_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:30:31
% EndTime: 2019-12-05 15:30:45
% DurationCPUTime: 4.91s
% Computational Cost: add. (1223->276), mult. (2654->368), div. (0->0), fcn. (1598->6), ass. (0->141)
t203 = Ifges(5,4) + Ifges(6,4);
t215 = Ifges(5,1) + Ifges(6,1);
t214 = Ifges(5,2) + Ifges(6,2);
t90 = cos(qJ(4));
t218 = t203 * t90;
t89 = sin(qJ(4));
t217 = t203 * t89;
t86 = sin(pkin(8));
t216 = -t86 / 0.2e1;
t202 = Ifges(5,5) + Ifges(6,5);
t201 = Ifges(5,6) + Ifges(6,6);
t213 = Ifges(5,3) + Ifges(6,3);
t212 = (-t214 * t89 + t218) * t86;
t211 = (t215 * t90 - t217) * t86;
t68 = qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t87 = cos(pkin(8));
t72 = -qJD(2) * t87 + qJD(4);
t210 = t90 * (t212 * qJD(2) + t201 * t72) + t89 * (t211 * qJD(2) + t202 * t72);
t144 = qJ(3) * qJD(2);
t80 = t87 * qJD(1);
t62 = t144 * t86 - t80;
t171 = t62 * t86;
t61 = (pkin(4) * t89 + qJ(3)) * t86;
t29 = qJD(2) * t61 + qJD(5) - t80;
t174 = t29 * t86;
t209 = (-t201 * t90 - t202 * t89) * t72 * t216 - (mrSges(5,1) * t90 - mrSges(5,2) * t89) * t171 - (mrSges(6,1) * t90 - mrSges(6,2) * t89) * t174;
t83 = t86 ^ 2;
t193 = mrSges(4,3) * (t87 ^ 2 + t83);
t185 = m(6) * pkin(4);
t206 = -mrSges(5,1) - mrSges(6,1);
t205 = mrSges(5,2) + mrSges(6,2);
t204 = -mrSges(4,3) + mrSges(3,2);
t110 = pkin(3) * t87 + pkin(6) * t86;
t64 = -pkin(2) - t110;
t162 = t87 * t90;
t71 = qJ(3) * t162;
t28 = t89 * t64 + t71;
t197 = t185 + mrSges(6,1);
t139 = qJD(2) * qJD(4);
t51 = (qJDD(2) * t90 - t139 * t89) * t86;
t52 = (-qJDD(2) * t89 - t139 * t90) * t86;
t142 = qJDD(2) * t87;
t70 = qJDD(4) - t142;
t196 = t201 * t52 + t202 * t51 + t213 * t70;
t53 = -t87 * qJDD(1) + t68 * t86;
t54 = qJDD(1) * t86 + t68 * t87;
t195 = t53 * t86 + t54 * t87;
t141 = m(4) + m(5) + m(6);
t194 = mrSges(5,1) + t197;
t108 = mrSges(4,1) * t87 - mrSges(4,2) * t86;
t190 = mrSges(3,1) + t108 + (mrSges(5,3) + mrSges(6,3)) * t86;
t149 = qJD(2) * t86;
t131 = t90 * t149;
t115 = mrSges(6,3) * t131;
t46 = mrSges(6,1) * t72 - t115;
t117 = mrSges(5,3) * t131;
t47 = mrSges(5,1) * t72 - t117;
t154 = t46 + t47;
t132 = t89 * t149;
t116 = mrSges(6,3) * t132;
t44 = -mrSges(6,2) * t72 - t116;
t118 = mrSges(5,3) * t132;
t45 = -mrSges(5,2) * t72 - t118;
t155 = t44 + t45;
t189 = t154 * t89 - t155 * t90;
t188 = ((-t215 * t89 - t218) * t90 / 0.2e1 - (-t214 * t90 - t217) * t89 / 0.2e1) * t83;
t145 = qJD(5) * t86;
t122 = qJD(2) * t145;
t146 = qJD(4) * t90;
t50 = qJDD(2) * t64 + qJDD(3);
t55 = qJD(2) * t64 + qJD(3);
t137 = t55 * t146 + t89 * t50 + t90 * t54;
t63 = qJD(1) * t86 + t144 * t87;
t2 = qJ(5) * t52 + (-qJD(4) * t63 - t122) * t89 + t137;
t147 = qJD(4) * t89;
t3 = -t147 * t63 + t137;
t14 = t55 * t89 + t63 * t90;
t4 = -qJD(4) * t14 + t90 * t50 - t54 * t89;
t187 = -t4 * mrSges(5,1) + t3 * mrSges(5,2) + t2 * mrSges(6,2);
t186 = qJD(2) ^ 2;
t184 = t51 / 0.2e1;
t183 = t52 / 0.2e1;
t182 = t70 / 0.2e1;
t180 = m(2) + m(3);
t179 = mrSges(6,2) * t90;
t85 = pkin(7) + qJ(2);
t81 = sin(t85);
t169 = t81 * t89;
t82 = cos(t85);
t168 = t82 * t89;
t165 = t86 * t89;
t164 = t86 * t90;
t163 = t87 * t89;
t19 = mrSges(6,1) * t70 - mrSges(6,3) * t51;
t20 = mrSges(5,1) * t70 - mrSges(5,3) * t51;
t157 = t19 + t20;
t21 = -mrSges(6,2) * t70 + mrSges(6,3) * t52;
t22 = -mrSges(5,2) * t70 + mrSges(5,3) * t52;
t156 = t21 + t22;
t148 = qJD(3) * t87;
t153 = t64 * t146 + t90 * t148;
t151 = qJ(3) * t89;
t150 = qJ(5) * t86;
t143 = qJDD(2) * t86;
t134 = t87 * t151;
t133 = t90 * t150;
t130 = t89 * t148;
t129 = t86 * t147;
t128 = t86 * t146;
t125 = qJ(5) * t149;
t10 = -t52 * mrSges(6,1) + t51 * mrSges(6,2);
t13 = t90 * t55 - t63 * t89;
t109 = -mrSges(4,1) * t142 + mrSges(4,2) * t143;
t107 = mrSges(5,1) * t89 + mrSges(5,2) * t90;
t106 = t63 * t87 + t171;
t105 = (pkin(4) * t90 + pkin(3)) * t87 - t86 * (-qJ(5) - pkin(6));
t37 = -t163 * t82 + t81 * t90;
t35 = t163 * t81 + t82 * t90;
t104 = (mrSges(6,1) * t89 + t179) * t86;
t8 = -t125 * t90 + t13;
t78 = -qJDD(2) * pkin(2) + qJDD(3);
t60 = t90 * t64;
t58 = (pkin(4) * t146 + qJD(3)) * t86;
t56 = t107 * t86;
t49 = t107 * t149;
t48 = qJD(2) * t104;
t38 = t162 * t82 + t169;
t36 = -t162 * t81 + t168;
t27 = t60 - t134;
t18 = -t150 * t89 + t28;
t17 = -qJD(4) * t28 - t130;
t16 = -qJD(4) * t134 + t153;
t15 = -t133 + t60 + (-pkin(4) - t151) * t87;
t12 = -pkin(4) * t52 + qJDD(5) + t53;
t11 = -mrSges(5,1) * t52 + mrSges(5,2) * t51;
t9 = -t125 * t89 + t14;
t7 = -t130 - t90 * t145 + (-t71 + (-t64 + t150) * t89) * qJD(4);
t6 = -t89 * t145 + (-t133 - t134) * qJD(4) + t153;
t5 = pkin(4) * t72 + t8;
t1 = pkin(4) * t70 - qJ(5) * t51 - t122 * t90 + t4;
t23 = [t180 * qJDD(1) + (-m(5) * t53 - m(6) * t12 - t10 - t11) * t87 + (-t141 - t180) * g(3) + (t156 * t90 - t157 * t89 + (-t154 * t90 - t155 * t89) * qJD(4) + m(5) * (-t13 * t146 - t14 * t147 + t3 * t90 - t4 * t89) + m(6) * (-t1 * t89 - t146 * t5 - t147 * t9 + t2 * t90)) * t86 + m(4) * (-t53 * t87 + t54 * t86); -pkin(2) * t109 - t78 * t108 + (-t169 * t185 + t204 * t81 + t206 * t38 - t205 * t37 - t141 * (t82 * pkin(2) + t81 * qJ(3)) + (-m(5) * t110 - m(6) * t105 - t190) * t82) * g(2) + (-t1 * t164 - t128 * t9 + t129 * t5 - t165 * t2) * mrSges(6,3) + (-t128 * t14 + t129 * t13 - t164 * t4 - t165 * t3) * mrSges(5,3) + t12 * t104 + m(5) * (t13 * t17 + t14 * t16 + t27 * t4 + t28 * t3) + (t210 * t216 - t209) * qJD(4) + (t202 * t70 + t203 * t52 + t215 * t51) * t164 / 0.2e1 + (-t1 * mrSges(6,1) + Ifges(4,4) * t143 + Ifges(4,2) * t142 - t202 * t184 - t201 * t183 - t213 * t182 + t187 - t196 / 0.2e1) * t87 - (t201 * t70 + t203 * t51 + t214 * t52) * t165 / 0.2e1 + m(6) * (t1 * t15 + t12 * t61 + t18 * t2 + t29 * t58 + t5 * t7 + t6 * t9) + t53 * t56 + t58 * t48 + t61 * t10 + t6 * t44 + t16 * t45 + t7 * t46 + t17 * t47 + t18 * t21 + t27 * t20 + t28 * t22 + t15 * t19 + (Ifges(4,4) * t142 + Ifges(4,1) * t143 + m(5) * (qJ(3) * t53 + qJD(3) * t62) + qJ(3) * t11 + qJD(3) * t49 + (-t201 * t89 + t202 * t90) * t182) * t86 + t195 * mrSges(4,3) + m(4) * (-pkin(2) * t78 + qJ(3) * t195 + t106 * qJD(3)) + t188 * t139 + (-t168 * t185 + t206 * t36 - t205 * t35 + (m(4) * pkin(2) - m(5) * t64 - m(6) * (-pkin(2) - t105) + t190) * t81 + (-t141 * qJ(3) + t204) * t82) * g(1) + t68 * t193 + t211 * t184 + t212 * t183 + Ifges(3,3) * qJDD(2); t157 * t90 + t156 * t89 - t186 * t193 - t189 * qJD(4) + m(6) * (t1 * t90 + t2 * t89 + (-t5 * t89 + t9 * t90) * qJD(4)) + m(5) * (t3 * t89 + t4 * t90 + (-t13 * t89 + t14 * t90) * qJD(4)) + m(4) * t78 + ((-t48 - t49) * t86 + t189 * t87 - m(5) * (-t13 * t163 + t14 * t162 + t171) - m(6) * (t162 * t9 - t163 * t5 + t174) - m(4) * t106) * qJD(2) + t109 + (-g(1) * t81 + g(2) * t82) * t141; -t188 * t186 + (-(-t197 * t89 - t179) * t86 + t56) * g(3) + (t117 + t47) * t14 + (t115 + t46 - m(6) * (-t5 + t8)) * t9 + (-t45 - t118) * t13 + t197 * t1 + (-t194 * t37 + t205 * t38) * g(1) + (t194 * t35 - t205 * t36) * g(2) - t8 * t44 - t187 - t5 * t116 + t196 + t210 * t149 / 0.2e1 + t209 * qJD(2) + ((-m(6) * t29 - t48) * t131 + t19) * pkin(4); (g(3) * t87 + t12) * m(6) + ((-g(1) * t82 - g(2) * t81) * m(6) + (t89 * t44 + t90 * t46 - m(6) * (-t5 * t90 - t89 * t9)) * qJD(2)) * t86 + t10;];
tau = t23;
