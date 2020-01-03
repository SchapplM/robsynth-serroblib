% Calculate vector of inverse dynamics joint torques for
% S4RPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRR8_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR8_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR8_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR8_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR8_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR8_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR8_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:55:05
% EndTime: 2019-12-31 16:55:09
% DurationCPUTime: 2.22s
% Computational Cost: add. (1190->234), mult. (2278->319), div. (0->0), fcn. (1237->8), ass. (0->113)
t89 = sin(qJ(3));
t161 = -t89 / 0.2e1;
t92 = cos(qJ(3));
t148 = t92 / 0.2e1;
t87 = qJ(3) + qJ(4);
t80 = sin(t87);
t81 = cos(t87);
t110 = mrSges(5,1) * t80 + mrSges(5,2) * t81;
t111 = mrSges(4,1) * t89 + mrSges(4,2) * t92;
t160 = t110 + t111;
t123 = qJD(1) * t92;
t124 = qJD(1) * t89;
t105 = -(qJD(3) * mrSges(4,1) - mrSges(4,3) * t123) * t89 + (-qJD(3) * mrSges(4,2) - mrSges(4,3) * t124) * t92;
t159 = t105 * qJD(3);
t93 = cos(qJ(1));
t158 = g(2) * t93;
t95 = -pkin(1) - pkin(5);
t67 = qJD(1) * t95 + qJD(2);
t38 = -pkin(6) * t124 + t67 * t89;
t91 = cos(qJ(4));
t133 = t38 * t91;
t39 = -pkin(6) * t123 + t92 * t67;
t37 = qJD(3) * pkin(3) + t39;
t88 = sin(qJ(4));
t18 = t37 * t88 + t133;
t157 = qJD(4) * t18;
t117 = qJD(1) * qJD(2);
t68 = qJDD(1) * qJ(2) + t117;
t116 = qJD(1) * qJD(3);
t55 = qJDD(1) * t92 - t116 * t89;
t56 = -qJDD(1) * t89 - t116 * t92;
t156 = t92 * (qJDD(3) * mrSges(4,1) - mrSges(4,3) * t55) + t89 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t56);
t122 = qJD(3) * t89;
t66 = qJDD(1) * t95 + qJDD(2);
t32 = -t122 * t67 + t92 * t66;
t121 = qJD(3) * t92;
t33 = t67 * t121 + t89 * t66;
t106 = t32 * t92 + t33 * t89;
t155 = -m(4) - m(3) - m(5);
t154 = mrSges(2,1) + mrSges(4,3) - mrSges(3,2);
t147 = pkin(3) * t89;
t153 = -m(5) * t147 + mrSges(2,2) - mrSges(3,3) - t160;
t127 = t91 * t92;
t104 = t88 * t89 - t127;
t134 = t38 * t88;
t17 = t37 * t91 - t134;
t21 = qJDD(3) * pkin(3) - pkin(6) * t55 + t32;
t25 = pkin(6) * t56 + t33;
t2 = qJD(4) * t17 + t21 * t88 + t25 * t91;
t119 = qJD(4) * t88;
t86 = qJD(3) + qJD(4);
t28 = -t119 * t89 - t122 * t88 + t127 * t86;
t50 = -t88 * t92 - t91 * t89;
t100 = t50 * qJD(4);
t29 = qJD(3) * t50 + t100;
t3 = t21 * t91 - t25 * t88 - t157;
t152 = t104 * t3 - t17 * t29 - t18 * t28 + t2 * t50;
t112 = mrSges(4,1) * t92 - mrSges(4,2) * t89;
t139 = Ifges(4,4) * t92;
t140 = Ifges(4,4) * t89;
t151 = (-Ifges(4,2) * t92 - t140) * t161 + (-Ifges(4,1) * t89 - t139) * t148 + qJ(2) * t112;
t46 = t123 * t91 - t124 * t88;
t149 = t46 / 0.2e1;
t144 = pkin(6) - t95;
t143 = mrSges(5,1) * t81;
t142 = mrSges(5,2) * t80;
t45 = t50 * qJD(1);
t141 = mrSges(5,3) * t45;
t132 = t46 * mrSges(5,3);
t131 = t46 * Ifges(5,4);
t118 = qJDD(1) * mrSges(3,2);
t59 = t144 * t92;
t73 = qJ(2) + t147;
t114 = (t117 + t68) * qJ(2);
t90 = sin(qJ(1));
t113 = -g(1) * t90 + t158;
t109 = Ifges(4,1) * t92 - t140;
t108 = -Ifges(4,2) * t89 + t139;
t107 = -Ifges(4,5) * t89 - Ifges(4,6) * t92;
t58 = t144 * t89;
t31 = -t58 * t91 - t59 * t88;
t30 = t58 * t88 - t59 * t91;
t96 = qJD(1) ^ 2;
t99 = -t96 * qJ(2) + t113;
t15 = qJD(1) * t100 + t55 * t91 + t56 * t88;
t16 = qJD(1) * qJD(4) * t104 - t55 * t88 + t56 * t91;
t23 = t45 * Ifges(5,2) + t86 * Ifges(5,6) + t131;
t42 = Ifges(5,4) * t45;
t24 = t46 * Ifges(5,1) + t86 * Ifges(5,5) + t42;
t62 = t73 * qJD(1);
t85 = qJDD(3) + qJDD(4);
t97 = t3 * mrSges(5,1) - t2 * mrSges(5,2) + t17 * t141 + t23 * t149 - t62 * (mrSges(5,1) * t46 + mrSges(5,2) * t45) - t46 * (Ifges(5,1) * t45 - t131) / 0.2e1 + Ifges(5,6) * t16 + Ifges(5,5) * t15 - t86 * (Ifges(5,5) * t45 - Ifges(5,6) * t46) / 0.2e1 + Ifges(5,3) * t85 - (-Ifges(5,2) * t46 + t24 + t42) * t45 / 0.2e1;
t94 = -pkin(6) - pkin(5);
t78 = -qJDD(1) * pkin(1) + qJDD(2);
t69 = pkin(3) * t121 + qJD(2);
t65 = t93 * t142;
t64 = t90 * t143;
t52 = t111 * qJD(1);
t48 = qJD(3) * t59;
t47 = t144 * t122;
t44 = Ifges(4,5) * qJD(3) + qJD(1) * t109;
t43 = Ifges(4,6) * qJD(3) + qJD(1) * t108;
t36 = -pkin(3) * t56 + t68;
t35 = mrSges(5,1) * t86 - t132;
t34 = -mrSges(5,2) * t86 + t141;
t27 = -mrSges(5,1) * t45 + mrSges(5,2) * t46;
t20 = t39 * t91 - t134;
t19 = -t39 * t88 - t133;
t10 = -mrSges(5,2) * t85 + mrSges(5,3) * t16;
t9 = mrSges(5,1) * t85 - mrSges(5,3) * t15;
t5 = -qJD(4) * t31 + t47 * t91 + t48 * t88;
t4 = qJD(4) * t30 + t47 * t88 - t48 * t91;
t1 = [(Ifges(4,4) * t55 + Ifges(4,2) * t56) * t161 + t95 * t159 + (0.2e1 * Ifges(4,5) * t148 - Ifges(4,6) * t89) * qJDD(3) + m(4) * t114 + m(5) * (t17 * t5 + t18 * t4 + t2 * t31 + t3 * t30 + t36 * t73 + t62 * t69) + (t152 - t158) * mrSges(5,3) + (t155 * (t93 * pkin(1) + t90 * qJ(2)) + (-m(4) * pkin(5) + m(5) * t94 - t154) * t93 + t153 * t90) * g(2) - t106 * mrSges(4,3) + (t106 * m(4) + t156) * t95 + t151 * t116 - t43 * t121 / 0.2e1 - t44 * t122 / 0.2e1 - pkin(1) * t118 + t86 * (Ifges(5,5) * t29 - Ifges(5,6) * t28) / 0.2e1 + t69 * t27 + t73 * (-mrSges(5,1) * t16 + mrSges(5,2) * t15) + t78 * mrSges(3,2) + qJ(2) * (-mrSges(4,1) * t56 + mrSges(4,2) * t55) + t62 * (mrSges(5,1) * t28 + mrSges(5,2) * t29) + qJD(2) * t52 + t45 * (Ifges(5,4) * t29 - Ifges(5,2) * t28) / 0.2e1 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + t4 * t34 + t5 * t35 - t28 * t23 / 0.2e1 + t29 * t24 / 0.2e1 + t30 * t9 + t31 * t10 + ((-m(5) * (-pkin(1) + t94) + mrSges(5,3) - m(4) * t95 + m(3) * pkin(1) + t154) * t90 + (t155 * qJ(2) + t153) * t93) * g(1) + (0.2e1 * mrSges(3,3) + t111) * t68 + m(3) * (-pkin(1) * t78 + t114) + t56 * t108 / 0.2e1 + t55 * t109 / 0.2e1 - (mrSges(5,2) * t36 + Ifges(5,1) * t15 + Ifges(5,4) * t16 + Ifges(5,5) * t85) * t104 + qJD(3) ^ 2 * t107 / 0.2e1 + (Ifges(4,1) * t55 + Ifges(4,4) * t56) * t148 + (Ifges(5,1) * t29 - Ifges(5,4) * t28) * t149 + (-mrSges(5,1) * t36 + Ifges(5,4) * t15 + Ifges(5,2) * t16 + Ifges(5,6) * t85) * t50; t118 - t96 * mrSges(3,3) - t50 * t10 + t28 * t34 + t29 * t35 - t104 * t9 + t159 + (-t27 - t52) * qJD(1) + (-t62 * qJD(1) + t113 - t152) * m(5) + (t106 + t99) * m(4) + (t78 + t99) * m(3) + t156; (-t64 + (-t112 + t142) * t90) * g(1) + (-t65 + (t112 + t143) * t93) * g(2) + t160 * g(3) - t105 * t67 + t18 * t132 + t97 + Ifges(4,5) * t55 + Ifges(4,6) * t56 + (t43 * t148 + t89 * t44 / 0.2e1 - qJD(3) * t107 / 0.2e1 - t151 * qJD(1) + (-m(5) * t62 - t27) * t92 * pkin(3)) * qJD(1) + ((g(3) * t89 + t113 * t92 - t119 * t17) * m(5) + (m(5) * t2 - qJD(4) * t35 + t10) * t88 + (t9 + t34 * qJD(4) + (t3 + t157) * m(5)) * t91) * pkin(3) - m(5) * (t17 * t19 + t18 * t20) - t20 * t34 - t19 * t35 + t32 * mrSges(4,1) - t33 * mrSges(4,2) + Ifges(4,3) * qJDD(3); (t35 + t132) * t18 + t97 - g(2) * (-t143 * t93 + t65) - g(1) * (-t142 * t90 + t64) + g(3) * t110 - t17 * t34;];
tau = t1;
