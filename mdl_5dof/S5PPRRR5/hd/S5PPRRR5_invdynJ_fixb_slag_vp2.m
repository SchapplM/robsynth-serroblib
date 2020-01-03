% Calculate vector of inverse dynamics joint torques for
% S5PPRRR5
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
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
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
% Datum: 2019-12-31 17:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRRR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR5_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR5_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR5_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR5_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR5_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:35:39
% EndTime: 2019-12-31 17:35:41
% DurationCPUTime: 1.28s
% Computational Cost: add. (1135->192), mult. (1988->269), div. (0->0), fcn. (1256->10), ass. (0->98)
t80 = sin(qJ(5));
t83 = cos(qJ(5));
t63 = -mrSges(6,1) * t83 + t80 * mrSges(6,2);
t160 = mrSges(5,1) - t63;
t159 = -mrSges(5,2) + mrSges(6,3);
t150 = t80 / 0.2e1;
t78 = sin(pkin(8));
t79 = cos(pkin(8));
t82 = sin(qJ(3));
t85 = cos(qJ(3));
t54 = -t78 * t85 + t79 * t82;
t118 = qJ(3) + qJ(4);
t105 = sin(t118);
t106 = cos(t118);
t43 = -t105 * t78 - t106 * t79;
t44 = t105 * t79 - t106 * t78;
t158 = -g(1) * t43 - g(2) * t44;
t103 = mrSges(6,1) * t80 + mrSges(6,2) * t83;
t122 = qJD(2) * t82;
t81 = sin(qJ(4));
t113 = t81 * t122;
t64 = qJD(3) * pkin(3) + qJD(2) * t85;
t84 = cos(qJ(4));
t30 = t64 * t84 - t113;
t77 = qJD(3) + qJD(4);
t22 = -pkin(4) * t77 - t30;
t157 = t22 * t103 + qJD(5) * (Ifges(6,5) * t83 - Ifges(6,6) * t80) / 0.2e1;
t156 = t159 * t44 - t160 * t43;
t155 = t159 * t43 + t160 * t44;
t45 = t63 * t77;
t154 = -t77 * mrSges(5,1) + t45;
t119 = qJD(5) * t83;
t120 = qJD(5) * t80;
t31 = t122 * t84 + t64 * t81;
t23 = pkin(7) * t77 + t31;
t14 = -qJD(1) * t83 - t23 * t80;
t97 = qJD(1) * t80 - t23 * t83;
t153 = -t14 * t119 + t120 * t97;
t76 = qJDD(3) + qJDD(4);
t46 = -t120 * t77 + t76 * t83;
t24 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t46;
t47 = t119 * t77 + t76 * t80;
t25 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t47;
t152 = t83 * t24 - t80 * t25;
t130 = t80 * mrSges(6,3);
t59 = qJD(5) * mrSges(6,1) - t130 * t77;
t134 = t77 * t83;
t60 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t134;
t151 = -t77 * mrSges(5,2) - t80 * t59 + t83 * t60;
t116 = qJD(2) * qJD(3);
t109 = t82 * t116;
t61 = qJDD(2) * t85 - t109;
t52 = qJDD(3) * pkin(3) + t61;
t108 = t85 * t116;
t62 = qJDD(2) * t82 + t108;
t9 = -qJD(4) * t31 + t52 * t84 - t62 * t81;
t148 = pkin(3) * t81;
t147 = pkin(3) * t84;
t121 = qJD(4) * t84;
t8 = -qJD(4) * t113 + t121 * t64 + t52 * t81 + t62 * t84;
t5 = pkin(7) * t76 + t8;
t2 = qJD(5) * t14 - qJDD(1) * t80 + t5 * t83;
t144 = t2 * t83;
t3 = qJD(5) * t97 - qJDD(1) * t83 - t5 * t80;
t143 = t3 * t80;
t141 = Ifges(6,1) * t80;
t140 = Ifges(6,4) * t80;
t139 = Ifges(6,4) * t83;
t138 = Ifges(6,2) * t83;
t137 = t76 * mrSges(5,2);
t123 = pkin(3) * qJD(4);
t117 = m(3) + m(4) + m(5);
t110 = m(6) + t117;
t104 = t54 * pkin(3);
t102 = t138 + t140;
t100 = t14 * t83 - t80 * t97;
t99 = t14 * t80 + t83 * t97;
t53 = -t78 * t82 - t79 * t85;
t58 = t81 * t85 + t82 * t84;
t57 = t81 * t82 - t84 * t85;
t94 = t53 * pkin(3);
t92 = t80 * (Ifges(6,1) * t83 - t140);
t90 = -qJD(5) * t100 - t143;
t89 = t90 + t144;
t88 = (-t59 * t83 - t60 * t80) * qJD(5) + t152;
t32 = Ifges(6,6) * qJD(5) + t102 * t77;
t68 = Ifges(6,4) * t134;
t33 = Ifges(6,5) * qJD(5) + t141 * t77 + t68;
t6 = -pkin(4) * t76 - t9;
t87 = t9 * mrSges(5,1) - t8 * mrSges(5,2) + mrSges(6,3) * t144 + t6 * t63 + (Ifges(6,1) * t47 + Ifges(6,4) * t46) * t150 + t83 * (Ifges(6,4) * t47 + Ifges(6,2) * t46) / 0.2e1 + t46 * t102 / 0.2e1 + t47 * (t139 + t141) / 0.2e1 - t32 * t120 / 0.2e1 + Ifges(5,3) * t76 + (t33 + t77 * (-Ifges(6,2) * t80 + t139)) * t119 / 0.2e1 + (0.2e1 * Ifges(6,5) * t150 + Ifges(6,6) * t83) * qJDD(5) + (t92 * t77 / 0.2e1 + t157) * qJD(5);
t86 = qJD(3) ^ 2;
t72 = -pkin(4) - t147;
t41 = t44 * pkin(4);
t40 = t43 * pkin(4);
t17 = t77 * t58;
t16 = t77 * t57;
t13 = -mrSges(6,1) * t46 + mrSges(6,2) * t47;
t1 = [-t60 * t119 + t59 * t120 - t80 * t24 - t83 * t25 + m(6) * (qJD(5) * t99 - t2 * t80 - t3 * t83) + (m(2) + t117) * qJDD(1) + (-m(2) - t110) * g(3); m(3) * qJDD(2) + t57 * t13 + t17 * t45 + (-qJDD(3) * t82 - t85 * t86) * mrSges(4,2) + (-t17 * t77 - t57 * t76) * mrSges(5,1) + (qJDD(3) * t85 - t82 * t86) * mrSges(4,1) - t151 * t16 + (t88 - t137) * t58 + m(5) * (-t16 * t31 - t17 * t30 - t57 * t9 + t58 * t8) + m(6) * (t16 * t99 + t17 * t22 + t57 * t6 + t58 * t89) + m(4) * (t61 * t85 + t62 * t82) + (-g(1) * t78 + g(2) * t79) * t110; -t137 * t148 + m(6) * (t6 * t72 + (t22 * t81 - t84 * t99) * t123) - t3 * t130 + t87 + t72 * t13 + t76 * mrSges(5,1) * t147 + Ifges(4,3) * qJDD(3) + t154 * t81 * t123 + t153 * mrSges(6,3) + (t108 - t62) * mrSges(4,2) + (t109 + t61) * mrSges(4,1) + (-m(5) * t94 - m(6) * (-pkin(7) * t44 + t40 + t94) - t53 * mrSges(4,1) - t54 * mrSges(4,2) + t156) * g(2) + (m(5) * t104 - m(6) * (-pkin(7) * t43 - t104 - t41) + t54 * mrSges(4,1) - t53 * mrSges(4,2) + t155) * g(1) + (m(6) * t89 - t119 * t59 - t120 * t60 + t152) * (pkin(7) + t148) + ((m(5) * t30 - m(6) * t22 - t154) * t58 - (-m(5) * t31 + m(6) * t99 - t151) * t57) * qJD(2) + (m(5) * (t8 * t81 + t84 * t9 + (-t30 * t81 + t31 * t84) * qJD(4)) + t151 * t121) * pkin(3); t88 * pkin(7) + t90 * mrSges(6,3) + t156 * g(2) + t155 * g(1) - t154 * t31 - t151 * t30 + t87 - pkin(4) * t13 + ((-t143 + t144 + t153 - t158) * pkin(7) - g(2) * t40 + g(1) * t41 - t22 * t31 + t99 * t30 - t6 * pkin(4)) * m(6); t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t47 + Ifges(6,6) * t46 + Ifges(6,3) * qJDD(5) - g(3) * t63 - t14 * t60 - t97 * t59 + (t32 * t150 + (-t92 / 0.2e1 + t138 * t150) * t77 + t100 * mrSges(6,3) - (t33 + t68) * t83 / 0.2e1 - t157) * t77 + t158 * t103;];
tau = t1;
