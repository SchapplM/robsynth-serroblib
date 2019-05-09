% Calculate vector of cutting forces with Newton-Euler
% S6RRPPRP4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
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
% f_new [3x7]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 09:23
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPPRP4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP4_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP4_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP4_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP4_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP4_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:18:42
% EndTime: 2019-05-06 09:18:47
% DurationCPUTime: 1.63s
% Computational Cost: add. (16115->196), mult. (35189->240), div. (0->0), fcn. (22319->8), ass. (0->92)
t153 = -2 * qJD(4);
t108 = cos(qJ(2));
t133 = qJD(1) * t108;
t104 = sin(pkin(9));
t106 = sin(qJ(2));
t134 = qJD(1) * t106;
t137 = cos(pkin(9));
t83 = -qJD(2) * t137 + t104 * t134;
t129 = t83 * t133;
t110 = qJD(2) ^ 2;
t111 = qJD(1) ^ 2;
t107 = sin(qJ(1));
t109 = cos(qJ(1));
t124 = -g(1) * t109 - g(2) * t107;
t79 = -pkin(1) * t111 + qJDD(1) * pkin(7) + t124;
t138 = -g(3) * t108 - t106 * t79;
t87 = (-pkin(2) * t108 - qJ(3) * t106) * qJD(1);
t44 = -qJDD(2) * pkin(2) - qJ(3) * t110 + t134 * t87 + qJDD(3) - t138;
t132 = qJD(1) * qJD(2);
t127 = t108 * t132;
t89 = qJDD(1) * t106 + t127;
t67 = -qJDD(2) * t137 + t104 * t89;
t68 = qJDD(2) * t104 + t137 * t89;
t84 = qJD(2) * t104 + t134 * t137;
t152 = t44 - (t129 + t68) * qJ(4) - (t133 * t84 - t67) * pkin(3) + t84 * t153;
t105 = sin(qJ(5));
t146 = cos(qJ(5));
t126 = g(1) * t107 - g(2) * t109;
t78 = -qJDD(1) * pkin(1) - pkin(7) * t111 - t126;
t96 = t106 * t132;
t90 = qJDD(1) * t108 - t96;
t40 = (-t89 - t127) * qJ(3) + (-t90 + t96) * pkin(2) + t78;
t128 = -g(3) * t106 + t108 * t79;
t45 = -pkin(2) * t110 + qJDD(2) * qJ(3) + t133 * t87 + t128;
t123 = -t104 * t45 + t137 * t40;
t135 = t108 ^ 2 * t111;
t56 = pkin(3) * t83 - qJ(4) * t84;
t23 = t90 * pkin(3) - qJ(4) * t135 + qJDD(4) - t123 + ((2 * qJD(3)) + t56) * t84;
t17 = (-t68 + t129) * pkin(8) + (t83 * t84 + t90) * pkin(4) + t23;
t149 = -2 * qJD(3);
t130 = t104 * t40 + t137 * t45 + t149 * t83;
t118 = -pkin(3) * t135 - t90 * qJ(4) + t133 * t153 - t83 * t56 + t130;
t69 = pkin(4) * t133 - pkin(8) * t84;
t81 = t83 ^ 2;
t19 = -pkin(4) * t81 + pkin(8) * t67 - t133 * t69 + t118;
t142 = t105 * t17 + t146 * t19;
t54 = t105 * t84 - t146 * t83;
t55 = t105 * t83 + t146 * t84;
t36 = pkin(5) * t54 - qJ(6) * t55;
t94 = qJD(5) + t133;
t48 = -mrSges(7,1) * t94 + mrSges(7,2) * t55;
t86 = qJDD(5) + t90;
t93 = t94 ^ 2;
t131 = m(7) * (-pkin(5) * t93 + qJ(6) * t86 + 0.2e1 * qJD(6) * t94 - t36 * t54 + t142) + t94 * t48 + t86 * mrSges(7,3);
t37 = mrSges(7,1) * t54 - mrSges(7,3) * t55;
t141 = -mrSges(6,1) * t54 - mrSges(6,2) * t55 - t37;
t143 = -mrSges(6,3) - mrSges(7,2);
t30 = qJD(5) * t55 + t105 * t68 - t146 * t67;
t47 = mrSges(6,1) * t94 - mrSges(6,3) * t55;
t10 = m(6) * t142 - t86 * mrSges(6,2) + t141 * t54 + t143 * t30 - t94 * t47 + t131;
t120 = -t105 * t19 + t146 * t17;
t147 = m(7) * (-pkin(5) * t86 - qJ(6) * t93 + t36 * t55 + qJDD(6) - t120);
t31 = -qJD(5) * t54 + t105 * t67 + t146 * t68;
t46 = -mrSges(6,2) * t94 - mrSges(6,3) * t54;
t49 = -mrSges(7,2) * t54 + mrSges(7,3) * t94;
t11 = m(6) * t120 - t147 + (t46 + t49) * t94 + (mrSges(6,1) + mrSges(7,1)) * t86 + t141 * t55 + t143 * t31;
t121 = m(5) * t118 - mrSges(5,3) * t90 + t10 * t146 - t105 * t11;
t66 = mrSges(5,1) * t133 + mrSges(5,2) * t84;
t139 = -mrSges(4,1) * t133 - mrSges(4,3) * t84 - t66;
t57 = mrSges(5,1) * t83 - mrSges(5,3) * t84;
t140 = -mrSges(4,1) * t83 - mrSges(4,2) * t84 - t57;
t144 = -mrSges(4,3) - mrSges(5,2);
t5 = m(4) * t130 + t90 * mrSges(4,2) + t133 * t139 + t140 * t83 + t144 * t67 + t121;
t119 = m(5) * t23 + t105 * t10 + t11 * t146;
t63 = -mrSges(5,2) * t83 - mrSges(5,3) * t133;
t64 = mrSges(4,2) * t133 - mrSges(4,3) * t83;
t6 = m(4) * t123 + (-mrSges(4,1) - mrSges(5,1)) * t90 + (m(4) * t149 + t140) * t84 + t144 * t68 + (-t63 - t64) * t133 - t119;
t91 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t134;
t92 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t133;
t151 = m(3) * t78 - t90 * mrSges(3,1) + t89 * mrSges(3,2) + (t106 * t91 - t108 * t92) * qJD(1) + t104 * t5 + t137 * t6;
t112 = -pkin(4) * t67 - pkin(8) * t81 + t69 * t84 - t152;
t125 = m(7) * (t112 + (t55 * t94 + t30) * pkin(5) + (t54 * t94 - t31) * qJ(6) - 0.2e1 * qJD(6) * t55) - t31 * mrSges(7,3) + t30 * mrSges(7,1) - t55 * t48 + t54 * t49;
t117 = m(6) * t112 + mrSges(6,1) * t30 + mrSges(6,2) * t31 + t46 * t54 + t47 * t55 + t125;
t116 = m(5) * t152 + mrSges(5,1) * t67 + t63 * t83 - t117;
t150 = m(4) * t44 + t67 * mrSges(4,1) + t139 * t84 + (mrSges(4,2) - mrSges(5,3)) * t68 + t83 * t64 + t116;
t88 = (-mrSges(3,1) * t108 + mrSges(3,2) * t106) * qJD(1);
t4 = m(3) * t128 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t90 - qJD(2) * t91 - t104 * t6 + t133 * t88 + t137 * t5;
t8 = m(3) * t138 + qJDD(2) * mrSges(3,1) - t89 * mrSges(3,3) + qJD(2) * t92 - t134 * t88 - t150;
t148 = t106 * t4 + t108 * t8;
t2 = m(2) * t126 + qJDD(1) * mrSges(2,1) - t111 * mrSges(2,2) - t151;
t1 = m(2) * t124 - mrSges(2,1) * t111 - qJDD(1) * mrSges(2,2) - t106 * t8 + t108 * t4;
t3 = [-m(1) * g(1) + t1 * t109 - t107 * t2, t1, t4, t5, -t67 * mrSges(5,2) - t133 * t66 - t83 * t57 + t121, t10, -t30 * mrSges(7,2) - t54 * t37 + t131; -m(1) * g(2) + t1 * t107 + t109 * t2, t2, t8, t6, -t68 * mrSges(5,3) - t84 * t66 + t116, t11, t125; (-m(1) - m(2)) * g(3) + t148, -m(2) * g(3) + t148, t151, t150, t90 * mrSges(5,1) + t68 * mrSges(5,2) + t133 * t63 + t84 * t57 + t119, t117, -t86 * mrSges(7,1) + t31 * mrSges(7,2) + t55 * t37 - t94 * t49 + t147;];
f_new  = t3;
