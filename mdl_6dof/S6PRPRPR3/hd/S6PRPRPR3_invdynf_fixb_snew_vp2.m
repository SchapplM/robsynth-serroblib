% Calculate vector of cutting forces with Newton-Euler
% S6PRPRPR3
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
% Datum: 2019-05-04 22:37
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRPRPR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR3_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR3_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR3_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR3_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR3_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:33:04
% EndTime: 2019-05-04 22:33:07
% DurationCPUTime: 1.26s
% Computational Cost: add. (14479->152), mult. (26995->195), div. (0->0), fcn. (17059->12), ass. (0->88)
t83 = sin(qJ(4));
t108 = t83 * qJD(2);
t64 = mrSges(6,1) * t108 + qJD(4) * mrSges(6,2);
t110 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t108 - t64;
t115 = mrSges(5,1) - mrSges(6,2);
t89 = qJD(2) ^ 2;
t121 = t89 * pkin(8);
t107 = qJD(2) * qJD(4);
t86 = cos(qJ(4));
t102 = t86 * t107;
t57 = qJDD(2) * t83 + t102;
t103 = t83 * t107;
t58 = qJDD(2) * t86 - t103;
t109 = qJD(2) * t86;
t62 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t109;
t122 = pkin(9) * t89;
t123 = -pkin(4) - pkin(9);
t77 = sin(pkin(10));
t80 = cos(pkin(10));
t59 = g(1) * t77 - g(2) * t80;
t75 = -g(3) + qJDD(1);
t78 = sin(pkin(6));
t81 = cos(pkin(6));
t100 = -t59 * t78 + t75 * t81;
t44 = qJDD(3) + t100;
t87 = cos(qJ(2));
t118 = t78 * t87;
t120 = t59 * t81;
t60 = -g(1) * t80 - g(2) * t77;
t84 = sin(qJ(2));
t99 = t118 * t75 + t120 * t87 - t60 * t84;
t33 = qJDD(2) * pkin(2) + t99;
t119 = t78 * t84;
t105 = t119 * t75 + t120 * t84 + t60 * t87;
t34 = -pkin(2) * t89 + t105;
t76 = sin(pkin(11));
t79 = cos(pkin(11));
t112 = t33 * t76 + t34 * t79;
t29 = -pkin(3) * t89 + qJDD(2) * pkin(8) + t112;
t26 = t83 * t29;
t54 = (-pkin(4) * t86 - qJ(5) * t83) * qJD(2);
t88 = qJD(4) ^ 2;
t98 = -t88 * qJ(5) + t108 * t54 + qJDD(5) + t26;
t20 = t57 * pkin(5) + t123 * qJDD(4) + (-pkin(5) * t107 - t122 * t83 - t44) * t86 + t98;
t67 = pkin(5) * t108 - qJD(4) * pkin(9);
t74 = t86 ^ 2;
t124 = -2 * qJD(5);
t101 = t79 * t33 - t34 * t76;
t97 = -qJDD(2) * pkin(3) - t101;
t90 = pkin(4) * t103 + t108 * t124 + (-t57 - t102) * qJ(5) + t97;
t21 = -t67 * t108 + (-pkin(5) * t74 - pkin(8)) * t89 + t123 * t58 + t90;
t82 = sin(qJ(6));
t85 = cos(qJ(6));
t52 = -qJD(4) * t82 - t109 * t85;
t38 = qJD(6) * t52 + qJDD(4) * t85 - t58 * t82;
t53 = qJD(4) * t85 - t109 * t82;
t39 = -mrSges(7,1) * t52 + mrSges(7,2) * t53;
t69 = qJD(6) + t108;
t42 = -mrSges(7,2) * t69 + mrSges(7,3) * t52;
t50 = qJDD(6) + t57;
t16 = m(7) * (t20 * t85 - t21 * t82) - t38 * mrSges(7,3) + t50 * mrSges(7,1) - t53 * t39 + t69 * t42;
t37 = -qJD(6) * t53 - qJDD(4) * t82 - t58 * t85;
t43 = mrSges(7,1) * t69 - mrSges(7,3) * t53;
t17 = m(7) * (t20 * t82 + t21 * t85) + t37 * mrSges(7,3) - t50 * mrSges(7,2) + t52 * t39 - t69 * t43;
t63 = -mrSges(6,1) * t109 - qJD(4) * mrSges(6,3);
t96 = t82 * t16 - t85 * t17 - m(6) * (-t58 * pkin(4) - t121 + t90) - t63 * t109 + t57 * mrSges(6,3);
t127 = (t110 * t83 - t86 * t62) * qJD(2) - t115 * t58 + m(5) * (t97 - t121) + t57 * mrSges(5,2) - t96;
t117 = t86 * t44;
t114 = mrSges(5,3) + mrSges(6,1);
t113 = t29 * t86 + t44 * t83;
t55 = (mrSges(6,2) * t86 - mrSges(6,3) * t83) * qJD(2);
t111 = t55 + (-mrSges(5,1) * t86 + mrSges(5,2) * t83) * qJD(2);
t95 = -m(6) * (-qJDD(4) * pkin(4) - t117 + t98) - t85 * t16 - t82 * t17;
t12 = m(5) * (-t26 + t117) - t114 * t57 + t115 * qJDD(4) + (t62 - t63) * qJD(4) - t111 * t108 + t95;
t92 = -t88 * pkin(4) + qJDD(4) * qJ(5) + t109 * t54 + t113;
t94 = -t37 * mrSges(7,1) - t52 * t42 + m(7) * (-t74 * t122 + t58 * pkin(5) + ((2 * qJD(5)) + t67) * qJD(4) + t92) + t38 * mrSges(7,2) + t53 * t43;
t93 = -m(6) * (qJD(4) * t124 - t92) + t94;
t14 = m(5) * t113 + t114 * t58 + (-mrSges(5,2) + mrSges(6,3)) * qJDD(4) - t110 * qJD(4) + t111 * t109 + t93;
t106 = m(4) * t44 + t12 * t86 + t14 * t83;
t10 = m(4) * t101 + qJDD(2) * mrSges(4,1) - t89 * mrSges(4,2) - t127;
t7 = m(4) * t112 - mrSges(4,1) * t89 - qJDD(2) * mrSges(4,2) - t12 * t83 + t14 * t86;
t5 = m(3) * t99 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t89 + t10 * t79 + t7 * t76;
t6 = m(3) * t105 - mrSges(3,1) * t89 - qJDD(2) * mrSges(3,2) - t10 * t76 + t7 * t79;
t9 = m(3) * t100 + t106;
t104 = m(2) * t75 + t118 * t5 + t119 * t6 + t81 * t9;
t2 = m(2) * t60 - t5 * t84 + t6 * t87;
t1 = m(2) * t59 - t78 * t9 + (t5 * t87 + t6 * t84) * t81;
t3 = [-m(1) * g(1) - t1 * t77 + t2 * t80, t2, t6, t7, t14, t58 * mrSges(6,2) - t108 * t64 - t96, t17; -m(1) * g(2) + t1 * t80 + t2 * t77, t1, t5, t10, t12, -t58 * mrSges(6,1) - qJDD(4) * mrSges(6,3) - qJD(4) * t64 - t109 * t55 - t93, t16; -m(1) * g(3) + t104, t104, t9, t106, t127, t57 * mrSges(6,1) + qJDD(4) * mrSges(6,2) + qJD(4) * t63 + t108 * t55 - t95, t94;];
f_new  = t3;
