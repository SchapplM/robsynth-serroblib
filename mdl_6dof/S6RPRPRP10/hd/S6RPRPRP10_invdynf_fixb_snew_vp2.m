% Calculate vector of cutting forces with Newton-Euler
% S6RPRPRP10
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
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
% Datum: 2019-05-05 18:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRPRP10_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP10_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP10_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP10_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP10_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP10_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP10_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:09:32
% EndTime: 2019-05-05 18:09:35
% DurationCPUTime: 0.86s
% Computational Cost: add. (5693->175), mult. (11124->204), div. (0->0), fcn. (5507->6), ass. (0->83)
t83 = qJD(1) ^ 2;
t128 = -pkin(1) - pkin(7);
t107 = t128 * t83;
t80 = cos(qJ(3));
t112 = t80 * qJD(1);
t78 = sin(qJ(3));
t113 = qJD(1) * t78;
t119 = mrSges(4,1) - mrSges(5,2);
t111 = qJD(1) * qJD(3);
t103 = t80 * t111;
t56 = t78 * qJDD(1) + t103;
t70 = t78 * t111;
t57 = t80 * qJDD(1) - t70;
t59 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t113;
t60 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t112;
t122 = cos(qJ(5));
t61 = mrSges(5,1) * t113 - qJD(3) * mrSges(5,3);
t62 = mrSges(5,1) * t112 + qJD(3) * mrSges(5,2);
t127 = pkin(3) + pkin(8);
t63 = pkin(4) * t112 - qJD(3) * pkin(8);
t76 = t78 ^ 2;
t132 = -2 * qJD(4);
t79 = sin(qJ(1));
t81 = cos(qJ(1));
t101 = -t81 * g(1) - t79 * g(2);
t96 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t101;
t89 = pkin(3) * t103 + t112 * t132 + t96 + (-t57 + t70) * qJ(4);
t14 = -t63 * t112 + t127 * t56 + (-pkin(4) * t76 + t128) * t83 + t89;
t124 = pkin(8) * t83;
t105 = t79 * g(1) - t81 * g(2);
t94 = -t83 * qJ(2) + qJDD(2) - t105;
t40 = t128 * qJDD(1) + t94;
t121 = t80 * t40;
t53 = (pkin(3) * t78 - qJ(4) * t80) * qJD(1);
t82 = qJD(3) ^ 2;
t93 = -t82 * qJ(4) + t53 * t112 + qJDD(4) - t121;
t18 = t57 * pkin(4) - t127 * qJDD(3) + (pkin(4) * t111 + t80 * t124 - g(3)) * t78 + t93;
t77 = sin(qJ(5));
t115 = t122 * t14 + t77 * t18;
t51 = t77 * qJD(3) - t122 * t113;
t52 = t122 * qJD(3) + t77 * t113;
t30 = t51 * pkin(5) - t52 * qJ(6);
t68 = qJD(5) + t112;
t35 = -t68 * mrSges(7,1) + t52 * mrSges(7,2);
t50 = qJDD(5) + t57;
t64 = t68 ^ 2;
t108 = t68 * t35 + t50 * mrSges(7,3) + m(7) * (-t64 * pkin(5) + t50 * qJ(6) + 0.2e1 * qJD(6) * t68 - t51 * t30 + t115);
t31 = t51 * mrSges(7,1) - t52 * mrSges(7,3);
t114 = -t51 * mrSges(6,1) - t52 * mrSges(6,2) - t31;
t116 = -mrSges(6,3) - mrSges(7,2);
t27 = t52 * qJD(5) + t77 * qJDD(3) - t122 * t56;
t34 = t68 * mrSges(6,1) - t52 * mrSges(6,3);
t7 = m(6) * t115 - t50 * mrSges(6,2) + t114 * t51 + t116 * t27 - t68 * t34 + t108;
t97 = t122 * t18 - t77 * t14;
t125 = m(7) * (-t50 * pkin(5) - t64 * qJ(6) + t52 * t30 + qJDD(6) - t97);
t28 = -t51 * qJD(5) + t122 * qJDD(3) + t77 * t56;
t33 = -t68 * mrSges(6,2) - t51 * mrSges(6,3);
t36 = -t51 * mrSges(7,2) + t68 * mrSges(7,3);
t8 = m(6) * t97 - t125 + (t33 + t36) * t68 + t114 * t52 + (mrSges(6,1) + mrSges(7,1)) * t50 + t116 * t28;
t88 = -(t78 * t61 + t80 * t62) * qJD(1) + t122 * t7 - t77 * t8 + m(5) * (t56 * pkin(3) + t107 + t89) - t57 * mrSges(5,3);
t84 = t119 * t56 + m(4) * (t107 + t96) + t59 * t113 + t60 * t112 + t57 * mrSges(4,2) + t88;
t133 = m(3) * (t83 * pkin(1) - t96) - t84;
t104 = -t80 * g(3) + t78 * t40;
t131 = t82 * pkin(3) - qJDD(3) * qJ(4) + qJD(3) * t132 + t53 * t113 - t104;
t129 = -m(2) - m(3);
t123 = t78 * g(3);
t120 = mrSges(2,1) - mrSges(3,2);
t118 = -mrSges(2,2) + mrSges(3,3);
t117 = -mrSges(4,3) - mrSges(5,1);
t54 = (-mrSges(5,2) * t78 - mrSges(5,3) * t80) * qJD(1);
t102 = qJD(1) * (-t54 - (mrSges(4,1) * t78 + mrSges(4,2) * t80) * qJD(1));
t91 = m(5) * (-qJDD(3) * pkin(3) - t123 + t93) + t122 * t8 + t77 * t7;
t3 = m(4) * (t121 + t123) + t117 * t57 + t119 * qJDD(3) + (t59 - t61) * qJD(3) + t80 * t102 - t91;
t86 = -t56 * pkin(4) + qJD(3) * t63 - t76 * t124 - t131;
t92 = -t28 * mrSges(7,3) - t52 * t35 + m(7) * (-0.2e1 * qJD(6) * t52 + (t51 * t68 - t28) * qJ(6) + (t52 * t68 + t27) * pkin(5) + t86) + t27 * mrSges(7,1) + t51 * t36;
t87 = m(6) * t86 + t27 * mrSges(6,1) + t28 * mrSges(6,2) + t51 * t33 + t52 * t34 + t92;
t85 = -m(5) * t131 + t87;
t5 = t85 + t117 * t56 + (-mrSges(4,2) + mrSges(5,3)) * qJDD(3) + (-t60 + t62) * qJD(3) + t78 * t102 + m(4) * t104;
t106 = -t78 * t3 + t80 * t5;
t95 = -m(3) * (-qJDD(1) * pkin(1) + t94) - t80 * t3 - t78 * t5;
t2 = m(2) * t101 + t118 * qJDD(1) - t120 * t83 - t133;
t1 = m(2) * t105 + t120 * qJDD(1) + t118 * t83 + t95;
t4 = [-m(1) * g(1) - t79 * t1 + t81 * t2, t2, -m(3) * g(3) + t106, t5, -t56 * mrSges(5,2) + t88, t7, -t27 * mrSges(7,2) - t51 * t31 + t108; -m(1) * g(2) + t81 * t1 + t79 * t2, t1, -t83 * mrSges(3,2) - qJDD(1) * mrSges(3,3) + t133, t3, t56 * mrSges(5,1) - qJDD(3) * mrSges(5,3) - qJD(3) * t62 + t54 * t113 - t85, t8, t92; (-m(1) + t129) * g(3) + t106, t129 * g(3) + t106, qJDD(1) * mrSges(3,2) - t83 * mrSges(3,3) - t95, t84, t57 * mrSges(5,1) + qJDD(3) * mrSges(5,2) + qJD(3) * t61 + t54 * t112 + t91, t87, -t50 * mrSges(7,1) + t28 * mrSges(7,2) + t52 * t31 - t68 * t36 + t125;];
f_new  = t4;
