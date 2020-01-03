% Calculate vector of cutting forces with Newton-Euler
% S5RRPRP10
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% f_new [3x6]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRPRP10_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP10_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP10_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP10_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP10_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP10_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP10_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:09:27
% EndTime: 2019-12-31 20:09:28
% DurationCPUTime: 0.71s
% Computational Cost: add. (4314->160), mult. (8804->191), div. (0->0), fcn. (4607->6), ass. (0->72)
t121 = -2 * qJD(3);
t78 = cos(qJ(2));
t103 = qJD(1) * t78;
t75 = sin(qJ(2));
t54 = (-pkin(2) * t78 - qJ(3) * t75) * qJD(1);
t80 = qJD(2) ^ 2;
t81 = qJD(1) ^ 2;
t76 = sin(qJ(1));
t79 = cos(qJ(1));
t91 = -t79 * g(1) - t76 * g(2);
t46 = -t81 * pkin(1) + qJDD(1) * pkin(6) + t91;
t95 = -t75 * g(3) + t78 * t46;
t118 = t80 * pkin(2) - qJDD(2) * qJ(3) + (qJD(2) * t121) - t54 * t103 - t95;
t74 = sin(qJ(4));
t77 = cos(qJ(4));
t53 = t77 * qJD(2) - t74 * t103;
t101 = qJD(1) * qJD(2);
t94 = t75 * t101;
t58 = t78 * qJDD(1) - t94;
t31 = -t53 * qJD(4) - t74 * qJDD(2) - t77 * t58;
t52 = -t74 * qJD(2) - t103 * t77;
t32 = t52 * qJD(4) + t77 * qJDD(2) - t74 * t58;
t102 = t75 * qJD(1);
t66 = qJD(4) + t102;
t36 = -t66 * mrSges(6,2) + t52 * mrSges(6,3);
t37 = -t66 * mrSges(5,2) + t52 * mrSges(5,3);
t40 = t66 * mrSges(5,1) - t53 * mrSges(5,3);
t64 = pkin(3) * t102 - (qJD(2) * pkin(7));
t73 = t78 ^ 2;
t83 = -t73 * t81 * pkin(7) + t58 * pkin(3) + qJD(2) * t64 - t118;
t38 = t66 * pkin(4) - t53 * qJ(5);
t39 = t66 * mrSges(6,1) - t53 * mrSges(6,3);
t50 = t52 ^ 2;
t97 = m(6) * (-t31 * pkin(4) - t50 * qJ(5) + t53 * t38 + qJDD(5) + t83) + t32 * mrSges(6,2) + t53 * t39;
t84 = m(5) * t83 + t32 * mrSges(5,2) + (-t37 - t36) * t52 - (mrSges(5,1) + mrSges(6,1)) * t31 + t53 * t40 + t97;
t120 = -m(4) * t118 + t84;
t63 = mrSges(4,1) * t102 + (qJD(2) * mrSges(4,2));
t104 = (qJD(2) * mrSges(3,1)) - mrSges(3,3) * t102 - t63;
t111 = mrSges(3,1) - mrSges(4,2);
t113 = t81 * pkin(6);
t93 = t78 * t101;
t57 = t75 * qJDD(1) + t93;
t61 = -(qJD(2) * mrSges(3,2)) + mrSges(3,3) * t103;
t96 = t76 * g(1) - t79 * g(2);
t88 = -qJDD(1) * pkin(1) - t96;
t34 = -t52 * mrSges(6,1) + t53 * mrSges(6,2);
t35 = -t52 * mrSges(5,1) + t53 * mrSges(5,2);
t51 = qJDD(4) + t57;
t82 = pkin(2) * t94 + t102 * t121 + (-t57 - t93) * qJ(3) + t88;
t15 = -t64 * t102 + (-pkin(3) * t73 - pkin(6)) * t81 + (-pkin(2) - pkin(7)) * t58 + t82;
t106 = -t78 * g(3) - t75 * t46;
t24 = -qJDD(2) * pkin(2) - t80 * qJ(3) + t54 * t102 + qJDD(3) - t106;
t20 = (-t75 * t78 * t81 - qJDD(2)) * pkin(7) + (t57 - t93) * pkin(3) + t24;
t92 = -t74 * t15 + t77 * t20;
t99 = m(6) * (-0.2e1 * qJD(5) * t53 + (t52 * t66 - t32) * qJ(5) + (t52 * t53 + t51) * pkin(4) + t92) + t66 * t36 + t51 * mrSges(6,1);
t5 = m(5) * t92 + t51 * mrSges(5,1) + t66 * t37 + (-t35 - t34) * t53 + (-mrSges(5,3) - mrSges(6,3)) * t32 + t99;
t62 = -mrSges(4,1) * t103 - (qJD(2) * mrSges(4,3));
t108 = t77 * t15 + t74 * t20;
t98 = m(6) * (-t50 * pkin(4) + t31 * qJ(5) + 0.2e1 * qJD(5) * t52 - t66 * t38 + t108) + t52 * t34 + t31 * mrSges(6,3);
t9 = m(5) * t108 + t31 * mrSges(5,3) + t52 * t35 + (-t40 - t39) * t66 + (-mrSges(5,2) - mrSges(6,2)) * t51 + t98;
t89 = t74 * t5 - t77 * t9 - m(4) * (-t58 * pkin(2) - t113 + t82) - t62 * t103 + t57 * mrSges(4,3);
t119 = (t104 * t75 - t78 * t61) * qJD(1) - t111 * t58 + m(3) * (t88 - t113) + t57 * mrSges(3,2) - t89;
t55 = (mrSges(4,2) * t78 - mrSges(4,3) * t75) * qJD(1);
t105 = t55 + (-mrSges(3,1) * t78 + mrSges(3,2) * t75) * qJD(1);
t109 = mrSges(3,3) + mrSges(4,1);
t87 = -m(4) * t24 - t77 * t5 - t74 * t9;
t4 = m(3) * t106 - t109 * t57 + t111 * qJDD(2) + (t61 - t62) * qJD(2) - t105 * t102 + t87;
t8 = m(3) * t95 + t109 * t58 + (-mrSges(3,2) + mrSges(4,3)) * qJDD(2) - t104 * qJD(2) + t105 * t103 + t120;
t115 = t78 * t4 + t75 * t8;
t2 = m(2) * t96 + qJDD(1) * mrSges(2,1) - t81 * mrSges(2,2) - t119;
t1 = m(2) * t91 - t81 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t75 * t4 + t78 * t8;
t3 = [-m(1) * g(1) + t79 * t1 - t76 * t2, t1, t8, t58 * mrSges(4,2) - t102 * t63 - t89, t9, -t51 * mrSges(6,2) - t66 * t39 + t98; -m(1) * g(2) + t76 * t1 + t79 * t2, t2, t4, -t58 * mrSges(4,1) - qJDD(2) * mrSges(4,3) - qJD(2) * t63 - t55 * t103 - t120, t5, -t32 * mrSges(6,3) - t53 * t34 + t99; (-m(1) - m(2)) * g(3) + t115, -m(2) * g(3) + t115, t119, t57 * mrSges(4,1) + qJDD(2) * mrSges(4,2) + qJD(2) * t62 + t55 * t102 - t87, t84, -t31 * mrSges(6,1) - t52 * t36 + t97;];
f_new = t3;
