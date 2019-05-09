% Calculate vector of cutting forces with Newton-Euler
% S6RPRPRP2
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-05-05 17:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRPRP2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP2_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP2_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP2_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP2_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP2_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:32:35
% EndTime: 2019-05-05 17:32:39
% DurationCPUTime: 1.80s
% Computational Cost: add. (20992->175), mult. (44933->225), div. (0->0), fcn. (29038->10), ass. (0->86)
t130 = -2 * qJD(4);
t91 = sin(qJ(1));
t93 = cos(qJ(1));
t109 = t91 * g(1) - t93 * g(2);
t69 = qJDD(1) * pkin(1) + t109;
t104 = -t93 * g(1) - t91 * g(2);
t95 = qJD(1) ^ 2;
t71 = -t95 * pkin(1) + t104;
t86 = sin(pkin(9));
t88 = cos(pkin(9));
t118 = t86 * t69 + t88 * t71;
t46 = -t95 * pkin(2) + qJDD(1) * pkin(7) + t118;
t84 = -g(3) + qJDD(2);
t90 = sin(qJ(3));
t92 = cos(qJ(3));
t106 = -t90 * t46 + t92 * t84;
t114 = qJD(1) * qJD(3);
t108 = t92 * t114;
t72 = t90 * qJDD(1) + t108;
t27 = (-t72 + t108) * qJ(4) + (t90 * t92 * t95 + qJDD(3)) * pkin(3) + t106;
t119 = t92 * t46 + t90 * t84;
t73 = t92 * qJDD(1) - t90 * t114;
t117 = qJD(1) * t90;
t74 = qJD(3) * pkin(3) - qJ(4) * t117;
t83 = t92 ^ 2;
t28 = -t83 * t95 * pkin(3) + t73 * qJ(4) - qJD(3) * t74 + t119;
t85 = sin(pkin(10));
t87 = cos(pkin(10));
t64 = (t85 * t92 + t87 * t90) * qJD(1);
t129 = t64 * t130 + t87 * t27 - t85 * t28;
t63 = (t85 * t90 - t87 * t92) * qJD(1);
t105 = t88 * t69 - t86 * t71;
t101 = -qJDD(1) * pkin(2) - t105;
t75 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t117;
t116 = qJD(1) * t92;
t76 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t116;
t125 = cos(qJ(5));
t110 = t63 * t130 + t85 * t27 + t87 * t28;
t49 = t63 * pkin(4) - t64 * pkin(8);
t94 = qJD(3) ^ 2;
t21 = -t94 * pkin(4) + qJDD(3) * pkin(8) - t63 * t49 + t110;
t53 = -t85 * t72 + t87 * t73;
t54 = t87 * t72 + t85 * t73;
t97 = -t73 * pkin(3) + qJDD(4) + t74 * t117 + (-qJ(4) * t83 - pkin(7)) * t95 + t101;
t23 = (qJD(3) * t63 - t54) * pkin(8) + (qJD(3) * t64 - t53) * pkin(4) + t97;
t89 = sin(qJ(5));
t122 = t125 * t21 + t89 * t23;
t55 = -t125 * qJD(3) + t89 * t64;
t56 = t89 * qJD(3) + t125 * t64;
t36 = t55 * pkin(5) - t56 * qJ(6);
t62 = qJD(5) + t63;
t43 = -t62 * mrSges(7,1) + t56 * mrSges(7,2);
t52 = qJDD(5) - t53;
t61 = t62 ^ 2;
t112 = m(7) * (-t61 * pkin(5) + t52 * qJ(6) + 0.2e1 * qJD(6) * t62 - t55 * t36 + t122) + t62 * t43 + t52 * mrSges(7,3);
t37 = t55 * mrSges(7,1) - t56 * mrSges(7,3);
t121 = -t55 * mrSges(6,1) - t56 * mrSges(6,2) - t37;
t123 = -mrSges(6,3) - mrSges(7,2);
t33 = t56 * qJD(5) - t125 * qJDD(3) + t89 * t54;
t42 = t62 * mrSges(6,1) - t56 * mrSges(6,3);
t12 = m(6) * t122 - t52 * mrSges(6,2) + t121 * t55 + t123 * t33 - t62 * t42 + t112;
t100 = t125 * t23 - t89 * t21;
t126 = m(7) * (-t52 * pkin(5) - t61 * qJ(6) + t56 * t36 + qJDD(6) - t100);
t34 = -t55 * qJD(5) + t89 * qJDD(3) + t125 * t54;
t40 = -t55 * mrSges(7,2) + t62 * mrSges(7,3);
t41 = -t62 * mrSges(6,2) - t55 * mrSges(6,3);
t14 = m(6) * t100 - t126 + (t41 + t40) * t62 + t121 * t56 + (mrSges(6,1) + mrSges(7,1)) * t52 + t123 * t34;
t57 = -qJD(3) * mrSges(5,2) - t63 * mrSges(5,3);
t58 = qJD(3) * mrSges(5,1) - t64 * mrSges(5,3);
t99 = -m(5) * t97 + t53 * mrSges(5,1) - t54 * mrSges(5,2) - t89 * t12 - t125 * t14 - t63 * t57 - t64 * t58;
t128 = (t90 * t75 - t92 * t76) * qJD(1) + m(4) * (-t95 * pkin(7) + t101) - t73 * mrSges(4,1) + t72 * mrSges(4,2) - t99;
t20 = -qJDD(3) * pkin(4) - t94 * pkin(8) + t64 * t49 - t129;
t111 = m(7) * (-0.2e1 * qJD(6) * t56 + (t55 * t62 - t34) * qJ(6) + (t56 * t62 + t33) * pkin(5) + t20) + t55 * t40 + t33 * mrSges(7,1);
t127 = m(6) * t20 + t33 * mrSges(6,1) + (t42 - t43) * t56 + (mrSges(6,2) - mrSges(7,3)) * t34 + t55 * t41 + t111;
t48 = t63 * mrSges(5,1) + t64 * mrSges(5,2);
t10 = m(5) * t129 + qJDD(3) * mrSges(5,1) - t54 * mrSges(5,3) + qJD(3) * t57 - t64 * t48 - t127;
t70 = (-mrSges(4,1) * t92 + mrSges(4,2) * t90) * qJD(1);
t9 = m(5) * t110 - qJDD(3) * mrSges(5,2) + t53 * mrSges(5,3) - qJD(3) * t58 + t125 * t12 - t89 * t14 - t63 * t48;
t6 = m(4) * t106 + qJDD(3) * mrSges(4,1) - t72 * mrSges(4,3) + qJD(3) * t76 + t87 * t10 - t70 * t117 + t85 * t9;
t7 = m(4) * t119 - qJDD(3) * mrSges(4,2) + t73 * mrSges(4,3) - qJD(3) * t75 - t85 * t10 + t70 * t116 + t87 * t9;
t113 = m(3) * t84 + t92 * t6 + t90 * t7;
t8 = m(3) * t105 + qJDD(1) * mrSges(3,1) - t95 * mrSges(3,2) - t128;
t3 = m(3) * t118 - t95 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t90 * t6 + t92 * t7;
t2 = m(2) * t104 - t95 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t88 * t3 - t86 * t8;
t1 = m(2) * t109 + qJDD(1) * mrSges(2,1) - t95 * mrSges(2,2) + t86 * t3 + t88 * t8;
t4 = [-m(1) * g(1) - t91 * t1 + t93 * t2, t2, t3, t7, t9, t12, -t33 * mrSges(7,2) - t55 * t37 + t112; -m(1) * g(2) + t93 * t1 + t91 * t2, t1, t8, t6, t10, t14, -t34 * mrSges(7,3) - t56 * t43 + t111; (-m(1) - m(2)) * g(3) + t113, -m(2) * g(3) + t113, t113, t128, -t99, t127, -t52 * mrSges(7,1) + t34 * mrSges(7,2) + t56 * t37 - t62 * t40 + t126;];
f_new  = t4;
