% Calculate vector of cutting forces with Newton-Euler
% S6RPPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-05-05 14:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPPRRP2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP2_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP2_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP2_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP2_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP2_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:46:59
% EndTime: 2019-05-05 14:47:02
% DurationCPUTime: 1.63s
% Computational Cost: add. (19242->164), mult. (42436->203), div. (0->0), fcn. (28956->10), ass. (0->86)
t90 = qJD(1) ^ 2;
t82 = cos(pkin(10));
t78 = t82 ^ 2;
t80 = sin(pkin(10));
t113 = t80 ^ 2 + t78;
t126 = t113 * mrSges(4,3);
t85 = sin(qJ(4));
t87 = cos(qJ(4));
t97 = t80 * t85 - t82 * t87;
t61 = t97 * qJD(1);
t98 = t80 * t87 + t82 * t85;
t62 = t98 * qJD(1);
t110 = t62 * qJD(4);
t53 = -t97 * qJDD(1) - t110;
t112 = pkin(7) * qJDD(1);
t109 = qJD(1) * qJD(3);
t79 = -g(3) + qJDD(2);
t114 = -0.2e1 * t80 * t109 + t82 * t79;
t123 = pkin(3) * t90;
t86 = sin(qJ(1));
t88 = cos(qJ(1));
t104 = t86 * g(1) - t88 * g(2);
t67 = qJDD(1) * pkin(1) + t104;
t101 = -t88 * g(1) - t86 * g(2);
t68 = -t90 * pkin(1) + t101;
t81 = sin(pkin(9));
t83 = cos(pkin(9));
t115 = t81 * t67 + t83 * t68;
t46 = -t90 * pkin(2) + qJDD(1) * qJ(3) + t115;
t29 = (t82 * t123 - t112 - t46) * t80 + t114;
t105 = t80 * t79 + (0.2e1 * t109 + t46) * t82;
t32 = t82 * t112 - t78 * t123 + t105;
t103 = t87 * t29 - t85 * t32;
t52 = t61 * pkin(4) - t62 * pkin(8);
t89 = qJD(4) ^ 2;
t20 = -qJDD(4) * pkin(4) - t89 * pkin(8) + t62 * t52 - t103;
t122 = cos(qJ(5));
t111 = t61 * qJD(4);
t54 = t98 * qJDD(1) - t111;
t84 = sin(qJ(5));
t56 = t84 * qJD(4) + t122 * t62;
t30 = t56 * qJD(5) - t122 * qJDD(4) + t84 * t54;
t55 = -t122 * qJD(4) + t84 * t62;
t31 = -t55 * qJD(5) + t84 * qJDD(4) + t122 * t54;
t60 = qJD(5) + t61;
t40 = -t55 * mrSges(7,2) + t60 * mrSges(7,3);
t106 = m(7) * (-0.2e1 * qJD(6) * t56 + (t55 * t60 - t31) * qJ(6) + (t56 * t60 + t30) * pkin(5) + t20) + t30 * mrSges(7,1) + t55 * t40;
t41 = -t60 * mrSges(6,2) - t55 * mrSges(6,3);
t42 = t60 * mrSges(6,1) - t56 * mrSges(6,3);
t43 = -t60 * mrSges(7,1) + t56 * mrSges(7,2);
t125 = m(6) * t20 + t30 * mrSges(6,1) + (t42 - t43) * t56 + (mrSges(6,2) - mrSges(7,3)) * t31 + t55 * t41 + t106;
t36 = t55 * pkin(5) - t56 * qJ(6);
t51 = qJDD(5) - t53;
t59 = t60 ^ 2;
t118 = t85 * t29 + t87 * t32;
t21 = -t89 * pkin(4) + qJDD(4) * pkin(8) - t61 * t52 + t118;
t102 = t83 * t67 - t81 * t68;
t100 = qJDD(3) - t102;
t91 = (-pkin(3) * t82 - pkin(2)) * qJDD(1) + (-t113 * pkin(7) - qJ(3)) * t90 + t100;
t23 = (-t54 + t111) * pkin(8) + (-t53 + t110) * pkin(4) + t91;
t95 = t122 * t23 - t84 * t21;
t124 = m(7) * (-t51 * pkin(5) - t59 * qJ(6) + t56 * t36 + qJDD(6) - t95);
t120 = -mrSges(6,3) - mrSges(7,2);
t119 = t122 * t21 + t84 * t23;
t37 = t55 * mrSges(7,1) - t56 * mrSges(7,3);
t117 = -t55 * mrSges(6,1) - t56 * mrSges(6,2) - t37;
t49 = t61 * mrSges(5,1) + t62 * mrSges(5,2);
t57 = -qJD(4) * mrSges(5,2) - t61 * mrSges(5,3);
t10 = m(5) * t103 + qJDD(4) * mrSges(5,1) - t54 * mrSges(5,3) + qJD(4) * t57 - t62 * t49 - t125;
t107 = m(7) * (-t59 * pkin(5) + t51 * qJ(6) + 0.2e1 * qJD(6) * t60 - t55 * t36 + t119) + t60 * t43 + t51 * mrSges(7,3);
t12 = m(6) * t119 - t51 * mrSges(6,2) + t117 * t55 + t120 * t30 - t60 * t42 + t107;
t14 = m(6) * t95 - t124 + (t41 + t40) * t60 + t117 * t56 + (mrSges(6,1) + mrSges(7,1)) * t51 + t120 * t31;
t58 = qJD(4) * mrSges(5,1) - t62 * mrSges(5,3);
t9 = m(5) * t118 - qJDD(4) * mrSges(5,2) + t53 * mrSges(5,3) - qJD(4) * t58 + t122 * t12 - t84 * t14 - t61 * t49;
t99 = -t82 * mrSges(4,1) + t80 * mrSges(4,2);
t96 = qJDD(1) * mrSges(4,3) + t90 * t99;
t6 = m(4) * t114 + t85 * t9 + t87 * t10 + (-m(4) * t46 - t96) * t80;
t7 = m(4) * t105 - t85 * t10 + t96 * t82 + t87 * t9;
t108 = m(3) * t79 + t82 * t6 + t80 * t7;
t94 = -m(5) * t91 + t53 * mrSges(5,1) - t54 * mrSges(5,2) - t84 * t12 - t122 * t14 - t61 * t57 - t62 * t58;
t92 = m(4) * (-qJDD(1) * pkin(2) - t90 * qJ(3) + t100) - t94;
t8 = m(3) * t102 + (-mrSges(3,2) + t126) * t90 + (mrSges(3,1) - t99) * qJDD(1) - t92;
t3 = m(3) * t115 - t90 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t80 * t6 + t82 * t7;
t2 = m(2) * t101 - t90 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t83 * t3 - t81 * t8;
t1 = m(2) * t104 + qJDD(1) * mrSges(2,1) - t90 * mrSges(2,2) + t81 * t3 + t83 * t8;
t4 = [-m(1) * g(1) - t86 * t1 + t88 * t2, t2, t3, t7, t9, t12, -t30 * mrSges(7,2) - t55 * t37 + t107; -m(1) * g(2) + t88 * t1 + t86 * t2, t1, t8, t6, t10, t14, -t31 * mrSges(7,3) - t56 * t43 + t106; (-m(1) - m(2)) * g(3) + t108, -m(2) * g(3) + t108, t108, t99 * qJDD(1) - t90 * t126 + t92, -t94, t125, -t51 * mrSges(7,1) + t31 * mrSges(7,2) + t56 * t37 - t60 * t40 + t124;];
f_new  = t4;
