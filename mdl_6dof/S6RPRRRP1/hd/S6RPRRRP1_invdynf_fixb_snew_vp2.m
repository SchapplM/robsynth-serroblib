% Calculate vector of cutting forces with Newton-Euler
% S6RPRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-05-06 01:11
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRRP1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:07:35
% EndTime: 2019-05-06 01:07:39
% DurationCPUTime: 1.88s
% Computational Cost: add. (24388->176), mult. (47639->224), div. (0->0), fcn. (31386->10), ass. (0->87)
t89 = sin(qJ(4));
t90 = sin(qJ(3));
t92 = cos(qJ(4));
t93 = cos(qJ(3));
t64 = (t89 * t90 - t92 * t93) * qJD(1);
t91 = sin(qJ(1));
t94 = cos(qJ(1));
t109 = t91 * g(1) - t94 * g(2);
t67 = qJDD(1) * pkin(1) + t109;
t104 = -t94 * g(1) - t91 * g(2);
t95 = qJD(1) ^ 2;
t69 = -t95 * pkin(1) + t104;
t86 = sin(pkin(10));
t87 = cos(pkin(10));
t105 = t87 * t67 - t86 * t69;
t101 = -qJDD(1) * pkin(2) - t105;
t113 = qJD(1) * qJD(3);
t108 = t93 * t113;
t70 = t90 * qJDD(1) + t108;
t71 = t93 * qJDD(1) - t90 * t113;
t115 = qJD(1) * t90;
t72 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t115;
t114 = qJD(1) * t93;
t73 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t114;
t124 = cos(qJ(5));
t116 = t86 * t67 + t87 * t69;
t52 = -t95 * pkin(2) + qJDD(1) * pkin(7) + t116;
t85 = -g(3) + qJDD(2);
t106 = -t90 * t52 + t93 * t85;
t30 = (-t70 + t108) * pkin(8) + (t90 * t93 * t95 + qJDD(3)) * pkin(3) + t106;
t117 = t93 * t52 + t90 * t85;
t74 = qJD(3) * pkin(3) - pkin(8) * t115;
t84 = t93 ^ 2;
t31 = -t84 * t95 * pkin(3) + t71 * pkin(8) - qJD(3) * t74 + t117;
t120 = t89 * t30 + t92 * t31;
t65 = (t89 * t93 + t90 * t92) * qJD(1);
t54 = t64 * pkin(4) - t65 * pkin(9);
t83 = qJD(3) + qJD(4);
t81 = t83 ^ 2;
t82 = qJDD(3) + qJDD(4);
t21 = -t81 * pkin(4) + t82 * pkin(9) - t64 * t54 + t120;
t43 = -t65 * qJD(4) - t89 * t70 + t92 * t71;
t44 = -t64 * qJD(4) + t92 * t70 + t89 * t71;
t97 = -t71 * pkin(3) + t74 * t115 + (-pkin(8) * t84 - pkin(7)) * t95 + t101;
t23 = (t64 * t83 - t44) * pkin(9) + (t65 * t83 - t43) * pkin(4) + t97;
t88 = sin(qJ(5));
t121 = t124 * t21 + t88 * t23;
t55 = -t124 * t83 + t88 * t65;
t56 = t124 * t65 + t88 * t83;
t36 = t55 * pkin(5) - t56 * qJ(6);
t42 = qJDD(5) - t43;
t60 = qJD(5) + t64;
t48 = -t60 * mrSges(7,1) + t56 * mrSges(7,2);
t59 = t60 ^ 2;
t111 = m(7) * (-t59 * pkin(5) + t42 * qJ(6) + 0.2e1 * qJD(6) * t60 - t55 * t36 + t121) + t60 * t48 + t42 * mrSges(7,3);
t37 = t55 * mrSges(7,1) - t56 * mrSges(7,3);
t119 = -t55 * mrSges(6,1) - t56 * mrSges(6,2) - t37;
t122 = -mrSges(6,3) - mrSges(7,2);
t25 = t56 * qJD(5) - t124 * t82 + t88 * t44;
t47 = t60 * mrSges(6,1) - t56 * mrSges(6,3);
t12 = m(6) * t121 - t42 * mrSges(6,2) + t119 * t55 + t122 * t25 - t60 * t47 + t111;
t100 = t124 * t23 - t88 * t21;
t125 = m(7) * (-t42 * pkin(5) - t59 * qJ(6) + t56 * t36 + qJDD(6) - t100);
t26 = -t55 * qJD(5) + t124 * t44 + t88 * t82;
t45 = -t55 * mrSges(7,2) + t60 * mrSges(7,3);
t46 = -t60 * mrSges(6,2) - t55 * mrSges(6,3);
t14 = m(6) * t100 - t125 + (t46 + t45) * t60 + t119 * t56 + (mrSges(6,1) + mrSges(7,1)) * t42 + t122 * t26;
t57 = -t83 * mrSges(5,2) - t64 * mrSges(5,3);
t58 = t83 * mrSges(5,1) - t65 * mrSges(5,3);
t99 = -m(5) * t97 + t43 * mrSges(5,1) - t44 * mrSges(5,2) - t88 * t12 - t124 * t14 - t64 * t57 - t65 * t58;
t127 = (t90 * t72 - t93 * t73) * qJD(1) + m(4) * (-t95 * pkin(7) + t101) - t71 * mrSges(4,1) + t70 * mrSges(4,2) - t99;
t107 = t92 * t30 - t89 * t31;
t20 = -t82 * pkin(4) - t81 * pkin(9) + t65 * t54 - t107;
t110 = m(7) * (-0.2e1 * qJD(6) * t56 + (t55 * t60 - t26) * qJ(6) + (t56 * t60 + t25) * pkin(5) + t20) + t25 * mrSges(7,1) + t55 * t45;
t126 = m(6) * t20 + t25 * mrSges(6,1) + (t47 - t48) * t56 + (mrSges(6,2) - mrSges(7,3)) * t26 + t55 * t46 + t110;
t53 = t64 * mrSges(5,1) + t65 * mrSges(5,2);
t10 = m(5) * t107 + t82 * mrSges(5,1) - t44 * mrSges(5,3) - t65 * t53 + t83 * t57 - t126;
t68 = (-mrSges(4,1) * t93 + mrSges(4,2) * t90) * qJD(1);
t9 = m(5) * t120 - t82 * mrSges(5,2) + t43 * mrSges(5,3) + t124 * t12 - t88 * t14 - t64 * t53 - t83 * t58;
t6 = m(4) * t106 + qJDD(3) * mrSges(4,1) - t70 * mrSges(4,3) + qJD(3) * t73 + t92 * t10 - t68 * t115 + t89 * t9;
t7 = m(4) * t117 - qJDD(3) * mrSges(4,2) + t71 * mrSges(4,3) - qJD(3) * t72 - t89 * t10 + t68 * t114 + t92 * t9;
t112 = m(3) * t85 + t93 * t6 + t90 * t7;
t8 = m(3) * t105 + qJDD(1) * mrSges(3,1) - t95 * mrSges(3,2) - t127;
t3 = m(3) * t116 - t95 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t90 * t6 + t93 * t7;
t2 = m(2) * t104 - t95 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t87 * t3 - t86 * t8;
t1 = m(2) * t109 + qJDD(1) * mrSges(2,1) - t95 * mrSges(2,2) + t86 * t3 + t87 * t8;
t4 = [-m(1) * g(1) - t91 * t1 + t94 * t2, t2, t3, t7, t9, t12, -t25 * mrSges(7,2) - t55 * t37 + t111; -m(1) * g(2) + t94 * t1 + t91 * t2, t1, t8, t6, t10, t14, -t26 * mrSges(7,3) - t56 * t48 + t110; (-m(1) - m(2)) * g(3) + t112, -m(2) * g(3) + t112, t112, t127, -t99, t126, -t42 * mrSges(7,1) + t26 * mrSges(7,2) + t56 * t37 - t60 * t45 + t125;];
f_new  = t4;
