% Calculate vector of cutting forces with Newton-Euler
% S6RPPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
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
% Datum: 2019-05-05 14:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPPRPR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR5_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR5_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR5_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR5_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR5_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:21:00
% EndTime: 2019-05-05 14:21:03
% DurationCPUTime: 0.96s
% Computational Cost: add. (10535->152), mult. (21169->187), div. (0->0), fcn. (11881->8), ass. (0->78)
t80 = sin(qJ(1));
t83 = cos(qJ(1));
t112 = t80 * g(1) - t83 * g(2);
t85 = qJD(1) ^ 2;
t48 = -qJDD(1) * pkin(1) - t85 * qJ(2) + qJDD(2) - t112;
t117 = -qJDD(1) * qJ(3) - (2 * qJD(3) * qJD(1)) + t48;
t97 = -t83 * g(1) - t80 * g(2);
t116 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t97;
t115 = -m(3) - m(4);
t114 = mrSges(3,2) - mrSges(4,3);
t113 = -mrSges(4,2) - mrSges(3,3);
t82 = cos(qJ(4));
t111 = qJD(1) * t82;
t79 = sin(qJ(4));
t110 = t79 * qJD(1);
t109 = -m(2) + t115;
t108 = qJD(1) * qJD(4);
t105 = mrSges(2,1) - t114;
t100 = t82 * t108;
t101 = t79 * t108;
t62 = t79 * qJDD(1) + t100;
t63 = t82 * qJDD(1) - t101;
t89 = -t85 * pkin(7) - t117;
t25 = (-t63 + t101) * qJ(5) + (t62 + t100) * pkin(4) + t89;
t87 = qJDD(3) + (-pkin(1) - qJ(3)) * t85 + t116;
t37 = -qJDD(1) * pkin(7) + t87;
t103 = -t82 * g(3) + t79 * t37;
t60 = (pkin(4) * t79 - qJ(5) * t82) * qJD(1);
t84 = qJD(4) ^ 2;
t28 = -t84 * pkin(4) + qJDD(4) * qJ(5) - t60 * t110 + t103;
t76 = sin(pkin(9));
t77 = cos(pkin(9));
t55 = t77 * qJD(4) - t76 * t111;
t104 = 0.2e1 * qJD(5) * t55 + t76 * t25 + t77 * t28;
t61 = (mrSges(5,1) * t79 + mrSges(5,2) * t82) * qJD(1);
t64 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t110;
t96 = t79 * g(3) + t82 * t37;
t27 = -qJDD(4) * pkin(4) - t84 * qJ(5) + t60 * t111 + qJDD(5) - t96;
t42 = -mrSges(6,2) * t110 + t55 * mrSges(6,3);
t56 = t76 * qJD(4) + t77 * t111;
t43 = mrSges(6,1) * t110 - t56 * mrSges(6,3);
t44 = t77 * qJDD(4) - t76 * t63;
t45 = t76 * qJDD(4) + t77 * t63;
t78 = sin(qJ(6));
t81 = cos(qJ(6));
t32 = t78 * t55 + t81 * t56;
t19 = -t32 * qJD(6) + t81 * t44 - t78 * t45;
t31 = t81 * t55 - t78 * t56;
t20 = t31 * qJD(6) + t78 * t44 + t81 * t45;
t66 = qJD(6) + t110;
t29 = -t66 * mrSges(7,2) + t31 * mrSges(7,3);
t30 = t66 * mrSges(7,1) - t32 * mrSges(7,3);
t46 = pkin(5) * t110 - t56 * pkin(8);
t54 = t55 ^ 2;
t90 = t19 * mrSges(7,1) + t31 * t29 - m(7) * (-t44 * pkin(5) - t54 * pkin(8) + t56 * t46 + t27) - t20 * mrSges(7,2) - t32 * t30;
t86 = m(6) * t27 - t44 * mrSges(6,1) + t45 * mrSges(6,2) - t55 * t42 + t56 * t43 - t90;
t11 = m(5) * t96 + qJDD(4) * mrSges(5,1) - t63 * mrSges(5,3) + qJD(4) * t64 - t61 * t111 - t86;
t65 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t111;
t98 = -0.2e1 * qJD(5) * t56 + t77 * t25 - t76 * t28;
t14 = (t55 * t110 - t45) * pkin(8) + (t55 * t56 + t62) * pkin(5) + t98;
t15 = -t54 * pkin(5) + t44 * pkin(8) - t46 * t110 + t104;
t22 = -t31 * mrSges(7,1) + t32 * mrSges(7,2);
t59 = qJDD(6) + t62;
t12 = m(7) * (t81 * t14 - t78 * t15) - t20 * mrSges(7,3) + t59 * mrSges(7,1) - t32 * t22 + t66 * t29;
t13 = m(7) * (t78 * t14 + t81 * t15) + t19 * mrSges(7,3) - t59 * mrSges(7,2) + t31 * t22 - t66 * t30;
t33 = -t55 * mrSges(6,1) + t56 * mrSges(6,2);
t8 = m(6) * t98 + t62 * mrSges(6,1) - t45 * mrSges(6,3) + t42 * t110 + t81 * t12 + t78 * t13 - t56 * t33;
t9 = m(6) * t104 - t62 * mrSges(6,2) + t44 * mrSges(6,3) - t43 * t110 - t78 * t12 + t81 * t13 + t55 * t33;
t5 = m(5) * t103 - qJDD(4) * mrSges(5,2) - t62 * mrSges(5,3) - qJD(4) * t65 - t61 * t110 - t76 * t8 + t77 * t9;
t102 = -t79 * t11 + t82 * t5;
t99 = m(4) * t87 + qJDD(1) * mrSges(4,2) + t82 * t11 + t79 * t5;
t94 = m(3) * (t85 * pkin(1) - t116) - t99;
t92 = m(5) * t89 + t62 * mrSges(5,1) + t63 * mrSges(5,2) + t64 * t110 + t65 * t111 + t76 * t9 + t77 * t8;
t91 = m(4) * t117 - t92;
t88 = m(3) * t48 + t91;
t2 = m(2) * t112 + (-mrSges(2,2) - t113) * t85 + t105 * qJDD(1) - t88;
t1 = m(2) * t97 + (-mrSges(2,2) + mrSges(3,3)) * qJDD(1) - t105 * t85 - t94;
t3 = [-m(1) * g(1) + t83 * t1 - t80 * t2, t1, t115 * g(3) + t102, -m(4) * g(3) + t102, t5, t9, t13; -m(1) * g(2) + t80 * t1 + t83 * t2, t2, -qJDD(1) * mrSges(3,3) - t114 * t85 + t94, -t85 * mrSges(4,2) - qJDD(1) * mrSges(4,3) + t91, t11, t8, t12; (-m(1) + t109) * g(3) + t102, t109 * g(3) + t102, t114 * qJDD(1) + t113 * t85 + t88, -t85 * mrSges(4,3) + t99, t92, t86, -t90;];
f_new  = t3;
