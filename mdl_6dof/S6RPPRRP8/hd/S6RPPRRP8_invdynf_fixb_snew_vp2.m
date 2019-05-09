% Calculate vector of cutting forces with Newton-Euler
% S6RPPRRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2019-05-05 15:08
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPPRRP8_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP8_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP8_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP8_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP8_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP8_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP8_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:06:13
% EndTime: 2019-05-05 15:06:16
% DurationCPUTime: 1.29s
% Computational Cost: add. (12854->163), mult. (28729->196), div. (0->0), fcn. (19425->8), ass. (0->84)
t119 = -pkin(1) - qJ(3);
t82 = sin(qJ(1));
t84 = cos(qJ(1));
t104 = t82 * g(1) - t84 * g(2);
t86 = qJD(1) ^ 2;
t93 = -t86 * qJ(2) + qJDD(2) - t104;
t131 = -(2 * qJD(1) * qJD(3)) + t119 * qJDD(1) + t93;
t78 = sin(pkin(9));
t75 = t78 ^ 2;
t79 = cos(pkin(9));
t114 = t79 ^ 2 + t75;
t103 = t114 * mrSges(4,3);
t100 = -t84 * g(1) - t82 * g(2);
t130 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t100;
t81 = sin(qJ(4));
t83 = cos(qJ(4));
t99 = t78 * t83 + t79 * t81;
t62 = t99 * qJD(1);
t98 = -t78 * t81 + t79 * t83;
t63 = t98 * qJD(1);
t112 = t63 * qJD(4);
t47 = -t99 * qJDD(1) - t112;
t107 = t78 * g(3) + t131 * t79;
t126 = pkin(3) * t86;
t33 = (-pkin(7) * qJDD(1) - t78 * t126) * t79 + t107;
t101 = -t79 * g(3) + t131 * t78;
t111 = qJDD(1) * t78;
t34 = -pkin(7) * t111 - t75 * t126 + t101;
t102 = t83 * t33 - t81 * t34;
t46 = t62 * pkin(4) - t63 * pkin(8);
t85 = qJD(4) ^ 2;
t17 = -qJDD(4) * pkin(4) - t85 * pkin(8) + t63 * t46 - t102;
t125 = cos(qJ(5));
t113 = t62 * qJD(4);
t48 = t98 * qJDD(1) - t113;
t80 = sin(qJ(5));
t51 = t80 * qJD(4) + t125 * t63;
t23 = t51 * qJD(5) - t125 * qJDD(4) + t80 * t48;
t50 = -t125 * qJD(4) + t80 * t63;
t24 = -t50 * qJD(5) + t80 * qJDD(4) + t125 * t48;
t60 = qJD(5) + t62;
t35 = -t50 * mrSges(7,2) + t60 * mrSges(7,3);
t108 = m(7) * (-0.2e1 * qJD(6) * t51 + (t50 * t60 - t24) * qJ(6) + (t51 * t60 + t23) * pkin(5) + t17) + t23 * mrSges(7,1) + t50 * t35;
t36 = -t60 * mrSges(6,2) - t50 * mrSges(6,3);
t37 = t60 * mrSges(6,1) - t51 * mrSges(6,3);
t38 = -t60 * mrSges(7,1) + t51 * mrSges(7,2);
t129 = m(6) * t17 + t23 * mrSges(6,1) + (t37 - t38) * t51 + (mrSges(6,2) - mrSges(7,3)) * t24 + t50 * t36 + t108;
t128 = -m(2) - m(3);
t27 = t50 * pkin(5) - t51 * qJ(6);
t45 = qJDD(5) - t47;
t59 = t60 ^ 2;
t116 = t81 * t33 + t83 * t34;
t18 = -t85 * pkin(4) + qJDD(4) * pkin(8) - t62 * t46 + t116;
t92 = qJDD(3) + t130;
t87 = pkin(3) * t111 + (-t114 * pkin(7) + t119) * t86 + t92;
t20 = (-t48 + t113) * pkin(8) + (-t47 + t112) * pkin(4) + t87;
t95 = t125 * t20 - t80 * t18;
t127 = m(7) * (-t45 * pkin(5) - t59 * qJ(6) + t51 * t27 + qJDD(6) - t95);
t124 = mrSges(4,2) * t79;
t123 = mrSges(2,1) - mrSges(3,2);
t122 = -mrSges(2,2) + mrSges(3,3);
t120 = -mrSges(6,3) - mrSges(7,2);
t118 = t125 * t18 + t80 * t20;
t28 = t50 * mrSges(7,1) - t51 * mrSges(7,3);
t117 = -t50 * mrSges(6,1) - t51 * mrSges(6,2) - t28;
t109 = m(7) * (-t59 * pkin(5) + t45 * qJ(6) + 0.2e1 * qJD(6) * t60 - t50 * t27 + t118) + t60 * t38 + t45 * mrSges(7,3);
t11 = m(6) * t95 - t127 + (t36 + t35) * t60 + t117 * t51 + (mrSges(6,1) + mrSges(7,1)) * t45 + t120 * t24;
t43 = t62 * mrSges(5,1) + t63 * mrSges(5,2);
t56 = qJD(4) * mrSges(5,1) - t63 * mrSges(5,3);
t9 = m(6) * t118 - t45 * mrSges(6,2) + t117 * t50 + t120 * t23 - t60 * t37 + t109;
t6 = m(5) * t116 - qJDD(4) * mrSges(5,2) + t47 * mrSges(5,3) - qJD(4) * t56 - t80 * t11 + t125 * t9 - t62 * t43;
t55 = -qJD(4) * mrSges(5,2) - t62 * mrSges(5,3);
t7 = m(5) * t102 + qJDD(4) * mrSges(5,1) - t48 * mrSges(5,3) + qJD(4) * t55 - t63 * t43 - t129;
t96 = -qJDD(1) * mrSges(4,3) - t86 * (mrSges(4,1) * t78 + t124);
t3 = m(4) * t107 + t81 * t6 + t83 * t7 + t96 * t79;
t4 = m(4) * t101 + t83 * t6 - t81 * t7 + t96 * t78;
t105 = -t78 * t3 + t79 * t4;
t94 = -m(3) * (-qJDD(1) * pkin(1) + t93) - t79 * t3 - t78 * t4;
t91 = -m(5) * t87 + t47 * mrSges(5,1) - t48 * mrSges(5,2) - t125 * t11 - t62 * t55 - t63 * t56 - t80 * t9;
t89 = m(4) * (t119 * t86 + t92) + mrSges(4,1) * t111 + qJDD(1) * t124 - t91;
t88 = m(3) * (t86 * pkin(1) - t130) - t89;
t5 = (-t103 - t123) * t86 + t122 * qJDD(1) - t88 + m(2) * t100;
t1 = m(2) * t104 + t123 * qJDD(1) + t122 * t86 + t94;
t2 = [-m(1) * g(1) - t82 * t1 + t84 * t5, t5, -m(3) * g(3) + t105, t4, t6, t9, -t23 * mrSges(7,2) - t50 * t28 + t109; -m(1) * g(2) + t84 * t1 + t82 * t5, t1, -qJDD(1) * mrSges(3,3) + (-mrSges(3,2) + t103) * t86 + t88, t3, t7, t11, -t24 * mrSges(7,3) - t51 * t38 + t108; (-m(1) + t128) * g(3) + t105, t128 * g(3) + t105, qJDD(1) * mrSges(3,2) - t86 * mrSges(3,3) - t94, -t86 * t103 + t89, -t91, t129, -t45 * mrSges(7,1) + t24 * mrSges(7,2) + t51 * t28 - t60 * t35 + t127;];
f_new  = t2;
