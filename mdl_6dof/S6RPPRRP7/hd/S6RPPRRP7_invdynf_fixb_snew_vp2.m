% Calculate vector of cutting forces with Newton-Euler
% S6RPPRRP7
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
% Datum: 2019-05-05 15:05
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPPRRP7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP7_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP7_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP7_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP7_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP7_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:03:14
% EndTime: 2019-05-05 15:03:18
% DurationCPUTime: 1.27s
% Computational Cost: add. (12983->163), mult. (29148->196), div. (0->0), fcn. (19795->8), ass. (0->83)
t121 = -pkin(1) - qJ(3);
t83 = sin(qJ(1));
t86 = cos(qJ(1));
t106 = t83 * g(1) - t86 * g(2);
t88 = qJD(1) ^ 2;
t95 = -t88 * qJ(2) + qJDD(2) - t106;
t130 = -(2 * qJD(1) * qJD(3)) + t121 * qJDD(1) + t95;
t79 = sin(pkin(9));
t76 = t79 ^ 2;
t80 = cos(pkin(9));
t117 = t80 ^ 2 + t76;
t105 = t117 * mrSges(4,3);
t101 = -t86 * g(1) - t83 * g(2);
t129 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t101;
t82 = sin(qJ(4));
t85 = cos(qJ(4));
t100 = t79 * t85 + t80 * t82;
t65 = t100 * qJD(1);
t99 = -t79 * t82 + t80 * t85;
t66 = t99 * qJD(1);
t115 = t66 * qJD(4);
t50 = -t100 * qJDD(1) - t115;
t109 = t79 * g(3) + t130 * t80;
t126 = pkin(3) * t88;
t35 = (-pkin(7) * qJDD(1) - t79 * t126) * t80 + t109;
t102 = -t80 * g(3) + t130 * t79;
t114 = qJDD(1) * t79;
t36 = -pkin(7) * t114 - t76 * t126 + t102;
t103 = t85 * t35 - t82 * t36;
t49 = t65 * pkin(4) - t66 * pkin(8);
t87 = qJD(4) ^ 2;
t17 = -qJDD(4) * pkin(4) - t87 * pkin(8) + t66 * t49 - t103;
t116 = t65 * qJD(4);
t51 = t99 * qJDD(1) - t116;
t81 = sin(qJ(5));
t84 = cos(qJ(5));
t55 = t81 * qJD(4) + t84 * t66;
t26 = -t55 * qJD(5) + t84 * qJDD(4) - t81 * t51;
t54 = t84 * qJD(4) - t81 * t66;
t27 = t54 * qJD(5) + t81 * qJDD(4) + t84 * t51;
t63 = qJD(5) + t65;
t39 = t63 * pkin(5) - t55 * qJ(6);
t40 = t63 * mrSges(7,1) - t55 * mrSges(7,3);
t53 = t54 ^ 2;
t110 = m(7) * (-t26 * pkin(5) - t53 * qJ(6) + t55 * t39 + qJDD(6) + t17) + t27 * mrSges(7,2) + t55 * t40;
t37 = -t63 * mrSges(7,2) + t54 * mrSges(7,3);
t38 = -t63 * mrSges(6,2) + t54 * mrSges(6,3);
t41 = t63 * mrSges(6,1) - t55 * mrSges(6,3);
t128 = m(6) * t17 + t27 * mrSges(6,2) - (t38 + t37) * t54 - (mrSges(6,1) + mrSges(7,1)) * t26 + t55 * t41 + t110;
t127 = -m(2) - m(3);
t125 = mrSges(4,2) * t80;
t124 = mrSges(2,1) - mrSges(3,2);
t122 = -mrSges(2,2) + mrSges(3,3);
t119 = t82 * t35 + t85 * t36;
t18 = -t87 * pkin(4) + qJDD(4) * pkin(8) - t65 * t49 + t119;
t94 = qJDD(3) + t129;
t89 = pkin(3) * t114 + (-t117 * pkin(7) + t121) * t88 + t94;
t21 = (-t51 + t116) * pkin(8) + (-t50 + t115) * pkin(4) + t89;
t120 = t84 * t18 + t81 * t21;
t104 = -t81 * t18 + t84 * t21;
t48 = qJDD(5) - t50;
t112 = m(7) * (-0.2e1 * qJD(6) * t55 + (t54 * t63 - t27) * qJ(6) + (t54 * t55 + t48) * pkin(5) + t104) + t63 * t37 + t48 * mrSges(7,1);
t30 = -t54 * mrSges(7,1) + t55 * mrSges(7,2);
t111 = m(7) * (-t53 * pkin(5) + t26 * qJ(6) + 0.2e1 * qJD(6) * t54 - t63 * t39 + t120) + t54 * t30 + t26 * mrSges(7,3);
t46 = t65 * mrSges(5,1) + t66 * mrSges(5,2);
t59 = -qJD(4) * mrSges(5,2) - t65 * mrSges(5,3);
t11 = m(5) * t103 + qJDD(4) * mrSges(5,1) - t51 * mrSges(5,3) + qJD(4) * t59 - t66 * t46 - t128;
t31 = -t54 * mrSges(6,1) + t55 * mrSges(6,2);
t10 = m(6) * t120 + t26 * mrSges(6,3) + t54 * t31 + (-t41 - t40) * t63 + (-mrSges(6,2) - mrSges(7,2)) * t48 + t111;
t60 = qJD(4) * mrSges(5,1) - t66 * mrSges(5,3);
t8 = m(6) * t104 + t48 * mrSges(6,1) + t63 * t38 + (-t31 - t30) * t55 + (-mrSges(6,3) - mrSges(7,3)) * t27 + t112;
t6 = m(5) * t119 - qJDD(4) * mrSges(5,2) + t50 * mrSges(5,3) - qJD(4) * t60 + t84 * t10 - t65 * t46 - t81 * t8;
t97 = -qJDD(1) * mrSges(4,3) - t88 * (mrSges(4,1) * t79 + t125);
t3 = m(4) * t109 + t85 * t11 + t82 * t6 + t97 * t80;
t4 = m(4) * t102 - t82 * t11 + t85 * t6 + t97 * t79;
t107 = -t79 * t3 + t80 * t4;
t96 = -m(3) * (-qJDD(1) * pkin(1) + t95) - t80 * t3 - t79 * t4;
t93 = -m(5) * t89 + t50 * mrSges(5,1) - t51 * mrSges(5,2) - t81 * t10 - t65 * t59 - t66 * t60 - t84 * t8;
t91 = m(4) * (t121 * t88 + t94) + mrSges(4,1) * t114 + qJDD(1) * t125 - t93;
t90 = m(3) * (t88 * pkin(1) - t129) - t91;
t5 = -t90 + (-t105 - t124) * t88 + t122 * qJDD(1) + m(2) * t101;
t1 = m(2) * t106 + t124 * qJDD(1) + t122 * t88 + t96;
t2 = [-m(1) * g(1) - t83 * t1 + t86 * t5, t5, -m(3) * g(3) + t107, t4, t6, t10, -t48 * mrSges(7,2) - t63 * t40 + t111; -m(1) * g(2) + t86 * t1 + t83 * t5, t1, -qJDD(1) * mrSges(3,3) + (-mrSges(3,2) + t105) * t88 + t90, t3, t11, t8, -t27 * mrSges(7,3) - t55 * t30 + t112; (-m(1) + t127) * g(3) + t107, t127 * g(3) + t107, qJDD(1) * mrSges(3,2) - t88 * mrSges(3,3) - t96, -t88 * t105 + t91, -t93, t128, -t26 * mrSges(7,1) - t54 * t37 + t110;];
f_new  = t2;
