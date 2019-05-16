% Calculate vector of cutting forces with Newton-Euler
% S6RPPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-05-05 16:07
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPPRRR7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR7_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR7_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR7_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR7_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR7_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:05:17
% EndTime: 2019-05-05 16:05:23
% DurationCPUTime: 2.48s
% Computational Cost: add. (30234->167), mult. (69431->211), div. (0->0), fcn. (51309->10), ass. (0->92)
t88 = sin(qJ(1));
t92 = cos(qJ(1));
t113 = t88 * g(1) - t92 * g(2);
t93 = qJD(1) ^ 2;
t102 = -t93 * qJ(2) + qJDD(2) - t113;
t123 = -pkin(1) - qJ(3);
t130 = -(2 * qJD(1) * qJD(3)) + t123 * qJDD(1) + t102;
t83 = sin(pkin(10));
t79 = t83 ^ 2;
t84 = cos(pkin(10));
t120 = t84 ^ 2 + t79;
t112 = t120 * mrSges(4,3);
t109 = -t92 * g(1) - t88 * g(2);
t129 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t109;
t128 = -m(2) - m(3);
t127 = pkin(3) * t93;
t126 = mrSges(4,2) * t84;
t125 = mrSges(2,1) - mrSges(3,2);
t124 = -mrSges(2,2) + mrSges(3,3);
t116 = t83 * g(3) + t130 * t84;
t38 = (-pkin(7) * qJDD(1) - t83 * t127) * t84 + t116;
t110 = -g(3) * t84 + t130 * t83;
t118 = qJDD(1) * t83;
t41 = -pkin(7) * t118 - t79 * t127 + t110;
t87 = sin(qJ(4));
t91 = cos(qJ(4));
t111 = t91 * t38 - t87 * t41;
t107 = -t83 * t91 - t84 * t87;
t66 = t107 * qJD(1);
t119 = t66 * qJD(4);
t106 = -t83 * t87 + t84 * t91;
t54 = t106 * qJDD(1) + t119;
t67 = t106 * qJD(1);
t18 = (-t54 + t119) * pkin(8) + (t66 * t67 + qJDD(4)) * pkin(4) + t111;
t121 = t87 * t38 + t91 * t41;
t53 = -t67 * qJD(4) + t107 * qJDD(1);
t61 = qJD(4) * pkin(4) - pkin(8) * t67;
t65 = t66 ^ 2;
t20 = -pkin(4) * t65 + pkin(8) * t53 - qJD(4) * t61 + t121;
t86 = sin(qJ(5));
t90 = cos(qJ(5));
t122 = t86 * t18 + t90 * t20;
t104 = -mrSges(4,3) * qJDD(1) - t93 * (mrSges(4,1) * t83 + t126);
t51 = -t66 * mrSges(5,1) + mrSges(5,2) * t67;
t59 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t66;
t47 = t66 * t90 - t67 * t86;
t48 = t66 * t86 + t67 * t90;
t34 = -pkin(5) * t47 - pkin(9) * t48;
t81 = qJD(4) + qJD(5);
t77 = t81 ^ 2;
t78 = qJDD(4) + qJDD(5);
t15 = -pkin(5) * t77 + pkin(9) * t78 + t47 * t34 + t122;
t27 = -t48 * qJD(5) + t53 * t90 - t54 * t86;
t28 = t47 * qJD(5) + t53 * t86 + t54 * t90;
t101 = qJDD(3) + t129;
t98 = pkin(3) * t118 + (-t120 * pkin(7) + t123) * t93 + t101;
t95 = -t53 * pkin(4) - t65 * pkin(8) + t67 * t61 + t98;
t16 = (-t47 * t81 - t28) * pkin(9) + (t48 * t81 - t27) * pkin(5) + t95;
t85 = sin(qJ(6));
t89 = cos(qJ(6));
t39 = -t48 * t85 + t81 * t89;
t22 = t39 * qJD(6) + t28 * t89 + t78 * t85;
t26 = qJDD(6) - t27;
t40 = t48 * t89 + t81 * t85;
t29 = -mrSges(7,1) * t39 + mrSges(7,2) * t40;
t46 = qJD(6) - t47;
t30 = -mrSges(7,2) * t46 + mrSges(7,3) * t39;
t12 = m(7) * (-t15 * t85 + t16 * t89) - t22 * mrSges(7,3) + t26 * mrSges(7,1) - t40 * t29 + t46 * t30;
t21 = -t40 * qJD(6) - t28 * t85 + t78 * t89;
t31 = mrSges(7,1) * t46 - mrSges(7,3) * t40;
t13 = m(7) * (t15 * t89 + t16 * t85) + t21 * mrSges(7,3) - t26 * mrSges(7,2) + t39 * t29 - t46 * t31;
t33 = -mrSges(6,1) * t47 + mrSges(6,2) * t48;
t43 = mrSges(6,1) * t81 - t48 * mrSges(6,3);
t8 = m(6) * t122 - t78 * mrSges(6,2) + t27 * mrSges(6,3) - t85 * t12 + t89 * t13 + t47 * t33 - t81 * t43;
t108 = t18 * t90 - t20 * t86;
t42 = -mrSges(6,2) * t81 + t47 * mrSges(6,3);
t99 = m(7) * (-pkin(5) * t78 - pkin(9) * t77 + t48 * t34 - t108) - t21 * mrSges(7,1) + t22 * mrSges(7,2) - t39 * t30 + t40 * t31;
t9 = m(6) * t108 + t78 * mrSges(6,1) - t28 * mrSges(6,3) - t48 * t33 + t81 * t42 - t99;
t5 = m(5) * t111 + qJDD(4) * mrSges(5,1) - t54 * mrSges(5,3) + qJD(4) * t59 - t67 * t51 + t86 * t8 + t90 * t9;
t60 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t67;
t6 = m(5) * t121 - qJDD(4) * mrSges(5,2) + t53 * mrSges(5,3) - qJD(4) * t60 + t66 * t51 + t90 * t8 - t86 * t9;
t3 = m(4) * t116 + t104 * t84 + t91 * t5 + t87 * t6;
t4 = m(4) * t110 + t104 * t83 - t87 * t5 + t91 * t6;
t114 = -t83 * t3 + t84 * t4;
t103 = -m(3) * (-qJDD(1) * pkin(1) + t102) - t84 * t3 - t83 * t4;
t100 = -m(6) * t95 + t27 * mrSges(6,1) - t28 * mrSges(6,2) - t89 * t12 - t85 * t13 + t47 * t42 - t48 * t43;
t97 = -m(5) * t98 + t53 * mrSges(5,1) - t54 * mrSges(5,2) + t66 * t59 - t67 * t60 + t100;
t96 = m(4) * (t123 * t93 + t101) + mrSges(4,1) * t118 + qJDD(1) * t126 - t97;
t94 = m(3) * (t93 * pkin(1) - t129) - t96;
t7 = (-t112 - t125) * t93 + t124 * qJDD(1) + m(2) * t109 - t94;
t1 = m(2) * t113 + t125 * qJDD(1) + t124 * t93 + t103;
t2 = [-m(1) * g(1) - t1 * t88 + t7 * t92, t7, -m(3) * g(3) + t114, t4, t6, t8, t13; -m(1) * g(2) + t1 * t92 + t7 * t88, t1, -qJDD(1) * mrSges(3,3) + t94 + (-mrSges(3,2) + t112) * t93, t3, t5, t9, t12; (-m(1) + t128) * g(3) + t114, t128 * g(3) + t114, qJDD(1) * mrSges(3,2) - t93 * mrSges(3,3) - t103, -t93 * t112 + t96, -t97, -t100, t99;];
f_new  = t2;
