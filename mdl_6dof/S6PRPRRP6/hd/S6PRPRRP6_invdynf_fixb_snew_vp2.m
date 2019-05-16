% Calculate vector of cutting forces with Newton-Euler
% S6PRPRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-05-05 00:03
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRPRRP6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP6_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP6_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP6_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP6_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP6_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:59:17
% EndTime: 2019-05-04 23:59:19
% DurationCPUTime: 0.95s
% Computational Cost: add. (10267->146), mult. (18726->182), div. (0->0), fcn. (11621->10), ass. (0->80)
t76 = sin(pkin(6));
t81 = sin(qJ(2));
t115 = t76 * t81;
t75 = sin(pkin(10));
t77 = cos(pkin(10));
t62 = t75 * g(1) - t77 * g(2);
t78 = cos(pkin(6));
t116 = t62 * t78;
t63 = -t77 * g(1) - t75 * g(2);
t74 = -g(3) + qJDD(1);
t83 = cos(qJ(2));
t100 = t74 * t115 + t81 * t116 + t83 * t63;
t123 = -qJDD(2) * qJ(3) - (2 * qJD(3) * qJD(2)) - t100;
t122 = -t81 * t63 + (t74 * t76 + t116) * t83;
t82 = cos(qJ(4));
t106 = qJD(2) * t82;
t80 = sin(qJ(4));
t59 = (pkin(4) * t80 - pkin(9) * t82) * qJD(2);
t84 = qJD(4) ^ 2;
t120 = -pkin(2) - pkin(8);
t85 = qJD(2) ^ 2;
t86 = -t85 * qJ(3) + qJDD(3) - t122;
t28 = t120 * qJDD(2) + t86;
t46 = -t76 * t62 + t78 * t74;
t96 = t82 * t28 - t80 * t46;
t21 = -qJDD(4) * pkin(4) - t84 * pkin(9) + t59 * t106 - t96;
t117 = cos(qJ(5));
t79 = sin(qJ(5));
t57 = t79 * qJD(4) + t117 * t106;
t104 = qJD(2) * qJD(4);
t98 = t80 * t104;
t61 = t82 * qJDD(2) - t98;
t33 = t57 * qJD(5) - t117 * qJDD(4) + t79 * t61;
t56 = -t117 * qJD(4) + t79 * t106;
t34 = -t56 * qJD(5) + t79 * qJDD(4) + t117 * t61;
t105 = t80 * qJD(2);
t68 = qJD(5) + t105;
t44 = -t56 * mrSges(7,2) + t68 * mrSges(7,3);
t101 = m(7) * (-0.2e1 * qJD(6) * t57 + (t56 * t68 - t34) * qJ(6) + (t57 * t68 + t33) * pkin(5) + t21) + t33 * mrSges(7,1) + t56 * t44;
t41 = -t68 * mrSges(6,2) - t56 * mrSges(6,3);
t42 = t68 * mrSges(6,1) - t57 * mrSges(6,3);
t43 = -t68 * mrSges(7,1) + t57 * mrSges(7,2);
t121 = m(6) * t21 + t33 * mrSges(6,1) + (t42 - t43) * t57 + (mrSges(6,2) - mrSges(7,3)) * t34 + t56 * t41 + t101;
t36 = t56 * pkin(5) - t57 * qJ(6);
t97 = t82 * t104;
t60 = -t80 * qJDD(2) - t97;
t53 = qJDD(5) - t60;
t67 = t68 ^ 2;
t108 = t80 * t28 + t82 * t46;
t22 = -t84 * pkin(4) + qJDD(4) * pkin(9) - t59 * t105 + t108;
t89 = t120 * t85 - t123;
t24 = (-t61 + t98) * pkin(9) + (-t60 + t97) * pkin(4) + t89;
t92 = t117 * t24 - t79 * t22;
t119 = m(7) * (-t53 * pkin(5) - t67 * qJ(6) + t57 * t36 + qJDD(6) - t92);
t113 = (-mrSges(3,2) + mrSges(4,3));
t114 = mrSges(3,1) - mrSges(4,2);
t110 = t117 * t22 + t79 * t24;
t102 = m(7) * (-t67 * pkin(5) + t53 * qJ(6) + 0.2e1 * qJD(6) * t68 - t56 * t36 + t110) + t68 * t43 + t53 * mrSges(7,3);
t37 = t56 * mrSges(7,1) - t57 * mrSges(7,3);
t109 = -t56 * mrSges(6,1) - t57 * mrSges(6,2) - t37;
t111 = -mrSges(6,3) - mrSges(7,2);
t13 = m(6) * t110 - t53 * mrSges(6,2) + t109 * t56 + t111 * t33 - t68 * t42 + t102;
t15 = m(6) * t92 - t119 + (t41 + t44) * t68 + t109 * t57 + (mrSges(6,1) + mrSges(7,1)) * t53 + t111 * t34;
t58 = (mrSges(5,1) * t80 + mrSges(5,2) * t82) * qJD(2);
t65 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t106;
t10 = m(5) * t108 - qJDD(4) * mrSges(5,2) + t60 * mrSges(5,3) - qJD(4) * t65 - t58 * t105 + t117 * t13 - t79 * t15;
t64 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t105;
t11 = m(5) * t96 + qJDD(4) * mrSges(5,1) - t61 * mrSges(5,3) + qJD(4) * t64 - t58 * t106 - t121;
t91 = -m(4) * (-qJDD(2) * pkin(2) + t86) - t80 * t10 - t82 * t11;
t4 = m(3) * t122 + t114 * qJDD(2) + (t113 * t85) + t91;
t118 = t4 * t83;
t94 = m(4) * t46 + t82 * t10 - t80 * t11;
t6 = m(3) * t46 + t94;
t90 = m(5) * t89 - t60 * mrSges(5,1) + t61 * mrSges(5,2) + t64 * t105 + t65 * t106 + t117 * t15 + t79 * t13;
t87 = -m(4) * (t85 * pkin(2) + t123) + t90;
t8 = m(3) * t100 + t113 * qJDD(2) - t114 * t85 + t87;
t99 = m(2) * t74 + t8 * t115 + t76 * t118 + t78 * t6;
t2 = m(2) * t63 - t81 * t4 + t83 * t8;
t1 = m(2) * t62 - t76 * t6 + (t8 * t81 + t118) * t78;
t3 = [-m(1) * g(1) - t75 * t1 + t77 * t2, t2, t8, t94, t10, t13, -t33 * mrSges(7,2) - t56 * t37 + t102; -m(1) * g(2) + t77 * t1 + t75 * t2, t1, t4, -(t85 * mrSges(4,2)) - qJDD(2) * mrSges(4,3) - t87, t11, t15, -t34 * mrSges(7,3) - t57 * t43 + t101; -m(1) * g(3) + t99, t99, t6, qJDD(2) * mrSges(4,2) - t85 * mrSges(4,3) - t91, t90, t121, -t53 * mrSges(7,1) + t34 * mrSges(7,2) + t57 * t37 - t68 * t44 + t119;];
f_new  = t3;
