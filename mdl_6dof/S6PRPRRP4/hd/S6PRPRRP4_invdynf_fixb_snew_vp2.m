% Calculate vector of cutting forces with Newton-Euler
% S6PRPRRP4
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-05-04 23:53
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRPRRP4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP4_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP4_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP4_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP4_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP4_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:47:21
% EndTime: 2019-05-04 23:47:26
% DurationCPUTime: 2.19s
% Computational Cost: add. (27766->163), mult. (60954->208), div. (0->0), fcn. (45126->12), ass. (0->91)
t95 = qJD(2) ^ 2;
t86 = cos(pkin(11));
t81 = t86 ^ 2;
t83 = sin(pkin(11));
t117 = t83 ^ 2 + t81;
t133 = t117 * mrSges(4,3);
t90 = sin(qJ(4));
t92 = cos(qJ(4));
t103 = t83 * t90 - t86 * t92;
t66 = t103 * qJD(2);
t84 = sin(pkin(10));
t87 = cos(pkin(10));
t72 = t84 * g(1) - t87 * g(2);
t88 = cos(pkin(6));
t126 = t72 * t88;
t73 = -t87 * g(1) - t84 * g(2);
t82 = -g(3) + qJDD(1);
t85 = sin(pkin(6));
t91 = sin(qJ(2));
t93 = cos(qJ(2));
t132 = (t82 * t85 + t126) * t93 - t91 * t73;
t104 = t83 * t92 + t86 * t90;
t67 = t104 * qJD(2);
t114 = t67 * qJD(4);
t55 = -t103 * qJDD(2) - t114;
t116 = pkin(8) * qJDD(2);
t113 = qJD(2) * qJD(3);
t63 = -t85 * t72 + t88 * t82;
t118 = -0.2e1 * t83 * t113 + t86 * t63;
t129 = pkin(3) * t95;
t125 = t85 * t91;
t109 = t82 * t125 + t91 * t126 + t93 * t73;
t48 = -t95 * pkin(2) + qJDD(2) * qJ(3) + t109;
t29 = (t86 * t129 - t116 - t48) * t83 + t118;
t110 = t83 * t63 + (0.2e1 * t113 + t48) * t86;
t30 = t86 * t116 - t81 * t129 + t110;
t107 = t92 * t29 - t90 * t30;
t54 = t66 * pkin(4) - t67 * pkin(9);
t94 = qJD(4) ^ 2;
t22 = -qJDD(4) * pkin(4) - t94 * pkin(9) + t67 * t54 - t107;
t128 = cos(qJ(5));
t115 = t66 * qJD(4);
t56 = t104 * qJDD(2) - t115;
t89 = sin(qJ(5));
t58 = t89 * qJD(4) + t128 * t67;
t33 = t58 * qJD(5) - t128 * qJDD(4) + t89 * t56;
t57 = -t128 * qJD(4) + t89 * t67;
t34 = -t57 * qJD(5) + t89 * qJDD(4) + t128 * t56;
t65 = qJD(5) + t66;
t44 = -t57 * mrSges(7,2) + t65 * mrSges(7,3);
t111 = m(7) * (-0.2e1 * qJD(6) * t58 + (t57 * t65 - t34) * qJ(6) + (t58 * t65 + t33) * pkin(5) + t22) + t33 * mrSges(7,1) + t57 * t44;
t45 = -t65 * mrSges(6,2) - t57 * mrSges(6,3);
t46 = t65 * mrSges(6,1) - t58 * mrSges(6,3);
t47 = -t65 * mrSges(7,1) + t58 * mrSges(7,2);
t131 = m(6) * t22 + t33 * mrSges(6,1) + (t46 - t47) * t58 + (mrSges(6,2) - mrSges(7,3)) * t34 + t57 * t45 + t111;
t121 = t90 * t29 + t92 * t30;
t23 = -t94 * pkin(4) + qJDD(4) * pkin(9) - t66 * t54 + t121;
t100 = qJDD(3) - t132;
t96 = (-pkin(3) * t86 - pkin(2)) * qJDD(2) + (-t117 * pkin(8) - qJ(3)) * t95 + t100;
t25 = (-t56 + t115) * pkin(9) + (-t55 + t114) * pkin(4) + t96;
t101 = t128 * t25 - t89 * t23;
t38 = t57 * pkin(5) - t58 * qJ(6);
t53 = qJDD(5) - t55;
t64 = t65 ^ 2;
t130 = m(7) * (-t53 * pkin(5) - t64 * qJ(6) + t58 * t38 + qJDD(6) - t101);
t106 = -t86 * mrSges(4,1) + t83 * mrSges(4,2);
t122 = t128 * t23 + t89 * t25;
t112 = m(7) * (-t64 * pkin(5) + t53 * qJ(6) + 0.2e1 * qJD(6) * t65 - t57 * t38 + t122) + t65 * t47 + t53 * mrSges(7,3);
t39 = t57 * mrSges(7,1) - t58 * mrSges(7,3);
t120 = -t57 * mrSges(6,1) - t58 * mrSges(6,2) - t39;
t123 = -mrSges(6,3) - mrSges(7,2);
t14 = m(6) * t122 - t53 * mrSges(6,2) + t120 * t57 + t123 * t33 - t65 * t46 + t112;
t16 = m(6) * t101 - t130 + (t45 + t44) * t65 + t120 * t58 + (mrSges(6,1) + mrSges(7,1)) * t53 + t123 * t34;
t61 = -qJD(4) * mrSges(5,2) - t66 * mrSges(5,3);
t62 = qJD(4) * mrSges(5,1) - t67 * mrSges(5,3);
t99 = -m(5) * t96 + t55 * mrSges(5,1) - t56 * mrSges(5,2) - t128 * t16 - t89 * t14 - t66 * t61 - t67 * t62;
t97 = m(4) * (-qJDD(2) * pkin(2) - t95 * qJ(3) + t100) - t99;
t10 = m(3) * t132 + (-mrSges(3,2) + t133) * t95 + (mrSges(3,1) - t106) * qJDD(2) - t97;
t127 = t10 * t93;
t102 = qJDD(2) * mrSges(4,3) + t95 * t106;
t51 = t66 * mrSges(5,1) + t67 * mrSges(5,2);
t11 = m(5) * t121 - qJDD(4) * mrSges(5,2) + t55 * mrSges(5,3) - qJD(4) * t62 + t128 * t14 - t89 * t16 - t66 * t51;
t12 = m(5) * t107 + qJDD(4) * mrSges(5,1) - t56 * mrSges(5,3) + qJD(4) * t61 - t67 * t51 - t131;
t7 = m(4) * t118 + t90 * t11 + t92 * t12 + (-m(4) * t48 - t102) * t83;
t8 = m(4) * t110 + t102 * t86 + t92 * t11 - t90 * t12;
t4 = m(3) * t109 - t95 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t83 * t7 + t86 * t8;
t6 = m(3) * t63 + t86 * t7 + t83 * t8;
t108 = m(2) * t82 + t4 * t125 + t85 * t127 + t88 * t6;
t2 = m(2) * t73 - t91 * t10 + t93 * t4;
t1 = m(2) * t72 - t85 * t6 + (t4 * t91 + t127) * t88;
t3 = [-m(1) * g(1) - t84 * t1 + t87 * t2, t2, t4, t8, t11, t14, -t33 * mrSges(7,2) - t57 * t39 + t112; -m(1) * g(2) + t87 * t1 + t84 * t2, t1, t10, t7, t12, t16, -t34 * mrSges(7,3) - t58 * t47 + t111; -m(1) * g(3) + t108, t108, t6, t106 * qJDD(2) - t95 * t133 + t97, -t99, t131, -t53 * mrSges(7,1) + t34 * mrSges(7,2) + t58 * t39 - t65 * t44 + t130;];
f_new  = t3;
