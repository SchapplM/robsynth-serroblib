% Calculate vector of cutting forces with Newton-Euler
% S6RPPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
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
% Datum: 2019-05-05 13:57
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPPRPR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:54:31
% EndTime: 2019-05-05 13:54:37
% DurationCPUTime: 3.00s
% Computational Cost: add. (39113->166), mult. (90445->219), div. (0->0), fcn. (64094->12), ass. (0->93)
t127 = cos(qJ(4));
t88 = sin(pkin(10));
t91 = cos(pkin(10));
t94 = sin(qJ(4));
t130 = -t91 * t127 + t88 * t94;
t99 = qJD(1) ^ 2;
t85 = t91 ^ 2;
t122 = t88 ^ 2 + t85;
t129 = t122 * mrSges(4,3);
t128 = pkin(3) * t99;
t121 = pkin(7) * qJDD(1);
t118 = qJD(1) * qJD(3);
t86 = -g(3) + qJDD(2);
t123 = -0.2e1 * t88 * t118 + t91 * t86;
t95 = sin(qJ(1));
t97 = cos(qJ(1));
t113 = t95 * g(1) - t97 * g(2);
t74 = qJDD(1) * pkin(1) + t113;
t110 = -t97 * g(1) - t95 * g(2);
t75 = -t99 * pkin(1) + t110;
t89 = sin(pkin(9));
t92 = cos(pkin(9));
t124 = t89 * t74 + t92 * t75;
t52 = -t99 * pkin(2) + qJDD(1) * qJ(3) + t124;
t37 = (t91 * t128 - t121 - t52) * t88 + t123;
t115 = t88 * t86 + (0.2e1 * t118 + t52) * t91;
t38 = t91 * t121 - t85 * t128 + t115;
t125 = t127 * t38 + t94 * t37;
t68 = t130 * qJD(1);
t120 = t68 * qJD(4);
t105 = t127 * t88 + t91 * t94;
t69 = t105 * qJD(1);
t119 = t69 * qJD(4);
t107 = -t91 * mrSges(4,1) + t88 * mrSges(4,2);
t106 = qJDD(1) * mrSges(4,3) + t99 * t107;
t109 = t127 * t37 - t94 * t38;
t54 = t68 * pkin(4) - t69 * qJ(5);
t98 = qJD(4) ^ 2;
t22 = -qJDD(4) * pkin(4) - t98 * qJ(5) + t69 * t54 + qJDD(5) - t109;
t87 = sin(pkin(11));
t90 = cos(pkin(11));
t63 = t90 * qJD(4) - t87 * t69;
t64 = t87 * qJD(4) + t90 * t69;
t93 = sin(qJ(6));
t96 = cos(qJ(6));
t42 = t93 * t63 + t96 * t64;
t59 = t105 * qJDD(1) - t120;
t48 = t90 * qJDD(4) - t87 * t59;
t49 = t87 * qJDD(4) + t90 * t59;
t28 = -t42 * qJD(6) + t96 * t48 - t93 * t49;
t41 = t96 * t63 - t93 * t64;
t29 = t41 * qJD(6) + t93 * t48 + t96 * t49;
t67 = qJD(6) + t68;
t32 = -t67 * mrSges(7,2) + t41 * mrSges(7,3);
t33 = t67 * mrSges(7,1) - t42 * mrSges(7,3);
t47 = t68 * pkin(5) - t64 * pkin(8);
t62 = t63 ^ 2;
t104 = t28 * mrSges(7,1) + t41 * t32 - m(7) * (-t48 * pkin(5) - t62 * pkin(8) + t64 * t47 + t22) - t29 * mrSges(7,2) - t42 * t33;
t45 = -t68 * mrSges(6,2) + t63 * mrSges(6,3);
t46 = t68 * mrSges(6,1) - t64 * mrSges(6,3);
t100 = m(6) * t22 - t48 * mrSges(6,1) + t49 * mrSges(6,2) - t63 * t45 + t64 * t46 - t104;
t55 = t68 * mrSges(5,1) + t69 * mrSges(5,2);
t65 = -qJD(4) * mrSges(5,2) - t68 * mrSges(5,3);
t14 = m(5) * t109 + qJDD(4) * mrSges(5,1) - t59 * mrSges(5,3) + qJD(4) * t65 - t69 * t55 - t100;
t23 = -t98 * pkin(4) + qJDD(4) * qJ(5) - t68 * t54 + t125;
t112 = t92 * t74 - t89 * t75;
t108 = qJDD(3) - t112;
t101 = (-pkin(3) * t91 - pkin(2)) * qJDD(1) + (-t122 * pkin(7) - qJ(3)) * t99 + t108;
t58 = t130 * qJDD(1) + t119;
t26 = (-t59 + t120) * qJ(5) + (t58 + t119) * pkin(4) + t101;
t111 = -0.2e1 * qJD(5) * t64 - t87 * t23 + t90 * t26;
t17 = (t63 * t68 - t49) * pkin(8) + (t63 * t64 + t58) * pkin(5) + t111;
t116 = 0.2e1 * qJD(5) * t63 + t90 * t23 + t87 * t26;
t18 = -t62 * pkin(5) + t48 * pkin(8) - t68 * t47 + t116;
t31 = -t41 * mrSges(7,1) + t42 * mrSges(7,2);
t57 = qJDD(6) + t58;
t15 = m(7) * (t96 * t17 - t93 * t18) - t29 * mrSges(7,3) + t57 * mrSges(7,1) - t42 * t31 + t67 * t32;
t16 = m(7) * (t93 * t17 + t96 * t18) + t28 * mrSges(7,3) - t57 * mrSges(7,2) + t41 * t31 - t67 * t33;
t43 = -t63 * mrSges(6,1) + t64 * mrSges(6,2);
t12 = m(6) * t111 + t58 * mrSges(6,1) - t49 * mrSges(6,3) + t96 * t15 + t93 * t16 - t64 * t43 + t68 * t45;
t13 = m(6) * t116 - t58 * mrSges(6,2) + t48 * mrSges(6,3) - t93 * t15 + t96 * t16 + t63 * t43 - t68 * t46;
t66 = qJD(4) * mrSges(5,1) - t69 * mrSges(5,3);
t9 = m(5) * t125 - qJDD(4) * mrSges(5,2) - t58 * mrSges(5,3) - qJD(4) * t66 - t87 * t12 + t90 * t13 - t68 * t55;
t6 = m(4) * t123 + t94 * t9 + t127 * t14 + (-m(4) * t52 - t106) * t88;
t7 = m(4) * t115 + t106 * t91 + t127 * t9 - t94 * t14;
t117 = m(3) * t86 + t91 * t6 + t88 * t7;
t103 = m(5) * t101 + t58 * mrSges(5,1) + t59 * mrSges(5,2) + t90 * t12 + t87 * t13 + t68 * t65 + t69 * t66;
t102 = m(4) * (-qJDD(1) * pkin(2) - t99 * qJ(3) + t108) + t103;
t8 = m(3) * t112 + (-mrSges(3,2) + t129) * t99 + (mrSges(3,1) - t107) * qJDD(1) - t102;
t3 = m(3) * t124 - t99 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t88 * t6 + t91 * t7;
t2 = m(2) * t110 - t99 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t92 * t3 - t89 * t8;
t1 = m(2) * t113 + qJDD(1) * mrSges(2,1) - t99 * mrSges(2,2) + t89 * t3 + t92 * t8;
t4 = [-m(1) * g(1) - t95 * t1 + t97 * t2, t2, t3, t7, t9, t13, t16; -m(1) * g(2) + t97 * t1 + t95 * t2, t1, t8, t6, t14, t12, t15; (-m(1) - m(2)) * g(3) + t117, -m(2) * g(3) + t117, t117, t107 * qJDD(1) - t129 * t99 + t102, t103, t100, -t104;];
f_new  = t4;
