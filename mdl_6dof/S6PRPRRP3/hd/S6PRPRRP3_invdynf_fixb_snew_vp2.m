% Calculate vector of cutting forces with Newton-Euler
% S6PRPRRP3
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
% Datum: 2019-05-04 23:46
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRPRRP3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP3_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP3_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP3_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP3_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP3_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:41:39
% EndTime: 2019-05-04 23:41:44
% DurationCPUTime: 2.28s
% Computational Cost: add. (28003->163), mult. (61724->208), div. (0->0), fcn. (45806->12), ass. (0->90)
t97 = qJD(2) ^ 2;
t87 = cos(pkin(11));
t82 = t87 ^ 2;
t84 = sin(pkin(11));
t120 = t84 ^ 2 + t82;
t132 = t120 * mrSges(4,3);
t91 = sin(qJ(4));
t94 = cos(qJ(4));
t104 = t84 * t91 - t87 * t94;
t69 = t104 * qJD(2);
t85 = sin(pkin(10));
t88 = cos(pkin(10));
t75 = t85 * g(1) - t88 * g(2);
t89 = cos(pkin(6));
t127 = t75 * t89;
t76 = -t88 * g(1) - t85 * g(2);
t83 = -g(3) + qJDD(1);
t86 = sin(pkin(6));
t92 = sin(qJ(2));
t95 = cos(qJ(2));
t131 = (t83 * t86 + t127) * t95 - t92 * t76;
t105 = t84 * t94 + t87 * t91;
t70 = t105 * qJD(2);
t117 = t70 * qJD(4);
t58 = -t104 * qJDD(2) - t117;
t119 = pkin(8) * qJDD(2);
t116 = qJD(2) * qJD(3);
t67 = -t86 * t75 + t89 * t83;
t121 = -0.2e1 * t84 * t116 + t87 * t67;
t129 = pkin(3) * t97;
t126 = t86 * t92;
t111 = t83 * t126 + t92 * t127 + t95 * t76;
t51 = -t97 * pkin(2) + qJDD(2) * qJ(3) + t111;
t30 = (t87 * t129 - t119 - t51) * t84 + t121;
t112 = t84 * t67 + (0.2e1 * t116 + t51) * t87;
t31 = t87 * t119 - t82 * t129 + t112;
t108 = t94 * t30 - t91 * t31;
t57 = t69 * pkin(4) - t70 * pkin(9);
t96 = qJD(4) ^ 2;
t22 = -qJDD(4) * pkin(4) - t96 * pkin(9) + t70 * t57 - t108;
t118 = t69 * qJD(4);
t59 = t105 * qJDD(2) - t118;
t90 = sin(qJ(5));
t93 = cos(qJ(5));
t62 = t90 * qJD(4) + t93 * t70;
t36 = -t62 * qJD(5) + t93 * qJDD(4) - t90 * t59;
t61 = t93 * qJD(4) - t90 * t70;
t37 = t61 * qJD(5) + t90 * qJDD(4) + t93 * t59;
t68 = qJD(5) + t69;
t48 = t68 * pkin(5) - t62 * qJ(6);
t49 = t68 * mrSges(7,1) - t62 * mrSges(7,3);
t60 = t61 ^ 2;
t113 = m(7) * (-t36 * pkin(5) - t60 * qJ(6) + t62 * t48 + qJDD(6) + t22) + t37 * mrSges(7,2) + t62 * t49;
t46 = -t68 * mrSges(7,2) + t61 * mrSges(7,3);
t47 = -t68 * mrSges(6,2) + t61 * mrSges(6,3);
t50 = t68 * mrSges(6,1) - t62 * mrSges(6,3);
t130 = m(6) * t22 + t37 * mrSges(6,2) - (t47 + t46) * t61 - (mrSges(6,1) + mrSges(7,1)) * t36 + t62 * t50 + t113;
t107 = -t87 * mrSges(4,1) + t84 * mrSges(4,2);
t123 = t91 * t30 + t94 * t31;
t23 = -t96 * pkin(4) + qJDD(4) * pkin(9) - t69 * t57 + t123;
t102 = qJDD(3) - t131;
t98 = (-pkin(3) * t87 - pkin(2)) * qJDD(2) + (-t120 * pkin(8) - qJ(3)) * t97 + t102;
t26 = (-t59 + t118) * pkin(9) + (-t58 + t117) * pkin(4) + t98;
t109 = -t90 * t23 + t93 * t26;
t56 = qJDD(5) - t58;
t115 = m(7) * (-0.2e1 * qJD(6) * t62 + (t61 * t68 - t37) * qJ(6) + (t61 * t62 + t56) * pkin(5) + t109) + t68 * t46 + t56 * mrSges(7,1);
t41 = -t61 * mrSges(7,1) + t62 * mrSges(7,2);
t42 = -t61 * mrSges(6,1) + t62 * mrSges(6,2);
t13 = m(6) * t109 + t56 * mrSges(6,1) + t68 * t47 + (-t42 - t41) * t62 + (-mrSges(6,3) - mrSges(7,3)) * t37 + t115;
t124 = t93 * t23 + t90 * t26;
t114 = m(7) * (-t60 * pkin(5) + t36 * qJ(6) + 0.2e1 * qJD(6) * t61 - t68 * t48 + t124) + t61 * t41 + t36 * mrSges(7,3);
t15 = m(6) * t124 + t36 * mrSges(6,3) + t61 * t42 + (-t50 - t49) * t68 + (-mrSges(6,2) - mrSges(7,2)) * t56 + t114;
t65 = -qJD(4) * mrSges(5,2) - t69 * mrSges(5,3);
t66 = qJD(4) * mrSges(5,1) - t70 * mrSges(5,3);
t101 = -m(5) * t98 + t58 * mrSges(5,1) - t59 * mrSges(5,2) - t93 * t13 - t90 * t15 - t69 * t65 - t70 * t66;
t99 = m(4) * (-qJDD(2) * pkin(2) - t97 * qJ(3) + t102) - t101;
t10 = m(3) * t131 + (-mrSges(3,2) + t132) * t97 + (mrSges(3,1) - t107) * qJDD(2) - t99;
t128 = t10 * t95;
t103 = qJDD(2) * mrSges(4,3) + t97 * t107;
t54 = t69 * mrSges(5,1) + t70 * mrSges(5,2);
t11 = m(5) * t123 - qJDD(4) * mrSges(5,2) + t58 * mrSges(5,3) - qJD(4) * t66 - t90 * t13 + t93 * t15 - t69 * t54;
t16 = m(5) * t108 + qJDD(4) * mrSges(5,1) - t59 * mrSges(5,3) + qJD(4) * t65 - t70 * t54 - t130;
t7 = m(4) * t121 + t91 * t11 + t94 * t16 + (-m(4) * t51 - t103) * t84;
t8 = m(4) * t112 + t103 * t87 + t94 * t11 - t91 * t16;
t4 = m(3) * t111 - t97 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t84 * t7 + t87 * t8;
t6 = m(3) * t67 + t87 * t7 + t84 * t8;
t110 = m(2) * t83 + t4 * t126 + t86 * t128 + t89 * t6;
t2 = m(2) * t76 - t92 * t10 + t95 * t4;
t1 = m(2) * t75 - t86 * t6 + (t4 * t92 + t128) * t89;
t3 = [-m(1) * g(1) - t85 * t1 + t88 * t2, t2, t4, t8, t11, t15, -t56 * mrSges(7,2) - t68 * t49 + t114; -m(1) * g(2) + t88 * t1 + t85 * t2, t1, t10, t7, t16, t13, -t37 * mrSges(7,3) - t62 * t41 + t115; -m(1) * g(3) + t110, t110, t6, t107 * qJDD(2) - t97 * t132 + t99, -t101, t130, -t36 * mrSges(7,1) - t61 * t46 + t113;];
f_new  = t3;
