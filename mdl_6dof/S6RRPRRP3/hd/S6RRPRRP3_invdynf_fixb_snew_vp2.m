% Calculate vector of cutting forces with Newton-Euler
% S6RRPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-05-06 17:34
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRRP3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP3_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP3_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP3_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP3_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP3_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 17:29:22
% EndTime: 2019-05-06 17:29:32
% DurationCPUTime: 3.88s
% Computational Cost: add. (45576->203), mult. (103954->261), div. (0->0), fcn. (74478->10), ass. (0->96)
t110 = qJD(2) ^ 2;
t100 = sin(pkin(10));
t101 = cos(pkin(10));
t104 = sin(qJ(2));
t108 = cos(qJ(2));
t132 = qJD(1) * qJD(2);
t111 = qJD(1) ^ 2;
t105 = sin(qJ(1));
t109 = cos(qJ(1));
t120 = -g(1) * t109 - g(2) * t105;
t89 = -pkin(1) * t111 + qJDD(1) * pkin(7) + t120;
t136 = t104 * t89;
t141 = pkin(2) * t111;
t92 = qJDD(1) * t104 + t108 * t132;
t57 = qJDD(2) * pkin(2) - t92 * qJ(3) - t136 + (qJ(3) * t132 + t104 * t141 - g(3)) * t108;
t124 = -g(3) * t104 + t108 * t89;
t93 = qJDD(1) * t108 - t104 * t132;
t134 = qJD(1) * t104;
t94 = qJD(2) * pkin(2) - qJ(3) * t134;
t99 = t108 ^ 2;
t58 = qJ(3) * t93 - qJD(2) * t94 - t99 * t141 + t124;
t86 = (t100 * t108 + t101 * t104) * qJD(1);
t145 = -0.2e1 * qJD(3) * t86 - t100 * t58 + t101 * t57;
t133 = qJD(1) * t108;
t85 = -t100 * t134 + t101 * t133;
t69 = -pkin(3) * t85 - pkin(8) * t86;
t32 = -qJDD(2) * pkin(3) - pkin(8) * t110 + t86 * t69 - t145;
t103 = sin(qJ(4));
t107 = cos(qJ(4));
t74 = t100 * t93 + t101 * t92;
t77 = qJD(2) * t103 + t107 * t86;
t49 = -qJD(4) * t77 + qJDD(2) * t107 - t103 * t74;
t84 = qJD(4) - t85;
t65 = pkin(4) * t84 - pkin(9) * t77;
t76 = qJD(2) * t107 - t103 * t86;
t75 = t76 ^ 2;
t115 = -pkin(4) * t49 - pkin(9) * t75 + t77 * t65 + t32;
t102 = sin(qJ(5));
t106 = cos(qJ(5));
t50 = qJD(4) * t76 + qJDD(2) * t103 + t107 * t74;
t56 = t102 * t76 + t106 * t77;
t29 = -qJD(5) * t56 - t102 * t50 + t106 * t49;
t55 = -t102 * t77 + t106 * t76;
t30 = qJD(5) * t55 + t102 * t49 + t106 * t50;
t80 = qJD(5) + t84;
t45 = pkin(5) * t80 - qJ(6) * t56;
t46 = mrSges(7,1) * t80 - mrSges(7,3) * t56;
t54 = t55 ^ 2;
t129 = m(7) * (-pkin(5) * t29 - qJ(6) * t54 + t45 * t56 + qJDD(6) + t115) + t30 * mrSges(7,2) + t56 * t46;
t43 = -mrSges(7,2) * t80 + mrSges(7,3) * t55;
t44 = -mrSges(6,2) * t80 + mrSges(6,3) * t55;
t47 = mrSges(6,1) * t80 - mrSges(6,3) * t56;
t147 = m(6) * t115 + t30 * mrSges(6,2) + t56 * t47 + t129 - (t44 + t43) * t55 - (mrSges(6,1) + mrSges(7,1)) * t29;
t63 = -mrSges(5,2) * t84 + mrSges(5,3) * t76;
t64 = mrSges(5,1) * t84 - mrSges(5,3) * t77;
t146 = m(5) * t32 - t49 * mrSges(5,1) + t50 * mrSges(5,2) - t76 * t63 + t77 * t64 + t147;
t128 = 0.2e1 * qJD(3) * t85 + t100 * t57 + t101 * t58;
t33 = -pkin(3) * t110 + qJDD(2) * pkin(8) + t69 * t85 + t128;
t125 = t105 * g(1) - t109 * g(2);
t117 = -qJDD(1) * pkin(1) - t125;
t114 = -t93 * pkin(2) + qJDD(3) + t94 * t134 + (-qJ(3) * t99 - pkin(7)) * t111 + t117;
t73 = -t100 * t92 + t101 * t93;
t37 = (-qJD(2) * t85 - t74) * pkin(8) + (qJD(2) * t86 - t73) * pkin(3) + t114;
t121 = -t103 * t33 + t107 * t37;
t72 = qJDD(4) - t73;
t21 = (t76 * t84 - t50) * pkin(9) + (t76 * t77 + t72) * pkin(4) + t121;
t138 = t103 * t37 + t107 * t33;
t23 = -pkin(4) * t75 + pkin(9) * t49 - t65 * t84 + t138;
t122 = -t102 * t23 + t106 * t21;
t70 = qJDD(5) + t72;
t131 = m(7) * (-0.2e1 * qJD(6) * t56 + (t55 * t80 - t30) * qJ(6) + (t55 * t56 + t70) * pkin(5) + t122) + t80 * t43 + t70 * mrSges(7,1);
t40 = -mrSges(7,1) * t55 + mrSges(7,2) * t56;
t41 = -mrSges(6,1) * t55 + mrSges(6,2) * t56;
t12 = m(6) * t122 + t70 * mrSges(6,1) + t80 * t44 + (-t41 - t40) * t56 + (-mrSges(6,3) - mrSges(7,3)) * t30 + t131;
t139 = t102 * t21 + t106 * t23;
t130 = m(7) * (-pkin(5) * t54 + qJ(6) * t29 + 0.2e1 * qJD(6) * t55 - t45 * t80 + t139) + t29 * mrSges(7,3) + t55 * t40;
t13 = m(6) * t139 + t29 * mrSges(6,3) + t55 * t41 + (-t47 - t46) * t80 + (-mrSges(6,2) - mrSges(7,2)) * t70 + t130;
t59 = -mrSges(5,1) * t76 + mrSges(5,2) * t77;
t10 = m(5) * t121 + t72 * mrSges(5,1) - t50 * mrSges(5,3) + t102 * t13 + t106 * t12 - t77 * t59 + t84 * t63;
t11 = m(5) * t138 - t72 * mrSges(5,2) + t49 * mrSges(5,3) - t102 * t12 + t106 * t13 + t76 * t59 - t84 * t64;
t78 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t85;
t79 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t86;
t116 = -m(4) * t114 + t73 * mrSges(4,1) - t74 * mrSges(4,2) - t107 * t10 - t103 * t11 + t85 * t78 - t86 * t79;
t95 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t134;
t96 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t133;
t143 = (t104 * t95 - t108 * t96) * qJD(1) + m(3) * (-t111 * pkin(7) + t117) - t93 * mrSges(3,1) + t92 * mrSges(3,2) - t116;
t67 = -mrSges(4,1) * t85 + mrSges(4,2) * t86;
t14 = m(4) * t145 + qJDD(2) * mrSges(4,1) - t74 * mrSges(4,3) + qJD(2) * t78 - t86 * t67 - t146;
t7 = m(4) * t128 - qJDD(2) * mrSges(4,2) + t73 * mrSges(4,3) - qJD(2) * t79 - t103 * t10 + t107 * t11 + t85 * t67;
t91 = (-mrSges(3,1) * t108 + mrSges(3,2) * t104) * qJD(1);
t4 = m(3) * (-t108 * g(3) - t136) - t92 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t91 * t134 + qJD(2) * t96 + t100 * t7 + t101 * t14;
t5 = m(3) * t124 - qJDD(2) * mrSges(3,2) + t93 * mrSges(3,3) - qJD(2) * t95 - t100 * t14 + t101 * t7 + t91 * t133;
t142 = t104 * t5 + t108 * t4;
t6 = m(2) * t125 + qJDD(1) * mrSges(2,1) - t111 * mrSges(2,2) - t143;
t1 = m(2) * t120 - t111 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t104 * t4 + t108 * t5;
t2 = [-m(1) * g(1) + t1 * t109 - t105 * t6, t1, t5, t7, t11, t13, -t70 * mrSges(7,2) - t80 * t46 + t130; -m(1) * g(2) + t1 * t105 + t109 * t6, t6, t4, t14, t10, t12, -t30 * mrSges(7,3) - t56 * t40 + t131; (-m(1) - m(2)) * g(3) + t142, -m(2) * g(3) + t142, t143, -t116, t146, t147, -t29 * mrSges(7,1) - t55 * t43 + t129;];
f_new  = t2;
