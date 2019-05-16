% Calculate vector of cutting forces with Newton-Euler
% S6RRPRRP1
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
% Datum: 2019-05-06 17:21
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRRP1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP1_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 17:16:28
% EndTime: 2019-05-06 17:16:38
% DurationCPUTime: 3.92s
% Computational Cost: add. (47470->204), mult. (109888->262), div. (0->0), fcn. (80761->10), ass. (0->96)
t103 = sin(qJ(2));
t107 = cos(qJ(2));
t109 = qJD(1) ^ 2;
t104 = sin(qJ(1));
t108 = cos(qJ(1));
t123 = t104 * g(1) - t108 * g(2);
t116 = -qJDD(1) * pkin(1) - t123;
t130 = qJD(1) * t103;
t128 = qJD(1) * qJD(2);
t90 = t107 * qJDD(1) - t103 * t128;
t91 = qJD(2) * pkin(2) - qJ(3) * t130;
t98 = t107 ^ 2;
t113 = -t90 * pkin(2) + qJDD(3) + t91 * t130 + (-qJ(3) * t98 - pkin(7)) * t109 + t116;
t101 = sin(qJ(5));
t105 = cos(qJ(5));
t102 = sin(qJ(4));
t106 = cos(qJ(4));
t100 = cos(pkin(10));
t118 = -t108 * g(1) - t104 * g(2);
t86 = -t109 * pkin(1) + qJDD(1) * pkin(7) + t118;
t131 = t103 * t86;
t136 = pkin(2) * t109;
t89 = t103 * qJDD(1) + t107 * t128;
t58 = qJDD(2) * pkin(2) - t89 * qJ(3) - t131 + (qJ(3) * t128 + t103 * t136 - g(3)) * t107;
t122 = -t103 * g(3) + t107 * t86;
t59 = t90 * qJ(3) - qJD(2) * t91 - t98 * t136 + t122;
t99 = sin(pkin(10));
t84 = (t100 * t103 + t107 * t99) * qJD(1);
t119 = -0.2e1 * qJD(3) * t84 + t100 * t58 - t99 * t59;
t75 = t100 * t89 + t99 * t90;
t83 = (t100 * t107 - t103 * t99) * qJD(1);
t26 = (qJD(2) * t83 - t75) * pkin(8) + (t83 * t84 + qJDD(2)) * pkin(3) + t119;
t124 = 0.2e1 * qJD(3) * t83 + t100 * t59 + t99 * t58;
t74 = t100 * t90 - t99 * t89;
t78 = qJD(2) * pkin(3) - t84 * pkin(8);
t82 = t83 ^ 2;
t29 = -t82 * pkin(3) + t74 * pkin(8) - qJD(2) * t78 + t124;
t133 = t102 * t26 + t106 * t29;
t68 = -t102 * t84 + t106 * t83;
t69 = t102 * t83 + t106 * t84;
t54 = -t68 * pkin(4) - t69 * pkin(9);
t97 = qJD(2) + qJD(4);
t95 = t97 ^ 2;
t96 = qJDD(2) + qJDD(4);
t21 = -t95 * pkin(4) + t96 * pkin(9) + t68 * t54 + t133;
t111 = -t74 * pkin(3) - t82 * pkin(8) + t84 * t78 + t113;
t43 = -t69 * qJD(4) - t102 * t75 + t106 * t74;
t44 = t68 * qJD(4) + t102 * t74 + t106 * t75;
t24 = (-t68 * t97 - t44) * pkin(9) + (t69 * t97 - t43) * pkin(4) + t111;
t121 = -t101 * t21 + t105 * t24;
t63 = -t101 * t69 + t105 * t97;
t33 = t63 * qJD(5) + t101 * t96 + t105 * t44;
t42 = qJDD(5) - t43;
t67 = qJD(5) - t68;
t47 = -t67 * mrSges(7,2) + t63 * mrSges(7,3);
t64 = t101 * t97 + t105 * t69;
t127 = m(7) * (-0.2e1 * qJD(6) * t64 + (t63 * t67 - t33) * qJ(6) + (t63 * t64 + t42) * pkin(5) + t121) + t67 * t47 + t42 * mrSges(7,1);
t45 = -t63 * mrSges(7,1) + t64 * mrSges(7,2);
t46 = -t63 * mrSges(6,1) + t64 * mrSges(6,2);
t48 = -t67 * mrSges(6,2) + t63 * mrSges(6,3);
t11 = m(6) * t121 + t42 * mrSges(6,1) + t67 * t48 + (-t46 - t45) * t64 + (-mrSges(6,3) - mrSges(7,3)) * t33 + t127;
t134 = t101 * t24 + t105 * t21;
t32 = -t64 * qJD(5) - t101 * t44 + t105 * t96;
t49 = t67 * pkin(5) - t64 * qJ(6);
t61 = t63 ^ 2;
t126 = m(7) * (-t61 * pkin(5) + t32 * qJ(6) + 0.2e1 * qJD(6) * t63 - t67 * t49 + t134) + t32 * mrSges(7,3) + t63 * t45;
t50 = t67 * mrSges(7,1) - t64 * mrSges(7,3);
t51 = t67 * mrSges(6,1) - t64 * mrSges(6,3);
t13 = m(6) * t134 + t32 * mrSges(6,3) + t63 * t46 + (-t51 - t50) * t67 + (-mrSges(6,2) - mrSges(7,2)) * t42 + t126;
t65 = -t97 * mrSges(5,2) + t68 * mrSges(5,3);
t66 = t97 * mrSges(5,1) - t69 * mrSges(5,3);
t115 = -m(5) * t111 + t43 * mrSges(5,1) - t44 * mrSges(5,2) - t101 * t13 - t105 * t11 + t68 * t65 - t69 * t66;
t76 = -qJD(2) * mrSges(4,2) + t83 * mrSges(4,3);
t77 = qJD(2) * mrSges(4,1) - t84 * mrSges(4,3);
t112 = -m(4) * t113 + t74 * mrSges(4,1) - t75 * mrSges(4,2) + t83 * t76 - t84 * t77 + t115;
t92 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t130;
t129 = qJD(1) * t107;
t93 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t129;
t139 = (t103 * t92 - t107 * t93) * qJD(1) + m(3) * (-t109 * pkin(7) + t116) - t90 * mrSges(3,1) + t89 * mrSges(3,2) - t112;
t120 = -t102 * t29 + t106 * t26;
t20 = -t96 * pkin(4) - t95 * pkin(9) + t69 * t54 - t120;
t125 = m(7) * (-t32 * pkin(5) - t61 * qJ(6) + t64 * t49 + qJDD(6) + t20) + t33 * mrSges(7,2) + t64 * t50;
t138 = m(6) * t20 + t33 * mrSges(6,2) - (t48 + t47) * t63 - (mrSges(6,1) + mrSges(7,1)) * t32 + t64 * t51 + t125;
t53 = -t68 * mrSges(5,1) + t69 * mrSges(5,2);
t14 = m(5) * t120 + t96 * mrSges(5,1) - t44 * mrSges(5,3) - t69 * t53 + t97 * t65 - t138;
t72 = -t83 * mrSges(4,1) + t84 * mrSges(4,2);
t9 = m(5) * t133 - t96 * mrSges(5,2) + t43 * mrSges(5,3) - t101 * t11 + t105 * t13 + t68 * t53 - t97 * t66;
t6 = m(4) * t119 + qJDD(2) * mrSges(4,1) - t75 * mrSges(4,3) + qJD(2) * t76 + t102 * t9 + t106 * t14 - t84 * t72;
t7 = m(4) * t124 - qJDD(2) * mrSges(4,2) + t74 * mrSges(4,3) - qJD(2) * t77 - t102 * t14 + t106 * t9 + t83 * t72;
t88 = (-mrSges(3,1) * t107 + mrSges(3,2) * t103) * qJD(1);
t4 = m(3) * (-t107 * g(3) - t131) - t89 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t88 * t130 + qJD(2) * t93 + t99 * t7 + t100 * t6;
t5 = m(3) * t122 - qJDD(2) * mrSges(3,2) + t90 * mrSges(3,3) - qJD(2) * t92 + t100 * t7 + t88 * t129 - t99 * t6;
t137 = t103 * t5 + t107 * t4;
t8 = m(2) * t123 + qJDD(1) * mrSges(2,1) - t109 * mrSges(2,2) - t139;
t1 = m(2) * t118 - t109 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t103 * t4 + t107 * t5;
t2 = [-m(1) * g(1) + t108 * t1 - t104 * t8, t1, t5, t7, t9, t13, -t42 * mrSges(7,2) - t67 * t50 + t126; -m(1) * g(2) + t104 * t1 + t108 * t8, t8, t4, t6, t14, t11, -t33 * mrSges(7,3) - t64 * t45 + t127; (-m(1) - m(2)) * g(3) + t137, -m(2) * g(3) + t137, t139, -t112, -t115, t138, -t32 * mrSges(7,1) - t63 * t47 + t125;];
f_new  = t2;
