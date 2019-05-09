% Calculate vector of cutting forces with Newton-Euler
% S6RRPRRP4
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
% Datum: 2019-05-06 17:44
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRRP4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP4_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP4_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP4_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP4_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP4_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP4_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP4_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 17:36:26
% EndTime: 2019-05-06 17:36:39
% DurationCPUTime: 3.50s
% Computational Cost: add. (44280->202), mult. (100466->261), div. (0->0), fcn. (71406->10), ass. (0->97)
t100 = cos(pkin(10));
t103 = sin(qJ(2));
t106 = cos(qJ(2));
t127 = qJD(1) * qJD(2);
t109 = qJD(1) ^ 2;
t104 = sin(qJ(1));
t107 = cos(qJ(1));
t120 = -t107 * g(1) - t104 * g(2);
t88 = -t109 * pkin(1) + qJDD(1) * pkin(7) + t120;
t131 = t103 * t88;
t136 = pkin(2) * t109;
t91 = t103 * qJDD(1) + t106 * t127;
t54 = qJDD(2) * pkin(2) - t91 * qJ(3) - t131 + (qJ(3) * t127 + t103 * t136 - g(3)) * t106;
t123 = -t103 * g(3) + t106 * t88;
t92 = t106 * qJDD(1) - t103 * t127;
t129 = qJD(1) * t103;
t93 = qJD(2) * pkin(2) - qJ(3) * t129;
t98 = t106 ^ 2;
t55 = t92 * qJ(3) - qJD(2) * t93 - t98 * t136 + t123;
t99 = sin(pkin(10));
t85 = (t100 * t103 + t106 * t99) * qJD(1);
t141 = -0.2e1 * qJD(3) * t85 + t100 * t54 - t99 * t55;
t101 = sin(qJ(5));
t137 = cos(qJ(5));
t102 = sin(qJ(4));
t105 = cos(qJ(4));
t108 = qJD(2) ^ 2;
t128 = qJD(1) * t106;
t84 = t100 * t128 - t99 * t129;
t125 = 0.2e1 * qJD(3) * t84 + t100 * t55 + t99 * t54;
t66 = -t84 * pkin(3) - t85 * pkin(8);
t31 = -t108 * pkin(3) + qJDD(2) * pkin(8) + t84 * t66 + t125;
t124 = t104 * g(1) - t107 * g(2);
t118 = -qJDD(1) * pkin(1) - t124;
t112 = -t92 * pkin(2) + qJDD(3) + t93 * t129 + (-qJ(3) * t98 - pkin(7)) * t109 + t118;
t70 = t100 * t92 - t99 * t91;
t71 = t100 * t91 + t99 * t92;
t34 = (-qJD(2) * t84 - t71) * pkin(8) + (qJD(2) * t85 - t70) * pkin(3) + t112;
t121 = -t102 * t31 + t105 * t34;
t74 = t105 * qJD(2) - t102 * t85;
t48 = t74 * qJD(4) + t102 * qJDD(2) + t105 * t71;
t69 = qJDD(4) - t70;
t75 = t102 * qJD(2) + t105 * t85;
t83 = qJD(4) - t84;
t20 = (t74 * t83 - t48) * pkin(9) + (t74 * t75 + t69) * pkin(4) + t121;
t133 = t102 * t34 + t105 * t31;
t47 = -t75 * qJD(4) + t105 * qJDD(2) - t102 * t71;
t62 = t83 * pkin(4) - t75 * pkin(9);
t73 = t74 ^ 2;
t22 = -t73 * pkin(4) + t47 * pkin(9) - t83 * t62 + t133;
t134 = t101 * t20 + t137 * t22;
t52 = t101 * t75 - t137 * t74;
t53 = t101 * t74 + t137 * t75;
t37 = t52 * pkin(5) - t53 * qJ(6);
t79 = qJD(5) + t83;
t44 = -t79 * mrSges(7,1) + t53 * mrSges(7,2);
t67 = qJDD(5) + t69;
t78 = t79 ^ 2;
t126 = m(7) * (-t78 * pkin(5) + t67 * qJ(6) + 0.2e1 * qJD(6) * t79 - t52 * t37 + t134) + t79 * t44 + t67 * mrSges(7,3);
t38 = t52 * mrSges(7,1) - t53 * mrSges(7,3);
t132 = -t52 * mrSges(6,1) - t53 * mrSges(6,2) - t38;
t135 = -mrSges(6,3) - mrSges(7,2);
t27 = t53 * qJD(5) + t101 * t48 - t137 * t47;
t43 = t79 * mrSges(6,1) - t53 * mrSges(6,3);
t12 = m(6) * t134 - t67 * mrSges(6,2) + t132 * t52 + t135 * t27 - t79 * t43 + t126;
t117 = -t101 * t22 + t137 * t20;
t138 = m(7) * (-t67 * pkin(5) - t78 * qJ(6) + t53 * t37 + qJDD(6) - t117);
t28 = -t52 * qJD(5) + t101 * t47 + t137 * t48;
t41 = -t52 * mrSges(7,2) + t79 * mrSges(7,3);
t42 = -t79 * mrSges(6,2) - t52 * mrSges(6,3);
t14 = m(6) * t117 - t138 + (t42 + t41) * t79 + (mrSges(6,1) + mrSges(7,1)) * t67 + t132 * t53 + t135 * t28;
t56 = -t74 * mrSges(5,1) + t75 * mrSges(5,2);
t60 = -t83 * mrSges(5,2) + t74 * mrSges(5,3);
t10 = m(5) * t121 + t69 * mrSges(5,1) - t48 * mrSges(5,3) + t101 * t12 + t137 * t14 - t75 * t56 + t83 * t60;
t61 = t83 * mrSges(5,1) - t75 * mrSges(5,3);
t11 = m(5) * t133 - t69 * mrSges(5,2) + t47 * mrSges(5,3) - t101 * t14 + t137 * t12 + t74 * t56 - t83 * t61;
t76 = -qJD(2) * mrSges(4,2) + t84 * mrSges(4,3);
t77 = qJD(2) * mrSges(4,1) - t85 * mrSges(4,3);
t115 = -m(4) * t112 + t70 * mrSges(4,1) - t71 * mrSges(4,2) - t105 * t10 - t102 * t11 + t84 * t76 - t85 * t77;
t94 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t129;
t95 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t128;
t140 = (t103 * t94 - t106 * t95) * qJD(1) + m(3) * (-t109 * pkin(7) + t118) - t92 * mrSges(3,1) + t91 * mrSges(3,2) - t115;
t30 = -qJDD(2) * pkin(3) - t108 * pkin(8) + t85 * t66 - t141;
t113 = -t47 * pkin(4) - t73 * pkin(9) + t75 * t62 + t30;
t116 = t28 * mrSges(7,3) + t53 * t44 - m(7) * (-0.2e1 * qJD(6) * t53 + (t52 * t79 - t28) * qJ(6) + (t53 * t79 + t27) * pkin(5) + t113) - t27 * mrSges(7,1) - t52 * t41;
t114 = m(6) * t113 + t27 * mrSges(6,1) + t28 * mrSges(6,2) + t52 * t42 + t53 * t43 - t116;
t110 = m(5) * t30 - t47 * mrSges(5,1) + t48 * mrSges(5,2) - t74 * t60 + t75 * t61 + t114;
t64 = -t84 * mrSges(4,1) + t85 * mrSges(4,2);
t13 = m(4) * t141 + qJDD(2) * mrSges(4,1) - t71 * mrSges(4,3) + qJD(2) * t76 - t85 * t64 - t110;
t7 = m(4) * t125 - qJDD(2) * mrSges(4,2) + t70 * mrSges(4,3) - qJD(2) * t77 - t102 * t10 + t105 * t11 + t84 * t64;
t90 = (-mrSges(3,1) * t106 + mrSges(3,2) * t103) * qJD(1);
t4 = m(3) * (-t106 * g(3) - t131) - t91 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t90 * t129 + qJD(2) * t95 + t99 * t7 + t100 * t13;
t5 = m(3) * t123 - qJDD(2) * mrSges(3,2) + t92 * mrSges(3,3) - qJD(2) * t94 + t100 * t7 + t90 * t128 - t99 * t13;
t139 = t103 * t5 + t106 * t4;
t6 = m(2) * t124 + qJDD(1) * mrSges(2,1) - t109 * mrSges(2,2) - t140;
t1 = m(2) * t120 - t109 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t103 * t4 + t106 * t5;
t2 = [-m(1) * g(1) + t107 * t1 - t104 * t6, t1, t5, t7, t11, t12, -t27 * mrSges(7,2) - t52 * t38 + t126; -m(1) * g(2) + t104 * t1 + t107 * t6, t6, t4, t13, t10, t14, -t116; (-m(1) - m(2)) * g(3) + t139, -m(2) * g(3) + t139, t140, -t115, t110, t114, -t67 * mrSges(7,1) + t28 * mrSges(7,2) + t53 * t38 - t79 * t41 + t138;];
f_new  = t2;
