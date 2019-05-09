% Calculate vector of cutting forces with Newton-Euler
% S6RPRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
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
% Datum: 2019-05-05 16:43
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRPPR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR3_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR3_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR3_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR3_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR3_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:41:29
% EndTime: 2019-05-05 16:41:31
% DurationCPUTime: 0.91s
% Computational Cost: add. (6371->183), mult. (12628->216), div. (0->0), fcn. (6010->8), ass. (0->88)
t101 = cos(qJ(3));
t123 = qJD(1) * qJD(3);
t118 = t101 * t123;
t104 = qJD(1) ^ 2;
t102 = cos(qJ(1));
t99 = sin(qJ(1));
t119 = t99 * g(1) - t102 * g(2);
t58 = qJDD(1) * pkin(1) + t119;
t116 = -t102 * g(1) - t99 * g(2);
t63 = -t104 * pkin(1) + t116;
t94 = sin(pkin(9));
t95 = cos(pkin(9));
t133 = t95 * t58 - t94 * t63;
t29 = -qJDD(1) * pkin(2) - t104 * pkin(7) - t133;
t98 = sin(qJ(3));
t65 = t98 * qJDD(1) + t118;
t120 = t98 * t123;
t66 = t101 * qJDD(1) - t120;
t114 = -t66 * pkin(3) + t29 + (-t118 - t65) * qJ(4);
t100 = cos(qJ(6));
t125 = t98 * qJD(1);
t128 = t104 * t101 ^ 2;
t140 = 2 * qJD(4);
t70 = -qJD(3) * pkin(4) - qJ(5) * t125;
t106 = -qJ(5) * t128 + qJDD(5) - t114 + (t140 + t70) * t125;
t138 = pkin(4) + pkin(8);
t139 = -pkin(3) - pkin(8);
t14 = t138 * t66 + (pkin(5) * t101 + t139 * t98) * t123 + t106 + t65 * pkin(5);
t103 = qJD(3) ^ 2;
t127 = t104 * t98;
t122 = qJD(1) * qJD(5);
t145 = -0.2e1 * t98 * t122 + (t118 - t65) * qJ(5);
t59 = (-pkin(3) * t101 - qJ(4) * t98) * qJD(1);
t147 = t59 * t125 + qJDD(4);
t132 = t94 * t58 + t95 * t63;
t30 = -t104 * pkin(2) + qJDD(1) * pkin(7) + t132;
t27 = t98 * t30;
t64 = (pkin(5) * t98 + pkin(8) * t101) * qJD(1);
t92 = -g(3) + qJDD(2);
t17 = -t64 * t125 + t27 + (-pkin(5) - qJ(4)) * t103 + (-pkin(4) * t127 - t92) * t101 + (-pkin(3) - t138) * qJDD(3) + t145 + t147;
t124 = qJD(1) * t101;
t97 = sin(qJ(6));
t56 = -t100 * qJD(3) + t124 * t97;
t34 = t56 * qJD(6) - t97 * qJDD(3) - t100 * t66;
t57 = -t97 * qJD(3) - t100 * t124;
t35 = -t56 * mrSges(7,1) + t57 * mrSges(7,2);
t79 = qJD(6) + t125;
t36 = -t79 * mrSges(7,2) + t56 * mrSges(7,3);
t54 = qJDD(6) + t65;
t12 = m(7) * (t100 * t14 - t97 * t17) - t34 * mrSges(7,3) + t54 * mrSges(7,1) - t57 * t35 + t79 * t36;
t33 = -t57 * qJD(6) - t100 * qJDD(3) + t97 * t66;
t37 = t79 * mrSges(7,1) - t57 * mrSges(7,3);
t13 = m(7) * (t100 * t17 + t97 * t14) + t33 * mrSges(7,3) - t54 * mrSges(7,2) + t56 * t35 - t79 * t37;
t71 = qJD(3) * mrSges(6,2) - mrSges(6,3) * t125;
t74 = -qJD(3) * mrSges(6,1) + mrSges(6,3) * t124;
t115 = t100 * t12 + t97 * t13 + m(6) * (-pkin(3) * t120 + t66 * pkin(4) + t106) + t71 * t125 - t74 * t124 + t65 * mrSges(6,1) - t66 * mrSges(6,2);
t141 = -2 * qJD(4);
t111 = m(5) * ((pkin(3) * qJD(3) + t141) * t125 + t114) - t66 * mrSges(5,1) - t115;
t76 = mrSges(5,2) * t124 + qJD(3) * mrSges(5,3);
t129 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t124 + t76;
t72 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t125;
t73 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t125;
t148 = -(t101 * t129 - (t72 - t73) * t98) * qJD(1) + m(4) * t29 - t66 * mrSges(4,1) + (mrSges(4,2) - mrSges(5,3)) * t65 + t111;
t134 = t101 * t30 + t98 * t92;
t146 = qJDD(3) * qJ(4) + t59 * t124 + t134;
t137 = t103 * pkin(3);
t60 = (-mrSges(5,1) * t101 - mrSges(5,3) * t98) * qJD(1);
t85 = qJD(3) * t140;
t144 = m(5) * (t85 - t137 + t146) + t60 * t124 + qJD(3) * t73 + qJDD(3) * mrSges(5,3);
t135 = mrSges(4,3) + mrSges(5,2);
t62 = (mrSges(6,1) * t98 - mrSges(6,2) * t101) * qJD(1);
t131 = (-mrSges(4,1) * t101 + mrSges(4,2) * t98) * qJD(1) - t62;
t117 = t101 * t92 - t27;
t26 = -qJDD(3) * pkin(3) - t103 * qJ(4) - t117 + t147;
t113 = -t100 * t13 + t97 * t12 - m(6) * ((-t101 * t127 - qJDD(3)) * pkin(4) + t26 + t145) + t65 * mrSges(6,3) - qJD(3) * t74 - qJDD(3) * mrSges(6,2);
t110 = m(5) * t26 - t113;
t6 = m(4) * t117 - t135 * t65 + (mrSges(4,1) + mrSges(5,1)) * qJDD(3) + t129 * qJD(3) + (-t60 - t131) * t125 - t110;
t109 = pkin(4) * t128 + t66 * qJ(5) - t146;
t112 = -t33 * mrSges(7,1) - t56 * t36 + m(7) * (qJDD(3) * pkin(5) + qJD(3) * t70 + t85 + t139 * t103 + (-0.2e1 * qJD(5) - t64) * t124 - t109) + t34 * mrSges(7,2) + t57 * t37;
t108 = -m(6) * (0.2e1 * t101 * t122 + t137 + (t141 - t70) * qJD(3) + t109) + t112 - t66 * mrSges(6,3) + qJD(3) * t71 + qJDD(3) * mrSges(6,1);
t8 = m(4) * t134 - qJDD(3) * mrSges(4,2) - qJD(3) * t72 + t131 * t124 + t135 * t66 + t108 + t144;
t121 = m(3) * t92 + t101 * t6 + t98 * t8;
t105 = t124 * t62 - t108;
t4 = m(3) * t133 + qJDD(1) * mrSges(3,1) - t104 * mrSges(3,2) - t148;
t3 = m(3) * t132 - t104 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t101 * t8 - t98 * t6;
t2 = m(2) * t116 - t104 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t95 * t3 - t94 * t4;
t1 = m(2) * t119 + qJDD(1) * mrSges(2,1) - t104 * mrSges(2,2) + t94 * t3 + t95 * t4;
t5 = [-m(1) * g(1) - t99 * t1 + t102 * t2, t2, t3, t8, t66 * mrSges(5,2) - t105 + t144, -t125 * t62 - t113, t13; -m(1) * g(2) + t102 * t1 + t99 * t2, t1, t4, t6, -t65 * mrSges(5,3) + (-t101 * t76 - t98 * t73) * qJD(1) + t111, t105, t12; (-m(1) - m(2)) * g(3) + t121, -m(2) * g(3) + t121, t121, t148, -qJDD(3) * mrSges(5,1) + t65 * mrSges(5,2) - qJD(3) * t76 + (t60 - t62) * t125 + t110, t115, t112;];
f_new  = t5;
