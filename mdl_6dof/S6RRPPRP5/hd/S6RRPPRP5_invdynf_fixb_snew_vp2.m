% Calculate vector of cutting forces with Newton-Euler
% S6RRPPRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
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
% Datum: 2019-05-06 09:29
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPPRP5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP5_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP5_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP5_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP5_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP5_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:25:08
% EndTime: 2019-05-06 09:25:12
% DurationCPUTime: 1.70s
% Computational Cost: add. (17136->200), mult. (37352->245), div. (0->0), fcn. (22144->8), ass. (0->92)
t147 = -2 * qJD(3);
t101 = cos(qJ(2));
t100 = sin(qJ(1));
t102 = cos(qJ(1));
t121 = t100 * g(1) - t102 * g(2);
t116 = -qJDD(1) * pkin(1) - t121;
t127 = qJD(1) * qJD(2);
t120 = t101 * t127;
t99 = sin(qJ(2));
t122 = t99 * t127;
t129 = t99 * qJD(1);
t78 = t99 * qJDD(1) + t120;
t110 = pkin(2) * t122 + t129 * t147 + (-t78 - t120) * qJ(3) + t116;
t128 = qJD(1) * t101;
t104 = qJD(1) ^ 2;
t139 = t104 * pkin(7);
t79 = t101 * qJDD(1) - t122;
t83 = pkin(3) * t129 - qJD(2) * qJ(4);
t95 = t101 ^ 2;
t23 = -t83 * t129 + (-pkin(2) - qJ(4)) * t79 + (-pkin(3) * t95 - pkin(7)) * t104 + t110;
t103 = qJD(2) ^ 2;
t118 = -t102 * g(1) - t100 * g(2);
t67 = -t104 * pkin(1) + qJDD(1) * pkin(7) + t118;
t133 = -t101 * g(3) - t99 * t67;
t75 = (-pkin(2) * t101 - qJ(3) * t99) * qJD(1);
t42 = -qJDD(2) * pkin(2) - t103 * qJ(3) + t75 * t129 + qJDD(3) - t133;
t32 = (-t101 * t104 * t99 - qJDD(2)) * qJ(4) + (t78 - t120) * pkin(3) + t42;
t96 = sin(pkin(9));
t97 = cos(pkin(9));
t72 = t97 * qJD(2) - t96 * t128;
t119 = -0.2e1 * qJD(4) * t72 - t96 * t23 + t97 * t32;
t140 = cos(qJ(5));
t71 = -t96 * qJD(2) - t97 * t128;
t51 = -t71 * mrSges(5,1) + t72 * mrSges(5,2);
t54 = -mrSges(5,2) * t129 + t71 * mrSges(5,3);
t57 = t97 * qJDD(2) - t96 * t79;
t17 = (t71 * t129 - t57) * pkin(8) + (t71 * t72 + t78) * pkin(4) + t119;
t124 = 0.2e1 * qJD(4) * t71 + t97 * t23 + t96 * t32;
t56 = -t96 * qJDD(2) - t97 * t79;
t58 = pkin(4) * t129 - t72 * pkin(8);
t70 = t71 ^ 2;
t19 = -t70 * pkin(4) + t56 * pkin(8) - t58 * t129 + t124;
t98 = sin(qJ(5));
t135 = t140 * t19 + t98 * t17;
t49 = -t140 * t71 + t98 * t72;
t50 = t140 * t72 + t98 * t71;
t35 = t49 * pkin(5) - t50 * qJ(6);
t88 = qJD(5) + t129;
t46 = -t88 * mrSges(7,1) + t50 * mrSges(7,2);
t74 = qJDD(5) + t78;
t86 = t88 ^ 2;
t125 = m(7) * (-t86 * pkin(5) + t74 * qJ(6) + 0.2e1 * qJD(6) * t88 - t49 * t35 + t135) + t88 * t46 + t74 * mrSges(7,3);
t36 = t49 * mrSges(7,1) - t50 * mrSges(7,3);
t134 = -t49 * mrSges(6,1) - t50 * mrSges(6,2) - t36;
t136 = -mrSges(6,3) - mrSges(7,2);
t26 = t50 * qJD(5) - t140 * t56 + t98 * t57;
t45 = t88 * mrSges(6,1) - t50 * mrSges(6,3);
t8 = m(6) * t135 - t74 * mrSges(6,2) + t134 * t49 + t136 * t26 - t88 * t45 + t125;
t115 = t140 * t17 - t98 * t19;
t142 = m(7) * (-t74 * pkin(5) - t86 * qJ(6) + t50 * t35 + qJDD(6) - t115);
t27 = -t49 * qJD(5) + t140 * t57 + t98 * t56;
t43 = -t88 * mrSges(6,2) - t49 * mrSges(6,3);
t44 = -t49 * mrSges(7,2) + t88 * mrSges(7,3);
t9 = m(6) * t115 - t142 + (t43 + t44) * t88 + (mrSges(6,1) + mrSges(7,1)) * t74 + t134 * t50 + t136 * t27;
t6 = m(5) * t119 + t78 * mrSges(5,1) - t57 * mrSges(5,3) + t54 * t129 + t140 * t9 - t72 * t51 + t98 * t8;
t55 = mrSges(5,1) * t129 - t72 * mrSges(5,3);
t7 = m(5) * t124 - t78 * mrSges(5,2) + t56 * mrSges(5,3) - t55 * t129 + t140 * t8 + t71 * t51 - t98 * t9;
t84 = -mrSges(4,1) * t128 - qJD(2) * mrSges(4,3);
t117 = t96 * t6 - t97 * t7 - m(4) * (-t79 * pkin(2) + t110 - t139) - t84 * t128 + t78 * mrSges(4,3);
t85 = mrSges(4,1) * t129 + qJD(2) * mrSges(4,2);
t131 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t129 - t85;
t138 = mrSges(3,1) - mrSges(4,2);
t82 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t128;
t146 = qJD(1) * (-t101 * t82 + t131 * t99) - t138 * t79 + m(3) * (t116 - t139) + t78 * mrSges(3,2) - t117;
t123 = -t99 * g(3) + t101 * t67;
t145 = t103 * pkin(2) - qJDD(2) * qJ(3) + qJD(2) * t147 - t75 * t128 - t123;
t108 = -t95 * t104 * qJ(4) + t79 * pkin(3) + qJD(2) * t83 + qJDD(4) - t145;
t107 = -t56 * pkin(4) - t70 * pkin(8) + t72 * t58 + t108;
t112 = t27 * mrSges(7,3) + t50 * t46 - m(7) * (-0.2e1 * qJD(6) * t50 + t107 + (t49 * t88 - t27) * qJ(6) + (t50 * t88 + t26) * pkin(5)) - t26 * mrSges(7,1) - t49 * t44;
t109 = m(6) * t107 + t26 * mrSges(6,1) + t27 * mrSges(6,2) + t49 * t43 + t50 * t45 - t112;
t106 = -m(5) * t108 + t56 * mrSges(5,1) - t57 * mrSges(5,2) + t71 * t54 - t72 * t55 - t109;
t105 = m(4) * t145 + t106;
t76 = (mrSges(4,2) * t101 - mrSges(4,3) * t99) * qJD(1);
t132 = t76 + (-mrSges(3,1) * t101 + mrSges(3,2) * t99) * qJD(1);
t137 = mrSges(3,3) + mrSges(4,1);
t11 = t132 * t128 - t105 + m(3) * t123 + t137 * t79 + (-mrSges(3,2) + mrSges(4,3)) * qJDD(2) - t131 * qJD(2);
t114 = -m(4) * t42 - t97 * t6 - t96 * t7;
t4 = m(3) * t133 - t137 * t78 + t138 * qJDD(2) + (t82 - t84) * qJD(2) - t132 * t129 + t114;
t141 = t101 * t4 + t99 * t11;
t2 = m(2) * t121 + qJDD(1) * mrSges(2,1) - t104 * mrSges(2,2) - t146;
t1 = m(2) * t118 - t104 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t101 * t11 - t99 * t4;
t3 = [-m(1) * g(1) + t102 * t1 - t100 * t2, t1, t11, t79 * mrSges(4,2) - t85 * t129 - t117, t7, t8, -t26 * mrSges(7,2) - t49 * t36 + t125; -m(1) * g(2) + t100 * t1 + t102 * t2, t2, t4, -t79 * mrSges(4,1) - qJDD(2) * mrSges(4,3) - qJD(2) * t85 - t76 * t128 + t105, t6, t9, -t112; (-m(1) - m(2)) * g(3) + t141, -m(2) * g(3) + t141, t146, t78 * mrSges(4,1) + qJDD(2) * mrSges(4,2) + qJD(2) * t84 + t76 * t129 - t114, -t106, t109, -t74 * mrSges(7,1) + t27 * mrSges(7,2) + t50 * t36 - t88 * t44 + t142;];
f_new  = t3;
