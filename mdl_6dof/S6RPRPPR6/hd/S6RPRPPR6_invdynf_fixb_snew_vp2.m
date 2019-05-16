% Calculate vector of cutting forces with Newton-Euler
% S6RPRPPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
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
% Datum: 2019-05-05 17:10
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRPPR6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR6_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR6_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR6_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR6_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR6_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:07:22
% EndTime: 2019-05-05 17:07:28
% DurationCPUTime: 2.40s
% Computational Cost: add. (28652->178), mult. (64887->233), div. (0->0), fcn. (43083->10), ass. (0->91)
t132 = -2 * qJD(4);
t96 = sin(qJ(1));
t99 = cos(qJ(1));
t112 = -t99 * g(1) - t96 * g(2);
t110 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t112;
t101 = qJD(1) ^ 2;
t122 = qJD(1) * qJD(3);
t95 = sin(qJ(3));
t115 = t95 * t122;
t117 = t96 * g(1) - t99 * g(2);
t108 = -t101 * qJ(2) + qJDD(2) - t117;
t129 = -pkin(1) - pkin(7);
t65 = t129 * qJDD(1) + t108;
t98 = cos(qJ(3));
t126 = t95 * g(3) + t98 * t65;
t80 = t98 * qJDD(1) - t115;
t34 = (-t80 - t115) * qJ(4) + (-t101 * t95 * t98 + qJDD(3)) * pkin(3) + t126;
t116 = -t98 * g(3) + t95 * t65;
t79 = -t95 * qJDD(1) - t98 * t122;
t124 = qJD(1) * t98;
t82 = qJD(3) * pkin(3) - qJ(4) * t124;
t89 = t95 ^ 2;
t35 = -t89 * t101 * pkin(3) + t79 * qJ(4) - qJD(3) * t82 + t116;
t125 = qJD(1) * t95;
t91 = sin(pkin(9));
t93 = cos(pkin(9));
t72 = t93 * t124 - t91 * t125;
t131 = t72 * t132 + t93 * t34 - t91 * t35;
t130 = -m(2) - m(3);
t128 = (mrSges(2,1) - mrSges(3,2));
t127 = -mrSges(2,2) + mrSges(3,3);
t100 = qJD(3) ^ 2;
t71 = (t91 * t98 + t93 * t95) * qJD(1);
t119 = t71 * t132 + t91 * t34 + t93 * t35;
t48 = t71 * pkin(4) - t72 * qJ(5);
t20 = -t100 * pkin(4) + qJDD(3) * qJ(5) - t71 * t48 + t119;
t104 = -t79 * pkin(3) + qJDD(4) + t82 * t124 + (-qJ(4) * t89 + t129) * t101 + t110;
t52 = -t93 * t79 + t91 * t80;
t53 = t91 * t79 + t93 * t80;
t23 = (qJD(3) * t71 - t53) * qJ(5) + (qJD(3) * t72 + t52) * pkin(4) + t104;
t90 = sin(pkin(10));
t92 = cos(pkin(10));
t57 = t92 * qJD(3) - t90 * t72;
t120 = 0.2e1 * qJD(5) * t57 + t92 * t20 + t90 * t23;
t19 = -qJDD(3) * pkin(4) - t100 * qJ(5) + t72 * t48 + qJDD(5) - t131;
t58 = t90 * qJD(3) + t92 * t72;
t94 = sin(qJ(6));
t97 = cos(qJ(6));
t38 = t94 * t57 + t97 * t58;
t45 = t92 * qJDD(3) - t90 * t53;
t46 = t90 * qJDD(3) + t92 * t53;
t25 = -t38 * qJD(6) + t97 * t45 - t94 * t46;
t37 = t97 * t57 - t94 * t58;
t26 = t37 * qJD(6) + t94 * t45 + t97 * t46;
t69 = qJD(6) + t71;
t29 = -t69 * mrSges(7,2) + t37 * mrSges(7,3);
t30 = t69 * mrSges(7,1) - t38 * mrSges(7,3);
t44 = t71 * pkin(5) - t58 * pkin(8);
t56 = t57 ^ 2;
t107 = t25 * mrSges(7,1) + t37 * t29 - m(7) * (-t45 * pkin(5) - t56 * pkin(8) + t58 * t44 + t19) - t26 * mrSges(7,2) - t38 * t30;
t42 = -t71 * mrSges(6,2) + t57 * mrSges(6,3);
t43 = t71 * mrSges(6,1) - t58 * mrSges(6,3);
t102 = m(6) * t19 - t45 * mrSges(6,1) + t46 * mrSges(6,2) - t57 * t42 + t58 * t43 - t107;
t49 = t71 * mrSges(5,1) + t72 * mrSges(5,2);
t63 = -qJD(3) * mrSges(5,2) - t71 * mrSges(5,3);
t11 = m(5) * t131 + qJDD(3) * mrSges(5,1) - t53 * mrSges(5,3) + qJD(3) * t63 - t72 * t49 - t102;
t113 = -0.2e1 * qJD(5) * t58 - t90 * t20 + t92 * t23;
t14 = (t57 * t71 - t46) * pkin(8) + (t57 * t58 + t52) * pkin(5) + t113;
t15 = -t56 * pkin(5) + t45 * pkin(8) - t71 * t44 + t120;
t28 = -t37 * mrSges(7,1) + t38 * mrSges(7,2);
t51 = qJDD(6) + t52;
t12 = m(7) * (t97 * t14 - t94 * t15) - t26 * mrSges(7,3) + t51 * mrSges(7,1) - t38 * t28 + t69 * t29;
t13 = m(7) * (t94 * t14 + t97 * t15) + t25 * mrSges(7,3) - t51 * mrSges(7,2) + t37 * t28 - t69 * t30;
t40 = -t57 * mrSges(6,1) + t58 * mrSges(6,2);
t10 = m(6) * t120 - t52 * mrSges(6,2) + t45 * mrSges(6,3) - t94 * t12 + t97 * t13 + t57 * t40 - t71 * t43;
t64 = qJD(3) * mrSges(5,1) - t72 * mrSges(5,3);
t9 = m(6) * t113 + t52 * mrSges(6,1) - t46 * mrSges(6,3) + t97 * t12 + t94 * t13 - t58 * t40 + t71 * t42;
t6 = m(5) * t119 - qJDD(3) * mrSges(5,2) - t52 * mrSges(5,3) - qJD(3) * t64 + t92 * t10 - t71 * t49 - t90 * t9;
t78 = (mrSges(4,1) * t95 + mrSges(4,2) * t98) * qJD(1);
t81 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t125;
t3 = m(4) * t126 + qJDD(3) * mrSges(4,1) - t80 * mrSges(4,3) + qJD(3) * t81 + t93 * t11 - t78 * t124 + t91 * t6;
t83 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t124;
t4 = m(4) * t116 - qJDD(3) * mrSges(4,2) + t79 * mrSges(4,3) - qJD(3) * t83 - t91 * t11 - t78 * t125 + t93 * t6;
t118 = -t95 * t3 + t98 * t4;
t109 = -m(3) * (-qJDD(1) * pkin(1) + t108) - t98 * t3 - t95 * t4;
t106 = m(5) * t104 + t52 * mrSges(5,1) + t53 * mrSges(5,2) + t90 * t10 + t71 * t63 + t72 * t64 + t92 * t9;
t105 = -t79 * mrSges(4,1) + t106 + m(4) * (t129 * t101 + t110) + t81 * t125 + t83 * t124 + t80 * mrSges(4,2);
t103 = -m(3) * (t101 * pkin(1) - t110) + t105;
t5 = m(2) * t112 + t127 * qJDD(1) - (t128 * t101) + t103;
t1 = m(2) * t117 + t128 * qJDD(1) + t127 * t101 + t109;
t2 = [-m(1) * g(1) - t96 * t1 + t99 * t5, t5, -m(3) * g(3) + t118, t4, t6, t10, t13; -m(1) * g(2) + t99 * t1 + t96 * t5, t1, -(t101 * mrSges(3,2)) - qJDD(1) * mrSges(3,3) - t103, t3, t11, t9, t12; (-m(1) + t130) * g(3) + t118, t130 * g(3) + t118, qJDD(1) * mrSges(3,2) - t101 * mrSges(3,3) - t109, t105, t106, t102, -t107;];
f_new  = t2;
