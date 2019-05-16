% Calculate vector of cutting forces with Newton-Euler
% S6RPRPRR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-05-05 20:04
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRPRR10_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR10_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR10_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR10_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR10_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR10_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR10_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 20:01:25
% EndTime: 2019-05-05 20:01:32
% DurationCPUTime: 2.84s
% Computational Cost: add. (38987->179), mult. (82387->228), div. (0->0), fcn. (54577->10), ass. (0->94)
t92 = sin(qJ(1));
t96 = cos(qJ(1));
t110 = -t96 * g(1) - t92 * g(2);
t127 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t110;
t126 = -m(2) - m(3);
t125 = (-pkin(1) - pkin(7));
t124 = (mrSges(2,1) - mrSges(3,2));
t123 = -mrSges(2,2) + mrSges(3,3);
t98 = qJD(1) ^ 2;
t102 = (t125 * t98) - t127;
t120 = qJD(1) * qJD(3);
t95 = cos(qJ(3));
t113 = t95 * t120;
t91 = sin(qJ(3));
t114 = t91 * t120;
t77 = qJDD(1) * t91 + t113;
t78 = qJDD(1) * t95 - t114;
t41 = (-t78 + t114) * qJ(4) + (t77 + t113) * pkin(3) + t102;
t117 = t92 * g(1) - g(2) * t96;
t106 = -t98 * qJ(2) + qJDD(2) - t117;
t61 = qJDD(1) * t125 + t106;
t116 = -g(3) * t95 + t61 * t91;
t75 = (pkin(3) * t91 - qJ(4) * t95) * qJD(1);
t83 = t91 * qJD(1);
t97 = qJD(3) ^ 2;
t44 = -pkin(3) * t97 + qJDD(3) * qJ(4) - t75 * t83 + t116;
t121 = qJD(1) * t95;
t87 = sin(pkin(10));
t88 = cos(pkin(10));
t72 = qJD(3) * t87 + t121 * t88;
t111 = -0.2e1 * qJD(4) * t72 + t41 * t88 - t87 * t44;
t59 = qJDD(3) * t87 + t78 * t88;
t71 = qJD(3) * t88 - t121 * t87;
t23 = (t71 * t83 - t59) * pkin(8) + (t71 * t72 + t77) * pkin(4) + t111;
t118 = 0.2e1 * qJD(4) * t71 + t41 * t87 + t44 * t88;
t58 = qJDD(3) * t88 - t78 * t87;
t60 = pkin(4) * t83 - pkin(8) * t72;
t70 = t71 ^ 2;
t25 = -pkin(4) * t70 + pkin(8) * t58 - t60 * t83 + t118;
t90 = sin(qJ(5));
t94 = cos(qJ(5));
t122 = t23 * t90 + t25 * t94;
t82 = t83 + qJD(5);
t109 = g(3) * t91 + t61 * t95;
t76 = (mrSges(4,1) * t91 + mrSges(4,2) * t95) * qJD(1);
t79 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t83;
t43 = -qJDD(3) * pkin(3) - qJ(4) * t97 + t121 * t75 + qJDD(4) - t109;
t100 = -pkin(4) * t58 - pkin(8) * t70 + t60 * t72 + t43;
t51 = t71 * t90 + t72 * t94;
t33 = -qJD(5) * t51 + t58 * t94 - t59 * t90;
t50 = t71 * t94 - t72 * t90;
t34 = qJD(5) * t50 + t58 * t90 + t59 * t94;
t89 = sin(qJ(6));
t93 = cos(qJ(6));
t37 = t50 * t89 + t51 * t93;
t19 = -qJD(6) * t37 + t33 * t93 - t34 * t89;
t36 = t50 * t93 - t51 * t89;
t20 = qJD(6) * t36 + t33 * t89 + t34 * t93;
t81 = qJD(6) + t82;
t30 = -mrSges(7,2) * t81 + mrSges(7,3) * t36;
t31 = mrSges(7,1) * t81 - mrSges(7,3) * t37;
t47 = pkin(5) * t82 - pkin(9) * t51;
t49 = t50 ^ 2;
t105 = t19 * mrSges(7,1) + t36 * t30 - m(7) * (-pkin(5) * t33 - pkin(9) * t49 + t47 * t51 + t100) - t20 * mrSges(7,2) - t37 * t31;
t45 = -mrSges(6,2) * t82 + mrSges(6,3) * t50;
t46 = mrSges(6,1) * t82 - mrSges(6,3) * t51;
t101 = -m(6) * t100 + t33 * mrSges(6,1) - mrSges(6,2) * t34 + t50 * t45 - t51 * t46 + t105;
t56 = -mrSges(5,2) * t83 + mrSges(5,3) * t71;
t57 = mrSges(5,1) * t83 - mrSges(5,3) * t72;
t99 = m(5) * t43 - t58 * mrSges(5,1) + t59 * mrSges(5,2) - t71 * t56 + t72 * t57 - t101;
t11 = m(4) * t109 + qJDD(3) * mrSges(4,1) - t78 * mrSges(4,3) + qJD(3) * t79 - t121 * t76 - t99;
t112 = t23 * t94 - t90 * t25;
t74 = qJDD(5) + t77;
t14 = (t50 * t82 - t34) * pkin(9) + (t50 * t51 + t74) * pkin(5) + t112;
t15 = -pkin(5) * t49 + pkin(9) * t33 - t47 * t82 + t122;
t27 = -mrSges(7,1) * t36 + mrSges(7,2) * t37;
t67 = qJDD(6) + t74;
t12 = m(7) * (t14 * t93 - t15 * t89) - t20 * mrSges(7,3) + t67 * mrSges(7,1) - t37 * t27 + t81 * t30;
t13 = m(7) * (t14 * t89 + t15 * t93) + t19 * mrSges(7,3) - t67 * mrSges(7,2) + t36 * t27 - t81 * t31;
t38 = -mrSges(6,1) * t50 + mrSges(6,2) * t51;
t10 = m(6) * t122 - mrSges(6,2) * t74 + mrSges(6,3) * t33 - t12 * t89 + t13 * t93 + t38 * t50 - t46 * t82;
t52 = -mrSges(5,1) * t71 + mrSges(5,2) * t72;
t9 = m(6) * t112 + mrSges(6,1) * t74 - mrSges(6,3) * t34 + t12 * t93 + t13 * t89 - t38 * t51 + t45 * t82;
t7 = m(5) * t111 + mrSges(5,1) * t77 - mrSges(5,3) * t59 + t10 * t90 - t52 * t72 + t56 * t83 + t9 * t94;
t8 = m(5) * t118 - mrSges(5,2) * t77 + mrSges(5,3) * t58 + t10 * t94 + t52 * t71 - t57 * t83 - t9 * t90;
t80 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t121;
t4 = m(4) * t116 - qJDD(3) * mrSges(4,2) - mrSges(4,3) * t77 - qJD(3) * t80 - t7 * t87 - t76 * t83 + t8 * t88;
t115 = -t91 * t11 + t4 * t95;
t107 = -m(3) * (-qJDD(1) * pkin(1) + t106) - t95 * t11 - t91 * t4;
t104 = m(4) * t102 + t77 * mrSges(4,1) + mrSges(4,2) * t78 + t121 * t80 + t7 * t88 + t79 * t83 + t8 * t87;
t103 = -m(3) * (t98 * pkin(1) + t127) + t104;
t2 = m(2) * t110 + qJDD(1) * t123 - (t124 * t98) + t103;
t1 = m(2) * t117 + qJDD(1) * t124 + t123 * t98 + t107;
t3 = [-m(1) * g(1) - t1 * t92 + t2 * t96, t2, -m(3) * g(3) + t115, t4, t8, t10, t13; -m(1) * g(2) + t1 * t96 + t2 * t92, t1, -(t98 * mrSges(3,2)) - qJDD(1) * mrSges(3,3) - t103, t11, t7, t9, t12; (-m(1) + t126) * g(3) + t115, g(3) * t126 + t115, qJDD(1) * mrSges(3,2) - t98 * mrSges(3,3) - t107, t104, t99, -t101, -t105;];
f_new  = t3;
