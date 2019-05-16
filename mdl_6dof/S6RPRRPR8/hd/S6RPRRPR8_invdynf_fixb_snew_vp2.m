% Calculate vector of cutting forces with Newton-Euler
% S6RPRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
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
% Datum: 2019-05-05 23:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRPR8_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR8_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR8_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR8_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR8_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR8_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR8_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR8_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR8_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 23:08:06
% EndTime: 2019-05-05 23:08:13
% DurationCPUTime: 2.98s
% Computational Cost: add. (41531->179), mult. (84262->228), div. (0->0), fcn. (55836->10), ass. (0->94)
t94 = sin(qJ(1));
t98 = cos(qJ(1));
t112 = -t98 * g(1) - t94 * g(2);
t128 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t112;
t127 = -m(2) - m(3);
t126 = (-pkin(1) - pkin(7));
t125 = (mrSges(2,1) - mrSges(3,2));
t124 = -mrSges(2,2) + mrSges(3,3);
t100 = qJD(1) ^ 2;
t104 = (t126 * t100) - t128;
t121 = qJD(1) * qJD(3);
t93 = sin(qJ(3));
t115 = t93 * t121;
t97 = cos(qJ(3));
t83 = t97 * t121;
t77 = -t93 * qJDD(1) - t83;
t78 = qJDD(1) * t97 - t115;
t41 = (-t78 + t115) * pkin(8) + (-t77 + t83) * pkin(3) + t104;
t118 = t94 * g(1) - t98 * g(2);
t108 = -t100 * qJ(2) + qJDD(2) - t118;
t63 = t126 * qJDD(1) + t108;
t117 = -g(3) * t97 + t93 * t63;
t76 = (pkin(3) * t93 - pkin(8) * t97) * qJD(1);
t85 = t93 * qJD(1);
t99 = qJD(3) ^ 2;
t44 = -pkin(3) * t99 + qJDD(3) * pkin(8) - t76 * t85 + t117;
t92 = sin(qJ(4));
t96 = cos(qJ(4));
t123 = t92 * t41 + t96 * t44;
t122 = qJD(1) * t97;
t82 = t85 + qJD(4);
t114 = t96 * t41 - t92 * t44;
t73 = qJD(3) * t96 - t92 * t122;
t52 = qJD(4) * t73 + qJDD(3) * t92 + t78 * t96;
t72 = qJDD(4) - t77;
t74 = qJD(3) * t92 + t96 * t122;
t23 = (t73 * t82 - t52) * qJ(5) + (t73 * t74 + t72) * pkin(4) + t114;
t51 = -qJD(4) * t74 + qJDD(3) * t96 - t78 * t92;
t58 = pkin(4) * t82 - qJ(5) * t74;
t71 = t73 ^ 2;
t25 = -pkin(4) * t71 + qJ(5) * t51 - t58 * t82 + t123;
t89 = sin(pkin(10));
t90 = cos(pkin(10));
t54 = t73 * t90 - t74 * t89;
t119 = 0.2e1 * qJD(5) * t54 + t89 * t23 + t90 * t25;
t111 = g(3) * t93 + t63 * t97;
t43 = -qJDD(3) * pkin(3) - pkin(8) * t99 + t76 * t122 - t111;
t102 = -pkin(4) * t51 - qJ(5) * t71 + t74 * t58 + qJDD(5) + t43;
t33 = t51 * t90 - t52 * t89;
t34 = t51 * t89 + t52 * t90;
t55 = t73 * t89 + t74 * t90;
t91 = sin(qJ(6));
t95 = cos(qJ(6));
t37 = t54 * t91 + t55 * t95;
t19 = -qJD(6) * t37 + t33 * t95 - t34 * t91;
t36 = t54 * t95 - t55 * t91;
t20 = qJD(6) * t36 + t33 * t91 + t34 * t95;
t81 = qJD(6) + t82;
t30 = -mrSges(7,2) * t81 + mrSges(7,3) * t36;
t31 = mrSges(7,1) * t81 - mrSges(7,3) * t37;
t47 = pkin(5) * t82 - pkin(9) * t55;
t53 = t54 ^ 2;
t107 = t19 * mrSges(7,1) + t36 * t30 - m(7) * (-pkin(5) * t33 - pkin(9) * t53 + t47 * t55 + t102) - t20 * mrSges(7,2) - t37 * t31;
t45 = -mrSges(6,2) * t82 + mrSges(6,3) * t54;
t46 = mrSges(6,1) * t82 - mrSges(6,3) * t55;
t103 = -m(6) * t102 + t33 * mrSges(6,1) - t34 * mrSges(6,2) + t54 * t45 - t55 * t46 + t107;
t57 = -mrSges(5,2) * t82 + mrSges(5,3) * t73;
t59 = mrSges(5,1) * t82 - mrSges(5,3) * t74;
t101 = m(5) * t43 - t51 * mrSges(5,1) + t52 * mrSges(5,2) - t73 * t57 + t74 * t59 - t103;
t75 = (mrSges(4,1) * t93 + mrSges(4,2) * t97) * qJD(1);
t79 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t85;
t13 = m(4) * t111 + qJDD(3) * mrSges(4,1) - t78 * mrSges(4,3) + qJD(3) * t79 - t75 * t122 - t101;
t113 = -0.2e1 * qJD(5) * t55 + t90 * t23 - t89 * t25;
t14 = (t54 * t82 - t34) * pkin(9) + (t54 * t55 + t72) * pkin(5) + t113;
t15 = -pkin(5) * t53 + pkin(9) * t33 - t47 * t82 + t119;
t27 = -mrSges(7,1) * t36 + mrSges(7,2) * t37;
t67 = qJDD(6) + t72;
t11 = m(7) * (t14 * t95 - t15 * t91) - t20 * mrSges(7,3) + t67 * mrSges(7,1) - t37 * t27 + t81 * t30;
t12 = m(7) * (t14 * t91 + t15 * t95) + t19 * mrSges(7,3) - t67 * mrSges(7,2) + t36 * t27 - t81 * t31;
t38 = -mrSges(6,1) * t54 + mrSges(6,2) * t55;
t10 = m(6) * t119 - t72 * mrSges(6,2) + t33 * mrSges(6,3) - t91 * t11 + t95 * t12 + t54 * t38 - t82 * t46;
t56 = -mrSges(5,1) * t73 + mrSges(5,2) * t74;
t9 = m(6) * t113 + t72 * mrSges(6,1) - t34 * mrSges(6,3) + t95 * t11 + t91 * t12 - t55 * t38 + t82 * t45;
t7 = m(5) * t114 + t72 * mrSges(5,1) - t52 * mrSges(5,3) + t89 * t10 - t74 * t56 + t82 * t57 + t90 * t9;
t8 = m(5) * t123 - t72 * mrSges(5,2) + t51 * mrSges(5,3) + t90 * t10 + t73 * t56 - t82 * t59 - t89 * t9;
t80 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t122;
t4 = m(4) * t117 - qJDD(3) * mrSges(4,2) + t77 * mrSges(4,3) - qJD(3) * t80 - t92 * t7 - t75 * t85 + t96 * t8;
t116 = -t93 * t13 + t97 * t4;
t109 = -m(3) * (-qJDD(1) * pkin(1) + t108) - t97 * t13 - t93 * t4;
t106 = m(4) * t104 - t77 * mrSges(4,1) + t78 * mrSges(4,2) + t80 * t122 + t96 * t7 + t79 * t85 + t92 * t8;
t105 = -m(3) * (t100 * pkin(1) + t128) + t106;
t2 = m(2) * t112 + t124 * qJDD(1) - (t125 * t100) + t105;
t1 = m(2) * t118 + t125 * qJDD(1) + t124 * t100 + t109;
t3 = [-m(1) * g(1) - t1 * t94 + t2 * t98, t2, -m(3) * g(3) + t116, t4, t8, t10, t12; -m(1) * g(2) + t1 * t98 + t2 * t94, t1, -(t100 * mrSges(3,2)) - qJDD(1) * mrSges(3,3) - t105, t13, t7, t9, t11; (-m(1) + t127) * g(3) + t116, t127 * g(3) + t116, qJDD(1) * mrSges(3,2) - t100 * mrSges(3,3) - t109, t106, t101, -t103, -t107;];
f_new  = t3;
