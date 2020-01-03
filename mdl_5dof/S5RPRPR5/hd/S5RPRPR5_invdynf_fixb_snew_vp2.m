% Calculate vector of cutting forces with Newton-Euler
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x6]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRPR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR5_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR5_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR5_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:41:37
% EndTime: 2020-01-03 11:41:44
% DurationCPUTime: 2.01s
% Computational Cost: add. (18677->153), mult. (47818->212), div. (0->0), fcn. (32490->10), ass. (0->88)
t80 = sin(pkin(8));
t76 = t80 ^ 2;
t82 = cos(pkin(8));
t118 = (t82 ^ 2 + t76) * mrSges(3,3);
t109 = qJDD(1) * mrSges(3,3);
t85 = sin(qJ(1));
t88 = cos(qJ(1));
t101 = -t85 * g(2) + t88 * g(3);
t89 = qJD(1) ^ 2;
t63 = -t89 * pkin(1) + qJDD(1) * qJ(2) + t101;
t114 = -t82 * g(1) - t80 * t63;
t97 = -t82 * mrSges(3,1) + t80 * mrSges(3,2);
t64 = t97 * qJD(1);
t103 = 0.2e1 * qJD(1) * qJD(2);
t111 = qJD(1) * t80;
t98 = -pkin(2) * t82 - pkin(6) * t80;
t65 = t98 * qJD(1);
t35 = t80 * t103 + t65 * t111 - t114;
t108 = qJD(1) * qJD(3);
t84 = sin(qJ(3));
t87 = cos(qJ(3));
t59 = (-qJDD(1) * t84 - t87 * t108) * t80;
t60 = (qJDD(1) * t87 - t84 * t108) * t80;
t79 = sin(pkin(9));
t81 = cos(pkin(9));
t39 = t81 * t59 - t79 * t60;
t40 = t79 * t59 + t81 * t60;
t51 = (-t79 * t87 - t81 * t84) * t111;
t110 = t82 * qJD(1);
t69 = qJD(3) - t110;
t41 = -t69 * mrSges(5,2) + t51 * mrSges(5,3);
t52 = (-t79 * t84 + t81 * t87) * t111;
t42 = t69 * mrSges(5,1) - t52 * mrSges(5,3);
t104 = t87 * t111;
t116 = t76 * t89;
t107 = t84 ^ 2 * t116;
t102 = qJ(4) * t111;
t56 = t69 * pkin(3) - t87 * t102;
t92 = -t59 * pkin(3) - qJ(4) * t107 + t56 * t104 + qJDD(4) + t35;
t83 = sin(qJ(5));
t86 = cos(qJ(5));
t34 = t83 * t51 + t86 * t52;
t18 = -t34 * qJD(5) + t86 * t39 - t83 * t40;
t33 = t86 * t51 - t83 * t52;
t19 = t33 * qJD(5) + t83 * t39 + t86 * t40;
t67 = qJD(5) + t69;
t29 = -t67 * mrSges(6,2) + t33 * mrSges(6,3);
t30 = t67 * mrSges(6,1) - t34 * mrSges(6,3);
t43 = t69 * pkin(4) - t52 * pkin(7);
t50 = t51 ^ 2;
t93 = t18 * mrSges(6,1) + t33 * t29 - m(6) * (-t39 * pkin(4) - t50 * pkin(7) + t52 * t43 + t92) - t19 * mrSges(6,2) - t34 * t30;
t91 = -m(5) * t92 + t39 * mrSges(5,1) - t40 * mrSges(5,2) + t51 * t41 - t52 * t42 + t93;
t90 = m(4) * t35 - t59 * mrSges(4,1) + t60 * mrSges(4,2) - t91;
t105 = t84 * t111;
t55 = -t69 * mrSges(4,2) - mrSges(4,3) * t105;
t57 = t69 * mrSges(4,1) - mrSges(4,3) * t104;
t96 = t55 * t84 + t57 * t87;
t10 = -t90 + (-t109 + (-0.2e1 * m(3) * qJD(2) - t64 - t96) * qJD(1)) * t80 + m(3) * t114;
t100 = -t80 * g(1) + (t103 + t63) * t82;
t36 = t65 * t110 + t100;
t113 = -t88 * g(2) - t85 * g(3);
t95 = -t89 * qJ(2) + qJDD(2) - t113;
t46 = (-pkin(1) + t98) * qJDD(1) + t95;
t45 = t87 * t46;
t58 = (mrSges(4,1) * t84 + mrSges(4,2) * t87) * t111;
t68 = -t82 * qJDD(1) + qJDD(3);
t23 = t68 * pkin(3) - t60 * qJ(4) + t45 + (-pkin(3) * t87 * t116 - t69 * t102 - t36) * t84;
t115 = t87 * t36 + t84 * t46;
t24 = -pkin(3) * t107 + t59 * qJ(4) - t69 * t56 + t115;
t99 = -0.2e1 * qJD(4) * t52 + t81 * t23 - t79 * t24;
t13 = (t51 * t69 - t40) * pkin(7) + (t51 * t52 + t68) * pkin(4) + t99;
t106 = 0.2e1 * qJD(4) * t51 + t79 * t23 + t81 * t24;
t14 = -t50 * pkin(4) + t39 * pkin(7) - t69 * t43 + t106;
t26 = -t33 * mrSges(6,1) + t34 * mrSges(6,2);
t66 = qJDD(5) + t68;
t11 = m(6) * (t86 * t13 - t83 * t14) - t19 * mrSges(6,3) + t66 * mrSges(6,1) - t34 * t26 + t67 * t29;
t12 = m(6) * (t83 * t13 + t86 * t14) + t18 * mrSges(6,3) - t66 * mrSges(6,2) + t33 * t26 - t67 * t30;
t37 = -t51 * mrSges(5,1) + t52 * mrSges(5,2);
t7 = m(5) * t99 + t68 * mrSges(5,1) - t40 * mrSges(5,3) + t86 * t11 + t83 * t12 - t52 * t37 + t69 * t41;
t8 = m(5) * t106 - t68 * mrSges(5,2) + t39 * mrSges(5,3) - t83 * t11 + t86 * t12 + t51 * t37 - t69 * t42;
t5 = m(4) * (-t84 * t36 + t45) - t60 * mrSges(4,3) + t68 * mrSges(4,1) - t58 * t104 + t69 * t55 + t79 * t8 + t81 * t7;
t6 = m(4) * t115 - t68 * mrSges(4,2) + t59 * mrSges(4,3) - t58 * t105 - t69 * t57 - t79 * t7 + t81 * t8;
t4 = m(3) * t100 + t87 * t6 - t84 * t5 + (qJD(1) * t64 + t109) * t82;
t117 = t82 * t10 + t80 * t4;
t94 = m(3) * (-qJDD(1) * pkin(1) + t95) + t87 * t5 + t84 * t6;
t2 = m(2) * t113 + (-mrSges(2,2) + t118) * t89 + (mrSges(2,1) - t97) * qJDD(1) - t94;
t1 = m(2) * t101 - t89 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t80 * t10 + t82 * t4;
t3 = [(-m(1) - m(2)) * g(1) + t117, t1, t4, t6, t8, t12; -m(1) * g(2) + t85 * t1 + t88 * t2, t2, t10, t5, t7, t11; -m(1) * g(3) - t88 * t1 + t85 * t2, -m(2) * g(1) + t117, t97 * qJDD(1) - t89 * t118 + t94, t96 * t111 + t90, -t91, -t93;];
f_new = t3;
