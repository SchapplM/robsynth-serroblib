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
% Datum: 2019-12-05 17:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:55:59
% EndTime: 2019-12-05 17:56:05
% DurationCPUTime: 1.87s
% Computational Cost: add. (18677->153), mult. (47818->212), div. (0->0), fcn. (32490->10), ass. (0->88)
t78 = sin(pkin(8));
t74 = t78 ^ 2;
t80 = cos(pkin(8));
t116 = (t80 ^ 2 + t74) * mrSges(3,3);
t108 = qJDD(1) * mrSges(3,3);
t83 = sin(qJ(1));
t86 = cos(qJ(1));
t100 = t83 * g(2) - t86 * g(3);
t87 = qJD(1) ^ 2;
t63 = -t87 * pkin(1) + qJDD(1) * qJ(2) + t100;
t112 = -t80 * g(1) - t78 * t63;
t95 = -t80 * mrSges(3,1) + t78 * mrSges(3,2);
t64 = t95 * qJD(1);
t102 = 0.2e1 * qJD(1) * qJD(2);
t110 = qJD(1) * t78;
t97 = -pkin(2) * t80 - pkin(6) * t78;
t65 = t97 * qJD(1);
t35 = t78 * t102 + t65 * t110 - t112;
t107 = qJD(1) * qJD(3);
t82 = sin(qJ(3));
t85 = cos(qJ(3));
t59 = (-qJDD(1) * t82 - t85 * t107) * t78;
t60 = (qJDD(1) * t85 - t82 * t107) * t78;
t77 = sin(pkin(9));
t79 = cos(pkin(9));
t39 = t79 * t59 - t77 * t60;
t40 = t77 * t59 + t79 * t60;
t51 = (-t77 * t85 - t79 * t82) * t110;
t109 = t80 * qJD(1);
t69 = qJD(3) - t109;
t41 = -t69 * mrSges(5,2) + t51 * mrSges(5,3);
t52 = (-t77 * t82 + t79 * t85) * t110;
t42 = t69 * mrSges(5,1) - t52 * mrSges(5,3);
t103 = t85 * t110;
t114 = t74 * t87;
t106 = t82 ^ 2 * t114;
t101 = qJ(4) * t110;
t56 = t69 * pkin(3) - t85 * t101;
t90 = -t59 * pkin(3) - qJ(4) * t106 + t56 * t103 + qJDD(4) + t35;
t81 = sin(qJ(5));
t84 = cos(qJ(5));
t34 = t81 * t51 + t84 * t52;
t18 = -t34 * qJD(5) + t84 * t39 - t81 * t40;
t33 = t84 * t51 - t81 * t52;
t19 = t33 * qJD(5) + t81 * t39 + t84 * t40;
t67 = qJD(5) + t69;
t29 = -t67 * mrSges(6,2) + t33 * mrSges(6,3);
t30 = t67 * mrSges(6,1) - t34 * mrSges(6,3);
t43 = t69 * pkin(4) - t52 * pkin(7);
t50 = t51 ^ 2;
t92 = t18 * mrSges(6,1) + t33 * t29 - m(6) * (-t39 * pkin(4) - t50 * pkin(7) + t52 * t43 + t90) - t19 * mrSges(6,2) - t34 * t30;
t89 = -m(5) * t90 + t39 * mrSges(5,1) - t40 * mrSges(5,2) + t51 * t41 - t52 * t42 + t92;
t88 = m(4) * t35 - t59 * mrSges(4,1) + t60 * mrSges(4,2) - t89;
t104 = t82 * t110;
t55 = -t69 * mrSges(4,2) - mrSges(4,3) * t104;
t57 = t69 * mrSges(4,1) - mrSges(4,3) * t103;
t94 = t55 * t82 + t57 * t85;
t10 = (-t108 + (-0.2e1 * m(3) * qJD(2) - t64 - t94) * qJD(1)) * t78 + m(3) * t112 - t88;
t99 = -t78 * g(1) + (t102 + t63) * t80;
t36 = t65 * t109 + t99;
t96 = t86 * g(2) + t83 * g(3);
t91 = -t87 * qJ(2) + qJDD(2) - t96;
t46 = (-pkin(1) + t97) * qJDD(1) + t91;
t45 = t85 * t46;
t58 = (mrSges(4,1) * t82 + mrSges(4,2) * t85) * t110;
t68 = -t80 * qJDD(1) + qJDD(3);
t23 = t68 * pkin(3) - t60 * qJ(4) + t45 + (-pkin(3) * t85 * t114 - t69 * t101 - t36) * t82;
t113 = t85 * t36 + t82 * t46;
t24 = -pkin(3) * t106 + t59 * qJ(4) - t69 * t56 + t113;
t98 = -0.2e1 * qJD(4) * t52 + t79 * t23 - t77 * t24;
t13 = (t51 * t69 - t40) * pkin(7) + (t51 * t52 + t68) * pkin(4) + t98;
t105 = 0.2e1 * qJD(4) * t51 + t77 * t23 + t79 * t24;
t14 = -t50 * pkin(4) + t39 * pkin(7) - t69 * t43 + t105;
t26 = -t33 * mrSges(6,1) + t34 * mrSges(6,2);
t66 = qJDD(5) + t68;
t11 = m(6) * (t84 * t13 - t81 * t14) - t19 * mrSges(6,3) + t66 * mrSges(6,1) - t34 * t26 + t67 * t29;
t12 = m(6) * (t81 * t13 + t84 * t14) + t18 * mrSges(6,3) - t66 * mrSges(6,2) + t33 * t26 - t67 * t30;
t37 = -t51 * mrSges(5,1) + t52 * mrSges(5,2);
t7 = m(5) * t98 + t68 * mrSges(5,1) - t40 * mrSges(5,3) + t84 * t11 + t81 * t12 - t52 * t37 + t69 * t41;
t8 = m(5) * t105 - t68 * mrSges(5,2) + t39 * mrSges(5,3) - t81 * t11 + t84 * t12 + t51 * t37 - t69 * t42;
t5 = m(4) * (-t82 * t36 + t45) - t60 * mrSges(4,3) + t68 * mrSges(4,1) - t58 * t103 + t69 * t55 + t77 * t8 + t79 * t7;
t6 = m(4) * t113 - t68 * mrSges(4,2) + t59 * mrSges(4,3) - t58 * t104 - t69 * t57 - t77 * t7 + t79 * t8;
t4 = m(3) * t99 + t85 * t6 - t82 * t5 + (qJD(1) * t64 + t108) * t80;
t115 = t80 * t10 + t78 * t4;
t93 = m(3) * (-qJDD(1) * pkin(1) + t91) + t85 * t5 + t82 * t6;
t2 = m(2) * t96 + (-mrSges(2,2) + t116) * t87 + (mrSges(2,1) - t95) * qJDD(1) - t93;
t1 = m(2) * t100 - t87 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t78 * t10 + t80 * t4;
t3 = [(-m(1) - m(2)) * g(1) + t115, t1, t4, t6, t8, t12; -m(1) * g(2) - t83 * t1 - t86 * t2, t2, t10, t5, t7, t11; -m(1) * g(3) + t86 * t1 - t83 * t2, -m(2) * g(1) + t115, t95 * qJDD(1) - t87 * t116 + t93, t94 * t110 + t88, -t89, -t92;];
f_new = t3;
