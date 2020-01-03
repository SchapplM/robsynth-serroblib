% Calculate vector of cutting forces with Newton-Euler
% S5RPRRP10
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRRP10_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP10_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP10_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP10_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP10_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP10_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP10_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:50:57
% EndTime: 2019-12-31 18:51:00
% DurationCPUTime: 1.01s
% Computational Cost: add. (9500->150), mult. (22386->186), div. (0->0), fcn. (15664->8), ass. (0->76)
t79 = qJD(1) ^ 2;
t71 = cos(pkin(8));
t69 = t71 ^ 2;
t70 = sin(pkin(8));
t102 = t70 ^ 2 + t69;
t110 = t102 * mrSges(3,3);
t73 = sin(qJ(3));
t76 = cos(qJ(3));
t85 = t70 * t73 - t71 * t76;
t61 = t85 * qJD(1);
t86 = t70 * t76 + t71 * t73;
t62 = t86 * qJD(1);
t99 = t62 * qJD(3);
t51 = -t85 * qJDD(1) - t99;
t49 = pkin(3) * t61 - pkin(7) * t62;
t78 = qJD(3) ^ 2;
t101 = pkin(6) * qJDD(1);
t107 = pkin(2) * t79;
t74 = sin(qJ(1));
t77 = cos(qJ(1));
t89 = -g(1) * t77 - g(2) * t74;
t63 = -pkin(1) * t79 + qJDD(1) * qJ(2) + t89;
t98 = qJD(1) * qJD(2);
t93 = -g(3) * t71 - 0.2e1 * t70 * t98;
t41 = (t71 * t107 - t101 - t63) * t70 + t93;
t90 = -g(3) * t70 + (t63 + 0.2e1 * t98) * t71;
t42 = t71 * t101 - t69 * t107 + t90;
t91 = t41 * t76 - t73 * t42;
t18 = -qJDD(3) * pkin(3) - pkin(7) * t78 + t62 * t49 - t91;
t100 = t61 * qJD(3);
t52 = t86 * qJDD(1) - t100;
t72 = sin(qJ(4));
t75 = cos(qJ(4));
t55 = qJD(3) * t72 + t62 * t75;
t27 = -t55 * qJD(4) + qJDD(3) * t75 - t52 * t72;
t54 = qJD(3) * t75 - t62 * t72;
t28 = t54 * qJD(4) + qJDD(3) * t72 + t52 * t75;
t59 = qJD(4) + t61;
t34 = -mrSges(6,2) * t59 + mrSges(6,3) * t54;
t35 = -mrSges(5,2) * t59 + mrSges(5,3) * t54;
t38 = mrSges(5,1) * t59 - mrSges(5,3) * t55;
t36 = pkin(4) * t59 - qJ(5) * t55;
t37 = mrSges(6,1) * t59 - mrSges(6,3) * t55;
t53 = t54 ^ 2;
t95 = m(6) * (-t27 * pkin(4) - t53 * qJ(5) + t55 * t36 + qJDD(5) + t18) + t28 * mrSges(6,2) + t55 * t37;
t109 = m(5) * t18 + t28 * mrSges(5,2) - (t35 + t34) * t54 - (mrSges(5,1) + mrSges(6,1)) * t27 + t55 * t38 + t95;
t46 = mrSges(4,1) * t61 + mrSges(4,2) * t62;
t56 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t61;
t12 = m(4) * t91 + qJDD(3) * mrSges(4,1) - t52 * mrSges(4,3) + qJD(3) * t56 - t62 * t46 - t109;
t103 = t73 * t41 + t76 * t42;
t19 = -pkin(3) * t78 + qJDD(3) * pkin(7) - t49 * t61 + t103;
t94 = g(1) * t74 - t77 * g(2);
t88 = qJDD(2) - t94;
t80 = (-pkin(2) * t71 - pkin(1)) * qJDD(1) + (-t102 * pkin(6) - qJ(2)) * t79 + t88;
t22 = (-t52 + t100) * pkin(7) + (-t51 + t99) * pkin(3) + t80;
t105 = t75 * t19 + t72 * t22;
t31 = -mrSges(5,1) * t54 + mrSges(5,2) * t55;
t48 = qJDD(4) - t51;
t30 = -mrSges(6,1) * t54 + mrSges(6,2) * t55;
t96 = m(6) * (-pkin(4) * t53 + qJ(5) * t27 + 0.2e1 * qJD(5) * t54 - t36 * t59 + t105) + t54 * t30 + t27 * mrSges(6,3);
t11 = m(5) * t105 + t27 * mrSges(5,3) + t54 * t31 + (-t38 - t37) * t59 + (-mrSges(5,2) - mrSges(6,2)) * t48 + t96;
t57 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t62;
t92 = -t19 * t72 + t75 * t22;
t97 = m(6) * (-0.2e1 * qJD(5) * t55 + (t54 * t59 - t28) * qJ(5) + (t54 * t55 + t48) * pkin(4) + t92) + t59 * t34 + t48 * mrSges(6,1);
t9 = m(5) * t92 + t48 * mrSges(5,1) + t59 * t35 + (-t31 - t30) * t55 + (-mrSges(5,3) - mrSges(6,3)) * t28 + t97;
t7 = m(4) * t103 - qJDD(3) * mrSges(4,2) + t51 * mrSges(4,3) - qJD(3) * t57 + t75 * t11 - t61 * t46 - t72 * t9;
t87 = -t71 * mrSges(3,1) + t70 * mrSges(3,2);
t84 = mrSges(3,3) * qJDD(1) + t79 * t87;
t4 = m(3) * t93 + t73 * t7 + t76 * t12 + (-m(3) * t63 - t84) * t70;
t5 = m(3) * t90 - t73 * t12 + t76 * t7 + t84 * t71;
t108 = t71 * t4 + t70 * t5;
t83 = -m(4) * t80 + t51 * mrSges(4,1) - t52 * mrSges(4,2) - t72 * t11 - t61 * t56 - t62 * t57 - t75 * t9;
t81 = m(3) * (-qJDD(1) * pkin(1) - qJ(2) * t79 + t88) - t83;
t6 = m(2) * t94 + (-mrSges(2,2) + t110) * t79 + (mrSges(2,1) - t87) * qJDD(1) - t81;
t1 = m(2) * t89 - t79 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t70 * t4 + t71 * t5;
t2 = [-m(1) * g(1) + t1 * t77 - t6 * t74, t1, t5, t7, t11, -t48 * mrSges(6,2) - t59 * t37 + t96; -m(1) * g(2) + t1 * t74 + t6 * t77, t6, t4, t12, t9, -t28 * mrSges(6,3) - t55 * t30 + t97; (-m(1) - m(2)) * g(3) + t108, -m(2) * g(3) + t108, t87 * qJDD(1) - t79 * t110 + t81, -t83, t109, -t27 * mrSges(6,1) - t54 * t34 + t95;];
f_new = t2;
