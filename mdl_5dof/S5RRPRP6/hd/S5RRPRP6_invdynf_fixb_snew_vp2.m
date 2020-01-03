% Calculate vector of cutting forces with Newton-Euler
% S5RRPRP6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRPRP6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP6_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP6_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP6_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP6_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP6_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:56:49
% EndTime: 2019-12-31 19:56:51
% DurationCPUTime: 1.17s
% Computational Cost: add. (10512->163), mult. (23822->209), div. (0->0), fcn. (15694->8), ass. (0->76)
t114 = -2 * qJD(3);
t101 = qJD(1) * qJD(2);
t84 = qJD(1) ^ 2;
t79 = sin(qJ(1));
t82 = cos(qJ(1));
t92 = -t82 * g(1) - t79 * g(2);
t65 = -t84 * pkin(1) + qJDD(1) * pkin(6) + t92;
t78 = sin(qJ(2));
t108 = t78 * t65;
t109 = pkin(2) * t84;
t81 = cos(qJ(2));
t68 = t78 * qJDD(1) + t81 * t101;
t33 = qJDD(2) * pkin(2) - t68 * qJ(3) - t108 + (qJ(3) * t101 + t78 * t109 - g(3)) * t81;
t69 = t81 * qJDD(1) - t78 * t101;
t104 = qJD(1) * t78;
t70 = qJD(2) * pkin(2) - qJ(3) * t104;
t74 = t81 ^ 2;
t95 = -t78 * g(3) + t81 * t65;
t34 = t69 * qJ(3) - qJD(2) * t70 - t74 * t109 + t95;
t75 = sin(pkin(8));
t76 = cos(pkin(8));
t63 = (t75 * t81 + t76 * t78) * qJD(1);
t113 = t63 * t114 + t76 * t33 - t75 * t34;
t62 = (t75 * t78 - t76 * t81) * qJD(1);
t71 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t104;
t103 = qJD(1) * t81;
t72 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t103;
t47 = t62 * pkin(3) - t63 * pkin(7);
t83 = qJD(2) ^ 2;
t97 = t62 * t114 + t75 * t33 + t76 * t34;
t19 = -t83 * pkin(3) + qJDD(2) * pkin(7) - t62 * t47 + t97;
t51 = -t75 * t68 + t76 * t69;
t52 = t76 * t68 + t75 * t69;
t96 = t79 * g(1) - t82 * g(2);
t89 = -qJDD(1) * pkin(1) - t96;
t86 = -t69 * pkin(2) + qJDD(3) + t70 * t104 + (-qJ(3) * t74 - pkin(6)) * t84 + t89;
t22 = (qJD(2) * t62 - t52) * pkin(7) + (qJD(2) * t63 - t51) * pkin(3) + t86;
t77 = sin(qJ(4));
t80 = cos(qJ(4));
t106 = t80 * t19 + t77 * t22;
t55 = t77 * qJD(2) + t80 * t63;
t27 = -t55 * qJD(4) + t80 * qJDD(2) - t77 * t52;
t54 = t80 * qJD(2) - t77 * t63;
t36 = -t54 * mrSges(5,1) + t55 * mrSges(5,2);
t61 = qJD(4) + t62;
t43 = t61 * mrSges(6,1) - t55 * mrSges(6,3);
t44 = t61 * mrSges(5,1) - t55 * mrSges(5,3);
t50 = qJDD(4) - t51;
t35 = -t54 * mrSges(6,1) + t55 * mrSges(6,2);
t42 = t61 * pkin(4) - t55 * qJ(5);
t53 = t54 ^ 2;
t99 = m(6) * (-t53 * pkin(4) + t27 * qJ(5) + 0.2e1 * qJD(5) * t54 - t61 * t42 + t106) + t54 * t35 + t27 * mrSges(6,3);
t11 = m(5) * t106 + t27 * mrSges(5,3) + t54 * t36 + (-t44 - t43) * t61 + (-mrSges(5,2) - mrSges(6,2)) * t50 + t99;
t56 = -qJD(2) * mrSges(4,2) - t62 * mrSges(4,3);
t57 = qJD(2) * mrSges(4,1) - t63 * mrSges(4,3);
t28 = t54 * qJD(4) + t77 * qJDD(2) + t80 * t52;
t40 = -t61 * mrSges(6,2) + t54 * mrSges(6,3);
t94 = -t77 * t19 + t80 * t22;
t100 = m(6) * (-0.2e1 * qJD(5) * t55 + (t54 * t61 - t28) * qJ(5) + (t54 * t55 + t50) * pkin(4) + t94) + t61 * t40 + t50 * mrSges(6,1);
t41 = -t61 * mrSges(5,2) + t54 * mrSges(5,3);
t9 = m(5) * t94 + t50 * mrSges(5,1) + t61 * t41 + (-t36 - t35) * t55 + (-mrSges(5,3) - mrSges(6,3)) * t28 + t100;
t88 = -m(4) * t86 + t51 * mrSges(4,1) - t52 * mrSges(4,2) - t77 * t11 - t62 * t56 - t63 * t57 - t80 * t9;
t112 = (t78 * t71 - t81 * t72) * qJD(1) + m(3) * (-t84 * pkin(6) + t89) - t69 * mrSges(3,1) + t68 * mrSges(3,2) - t88;
t18 = -qJDD(2) * pkin(3) - t83 * pkin(7) + t63 * t47 - t113;
t98 = m(6) * (-t27 * pkin(4) - t53 * qJ(5) + t55 * t42 + qJDD(5) + t18) + t55 * t43 + t28 * mrSges(6,2);
t111 = -m(5) * t18 - t28 * mrSges(5,2) + (t41 + t40) * t54 + (mrSges(5,1) + mrSges(6,1)) * t27 - t55 * t44 - t98;
t46 = t62 * mrSges(4,1) + t63 * mrSges(4,2);
t12 = m(4) * t113 + qJDD(2) * mrSges(4,1) - t52 * mrSges(4,3) + qJD(2) * t56 - t63 * t46 + t111;
t67 = (-mrSges(3,1) * t81 + mrSges(3,2) * t78) * qJD(1);
t7 = m(4) * t97 - qJDD(2) * mrSges(4,2) + t51 * mrSges(4,3) - qJD(2) * t57 + t80 * t11 - t62 * t46 - t77 * t9;
t4 = m(3) * (-t81 * g(3) - t108) - t68 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t67 * t104 + qJD(2) * t72 + t75 * t7 + t76 * t12;
t5 = m(3) * t95 - qJDD(2) * mrSges(3,2) + t69 * mrSges(3,3) - qJD(2) * t71 + t67 * t103 - t75 * t12 + t76 * t7;
t110 = t81 * t4 + t78 * t5;
t6 = m(2) * t96 + qJDD(1) * mrSges(2,1) - t84 * mrSges(2,2) - t112;
t1 = m(2) * t92 - t84 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t78 * t4 + t81 * t5;
t2 = [-m(1) * g(1) + t82 * t1 - t79 * t6, t1, t5, t7, t11, -t50 * mrSges(6,2) - t61 * t43 + t99; -m(1) * g(2) + t79 * t1 + t82 * t6, t6, t4, t12, t9, -t28 * mrSges(6,3) - t55 * t35 + t100; (-m(1) - m(2)) * g(3) + t110, -m(2) * g(3) + t110, t112, -t88, -t111, -t27 * mrSges(6,1) - t54 * t40 + t98;];
f_new = t2;
