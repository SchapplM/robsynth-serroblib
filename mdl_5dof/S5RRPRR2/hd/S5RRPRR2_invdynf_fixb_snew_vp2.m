% Calculate vector of cutting forces with Newton-Euler
% S5RRPRR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRPRR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR2_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR2_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR2_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR2_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR2_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:27:13
% EndTime: 2019-12-05 18:27:22
% DurationCPUTime: 3.00s
% Computational Cost: add. (34820->167), mult. (82810->224), div. (0->0), fcn. (58546->10), ass. (0->85)
t106 = qJD(1) * qJD(2);
t85 = sin(qJ(2));
t89 = cos(qJ(2));
t70 = t85 * qJDD(1) + t89 * t106;
t71 = t89 * qJDD(1) - t85 * t106;
t108 = qJD(1) * t85;
t73 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t108;
t107 = qJD(1) * t89;
t74 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t107;
t91 = qJD(1) ^ 2;
t81 = sin(pkin(9));
t82 = cos(pkin(9));
t55 = -t81 * t70 + t82 * t71;
t56 = t82 * t70 + t81 * t71;
t64 = (-t81 * t85 + t82 * t89) * qJD(1);
t57 = -qJD(2) * mrSges(4,2) + t64 * mrSges(4,3);
t65 = (t81 * t89 + t82 * t85) * qJD(1);
t58 = qJD(2) * mrSges(4,1) - t65 * mrSges(4,3);
t72 = qJD(2) * pkin(2) - qJ(3) * t108;
t80 = t89 ^ 2;
t86 = sin(qJ(1));
t90 = cos(qJ(1));
t104 = t86 * g(1) - t90 * g(2);
t98 = -qJDD(1) * pkin(1) - t104;
t95 = -t71 * pkin(2) + qJDD(3) + t72 * t108 + (-qJ(3) * t80 - pkin(6)) * t91 + t98;
t84 = sin(qJ(4));
t88 = cos(qJ(4));
t50 = t84 * t64 + t88 * t65;
t30 = -t50 * qJD(4) + t88 * t55 - t84 * t56;
t49 = t88 * t64 - t84 * t65;
t31 = t49 * qJD(4) + t84 * t55 + t88 * t56;
t79 = qJD(2) + qJD(4);
t45 = -t79 * mrSges(5,2) + t49 * mrSges(5,3);
t46 = t79 * mrSges(5,1) - t50 * mrSges(5,3);
t59 = qJD(2) * pkin(3) - t65 * pkin(7);
t63 = t64 ^ 2;
t94 = -t55 * pkin(3) - t63 * pkin(7) + t65 * t59 + t95;
t83 = sin(qJ(5));
t87 = cos(qJ(5));
t36 = t83 * t49 + t87 * t50;
t18 = -t36 * qJD(5) + t87 * t30 - t83 * t31;
t35 = t87 * t49 - t83 * t50;
t19 = t35 * qJD(5) + t83 * t30 + t87 * t31;
t76 = qJD(5) + t79;
t32 = -t76 * mrSges(6,2) + t35 * mrSges(6,3);
t33 = t76 * mrSges(6,1) - t36 * mrSges(6,3);
t47 = t79 * pkin(4) - t50 * pkin(8);
t48 = t49 ^ 2;
t97 = t18 * mrSges(6,1) + t35 * t32 - m(6) * (-t30 * pkin(4) - t48 * pkin(8) + t50 * t47 + t94) - t19 * mrSges(6,2) - t36 * t33;
t96 = -m(5) * t94 + t30 * mrSges(5,1) - t31 * mrSges(5,2) + t49 * t45 - t50 * t46 + t97;
t93 = -m(4) * t95 + t55 * mrSges(4,1) - t56 * mrSges(4,2) + t64 * t57 - t65 * t58 + t96;
t113 = (t85 * t73 - t89 * t74) * qJD(1) + m(3) * (-t91 * pkin(6) + t98) - t71 * mrSges(3,1) + t70 * mrSges(3,2) - t93;
t100 = -t90 * g(1) - t86 * g(2);
t67 = -t91 * pkin(1) + qJDD(1) * pkin(6) + t100;
t110 = t85 * t67;
t111 = pkin(2) * t91;
t41 = qJDD(2) * pkin(2) - t70 * qJ(3) - t110 + (qJ(3) * t106 + t85 * t111 - g(3)) * t89;
t103 = -t85 * g(3) + t89 * t67;
t42 = t71 * qJ(3) - qJD(2) * t72 - t80 * t111 + t103;
t101 = -0.2e1 * qJD(3) * t65 + t82 * t41 - t81 * t42;
t53 = -t64 * mrSges(4,1) + t65 * mrSges(4,2);
t22 = (qJD(2) * t64 - t56) * pkin(7) + (t64 * t65 + qJDD(2)) * pkin(3) + t101;
t105 = 0.2e1 * qJD(3) * t64 + t81 * t41 + t82 * t42;
t24 = -t63 * pkin(3) + t55 * pkin(7) - qJD(2) * t59 + t105;
t102 = t88 * t22 - t84 * t24;
t78 = qJDD(2) + qJDD(4);
t13 = (t49 * t79 - t31) * pkin(8) + (t49 * t50 + t78) * pkin(4) + t102;
t109 = t84 * t22 + t88 * t24;
t14 = -t48 * pkin(4) + t30 * pkin(8) - t79 * t47 + t109;
t26 = -t35 * mrSges(6,1) + t36 * mrSges(6,2);
t75 = qJDD(5) + t78;
t11 = m(6) * (t87 * t13 - t83 * t14) - t19 * mrSges(6,3) + t75 * mrSges(6,1) - t36 * t26 + t76 * t32;
t12 = m(6) * (t83 * t13 + t87 * t14) + t18 * mrSges(6,3) - t75 * mrSges(6,2) + t35 * t26 - t76 * t33;
t37 = -t49 * mrSges(5,1) + t50 * mrSges(5,2);
t8 = m(5) * t102 + t78 * mrSges(5,1) - t31 * mrSges(5,3) + t87 * t11 + t83 * t12 - t50 * t37 + t79 * t45;
t9 = m(5) * t109 - t78 * mrSges(5,2) + t30 * mrSges(5,3) - t83 * t11 + t87 * t12 + t49 * t37 - t79 * t46;
t6 = m(4) * t101 + qJDD(2) * mrSges(4,1) - t56 * mrSges(4,3) + qJD(2) * t57 - t65 * t53 + t88 * t8 + t84 * t9;
t69 = (-mrSges(3,1) * t89 + mrSges(3,2) * t85) * qJD(1);
t7 = m(4) * t105 - qJDD(2) * mrSges(4,2) + t55 * mrSges(4,3) - qJD(2) * t58 + t64 * t53 - t84 * t8 + t88 * t9;
t4 = m(3) * (-t89 * g(3) - t110) - t70 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t69 * t108 + qJD(2) * t74 + t81 * t7 + t82 * t6;
t5 = m(3) * t103 - qJDD(2) * mrSges(3,2) + t71 * mrSges(3,3) - qJD(2) * t73 + t69 * t107 - t81 * t6 + t82 * t7;
t112 = t89 * t4 + t85 * t5;
t10 = m(2) * t104 + qJDD(1) * mrSges(2,1) - t91 * mrSges(2,2) - t113;
t1 = m(2) * t100 - t91 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t85 * t4 + t89 * t5;
t2 = [-m(1) * g(1) + t90 * t1 - t86 * t10, t1, t5, t7, t9, t12; -m(1) * g(2) + t86 * t1 + t90 * t10, t10, t4, t6, t8, t11; (-m(1) - m(2)) * g(3) + t112, -m(2) * g(3) + t112, t113, -t93, -t96, -t97;];
f_new = t2;
