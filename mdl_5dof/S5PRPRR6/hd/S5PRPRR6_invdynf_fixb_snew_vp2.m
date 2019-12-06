% Calculate vector of cutting forces with Newton-Euler
% S5PRPRR6
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRPRR6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR6_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR6_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR6_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR6_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR6_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:56:37
% EndTime: 2019-12-05 15:56:39
% DurationCPUTime: 1.29s
% Computational Cost: add. (14403->125), mult. (31323->170), div. (0->0), fcn. (22787->12), ass. (0->79)
t77 = qJD(2) ^ 2;
t67 = cos(pkin(10));
t62 = t67 ^ 2;
t64 = sin(pkin(10));
t96 = t64 ^ 2 + t62;
t104 = t96 * mrSges(4,3);
t71 = sin(qJ(4));
t74 = cos(qJ(4));
t84 = t64 * t71 - t67 * t74;
t49 = t84 * qJD(2);
t65 = sin(pkin(9));
t68 = cos(pkin(9));
t55 = g(1) * t65 - g(2) * t68;
t69 = cos(pkin(5));
t100 = t55 * t69;
t56 = -g(1) * t68 - g(2) * t65;
t63 = -g(3) + qJDD(1);
t66 = sin(pkin(5));
t72 = sin(qJ(2));
t75 = cos(qJ(2));
t103 = -t72 * t56 + (t63 * t66 + t100) * t75;
t85 = t64 * t74 + t67 * t71;
t50 = t85 * qJD(2);
t93 = t50 * qJD(4);
t39 = -qJDD(2) * t84 - t93;
t102 = pkin(3) * t77;
t38 = pkin(4) * t49 - pkin(8) * t50;
t76 = qJD(4) ^ 2;
t99 = t66 * t72;
t90 = t72 * t100 + t75 * t56 + t63 * t99;
t33 = -pkin(2) * t77 + qJDD(2) * qJ(3) + t90;
t95 = pkin(7) * qJDD(2);
t47 = -t55 * t66 + t63 * t69;
t92 = qJD(2) * qJD(3);
t97 = t67 * t47 - 0.2e1 * t64 * t92;
t22 = (t102 * t67 - t33 - t95) * t64 + t97;
t91 = t64 * t47 + (t33 + 0.2e1 * t92) * t67;
t23 = -t102 * t62 + t67 * t95 + t91;
t98 = t71 * t22 + t74 * t23;
t18 = -pkin(4) * t76 + qJDD(4) * pkin(8) - t38 * t49 + t98;
t94 = t49 * qJD(4);
t40 = qJDD(2) * t85 - t94;
t82 = qJDD(3) - t103;
t78 = (-pkin(3) * t67 - pkin(2)) * qJDD(2) + (-pkin(7) * t96 - qJ(3)) * t77 + t82;
t19 = (-t40 + t94) * pkin(8) + (-t39 + t93) * pkin(4) + t78;
t70 = sin(qJ(5));
t73 = cos(qJ(5));
t41 = qJD(4) * t73 - t50 * t70;
t25 = qJD(5) * t41 + qJDD(4) * t70 + t40 * t73;
t42 = qJD(4) * t70 + t50 * t73;
t28 = -mrSges(6,1) * t41 + mrSges(6,2) * t42;
t48 = qJD(5) + t49;
t31 = -mrSges(6,2) * t48 + mrSges(6,3) * t41;
t37 = qJDD(5) - t39;
t15 = m(6) * (-t18 * t70 + t19 * t73) - t25 * mrSges(6,3) + t37 * mrSges(6,1) - t42 * t28 + t48 * t31;
t24 = -qJD(5) * t42 + qJDD(4) * t73 - t40 * t70;
t32 = mrSges(6,1) * t48 - mrSges(6,3) * t42;
t16 = m(6) * (t18 * t73 + t19 * t70) + t24 * mrSges(6,3) - t37 * mrSges(6,2) + t41 * t28 - t48 * t32;
t45 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t49;
t46 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t50;
t81 = -m(5) * t78 + t39 * mrSges(5,1) - t40 * mrSges(5,2) - t15 * t73 - t16 * t70 - t49 * t45 - t50 * t46;
t80 = m(4) * (-qJDD(2) * pkin(2) - t77 * qJ(3) + t82) - t81;
t88 = -mrSges(4,1) * t67 + mrSges(4,2) * t64;
t10 = m(3) * t103 + (-mrSges(3,2) + t104) * t77 + (mrSges(3,1) - t88) * qJDD(2) - t80;
t101 = t10 * t75;
t35 = mrSges(5,1) * t49 + mrSges(5,2) * t50;
t11 = m(5) * t98 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t39 - qJD(4) * t46 - t15 * t70 + t16 * t73 - t35 * t49;
t87 = t22 * t74 - t23 * t71;
t79 = m(6) * (-qJDD(4) * pkin(4) - pkin(8) * t76 + t38 * t50 - t87) - t24 * mrSges(6,1) + t25 * mrSges(6,2) - t41 * t31 + t42 * t32;
t12 = m(5) * t87 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t40 + qJD(4) * t45 - t35 * t50 - t79;
t83 = mrSges(4,3) * qJDD(2) + t77 * t88;
t7 = m(4) * t97 + t71 * t11 + t74 * t12 + (-m(4) * t33 - t83) * t64;
t8 = m(4) * t91 + t74 * t11 - t71 * t12 + t67 * t83;
t4 = m(3) * t90 - mrSges(3,1) * t77 - qJDD(2) * mrSges(3,2) - t64 * t7 + t67 * t8;
t6 = m(3) * t47 + t64 * t8 + t67 * t7;
t89 = m(2) * t63 + t101 * t66 + t4 * t99 + t6 * t69;
t2 = m(2) * t56 - t10 * t72 + t4 * t75;
t1 = m(2) * t55 - t66 * t6 + (t4 * t72 + t101) * t69;
t3 = [-m(1) * g(1) - t1 * t65 + t2 * t68, t2, t4, t8, t11, t16; -m(1) * g(2) + t1 * t68 + t2 * t65, t1, t10, t7, t12, t15; -m(1) * g(3) + t89, t89, t6, qJDD(2) * t88 - t104 * t77 + t80, -t81, t79;];
f_new = t3;
