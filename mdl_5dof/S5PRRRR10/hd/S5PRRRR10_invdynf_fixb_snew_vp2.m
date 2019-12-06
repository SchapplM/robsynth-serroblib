% Calculate vector of cutting forces with Newton-Euler
% S5PRRRR10
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRRRR10_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR10_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR10_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR10_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR10_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR10_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR10_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:23:37
% EndTime: 2019-12-05 17:23:44
% DurationCPUTime: 2.83s
% Computational Cost: add. (35778->142), mult. (73570->206), div. (0->0), fcn. (57600->14), ass. (0->88)
t72 = sin(pkin(6));
t106 = pkin(8) * t72;
t85 = qJD(2) ^ 2;
t73 = sin(pkin(5));
t84 = cos(qJ(2));
t101 = t73 * t84;
t71 = sin(pkin(11));
t74 = cos(pkin(11));
t62 = t71 * g(1) - t74 * g(2);
t76 = cos(pkin(5));
t103 = t62 * t76;
t63 = -t74 * g(1) - t71 * g(2);
t70 = -g(3) + qJDD(1);
t80 = sin(qJ(2));
t92 = t70 * t101 + t84 * t103 - t80 * t63;
t35 = qJDD(2) * pkin(2) + t85 * t106 + t92;
t49 = -t73 * t62 + t76 * t70;
t75 = cos(pkin(6));
t108 = t35 * t75 + t49 * t72;
t79 = sin(qJ(3));
t83 = cos(qJ(3));
t98 = qJD(2) * qJD(3);
t56 = (-qJDD(2) * t83 + t79 * t98) * t72;
t102 = t73 * t80;
t96 = t70 * t102 + t80 * t103 + t84 * t63;
t36 = -t85 * pkin(2) + qJDD(2) * t106 + t96;
t107 = t108 * t83 - t79 * t36;
t99 = qJD(2) * t72;
t54 = (-pkin(3) * t83 - pkin(9) * t79) * t99;
t68 = t75 * qJD(2) + qJD(3);
t66 = t68 ^ 2;
t67 = t75 * qJDD(2) + qJDD(3);
t94 = t83 * t99;
t97 = t108 * t79 + t83 * t36;
t21 = -t66 * pkin(3) + t67 * pkin(9) + t54 * t94 + t97;
t45 = t75 * t49;
t55 = (qJDD(2) * t79 + t83 * t98) * t72;
t23 = t56 * pkin(3) - t55 * pkin(9) + t45 + (-t35 + (pkin(3) * t79 - pkin(9) * t83) * t68 * qJD(2)) * t72;
t78 = sin(qJ(4));
t82 = cos(qJ(4));
t100 = t82 * t21 + t78 * t23;
t95 = t79 * t99;
t47 = t82 * t68 - t78 * t95;
t48 = t78 * t68 + t82 * t95;
t38 = -t47 * pkin(4) - t48 * pkin(10);
t50 = qJDD(4) + t56;
t61 = qJD(4) - t94;
t60 = t61 ^ 2;
t17 = -t60 * pkin(4) + t50 * pkin(10) + t47 * t38 + t100;
t20 = -t67 * pkin(3) - t66 * pkin(9) + t54 * t95 - t107;
t30 = -t48 * qJD(4) - t78 * t55 + t82 * t67;
t31 = t47 * qJD(4) + t82 * t55 + t78 * t67;
t18 = (-t47 * t61 - t31) * pkin(10) + (t48 * t61 - t30) * pkin(4) + t20;
t77 = sin(qJ(5));
t81 = cos(qJ(5));
t39 = -t77 * t48 + t81 * t61;
t25 = t39 * qJD(5) + t81 * t31 + t77 * t50;
t40 = t81 * t48 + t77 * t61;
t26 = -t39 * mrSges(6,1) + t40 * mrSges(6,2);
t46 = qJD(5) - t47;
t27 = -t46 * mrSges(6,2) + t39 * mrSges(6,3);
t29 = qJDD(5) - t30;
t14 = m(6) * (-t77 * t17 + t81 * t18) - t25 * mrSges(6,3) + t29 * mrSges(6,1) - t40 * t26 + t46 * t27;
t24 = -t40 * qJD(5) - t77 * t31 + t81 * t50;
t28 = t46 * mrSges(6,1) - t40 * mrSges(6,3);
t15 = m(6) * (t81 * t17 + t77 * t18) + t24 * mrSges(6,3) - t29 * mrSges(6,2) + t39 * t26 - t46 * t28;
t37 = -t47 * mrSges(5,1) + t48 * mrSges(5,2);
t42 = t61 * mrSges(5,1) - t48 * mrSges(5,3);
t12 = m(5) * t100 - t50 * mrSges(5,2) + t30 * mrSges(5,3) - t77 * t14 + t81 * t15 + t47 * t37 - t61 * t42;
t41 = -t61 * mrSges(5,2) + t47 * mrSges(5,3);
t90 = -t78 * t21 + t82 * t23;
t87 = m(6) * (-t50 * pkin(4) - t60 * pkin(10) + t48 * t38 - t90) - t24 * mrSges(6,1) + t25 * mrSges(6,2) - t39 * t27 + t40 * t28;
t13 = m(5) * t90 + t50 * mrSges(5,1) - t31 * mrSges(5,3) - t48 * t37 + t61 * t41 - t87;
t51 = t68 * mrSges(4,1) - mrSges(4,3) * t95;
t52 = -t68 * mrSges(4,2) + mrSges(4,3) * t94;
t10 = m(4) * (-t72 * t35 + t45) + t55 * mrSges(4,2) + t56 * mrSges(4,1) + t78 * t12 + t82 * t13 + (t51 * t79 - t52 * t83) * t99;
t53 = (-mrSges(4,1) * t83 + mrSges(4,2) * t79) * t99;
t86 = m(5) * t20 - t30 * mrSges(5,1) + t31 * mrSges(5,2) + t81 * t14 + t77 * t15 - t47 * t41 + t48 * t42;
t11 = m(4) * t107 + t67 * mrSges(4,1) - t55 * mrSges(4,3) + t68 * t52 - t53 * t95 - t86;
t9 = m(4) * t97 - t67 * mrSges(4,2) - t56 * mrSges(4,3) + t82 * t12 - t78 * t13 - t68 * t51 + t53 * t94;
t91 = t83 * t11 + t79 * t9;
t4 = m(3) * t92 + qJDD(2) * mrSges(3,1) - t85 * mrSges(3,2) - t72 * t10 + t91 * t75;
t6 = m(3) * t49 + t75 * t10 + t91 * t72;
t8 = m(3) * t96 - t85 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t79 * t11 + t83 * t9;
t93 = m(2) * t70 + t4 * t101 + t8 * t102 + t76 * t6;
t2 = m(2) * t63 - t80 * t4 + t84 * t8;
t1 = m(2) * t62 - t73 * t6 + (t4 * t84 + t8 * t80) * t76;
t3 = [-m(1) * g(1) - t71 * t1 + t74 * t2, t2, t8, t9, t12, t15; -m(1) * g(2) + t74 * t1 + t71 * t2, t1, t4, t11, t13, t14; -m(1) * g(3) + t93, t93, t6, t10, t86, t87;];
f_new = t3;
