% Calculate vector of cutting forces with Newton-Euler
% S5PRRRR8
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRRRR8_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR8_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR8_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR8_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR8_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR8_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR8_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR8_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR8_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:15:09
% EndTime: 2019-12-05 17:15:15
% DurationCPUTime: 1.53s
% Computational Cost: add. (17882->137), mult. (35018->191), div. (0->0), fcn. (24473->12), ass. (0->80)
t75 = sin(qJ(4));
t76 = sin(qJ(3));
t79 = cos(qJ(4));
t80 = cos(qJ(3));
t50 = (t75 * t76 - t79 * t80) * qJD(2);
t70 = sin(pkin(10));
t72 = cos(pkin(10));
t58 = t70 * g(1) - t72 * g(2);
t73 = cos(pkin(5));
t102 = t58 * t73;
t59 = -t72 * g(1) - t70 * g(2);
t69 = -g(3) + qJDD(1);
t71 = sin(pkin(5));
t77 = sin(qJ(2));
t81 = cos(qJ(2));
t105 = -t77 * t59 + (t69 * t71 + t102) * t81;
t96 = qJD(2) * qJD(3);
t93 = t80 * t96;
t56 = t76 * qJDD(2) + t93;
t57 = t80 * qJDD(2) - t76 * t96;
t98 = qJD(2) * t76;
t60 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t98;
t97 = qJD(2) * t80;
t61 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t97;
t82 = qJD(2) ^ 2;
t86 = -qJDD(2) * pkin(2) - t105;
t101 = t71 * t77;
t95 = t69 * t101 + t77 * t102 + t81 * t59;
t37 = -t82 * pkin(2) + qJDD(2) * pkin(7) + t95;
t47 = -t71 * t58 + t73 * t69;
t92 = -t76 * t37 + t80 * t47;
t22 = (-t56 + t93) * pkin(8) + (t76 * t80 * t82 + qJDD(3)) * pkin(3) + t92;
t63 = qJD(3) * pkin(3) - pkin(8) * t98;
t68 = t80 ^ 2;
t99 = t80 * t37 + t76 * t47;
t23 = -t68 * t82 * pkin(3) + t57 * pkin(8) - qJD(3) * t63 + t99;
t100 = t75 * t22 + t79 * t23;
t51 = (t75 * t80 + t76 * t79) * qJD(2);
t40 = t50 * pkin(4) - t51 * pkin(9);
t67 = qJD(3) + qJD(4);
t65 = t67 ^ 2;
t66 = qJDD(3) + qJDD(4);
t18 = -t65 * pkin(4) + t66 * pkin(9) - t50 * t40 + t100;
t32 = -t51 * qJD(4) - t75 * t56 + t79 * t57;
t33 = -t50 * qJD(4) + t79 * t56 + t75 * t57;
t83 = -t57 * pkin(3) + t63 * t98 + (-pkin(8) * t68 - pkin(7)) * t82 + t86;
t19 = (t50 * t67 - t33) * pkin(9) + (t51 * t67 - t32) * pkin(4) + t83;
t74 = sin(qJ(5));
t78 = cos(qJ(5));
t41 = -t74 * t51 + t78 * t67;
t25 = t41 * qJD(5) + t78 * t33 + t74 * t66;
t42 = t78 * t51 + t74 * t67;
t28 = -t41 * mrSges(6,1) + t42 * mrSges(6,2);
t30 = qJDD(5) - t32;
t48 = qJD(5) + t50;
t34 = -t48 * mrSges(6,2) + t41 * mrSges(6,3);
t15 = m(6) * (-t74 * t18 + t78 * t19) - t25 * mrSges(6,3) + t30 * mrSges(6,1) - t42 * t28 + t48 * t34;
t24 = -t42 * qJD(5) - t74 * t33 + t78 * t66;
t35 = t48 * mrSges(6,1) - t42 * mrSges(6,3);
t16 = m(6) * (t78 * t18 + t74 * t19) + t24 * mrSges(6,3) - t30 * mrSges(6,2) + t41 * t28 - t48 * t35;
t45 = -t67 * mrSges(5,2) - t50 * mrSges(5,3);
t46 = t67 * mrSges(5,1) - t51 * mrSges(5,3);
t87 = -m(5) * t83 + t32 * mrSges(5,1) - t33 * mrSges(5,2) - t78 * t15 - t74 * t16 - t50 * t45 - t51 * t46;
t104 = (t76 * t60 - t80 * t61) * qJD(2) + m(4) * (-t82 * pkin(7) + t86) - t57 * mrSges(4,1) + t56 * mrSges(4,2) - t87;
t10 = m(3) * t105 + qJDD(2) * mrSges(3,1) - t82 * mrSges(3,2) - t104;
t103 = t10 * t81;
t39 = t50 * mrSges(5,1) + t51 * mrSges(5,2);
t11 = m(5) * t100 - t66 * mrSges(5,2) + t32 * mrSges(5,3) - t74 * t15 + t78 * t16 - t50 * t39 - t67 * t46;
t91 = t79 * t22 - t75 * t23;
t85 = m(6) * (-t66 * pkin(4) - t65 * pkin(9) + t51 * t40 - t91) - t24 * mrSges(6,1) + t25 * mrSges(6,2) - t41 * t34 + t42 * t35;
t12 = m(5) * t91 + t66 * mrSges(5,1) - t33 * mrSges(5,3) - t51 * t39 + t67 * t45 - t85;
t55 = (-mrSges(4,1) * t80 + mrSges(4,2) * t76) * qJD(2);
t7 = m(4) * t92 + qJDD(3) * mrSges(4,1) - t56 * mrSges(4,3) + qJD(3) * t61 + t75 * t11 + t79 * t12 - t55 * t98;
t8 = m(4) * t99 - qJDD(3) * mrSges(4,2) + t57 * mrSges(4,3) - qJD(3) * t60 + t79 * t11 - t75 * t12 + t55 * t97;
t4 = m(3) * t95 - t82 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t76 * t7 + t80 * t8;
t6 = m(3) * t47 + t80 * t7 + t76 * t8;
t94 = m(2) * t69 + t4 * t101 + t71 * t103 + t73 * t6;
t2 = m(2) * t59 - t77 * t10 + t81 * t4;
t1 = m(2) * t58 - t71 * t6 + (t4 * t77 + t103) * t73;
t3 = [-m(1) * g(1) - t70 * t1 + t72 * t2, t2, t4, t8, t11, t16; -m(1) * g(2) + t72 * t1 + t70 * t2, t1, t10, t7, t12, t15; -m(1) * g(3) + t94, t94, t6, t104, -t87, t85;];
f_new = t3;
