% Calculate vector of cutting forces with Newton-Euler
% S5PRPRR7
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
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
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRPRR7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR7_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR7_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR7_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR7_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR7_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:59:46
% EndTime: 2019-12-05 15:59:48
% DurationCPUTime: 0.47s
% Computational Cost: add. (3694->106), mult. (6799->139), div. (0->0), fcn. (3951->8), ass. (0->59)
t56 = sin(pkin(8));
t78 = cos(pkin(8));
t43 = -t78 * g(1) - t56 * g(2);
t55 = -g(3) + qJDD(1);
t59 = sin(qJ(2));
t62 = cos(qJ(2));
t79 = t62 * t43 + t59 * t55;
t70 = qJDD(2) * qJ(3) + (2 * qJD(3) * qJD(2)) + t79;
t84 = -pkin(2) - pkin(6);
t63 = qJD(2) ^ 2;
t83 = pkin(4) * t63;
t82 = mrSges(3,1) - mrSges(4,2);
t81 = (-mrSges(3,2) + mrSges(4,3));
t71 = -t59 * t43 + t62 * t55;
t67 = -t63 * qJ(3) + qJDD(3) - t71;
t25 = t84 * qJDD(2) + t67;
t42 = t56 * g(1) - t78 * g(2);
t58 = sin(qJ(4));
t61 = cos(qJ(4));
t80 = t58 * t25 - t61 * t42;
t77 = qJD(2) * t58;
t76 = qJD(2) * t61;
t75 = qJD(2) * qJD(4);
t22 = t61 * t25;
t40 = t61 * qJDD(2) - t58 * t75;
t12 = (qJDD(4) * pkin(4)) - t40 * pkin(7) + t22 + (-pkin(7) * t75 - t61 * t83 + t42) * t58;
t39 = -t58 * qJDD(2) - t61 * t75;
t46 = (qJD(4) * pkin(4)) - pkin(7) * t76;
t54 = t58 ^ 2;
t13 = t39 * pkin(7) - qJD(4) * t46 - t54 * t83 + t80;
t57 = sin(qJ(5));
t60 = cos(qJ(5));
t30 = (-t57 * t61 - t58 * t60) * qJD(2);
t18 = t30 * qJD(5) + t57 * t39 + t60 * t40;
t31 = (-t57 * t58 + t60 * t61) * qJD(2);
t20 = -t30 * mrSges(6,1) + t31 * mrSges(6,2);
t52 = qJD(4) + qJD(5);
t28 = -t52 * mrSges(6,2) + t30 * mrSges(6,3);
t51 = qJDD(4) + qJDD(5);
t10 = m(6) * (t60 * t12 - t57 * t13) - t18 * mrSges(6,3) + t51 * mrSges(6,1) - t31 * t20 + t52 * t28;
t17 = -t31 * qJD(5) + t60 * t39 - t57 * t40;
t29 = t52 * mrSges(6,1) - t31 * mrSges(6,3);
t11 = m(6) * (t57 * t12 + t60 * t13) + t17 * mrSges(6,3) - t51 * mrSges(6,2) + t30 * t20 - t52 * t29;
t38 = (mrSges(5,1) * t58 + mrSges(5,2) * t61) * qJD(2);
t44 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t77;
t6 = m(5) * (t58 * t42 + t22) - t40 * mrSges(5,3) + qJDD(4) * mrSges(5,1) - t38 * t76 + qJD(4) * t44 + t57 * t11 + t60 * t10;
t45 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t76;
t7 = m(5) * t80 - qJDD(4) * mrSges(5,2) + t39 * mrSges(5,3) - qJD(4) * t45 - t57 * t10 + t60 * t11 - t38 * t77;
t68 = -m(4) * (-qJDD(2) * pkin(2) + t67) - t58 * t7 - t61 * t6;
t3 = m(3) * t71 + t82 * qJDD(2) + (t81 * t63) + t68;
t66 = -t17 * mrSges(6,1) - t30 * t28 + m(6) * (t46 * t76 - t39 * pkin(4) + (-pkin(7) * t54 + t84) * t63 + t70) + t18 * mrSges(6,2) + t31 * t29;
t65 = -t39 * mrSges(5,1) + m(5) * (t84 * t63 + t70) + t44 * t77 + t45 * t76 + t40 * mrSges(5,2) + t66;
t64 = -m(4) * (t63 * pkin(2) - t70) + t65;
t9 = m(3) * t79 + t81 * qJDD(2) - t82 * t63 + t64;
t73 = m(2) * t55 + t62 * t3 + t59 * t9;
t69 = m(4) * t42 + t58 * t6 - t61 * t7;
t4 = (m(2) + m(3)) * t42 + t69;
t1 = m(2) * t43 - t59 * t3 + t62 * t9;
t2 = [-m(1) * g(1) + t78 * t1 - t56 * t4, t1, t9, -t69, t7, t11; -m(1) * g(2) + t56 * t1 + t78 * t4, t4, t3, -(t63 * mrSges(4,2)) - qJDD(2) * mrSges(4,3) - t64, t6, t10; -m(1) * g(3) + t73, t73, -m(3) * t42 - t69, qJDD(2) * mrSges(4,2) - t63 * mrSges(4,3) - t68, t65, t66;];
f_new = t2;
