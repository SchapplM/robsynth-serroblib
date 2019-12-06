% Calculate vector of cutting forces with Newton-Euler
% S5PRPRR4
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
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRPRR4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR4_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR4_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR4_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR4_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR4_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:49:52
% EndTime: 2019-12-05 15:49:56
% DurationCPUTime: 0.92s
% Computational Cost: add. (10304->109), mult. (18681->153), div. (0->0), fcn. (12767->12), ass. (0->70)
t63 = sin(qJ(4));
t66 = cos(qJ(4));
t43 = (-pkin(4) * t66 - pkin(8) * t63) * qJD(2);
t68 = qJD(4) ^ 2;
t82 = qJD(2) * t66;
t69 = qJD(2) ^ 2;
t57 = sin(pkin(9));
t60 = cos(pkin(9));
t47 = -g(1) * t60 - g(2) * t57;
t55 = -g(3) + qJDD(1);
t64 = sin(qJ(2));
t67 = cos(qJ(2));
t58 = sin(pkin(5));
t87 = t58 * t67;
t46 = g(1) * t57 - g(2) * t60;
t61 = cos(pkin(5));
t89 = t46 * t61;
t73 = -t47 * t64 + t55 * t87 + t67 * t89;
t26 = qJDD(2) * pkin(2) + t73;
t88 = t58 * t64;
t79 = t67 * t47 + t55 * t88 + t64 * t89;
t27 = -pkin(2) * t69 + t79;
t56 = sin(pkin(10));
t59 = cos(pkin(10));
t84 = t56 * t26 + t59 * t27;
t22 = -pkin(3) * t69 + qJDD(2) * pkin(7) + t84;
t74 = -t46 * t58 + t61 * t55;
t35 = qJDD(3) + t74;
t85 = t66 * t22 + t63 * t35;
t18 = -pkin(4) * t68 + qJDD(4) * pkin(8) + t43 * t82 + t85;
t75 = t59 * t26 - t56 * t27;
t21 = -qJDD(2) * pkin(3) - t69 * pkin(7) - t75;
t81 = qJD(2) * qJD(4);
t76 = t66 * t81;
t44 = qJDD(2) * t63 + t76;
t77 = t63 * t81;
t45 = qJDD(2) * t66 - t77;
t19 = (-t44 - t76) * pkin(8) + (-t45 + t77) * pkin(4) + t21;
t62 = sin(qJ(5));
t65 = cos(qJ(5));
t83 = qJD(2) * t63;
t40 = qJD(4) * t65 - t62 * t83;
t29 = qJD(5) * t40 + qJDD(4) * t62 + t44 * t65;
t41 = qJD(4) * t62 + t65 * t83;
t30 = -mrSges(6,1) * t40 + mrSges(6,2) * t41;
t52 = qJD(5) - t82;
t33 = -mrSges(6,2) * t52 + mrSges(6,3) * t40;
t38 = qJDD(5) - t45;
t15 = m(6) * (-t18 * t62 + t19 * t65) - t29 * mrSges(6,3) + t38 * mrSges(6,1) - t41 * t30 + t52 * t33;
t28 = -qJD(5) * t41 + qJDD(4) * t65 - t44 * t62;
t34 = mrSges(6,1) * t52 - mrSges(6,3) * t41;
t16 = m(6) * (t18 * t65 + t19 * t62) + t28 * mrSges(6,3) - t38 * mrSges(6,2) + t40 * t30 - t52 * t34;
t48 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t83;
t49 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t82;
t90 = m(5) * t21 - t45 * mrSges(5,1) + t44 * mrSges(5,2) + t65 * t15 + t62 * t16 + (t48 * t63 - t49 * t66) * qJD(2);
t86 = t66 * t35;
t42 = (-mrSges(5,1) * t66 + mrSges(5,2) * t63) * qJD(2);
t12 = m(5) * t85 - qJDD(4) * mrSges(5,2) + t45 * mrSges(5,3) - qJD(4) * t48 - t62 * t15 + t65 * t16 + t42 * t82;
t70 = m(6) * (-qJDD(4) * pkin(4) - t68 * pkin(8) - t86 + (qJD(2) * t43 + t22) * t63) - t28 * mrSges(6,1) + t29 * mrSges(6,2) - t40 * t33 + t41 * t34;
t14 = m(5) * (-t22 * t63 + t86) - t44 * mrSges(5,3) + qJDD(4) * mrSges(5,1) - t42 * t83 + qJD(4) * t49 - t70;
t80 = m(4) * t35 + t63 * t12 + t66 * t14;
t10 = m(4) * t75 + qJDD(2) * mrSges(4,1) - t69 * mrSges(4,2) - t90;
t7 = m(4) * t84 - t69 * mrSges(4,1) - qJDD(2) * mrSges(4,2) + t66 * t12 - t63 * t14;
t5 = m(3) * t73 + qJDD(2) * mrSges(3,1) - t69 * mrSges(3,2) + t59 * t10 + t56 * t7;
t6 = m(3) * t79 - t69 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t56 * t10 + t59 * t7;
t9 = m(3) * t74 + t80;
t78 = m(2) * t55 + t5 * t87 + t6 * t88 + t61 * t9;
t2 = m(2) * t47 - t5 * t64 + t6 * t67;
t1 = m(2) * t46 - t58 * t9 + (t5 * t67 + t6 * t64) * t61;
t3 = [-m(1) * g(1) - t1 * t57 + t2 * t60, t2, t6, t7, t12, t16; -m(1) * g(2) + t1 * t60 + t2 * t57, t1, t5, t10, t14, t15; -m(1) * g(3) + t78, t78, t9, t80, t90, t70;];
f_new = t3;
