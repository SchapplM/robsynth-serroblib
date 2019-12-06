% Calculate vector of cutting forces with Newton-Euler
% S5PRPRR3
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRPRR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR3_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR3_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR3_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR3_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR3_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR3_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR3_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:46:42
% EndTime: 2019-12-05 15:46:44
% DurationCPUTime: 0.63s
% Computational Cost: add. (6633->105), mult. (11908->145), div. (0->0), fcn. (7335->10), ass. (0->62)
t62 = sin(qJ(4));
t65 = cos(qJ(4));
t78 = qJD(2) * qJD(4);
t75 = t65 * t78;
t41 = t62 * qJDD(2) + t75;
t42 = t65 * qJDD(2) - t62 * t78;
t80 = qJD(2) * t62;
t46 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t80;
t79 = qJD(2) * t65;
t47 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t79;
t67 = qJD(2) ^ 2;
t61 = sin(qJ(5));
t64 = cos(qJ(5));
t35 = (t61 * t65 + t62 * t64) * qJD(2);
t23 = -t35 * qJD(5) - t61 * t41 + t64 * t42;
t34 = (-t61 * t62 + t64 * t65) * qJD(2);
t24 = t34 * qJD(5) + t64 * t41 + t61 * t42;
t54 = qJD(4) + qJD(5);
t32 = -t54 * mrSges(6,2) + t34 * mrSges(6,3);
t33 = t54 * mrSges(6,1) - t35 * mrSges(6,3);
t48 = qJD(4) * pkin(4) - pkin(7) * t80;
t55 = t65 ^ 2;
t58 = sin(pkin(8));
t60 = cos(pkin(8));
t45 = -t60 * g(1) - t58 * g(2);
t56 = -g(3) + qJDD(1);
t63 = sin(qJ(2));
t66 = cos(qJ(2));
t72 = -t63 * t45 + t66 * t56;
t30 = qJDD(2) * pkin(2) + t72;
t81 = t66 * t45 + t63 * t56;
t31 = -t67 * pkin(2) + t81;
t57 = sin(pkin(9));
t59 = cos(pkin(9));
t73 = t59 * t30 - t57 * t31;
t70 = -qJDD(2) * pkin(3) - t73;
t69 = t23 * mrSges(6,1) + t34 * t32 - m(6) * (t48 * t80 - t42 * pkin(4) + (-pkin(7) * t55 - pkin(6)) * t67 + t70) - t24 * mrSges(6,2) - t35 * t33;
t84 = (t62 * t46 - t65 * t47) * qJD(2) + m(5) * (-t67 * pkin(6) + t70) - t42 * mrSges(5,1) + t41 * mrSges(5,2) - t69;
t82 = t57 * t30 + t59 * t31;
t21 = -t67 * pkin(3) + qJDD(2) * pkin(6) + t82;
t44 = t58 * g(1) - t60 * g(2);
t43 = qJDD(3) - t44;
t83 = t65 * t21 + t62 * t43;
t74 = -t62 * t21 + t65 * t43;
t15 = (-t41 + t75) * pkin(7) + (t62 * t65 * t67 + qJDD(4)) * pkin(4) + t74;
t16 = -t55 * t67 * pkin(4) + t42 * pkin(7) - qJD(4) * t48 + t83;
t26 = -t34 * mrSges(6,1) + t35 * mrSges(6,2);
t53 = qJDD(4) + qJDD(5);
t13 = m(6) * (t64 * t15 - t61 * t16) - t24 * mrSges(6,3) + t53 * mrSges(6,1) - t35 * t26 + t54 * t32;
t14 = m(6) * (t61 * t15 + t64 * t16) + t23 * mrSges(6,3) - t53 * mrSges(6,2) + t34 * t26 - t54 * t33;
t40 = (-mrSges(5,1) * t65 + mrSges(5,2) * t62) * qJD(2);
t10 = m(5) * t74 + qJDD(4) * mrSges(5,1) - t41 * mrSges(5,3) + qJD(4) * t47 + t64 * t13 + t61 * t14 - t40 * t80;
t11 = m(5) * t83 - qJDD(4) * mrSges(5,2) + t42 * mrSges(5,3) - qJD(4) * t46 - t61 * t13 + t64 * t14 + t40 * t79;
t77 = m(4) * t43 + t65 * t10 + t62 * t11;
t12 = m(4) * t73 + qJDD(2) * mrSges(4,1) - t67 * mrSges(4,2) - t84;
t6 = m(4) * t82 - t67 * mrSges(4,1) - qJDD(2) * mrSges(4,2) - t62 * t10 + t65 * t11;
t4 = m(3) * t72 + qJDD(2) * mrSges(3,1) - t67 * mrSges(3,2) + t59 * t12 + t57 * t6;
t5 = m(3) * t81 - t67 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t57 * t12 + t59 * t6;
t76 = m(2) * t56 + t66 * t4 + t63 * t5;
t7 = (m(2) + m(3)) * t44 - t77;
t1 = m(2) * t45 - t63 * t4 + t66 * t5;
t2 = [-m(1) * g(1) + t60 * t1 - t58 * t7, t1, t5, t6, t11, t14; -m(1) * g(2) + t58 * t1 + t60 * t7, t7, t4, t12, t10, t13; -m(1) * g(3) + t76, t76, -m(3) * t44 + t77, t77, t84, -t69;];
f_new = t2;
