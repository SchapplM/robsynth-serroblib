% Calculate vector of cutting forces with Newton-Euler
% S5PRPPR2
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
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
% Datum: 2019-12-05 15:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRPPR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR2_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR2_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR2_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR2_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR2_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR2_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR2_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:24:03
% EndTime: 2019-12-05 15:24:04
% DurationCPUTime: 0.57s
% Computational Cost: add. (5670->93), mult. (10727->126), div. (0->0), fcn. (6839->10), ass. (0->59)
t63 = qJD(2) ^ 2;
t56 = cos(pkin(9));
t51 = t56 ^ 2;
t53 = sin(pkin(9));
t78 = t53 ^ 2 + t51;
t83 = t78 * mrSges(5,3);
t82 = pkin(4) * t63;
t55 = sin(pkin(7));
t58 = cos(pkin(7));
t43 = -t58 * g(1) - t55 * g(2);
t52 = -g(3) + qJDD(1);
t60 = sin(qJ(2));
t62 = cos(qJ(2));
t71 = -t60 * t43 + t62 * t52;
t30 = qJDD(2) * pkin(2) + t71;
t79 = t62 * t43 + t60 * t52;
t31 = -t63 * pkin(2) + t79;
t54 = sin(pkin(8));
t57 = cos(pkin(8));
t81 = t54 * t30 + t57 * t31;
t42 = t55 * g(1) - t58 * g(2);
t41 = qJDD(3) - t42;
t76 = qJD(2) * qJD(4);
t80 = t56 * t41 - 0.2e1 * t53 * t76;
t77 = pkin(6) * qJDD(2);
t21 = -t63 * pkin(3) + qJDD(2) * qJ(4) + t81;
t15 = (t56 * t82 - t21 - t77) * t53 + t80;
t73 = t53 * t41 + (t21 + 0.2e1 * t76) * t56;
t16 = -t51 * t82 + t56 * t77 + t73;
t59 = sin(qJ(5));
t61 = cos(qJ(5));
t67 = -t53 * t59 + t56 * t61;
t34 = t67 * qJD(2);
t68 = t53 * t61 + t56 * t59;
t35 = t68 * qJD(2);
t23 = -t34 * mrSges(6,1) + t35 * mrSges(6,2);
t26 = t34 * qJD(5) + qJDD(2) * t68;
t32 = -qJD(5) * mrSges(6,2) + t34 * mrSges(6,3);
t13 = m(6) * (t61 * t15 - t59 * t16) - t26 * mrSges(6,3) + qJDD(5) * mrSges(6,1) - t35 * t23 + qJD(5) * t32;
t25 = -t35 * qJD(5) + qJDD(2) * t67;
t33 = qJD(5) * mrSges(6,1) - t35 * mrSges(6,3);
t14 = m(6) * (t59 * t15 + t61 * t16) + t25 * mrSges(6,3) - qJDD(5) * mrSges(6,2) + t34 * t23 - qJD(5) * t33;
t69 = -t56 * mrSges(5,1) + t53 * mrSges(5,2);
t66 = qJDD(2) * mrSges(5,3) + t63 * t69;
t10 = m(5) * t80 + t59 * t14 + t61 * t13 + (-m(5) * t21 - t66) * t53;
t11 = m(5) * t73 - t59 * t13 + t61 * t14 + t66 * t56;
t75 = m(4) * t41 + t56 * t10 + t53 * t11;
t72 = t57 * t30 - t54 * t31;
t70 = qJDD(4) - t72;
t65 = t25 * mrSges(6,1) + t34 * t32 - m(6) * ((-pkin(4) * t56 - pkin(3)) * qJDD(2) + (-t78 * pkin(6) - qJ(4)) * t63 + t70) - t35 * t33 - t26 * mrSges(6,2);
t64 = m(5) * (-qJDD(2) * pkin(3) - t63 * qJ(4) + t70) - t65;
t12 = m(4) * t72 + (-mrSges(4,2) + t83) * t63 + (mrSges(4,1) - t69) * qJDD(2) - t64;
t6 = m(4) * t81 - t63 * mrSges(4,1) - qJDD(2) * mrSges(4,2) - t53 * t10 + t56 * t11;
t4 = m(3) * t71 + qJDD(2) * mrSges(3,1) - t63 * mrSges(3,2) + t57 * t12 + t54 * t6;
t5 = m(3) * t79 - t63 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t54 * t12 + t57 * t6;
t74 = m(2) * t52 + t62 * t4 + t60 * t5;
t7 = (m(2) + m(3)) * t42 - t75;
t1 = m(2) * t43 - t60 * t4 + t62 * t5;
t2 = [-m(1) * g(1) + t58 * t1 - t55 * t7, t1, t5, t6, t11, t14; -m(1) * g(2) + t55 * t1 + t58 * t7, t7, t4, t12, t10, t13; -m(1) * g(3) + t74, t74, -m(3) * t42 + t75, t75, t69 * qJDD(2) - t63 * t83 + t64, -t65;];
f_new = t2;
