% Calculate vector of cutting forces with Newton-Euler
% S5PRRPR1
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRRPR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR1_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR1_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR1_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR1_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR1_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR1_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR1_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:15:39
% EndTime: 2019-12-05 16:15:40
% DurationCPUTime: 0.63s
% Computational Cost: add. (8446->93), mult. (12127->127), div. (0->0), fcn. (7969->10), ass. (0->62)
t51 = qJD(2) + qJD(3);
t47 = t51 ^ 2;
t55 = cos(pkin(9));
t50 = t55 ^ 2;
t53 = sin(pkin(9));
t78 = t53 ^ 2 + t50;
t84 = t78 * mrSges(5,3);
t83 = pkin(4) * t55;
t48 = qJDD(2) + qJDD(3);
t82 = pkin(7) * t48;
t54 = sin(pkin(8));
t56 = cos(pkin(8));
t38 = t54 * g(1) - t56 * g(2);
t39 = -t56 * g(1) - t54 * g(2);
t59 = sin(qJ(2));
t62 = cos(qJ(2));
t72 = t62 * t38 - t59 * t39;
t28 = qJDD(2) * pkin(2) + t72;
t63 = qJD(2) ^ 2;
t80 = t59 * t38 + t62 * t39;
t29 = -t63 * pkin(2) + t80;
t58 = sin(qJ(3));
t61 = cos(qJ(3));
t81 = t58 * t28 + t61 * t29;
t52 = -g(3) + qJDD(1);
t77 = qJD(4) * t51;
t79 = t55 * t52 - 0.2e1 * t53 * t77;
t19 = -t47 * pkin(3) + t48 * qJ(4) + t81;
t13 = (t47 * t83 - t19 - t82) * t53 + t79;
t75 = t53 * t52 + (t19 + 0.2e1 * t77) * t55;
t14 = -t50 * t47 * pkin(4) + t55 * t82 + t75;
t57 = sin(qJ(5));
t60 = cos(qJ(5));
t66 = -t53 * t57 + t55 * t60;
t32 = t66 * t51;
t67 = t53 * t60 + t55 * t57;
t33 = t67 * t51;
t22 = -t32 * mrSges(6,1) + t33 * mrSges(6,2);
t24 = t32 * qJD(5) + t67 * t48;
t30 = -qJD(5) * mrSges(6,2) + t32 * mrSges(6,3);
t11 = m(6) * (t60 * t13 - t57 * t14) - t24 * mrSges(6,3) + qJDD(5) * mrSges(6,1) - t33 * t22 + qJD(5) * t30;
t23 = -t33 * qJD(5) + t66 * t48;
t31 = qJD(5) * mrSges(6,1) - t33 * mrSges(6,3);
t12 = m(6) * (t57 * t13 + t60 * t14) + t23 * mrSges(6,3) - qJDD(5) * mrSges(6,2) + t32 * t22 - qJD(5) * t31;
t69 = -t55 * mrSges(5,1) + t53 * mrSges(5,2);
t68 = t48 * mrSges(5,3) + t47 * t69;
t8 = m(5) * t79 + t57 * t12 + t60 * t11 + (-m(5) * t19 - t68) * t53;
t9 = m(5) * t75 - t57 * t11 + t60 * t12 + t68 * t55;
t76 = m(4) * t52 + t53 * t9 + t55 * t8;
t74 = m(3) * t52 + t76;
t73 = t61 * t28 - t58 * t29;
t71 = m(2) * t52 + t74;
t70 = qJDD(4) - t73;
t65 = t23 * mrSges(6,1) + t32 * t30 - m(6) * ((-pkin(3) - t83) * t48 + (-t78 * pkin(7) - qJ(4)) * t47 + t70) - t33 * t31 - t24 * mrSges(6,2);
t64 = m(5) * (-t48 * pkin(3) - t47 * qJ(4) + t70) - t65;
t10 = m(4) * t73 + (mrSges(4,1) - t69) * t48 + (-mrSges(4,2) + t84) * t47 - t64;
t5 = m(4) * t81 - t47 * mrSges(4,1) - t48 * mrSges(4,2) - t53 * t8 + t55 * t9;
t4 = m(3) * t80 - t63 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t58 * t10 + t61 * t5;
t3 = m(3) * t72 + qJDD(2) * mrSges(3,1) - t63 * mrSges(3,2) + t61 * t10 + t58 * t5;
t2 = m(2) * t39 - t59 * t3 + t62 * t4;
t1 = m(2) * t38 + t62 * t3 + t59 * t4;
t6 = [-m(1) * g(1) - t54 * t1 + t56 * t2, t2, t4, t5, t9, t12; -m(1) * g(2) + t56 * t1 + t54 * t2, t1, t3, t10, t8, t11; -m(1) * g(3) + t71, t71, t74, t76, -t47 * t84 + t69 * t48 + t64, -t65;];
f_new = t6;
