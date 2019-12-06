% Calculate vector of cutting forces with Newton-Euler
% S5RRPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-05 18:19
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRPPR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR1_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR1_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR1_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR1_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:18:14
% EndTime: 2019-12-05 18:18:16
% DurationCPUTime: 0.70s
% Computational Cost: add. (10111->101), mult. (13789->133), div. (0->0), fcn. (7969->10), ass. (0->63)
t52 = qJD(1) + qJD(2);
t48 = t52 ^ 2;
t56 = cos(pkin(9));
t51 = t56 ^ 2;
t54 = sin(pkin(9));
t78 = t54 ^ 2 + t51;
t86 = t78 * mrSges(5,3);
t85 = -m(2) - m(3);
t84 = pkin(4) * t56;
t49 = qJDD(1) + qJDD(2);
t83 = pkin(7) * t49;
t60 = sin(qJ(1));
t63 = cos(qJ(1));
t79 = t63 * g(2) + t60 * g(3);
t38 = qJDD(1) * pkin(1) + t79;
t64 = qJD(1) ^ 2;
t74 = t60 * g(2) - t63 * g(3);
t39 = -t64 * pkin(1) + t74;
t59 = sin(qJ(2));
t62 = cos(qJ(2));
t72 = t62 * t38 - t59 * t39;
t28 = t49 * pkin(2) + t72;
t81 = t59 * t38 + t62 * t39;
t29 = -t48 * pkin(2) + t81;
t55 = sin(pkin(8));
t57 = cos(pkin(8));
t82 = t55 * t28 + t57 * t29;
t53 = -g(1) + qJDD(3);
t77 = qJD(4) * t52;
t80 = t56 * t53 - 0.2e1 * t54 * t77;
t19 = -t48 * pkin(3) + t49 * qJ(4) + t82;
t13 = (t48 * t84 - t19 - t83) * t54 + t80;
t75 = t54 * t53 + (t19 + 0.2e1 * t77) * t56;
t14 = -t51 * t48 * pkin(4) + t56 * t83 + t75;
t58 = sin(qJ(5));
t61 = cos(qJ(5));
t67 = -t54 * t58 + t56 * t61;
t32 = t67 * t52;
t68 = t54 * t61 + t56 * t58;
t33 = t68 * t52;
t22 = -t32 * mrSges(6,1) + t33 * mrSges(6,2);
t24 = t32 * qJD(5) + t68 * t49;
t30 = -qJD(5) * mrSges(6,2) + t32 * mrSges(6,3);
t11 = m(6) * (t61 * t13 - t58 * t14) - t24 * mrSges(6,3) + qJDD(5) * mrSges(6,1) - t33 * t22 + qJD(5) * t30;
t23 = -t33 * qJD(5) + t67 * t49;
t31 = qJD(5) * mrSges(6,1) - t33 * mrSges(6,3);
t12 = m(6) * (t58 * t13 + t61 * t14) + t23 * mrSges(6,3) - qJDD(5) * mrSges(6,2) + t32 * t22 - qJD(5) * t31;
t70 = -t56 * mrSges(5,1) + t54 * mrSges(5,2);
t69 = t49 * mrSges(5,3) + t48 * t70;
t8 = m(5) * t80 + t58 * t12 + t61 * t11 + (-m(5) * t19 - t69) * t54;
t9 = m(5) * t75 - t58 * t11 + t61 * t12 + t69 * t56;
t76 = m(4) * t53 + t54 * t9 + t56 * t8;
t73 = t57 * t28 - t55 * t29;
t71 = qJDD(4) - t73;
t66 = t23 * mrSges(6,1) + t32 * t30 - m(6) * ((-pkin(3) - t84) * t49 + (-t78 * pkin(7) - qJ(4)) * t48 + t71) - t33 * t31 - t24 * mrSges(6,2);
t65 = m(5) * (-t49 * pkin(3) - t48 * qJ(4) + t71) - t66;
t10 = m(4) * t73 + (mrSges(4,1) - t70) * t49 + (-mrSges(4,2) + t86) * t48 - t65;
t5 = m(4) * t82 - t48 * mrSges(4,1) - t49 * mrSges(4,2) - t54 * t8 + t56 * t9;
t4 = m(3) * t81 - t48 * mrSges(3,1) - t49 * mrSges(3,2) - t55 * t10 + t57 * t5;
t3 = m(3) * t72 + t49 * mrSges(3,1) - t48 * mrSges(3,2) + t57 * t10 + t55 * t5;
t2 = m(2) * t79 + qJDD(1) * mrSges(2,1) - t64 * mrSges(2,2) + t62 * t3 + t59 * t4;
t1 = m(2) * t74 - t64 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t59 * t3 + t62 * t4;
t6 = [(-m(1) + t85) * g(1) + t76, t1, t4, t5, t9, t12; -m(1) * g(2) - t60 * t1 - t63 * t2, t2, t3, t10, t8, t11; -m(1) * g(3) + t63 * t1 - t60 * t2, t85 * g(1) + t76, -m(3) * g(1) + t76, t76, -t48 * t86 + t70 * t49 + t65, -t66;];
f_new = t6;
