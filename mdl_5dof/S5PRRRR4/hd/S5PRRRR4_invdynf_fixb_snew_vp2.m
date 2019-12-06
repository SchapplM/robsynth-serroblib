% Calculate vector of cutting forces with Newton-Euler
% S5PRRRR4
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRRRR4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR4_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR4_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR4_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR4_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR4_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR4_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR4_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:07:46
% EndTime: 2019-12-05 17:07:47
% DurationCPUTime: 0.72s
% Computational Cost: add. (10005->105), mult. (13390->146), div. (0->0), fcn. (8503->10), ass. (0->65)
t52 = qJDD(2) + qJDD(3);
t60 = sin(qJ(4));
t64 = cos(qJ(4));
t54 = qJD(2) + qJD(3);
t79 = qJD(4) * t54;
t77 = t64 * t79;
t35 = t60 * t52 + t77;
t36 = t64 * t52 - t60 * t79;
t84 = t54 * t60;
t40 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t84;
t83 = t54 * t64;
t41 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t83;
t50 = t54 ^ 2;
t59 = sin(qJ(5));
t63 = cos(qJ(5));
t33 = (t59 * t64 + t60 * t63) * t54;
t21 = -t33 * qJD(5) - t59 * t35 + t63 * t36;
t32 = (-t59 * t60 + t63 * t64) * t54;
t22 = t32 * qJD(5) + t63 * t35 + t59 * t36;
t53 = qJD(4) + qJD(5);
t30 = -t53 * mrSges(6,2) + t32 * mrSges(6,3);
t31 = t53 * mrSges(6,1) - t33 * mrSges(6,3);
t42 = qJD(4) * pkin(4) - pkin(8) * t84;
t55 = t64 ^ 2;
t57 = sin(pkin(9));
t58 = cos(pkin(9));
t43 = t57 * g(1) - t58 * g(2);
t44 = -t58 * g(1) - t57 * g(2);
t62 = sin(qJ(2));
t66 = cos(qJ(2));
t73 = t66 * t43 - t62 * t44;
t28 = qJDD(2) * pkin(2) + t73;
t67 = qJD(2) ^ 2;
t80 = t62 * t43 + t66 * t44;
t29 = -t67 * pkin(2) + t80;
t61 = sin(qJ(3));
t65 = cos(qJ(3));
t74 = t65 * t28 - t61 * t29;
t70 = -t52 * pkin(3) - t74;
t69 = t21 * mrSges(6,1) + t32 * t30 - m(6) * (t42 * t84 - t36 * pkin(4) + (-pkin(8) * t55 - pkin(7)) * t50 + t70) - t22 * mrSges(6,2) - t33 * t31;
t85 = (t60 * t40 - t64 * t41) * t54 + m(5) * (-t50 * pkin(7) + t70) - t36 * mrSges(5,1) + t35 * mrSges(5,2) - t69;
t81 = t61 * t28 + t65 * t29;
t19 = -t50 * pkin(3) + t52 * pkin(7) + t81;
t56 = -g(3) + qJDD(1);
t82 = t64 * t19 + t60 * t56;
t75 = -t60 * t19 + t64 * t56;
t13 = (-t35 + t77) * pkin(8) + (t50 * t60 * t64 + qJDD(4)) * pkin(4) + t75;
t14 = -t55 * t50 * pkin(4) + t36 * pkin(8) - qJD(4) * t42 + t82;
t24 = -t32 * mrSges(6,1) + t33 * mrSges(6,2);
t51 = qJDD(4) + qJDD(5);
t11 = m(6) * (t63 * t13 - t59 * t14) - t22 * mrSges(6,3) + t51 * mrSges(6,1) - t33 * t24 + t53 * t30;
t12 = m(6) * (t59 * t13 + t63 * t14) + t21 * mrSges(6,3) - t51 * mrSges(6,2) + t32 * t24 - t53 * t31;
t34 = (-mrSges(5,1) * t64 + mrSges(5,2) * t60) * t54;
t8 = m(5) * t75 + qJDD(4) * mrSges(5,1) - t35 * mrSges(5,3) + qJD(4) * t41 + t63 * t11 + t59 * t12 - t34 * t84;
t9 = m(5) * t82 - qJDD(4) * mrSges(5,2) + t36 * mrSges(5,3) - qJD(4) * t40 - t59 * t11 + t63 * t12 + t34 * t83;
t78 = m(4) * t56 + t60 * t9 + t64 * t8;
t76 = m(3) * t56 + t78;
t72 = m(2) * t56 + t76;
t10 = m(4) * t74 + t52 * mrSges(4,1) - t50 * mrSges(4,2) - t85;
t5 = m(4) * t81 - t50 * mrSges(4,1) - t52 * mrSges(4,2) - t60 * t8 + t64 * t9;
t4 = m(3) * t80 - t67 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t61 * t10 + t65 * t5;
t3 = m(3) * t73 + qJDD(2) * mrSges(3,1) - t67 * mrSges(3,2) + t65 * t10 + t61 * t5;
t2 = m(2) * t44 - t62 * t3 + t66 * t4;
t1 = m(2) * t43 + t66 * t3 + t62 * t4;
t6 = [-m(1) * g(1) - t57 * t1 + t58 * t2, t2, t4, t5, t9, t12; -m(1) * g(2) + t58 * t1 + t57 * t2, t1, t3, t10, t8, t11; -m(1) * g(3) + t72, t72, t76, t78, t85, -t69;];
f_new = t6;
