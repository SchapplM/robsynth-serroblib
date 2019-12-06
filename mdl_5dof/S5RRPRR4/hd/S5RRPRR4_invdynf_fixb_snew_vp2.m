% Calculate vector of cutting forces with Newton-Euler
% S5RRPRR4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRPRR4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR4_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR4_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR4_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR4_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:31:54
% EndTime: 2019-12-05 18:31:56
% DurationCPUTime: 0.77s
% Computational Cost: add. (11670->113), mult. (15052->152), div. (0->0), fcn. (8503->10), ass. (0->66)
t53 = qJDD(1) + qJDD(2);
t61 = sin(qJ(4));
t65 = cos(qJ(4));
t55 = qJD(1) + qJD(2);
t79 = qJD(4) * t55;
t77 = t65 * t79;
t35 = t61 * t53 + t77;
t36 = t65 * t53 - t61 * t79;
t85 = t55 * t61;
t42 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t85;
t84 = t55 * t65;
t43 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t84;
t51 = t55 ^ 2;
t60 = sin(qJ(5));
t64 = cos(qJ(5));
t33 = (t60 * t65 + t61 * t64) * t55;
t21 = -t33 * qJD(5) - t60 * t35 + t64 * t36;
t32 = (-t60 * t61 + t64 * t65) * t55;
t22 = t32 * qJD(5) + t64 * t35 + t60 * t36;
t54 = qJD(4) + qJD(5);
t30 = -t54 * mrSges(6,2) + t32 * mrSges(6,3);
t31 = t54 * mrSges(6,1) - t33 * mrSges(6,3);
t44 = qJD(4) * pkin(4) - pkin(8) * t85;
t56 = t65 ^ 2;
t63 = sin(qJ(1));
t67 = cos(qJ(1));
t80 = t67 * g(2) + t63 * g(3);
t40 = qJDD(1) * pkin(1) + t80;
t68 = qJD(1) ^ 2;
t76 = t63 * g(2) - t67 * g(3);
t41 = -t68 * pkin(1) + t76;
t62 = sin(qJ(2));
t66 = cos(qJ(2));
t73 = t66 * t40 - t62 * t41;
t28 = t53 * pkin(2) + t73;
t81 = t62 * t40 + t66 * t41;
t29 = -t51 * pkin(2) + t81;
t58 = sin(pkin(9));
t59 = cos(pkin(9));
t74 = t59 * t28 - t58 * t29;
t71 = -t53 * pkin(3) - t74;
t70 = t21 * mrSges(6,1) + t32 * t30 - m(6) * (t44 * t85 - t36 * pkin(4) + (-pkin(8) * t56 - pkin(7)) * t51 + t71) - t22 * mrSges(6,2) - t33 * t31;
t87 = (t61 * t42 - t65 * t43) * t55 + m(5) * (-t51 * pkin(7) + t71) - t36 * mrSges(5,1) + t35 * mrSges(5,2) - t70;
t86 = -m(2) - m(3);
t82 = t58 * t28 + t59 * t29;
t19 = -t51 * pkin(3) + t53 * pkin(7) + t82;
t57 = -g(1) + qJDD(3);
t83 = t65 * t19 + t61 * t57;
t75 = -t61 * t19 + t65 * t57;
t13 = (-t35 + t77) * pkin(8) + (t51 * t61 * t65 + qJDD(4)) * pkin(4) + t75;
t14 = -t56 * t51 * pkin(4) + t36 * pkin(8) - qJD(4) * t44 + t83;
t24 = -t32 * mrSges(6,1) + t33 * mrSges(6,2);
t52 = qJDD(4) + qJDD(5);
t11 = m(6) * (t64 * t13 - t60 * t14) - t22 * mrSges(6,3) + t52 * mrSges(6,1) - t33 * t24 + t54 * t30;
t12 = m(6) * (t60 * t13 + t64 * t14) + t21 * mrSges(6,3) - t52 * mrSges(6,2) + t32 * t24 - t54 * t31;
t34 = (-mrSges(5,1) * t65 + mrSges(5,2) * t61) * t55;
t8 = m(5) * t75 + qJDD(4) * mrSges(5,1) - t35 * mrSges(5,3) + qJD(4) * t43 + t64 * t11 + t60 * t12 - t34 * t85;
t9 = m(5) * t83 - qJDD(4) * mrSges(5,2) + t36 * mrSges(5,3) - qJD(4) * t42 - t60 * t11 + t64 * t12 + t34 * t84;
t78 = m(4) * t57 + t61 * t9 + t65 * t8;
t10 = m(4) * t74 + t53 * mrSges(4,1) - t51 * mrSges(4,2) - t87;
t5 = m(4) * t82 - t51 * mrSges(4,1) - t53 * mrSges(4,2) - t61 * t8 + t65 * t9;
t4 = m(3) * t81 - t51 * mrSges(3,1) - t53 * mrSges(3,2) - t58 * t10 + t59 * t5;
t3 = m(3) * t73 + t53 * mrSges(3,1) - t51 * mrSges(3,2) + t59 * t10 + t58 * t5;
t2 = m(2) * t80 + qJDD(1) * mrSges(2,1) - t68 * mrSges(2,2) + t66 * t3 + t62 * t4;
t1 = m(2) * t76 - t68 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t62 * t3 + t66 * t4;
t6 = [(-m(1) + t86) * g(1) + t78, t1, t4, t5, t9, t12; -m(1) * g(2) - t63 * t1 - t67 * t2, t2, t3, t10, t8, t11; -m(1) * g(3) + t67 * t1 - t63 * t2, t86 * g(1) + t78, -m(3) * g(1) + t78, t78, t87, -t70;];
f_new = t6;
