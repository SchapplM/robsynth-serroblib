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
% m [6x1]
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
% Datum: 2022-01-20 10:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 10:48:10
% EndTime: 2022-01-20 10:48:13
% DurationCPUTime: 0.85s
% Computational Cost: add. (11670->113), mult. (15052->152), div. (0->0), fcn. (8503->10), ass. (0->66)
t51 = qJDD(1) + qJDD(2);
t59 = sin(qJ(4));
t63 = cos(qJ(4));
t53 = qJD(1) + qJD(2);
t78 = qJD(4) * t53;
t76 = t63 * t78;
t35 = t59 * t51 + t76;
t36 = t63 * t51 - t59 * t78;
t83 = t53 * t59;
t42 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t83;
t82 = t53 * t63;
t43 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t82;
t49 = t53 ^ 2;
t58 = sin(qJ(5));
t62 = cos(qJ(5));
t33 = (t58 * t63 + t59 * t62) * t53;
t21 = -t33 * qJD(5) - t58 * t35 + t62 * t36;
t32 = (-t58 * t59 + t62 * t63) * t53;
t22 = t32 * qJD(5) + t62 * t35 + t58 * t36;
t52 = qJD(4) + qJD(5);
t30 = -t52 * mrSges(6,2) + t32 * mrSges(6,3);
t31 = t52 * mrSges(6,1) - t33 * mrSges(6,3);
t44 = qJD(4) * pkin(4) - pkin(8) * t83;
t54 = t63 ^ 2;
t61 = sin(qJ(1));
t65 = cos(qJ(1));
t75 = t61 * g(1) - t65 * g(2);
t40 = qJDD(1) * pkin(1) + t75;
t66 = qJD(1) ^ 2;
t71 = -t65 * g(1) - t61 * g(2);
t41 = -t66 * pkin(1) + t71;
t60 = sin(qJ(2));
t64 = cos(qJ(2));
t72 = t64 * t40 - t60 * t41;
t28 = t51 * pkin(2) + t72;
t79 = t60 * t40 + t64 * t41;
t29 = -t49 * pkin(2) + t79;
t56 = sin(pkin(9));
t57 = cos(pkin(9));
t73 = t57 * t28 - t56 * t29;
t69 = -t51 * pkin(3) - t73;
t68 = t21 * mrSges(6,1) + t32 * t30 - m(6) * (t44 * t83 - t36 * pkin(4) + (-pkin(8) * t54 - pkin(7)) * t49 + t69) - t22 * mrSges(6,2) - t33 * t31;
t85 = (t59 * t42 - t63 * t43) * t53 + m(5) * (-t49 * pkin(7) + t69) - t36 * mrSges(5,1) + t35 * mrSges(5,2) - t68;
t84 = -m(2) - m(3);
t80 = t56 * t28 + t57 * t29;
t19 = -t49 * pkin(3) + t51 * pkin(7) + t80;
t55 = -g(3) + qJDD(3);
t81 = t63 * t19 + t59 * t55;
t74 = -t59 * t19 + t63 * t55;
t13 = (-t35 + t76) * pkin(8) + (t49 * t59 * t63 + qJDD(4)) * pkin(4) + t74;
t14 = -t54 * t49 * pkin(4) + t36 * pkin(8) - qJD(4) * t44 + t81;
t24 = -t32 * mrSges(6,1) + t33 * mrSges(6,2);
t50 = qJDD(4) + qJDD(5);
t11 = m(6) * (t62 * t13 - t58 * t14) - t22 * mrSges(6,3) + t50 * mrSges(6,1) - t33 * t24 + t52 * t30;
t12 = m(6) * (t58 * t13 + t62 * t14) + t21 * mrSges(6,3) - t50 * mrSges(6,2) + t32 * t24 - t52 * t31;
t34 = (-mrSges(5,1) * t63 + mrSges(5,2) * t59) * t53;
t8 = m(5) * t74 + qJDD(4) * mrSges(5,1) - t35 * mrSges(5,3) + qJD(4) * t43 + t62 * t11 + t58 * t12 - t34 * t83;
t9 = m(5) * t81 - qJDD(4) * mrSges(5,2) + t36 * mrSges(5,3) - qJD(4) * t42 - t58 * t11 + t62 * t12 + t34 * t82;
t77 = m(4) * t55 + t59 * t9 + t63 * t8;
t10 = m(4) * t73 + t51 * mrSges(4,1) - t49 * mrSges(4,2) - t85;
t5 = m(4) * t80 - t49 * mrSges(4,1) - t51 * mrSges(4,2) - t59 * t8 + t63 * t9;
t4 = m(3) * t79 - t49 * mrSges(3,1) - t51 * mrSges(3,2) - t56 * t10 + t57 * t5;
t3 = m(3) * t72 + t51 * mrSges(3,1) - t49 * mrSges(3,2) + t57 * t10 + t56 * t5;
t2 = m(2) * t71 - t66 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t60 * t3 + t64 * t4;
t1 = m(2) * t75 + qJDD(1) * mrSges(2,1) - t66 * mrSges(2,2) + t64 * t3 + t60 * t4;
t6 = [-m(1) * g(1) - t61 * t1 + t65 * t2, t2, t4, t5, t9, t12; -m(1) * g(2) + t65 * t1 + t61 * t2, t1, t3, t10, t8, t11; (-m(1) + t84) * g(3) + t77, t84 * g(3) + t77, -m(3) * g(3) + t77, t77, t85, -t68;];
f_new = t6;
