% Calculate vector of cutting forces with Newton-Euler
% S5RRPRR3
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
% Datum: 2022-01-20 10:34
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRPRR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR3_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR3_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR3_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR3_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:34:05
% EndTime: 2022-01-20 10:34:06
% DurationCPUTime: 0.66s
% Computational Cost: add. (10150->85), mult. (12295->114), div. (0->0), fcn. (6580->10), ass. (0->56)
t41 = qJDD(1) + qJDD(2);
t35 = qJDD(4) + t41;
t46 = sin(qJ(5));
t50 = cos(qJ(5));
t42 = qJD(1) + qJD(2);
t36 = qJD(4) + t42;
t64 = qJD(5) * t36;
t25 = t46 * t35 + t50 * t64;
t26 = t50 * t35 - t46 * t64;
t69 = t36 * t46;
t30 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t69;
t68 = t36 * t50;
t31 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t68;
t34 = t36 ^ 2;
t49 = sin(qJ(1));
t53 = cos(qJ(1));
t62 = t49 * g(1) - t53 * g(2);
t32 = qJDD(1) * pkin(1) + t62;
t54 = qJD(1) ^ 2;
t58 = -t53 * g(1) - t49 * g(2);
t33 = -t54 * pkin(1) + t58;
t48 = sin(qJ(2));
t52 = cos(qJ(2));
t60 = t52 * t32 - t48 * t33;
t22 = t41 * pkin(2) + t60;
t40 = t42 ^ 2;
t65 = t48 * t32 + t52 * t33;
t23 = -t40 * pkin(2) + t65;
t44 = sin(pkin(9));
t45 = cos(pkin(9));
t61 = t45 * t22 - t44 * t23;
t17 = t41 * pkin(3) + t61;
t66 = t44 * t22 + t45 * t23;
t18 = -t40 * pkin(3) + t66;
t47 = sin(qJ(4));
t51 = cos(qJ(4));
t57 = t51 * t17 - t47 * t18;
t71 = (t46 * t30 - t50 * t31) * t36 + m(6) * (-t35 * pkin(4) - t34 * pkin(8) - t57) - t26 * mrSges(6,1) + t25 * mrSges(6,2);
t70 = -m(2) - m(3);
t67 = t47 * t17 + t51 * t18;
t14 = -t34 * pkin(4) + t35 * pkin(8) + t67;
t24 = (-mrSges(6,1) * t50 + mrSges(6,2) * t46) * t36;
t43 = -g(3) + qJDD(3);
t11 = m(6) * (-t46 * t14 + t50 * t43) - t25 * mrSges(6,3) + qJDD(5) * mrSges(6,1) - t24 * t69 + qJD(5) * t31;
t12 = m(6) * (t50 * t14 + t46 * t43) + t26 * mrSges(6,3) - qJDD(5) * mrSges(6,2) + t24 * t68 - qJD(5) * t30;
t63 = m(5) * t43 + t50 * t11 + t46 * t12;
t59 = m(4) * t43 + t63;
t8 = m(5) * t57 + t35 * mrSges(5,1) - t34 * mrSges(5,2) - t71;
t7 = m(5) * t67 - t34 * mrSges(5,1) - t35 * mrSges(5,2) - t46 * t11 + t50 * t12;
t6 = m(4) * t66 - t40 * mrSges(4,1) - t41 * mrSges(4,2) - t47 * t8 + t51 * t7;
t5 = m(4) * t61 + t41 * mrSges(4,1) - t40 * mrSges(4,2) + t47 * t7 + t51 * t8;
t4 = m(3) * t65 - t40 * mrSges(3,1) - t41 * mrSges(3,2) - t44 * t5 + t45 * t6;
t3 = m(3) * t60 + t41 * mrSges(3,1) - t40 * mrSges(3,2) + t44 * t6 + t45 * t5;
t2 = m(2) * t58 - t54 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t48 * t3 + t52 * t4;
t1 = m(2) * t62 + qJDD(1) * mrSges(2,1) - t54 * mrSges(2,2) + t52 * t3 + t48 * t4;
t9 = [-m(1) * g(1) - t49 * t1 + t53 * t2, t2, t4, t6, t7, t12; -m(1) * g(2) + t53 * t1 + t49 * t2, t1, t3, t5, t8, t11; (-m(1) + t70) * g(3) + t59, t70 * g(3) + t59, -m(3) * g(3) + t59, t59, t63, t71;];
f_new = t9;
