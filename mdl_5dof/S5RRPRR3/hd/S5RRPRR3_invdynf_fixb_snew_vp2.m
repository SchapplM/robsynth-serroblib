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
% Datum: 2019-12-05 18:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:30:23
% EndTime: 2019-12-05 18:30:23
% DurationCPUTime: 0.59s
% Computational Cost: add. (10150->85), mult. (12295->114), div. (0->0), fcn. (6580->10), ass. (0->56)
t43 = qJDD(1) + qJDD(2);
t35 = qJDD(4) + t43;
t48 = sin(qJ(5));
t52 = cos(qJ(5));
t44 = qJD(1) + qJD(2);
t36 = qJD(4) + t44;
t65 = qJD(5) * t36;
t25 = t48 * t35 + t52 * t65;
t26 = t52 * t35 - t48 * t65;
t71 = t36 * t48;
t30 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t71;
t70 = t36 * t52;
t31 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t70;
t34 = t36 ^ 2;
t51 = sin(qJ(1));
t55 = cos(qJ(1));
t66 = t55 * g(2) + t51 * g(3);
t32 = qJDD(1) * pkin(1) + t66;
t56 = qJD(1) ^ 2;
t63 = t51 * g(2) - t55 * g(3);
t33 = -t56 * pkin(1) + t63;
t50 = sin(qJ(2));
t54 = cos(qJ(2));
t61 = t54 * t32 - t50 * t33;
t22 = t43 * pkin(2) + t61;
t42 = t44 ^ 2;
t67 = t50 * t32 + t54 * t33;
t23 = -t42 * pkin(2) + t67;
t46 = sin(pkin(9));
t47 = cos(pkin(9));
t62 = t47 * t22 - t46 * t23;
t17 = t43 * pkin(3) + t62;
t68 = t46 * t22 + t47 * t23;
t18 = -t42 * pkin(3) + t68;
t49 = sin(qJ(4));
t53 = cos(qJ(4));
t59 = t53 * t17 - t49 * t18;
t73 = (t48 * t30 - t52 * t31) * t36 + m(6) * (-t35 * pkin(4) - t34 * pkin(8) - t59) - t26 * mrSges(6,1) + t25 * mrSges(6,2);
t72 = -m(2) - m(3);
t69 = t49 * t17 + t53 * t18;
t14 = -t34 * pkin(4) + t35 * pkin(8) + t69;
t24 = (-mrSges(6,1) * t52 + mrSges(6,2) * t48) * t36;
t45 = -g(1) + qJDD(3);
t11 = m(6) * (-t48 * t14 + t52 * t45) - t25 * mrSges(6,3) + qJDD(5) * mrSges(6,1) - t24 * t71 + qJD(5) * t31;
t12 = m(6) * (t52 * t14 + t48 * t45) + t26 * mrSges(6,3) - qJDD(5) * mrSges(6,2) + t24 * t70 - qJD(5) * t30;
t64 = m(5) * t45 + t52 * t11 + t48 * t12;
t60 = m(4) * t45 + t64;
t8 = m(5) * t59 + t35 * mrSges(5,1) - t34 * mrSges(5,2) - t73;
t7 = m(5) * t69 - t34 * mrSges(5,1) - t35 * mrSges(5,2) - t48 * t11 + t52 * t12;
t6 = m(4) * t68 - t42 * mrSges(4,1) - t43 * mrSges(4,2) - t49 * t8 + t53 * t7;
t5 = m(4) * t62 + t43 * mrSges(4,1) - t42 * mrSges(4,2) + t49 * t7 + t53 * t8;
t4 = m(3) * t67 - t42 * mrSges(3,1) - t43 * mrSges(3,2) - t46 * t5 + t47 * t6;
t3 = m(3) * t61 + t43 * mrSges(3,1) - t42 * mrSges(3,2) + t46 * t6 + t47 * t5;
t2 = m(2) * t66 + qJDD(1) * mrSges(2,1) - t56 * mrSges(2,2) + t54 * t3 + t50 * t4;
t1 = m(2) * t63 - t56 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t50 * t3 + t54 * t4;
t9 = [(-m(1) + t72) * g(1) + t60, t1, t4, t6, t7, t12; -m(1) * g(2) - t51 * t1 - t55 * t2, t2, t3, t5, t8, t11; -m(1) * g(3) + t55 * t1 - t51 * t2, t72 * g(1) + t60, -m(3) * g(1) + t60, t60, t64, t73;];
f_new = t9;
