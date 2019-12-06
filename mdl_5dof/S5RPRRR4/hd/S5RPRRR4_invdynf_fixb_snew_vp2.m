% Calculate vector of cutting forces with Newton-Euler
% S5RPRRR4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-05 18:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRRR4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR4_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR4_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR4_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR4_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:14:25
% EndTime: 2019-12-05 18:14:26
% DurationCPUTime: 0.59s
% Computational Cost: add. (9393->84), mult. (12295->114), div. (0->0), fcn. (6580->10), ass. (0->56)
t44 = qJDD(1) + qJDD(3);
t35 = qJDD(4) + t44;
t49 = sin(qJ(5));
t53 = cos(qJ(5));
t45 = qJD(1) + qJD(3);
t36 = qJD(4) + t45;
t67 = qJD(5) * t36;
t25 = t49 * t35 + t53 * t67;
t26 = t53 * t35 - t49 * t67;
t73 = t36 * t49;
t30 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t73;
t72 = t36 * t53;
t31 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t72;
t34 = t36 ^ 2;
t52 = sin(qJ(1));
t56 = cos(qJ(1));
t68 = t56 * g(2) + t52 * g(3);
t32 = qJDD(1) * pkin(1) + t68;
t57 = qJD(1) ^ 2;
t65 = t52 * g(2) - t56 * g(3);
t33 = -t57 * pkin(1) + t65;
t47 = sin(pkin(9));
t48 = cos(pkin(9));
t63 = t48 * t32 - t47 * t33;
t22 = qJDD(1) * pkin(2) + t63;
t69 = t47 * t32 + t48 * t33;
t23 = -t57 * pkin(2) + t69;
t51 = sin(qJ(3));
t55 = cos(qJ(3));
t64 = t55 * t22 - t51 * t23;
t17 = t44 * pkin(3) + t64;
t43 = t45 ^ 2;
t70 = t51 * t22 + t55 * t23;
t18 = -t43 * pkin(3) + t70;
t50 = sin(qJ(4));
t54 = cos(qJ(4));
t60 = t54 * t17 - t50 * t18;
t74 = (t49 * t30 - t53 * t31) * t36 + m(6) * (-t35 * pkin(4) - t34 * pkin(8) - t60) - t26 * mrSges(6,1) + t25 * mrSges(6,2);
t71 = t50 * t17 + t54 * t18;
t14 = -t34 * pkin(4) + t35 * pkin(8) + t71;
t24 = (-mrSges(6,1) * t53 + mrSges(6,2) * t49) * t36;
t46 = -g(1) + qJDD(2);
t11 = m(6) * (-t49 * t14 + t53 * t46) - t25 * mrSges(6,3) + qJDD(5) * mrSges(6,1) - t24 * t73 + qJD(5) * t31;
t12 = m(6) * (t53 * t14 + t49 * t46) + t26 * mrSges(6,3) - qJDD(5) * mrSges(6,2) + t24 * t72 - qJD(5) * t30;
t66 = m(5) * t46 + t53 * t11 + t49 * t12;
t62 = m(4) * t46 + t66;
t61 = m(3) * t46 + t62;
t8 = m(5) * t60 + t35 * mrSges(5,1) - t34 * mrSges(5,2) - t74;
t7 = m(5) * t71 - t34 * mrSges(5,1) - t35 * mrSges(5,2) - t49 * t11 + t53 * t12;
t6 = m(4) * t70 - t43 * mrSges(4,1) - t44 * mrSges(4,2) - t50 * t8 + t54 * t7;
t5 = m(4) * t64 + t44 * mrSges(4,1) - t43 * mrSges(4,2) + t50 * t7 + t54 * t8;
t4 = m(3) * t69 - t57 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t51 * t5 + t55 * t6;
t3 = m(3) * t63 + qJDD(1) * mrSges(3,1) - t57 * mrSges(3,2) + t55 * t5 + t51 * t6;
t2 = m(2) * t68 + qJDD(1) * mrSges(2,1) - t57 * mrSges(2,2) + t48 * t3 + t47 * t4;
t1 = m(2) * t65 - t57 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t47 * t3 + t48 * t4;
t9 = [(-m(1) - m(2)) * g(1) + t61, t1, t4, t6, t7, t12; -m(1) * g(2) - t52 * t1 - t56 * t2, t2, t3, t5, t8, t11; -m(1) * g(3) + t56 * t1 - t52 * t2, -m(2) * g(1) + t61, t61, t62, t66, t74;];
f_new = t9;
