% Calculate vector of cutting forces with Newton-Euler
% S5PRPRR2
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
% Datum: 2019-12-05 15:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRPRR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR2_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR2_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR2_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR2_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR2_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:44:51
% EndTime: 2019-12-05 15:44:52
% DurationCPUTime: 0.48s
% Computational Cost: add. (5862->77), mult. (7808->107), div. (0->0), fcn. (4710->10), ass. (0->52)
t42 = qJDD(2) + qJDD(4);
t49 = sin(qJ(5));
t52 = cos(qJ(5));
t43 = qJD(2) + qJD(4);
t64 = qJD(5) * t43;
t27 = t49 * t42 + t52 * t64;
t28 = t52 * t42 - t49 * t64;
t69 = t43 * t49;
t32 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t69;
t68 = t43 * t52;
t33 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t68;
t41 = t43 ^ 2;
t46 = sin(pkin(8));
t48 = cos(pkin(8));
t36 = -t48 * g(1) - t46 * g(2);
t44 = -g(3) + qJDD(1);
t51 = sin(qJ(2));
t54 = cos(qJ(2));
t60 = -t51 * t36 + t54 * t44;
t24 = qJDD(2) * pkin(2) + t60;
t55 = qJD(2) ^ 2;
t65 = t54 * t36 + t51 * t44;
t25 = -t55 * pkin(2) + t65;
t45 = sin(pkin(9));
t47 = cos(pkin(9));
t61 = t47 * t24 - t45 * t25;
t19 = qJDD(2) * pkin(3) + t61;
t66 = t45 * t24 + t47 * t25;
t20 = -t55 * pkin(3) + t66;
t50 = sin(qJ(4));
t53 = cos(qJ(4));
t58 = t53 * t19 - t50 * t20;
t70 = (t49 * t32 - t52 * t33) * t43 + m(6) * (-t42 * pkin(4) - t41 * pkin(7) - t58) - t28 * mrSges(6,1) + t27 * mrSges(6,2);
t67 = t50 * t19 + t53 * t20;
t10 = m(5) * t58 + t42 * mrSges(5,1) - t41 * mrSges(5,2) - t70;
t16 = -t41 * pkin(4) + t42 * pkin(7) + t67;
t26 = (-mrSges(6,1) * t52 + mrSges(6,2) * t49) * t43;
t35 = t46 * g(1) - t48 * g(2);
t34 = qJDD(3) - t35;
t13 = m(6) * (-t49 * t16 + t52 * t34) - t27 * mrSges(6,3) + qJDD(5) * mrSges(6,1) - t26 * t69 + qJD(5) * t33;
t14 = m(6) * (t52 * t16 + t49 * t34) + t28 * mrSges(6,3) - qJDD(5) * mrSges(6,2) + t26 * t68 - qJD(5) * t32;
t8 = m(5) * t67 - t41 * mrSges(5,1) - t42 * mrSges(5,2) - t49 * t13 + t52 * t14;
t6 = m(4) * t61 + qJDD(2) * mrSges(4,1) - t55 * mrSges(4,2) + t53 * t10 + t50 * t8;
t7 = m(4) * t66 - t55 * mrSges(4,1) - qJDD(2) * mrSges(4,2) - t50 * t10 + t53 * t8;
t4 = m(3) * t60 + qJDD(2) * mrSges(3,1) - t55 * mrSges(3,2) + t45 * t7 + t47 * t6;
t5 = m(3) * t65 - t55 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t45 * t6 + t47 * t7;
t63 = m(2) * t44 + t54 * t4 + t51 * t5;
t62 = m(5) * t34 + t52 * t13 + t49 * t14;
t59 = m(4) * t34 + t62;
t9 = (m(2) + m(3)) * t35 - t59;
t1 = m(2) * t36 - t51 * t4 + t54 * t5;
t2 = [-m(1) * g(1) + t48 * t1 - t46 * t9, t1, t5, t7, t8, t14; -m(1) * g(2) + t46 * t1 + t48 * t9, t9, t4, t6, t10, t13; -m(1) * g(3) + t63, t63, -m(3) * t35 + t59, t59, t62, t70;];
f_new = t2;
