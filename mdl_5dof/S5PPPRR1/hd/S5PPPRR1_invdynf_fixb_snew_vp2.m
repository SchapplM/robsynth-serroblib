% Calculate vector of cutting forces with Newton-Euler
% S5PPPRR1
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
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
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
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PPPRR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR1_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPPRR1_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR1_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPPRR1_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPPRR1_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:57:55
% EndTime: 2019-12-05 14:57:56
% DurationCPUTime: 0.35s
% Computational Cost: add. (3394->63), mult. (5232->94), div. (0->0), fcn. (3816->10), ass. (0->47)
t39 = sin(qJ(5));
t41 = cos(qJ(5));
t53 = qJD(4) * qJD(5);
t25 = t39 * qJDD(4) + t41 * t53;
t26 = t41 * qJDD(4) - t39 * t53;
t55 = qJD(4) * t39;
t29 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t55;
t54 = qJD(4) * t41;
t30 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t54;
t43 = qJD(4) ^ 2;
t35 = sin(pkin(7));
t38 = cos(pkin(7));
t28 = -t38 * g(1) - t35 * g(2);
t32 = -g(3) + qJDD(1);
t34 = sin(pkin(8));
t37 = cos(pkin(8));
t23 = t37 * t28 + t34 * t32;
t50 = t35 * g(1) - t38 * g(2);
t27 = qJDD(2) - t50;
t33 = sin(pkin(9));
t36 = cos(pkin(9));
t19 = -t33 * t23 + t36 * t27;
t20 = t36 * t23 + t33 * t27;
t40 = sin(qJ(4));
t42 = cos(qJ(4));
t49 = t42 * t19 - t40 * t20;
t57 = (t39 * t29 - t41 * t30) * qJD(4) + m(6) * (-qJDD(4) * pkin(4) - t43 * pkin(6) - t49) - t26 * mrSges(6,1) + t25 * mrSges(6,2);
t56 = t40 * t19 + t42 * t20;
t10 = m(5) * t49 + qJDD(4) * mrSges(5,1) - t43 * mrSges(5,2) - t57;
t16 = -t43 * pkin(4) + qJDD(4) * pkin(6) + t56;
t48 = -t34 * t28 + t37 * t32;
t22 = qJDD(3) - t48;
t24 = (-mrSges(6,1) * t41 + mrSges(6,2) * t39) * qJD(4);
t13 = m(6) * (-t39 * t16 + t41 * t22) - t25 * mrSges(6,3) + qJDD(5) * mrSges(6,1) - t24 * t55 + qJD(5) * t30;
t14 = m(6) * (t41 * t16 + t39 * t22) + t26 * mrSges(6,3) - qJDD(5) * mrSges(6,2) + t24 * t54 - qJD(5) * t29;
t7 = m(5) * t56 - t43 * mrSges(5,1) - qJDD(4) * mrSges(5,2) - t39 * t13 + t41 * t14;
t5 = m(4) * t19 + t42 * t10 + t40 * t7;
t6 = m(4) * t20 - t40 * t10 + t42 * t7;
t4 = m(3) * t23 - t33 * t5 + t36 * t6;
t51 = m(5) * t22 + t41 * t13 + t39 * t14;
t46 = m(4) * t22 + t51;
t9 = m(3) * t48 - t46;
t52 = m(2) * t32 + t34 * t4 + t37 * t9;
t45 = m(3) * t27 + t33 * t6 + t36 * t5;
t3 = m(2) * t50 - t45;
t1 = m(2) * t28 - t34 * t9 + t37 * t4;
t2 = [-m(1) * g(1) + t38 * t1 - t35 * t3, t1, t4, t6, t7, t14; -m(1) * g(2) + t35 * t1 + t38 * t3, t3, t9, t5, t10, t13; -m(1) * g(3) + t52, t52, t45, t46, t51, t57;];
f_new = t2;
