% Calculate vector of cutting forces with Newton-Euler
% S5PPPRR2
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
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PPPRR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR2_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPPRR2_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR2_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPPRR2_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPPRR2_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:28
% EndTime: 2019-12-05 14:59:29
% DurationCPUTime: 0.33s
% Computational Cost: add. (3098->64), mult. (4724->94), div. (0->0), fcn. (3392->10), ass. (0->47)
t36 = sin(qJ(5));
t38 = cos(qJ(5));
t50 = qJD(4) * qJD(5);
t22 = t36 * qJDD(4) + t38 * t50;
t23 = t38 * qJDD(4) - t36 * t50;
t52 = qJD(4) * t36;
t26 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t52;
t51 = qJD(4) * t38;
t27 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t51;
t40 = qJD(4) ^ 2;
t32 = sin(pkin(7));
t35 = cos(pkin(7));
t25 = -t35 * g(1) - t32 * g(2);
t29 = -g(3) + qJDD(1);
t31 = sin(pkin(8));
t34 = cos(pkin(8));
t20 = t34 * t25 + t31 * t29;
t48 = t32 * g(1) - t35 * g(2);
t24 = qJDD(2) - t48;
t30 = sin(pkin(9));
t33 = cos(pkin(9));
t17 = t33 * t20 + t30 * t24;
t45 = -t31 * t25 + t34 * t29;
t19 = qJDD(3) - t45;
t37 = sin(qJ(4));
t39 = cos(qJ(4));
t46 = -t37 * t17 + t39 * t19;
t54 = (t36 * t26 - t38 * t27) * qJD(4) + m(6) * (-qJDD(4) * pkin(4) - t40 * pkin(6) - t46) - t23 * mrSges(6,1) + t22 * mrSges(6,2);
t53 = t39 * t17 + t37 * t19;
t10 = m(5) * t46 + qJDD(4) * mrSges(5,1) - t40 * mrSges(5,2) - t54;
t14 = -t40 * pkin(4) + qJDD(4) * pkin(6) + t53;
t16 = t30 * t20 - t33 * t24;
t21 = (-mrSges(6,1) * t38 + mrSges(6,2) * t36) * qJD(4);
t11 = m(6) * (-t36 * t14 + t38 * t16) - t22 * mrSges(6,3) + qJDD(5) * mrSges(6,1) - t21 * t52 + qJD(5) * t27;
t12 = m(6) * (t38 * t14 + t36 * t16) + t23 * mrSges(6,3) - qJDD(5) * mrSges(6,2) + t21 * t51 - qJD(5) * t26;
t9 = m(5) * t53 - t40 * mrSges(5,1) - qJDD(4) * mrSges(5,2) - t36 * t11 + t38 * t12;
t7 = m(4) * t17 - t37 * t10 + t39 * t9;
t47 = -t38 * t11 - t36 * t12;
t8 = (-m(4) - m(5)) * t16 + t47;
t4 = m(3) * t20 - t30 * t8 + t33 * t7;
t42 = m(4) * t19 + t39 * t10 + t37 * t9;
t6 = m(3) * t45 - t42;
t49 = m(2) * t29 + t31 * t4 + t34 * t6;
t43 = m(3) * t24 + t30 * t7 + t33 * t8;
t3 = m(2) * t48 - t43;
t1 = m(2) * t25 - t31 * t6 + t34 * t4;
t2 = [-m(1) * g(1) + t35 * t1 - t32 * t3, t1, t4, t7, t9, t12; -m(1) * g(2) + t32 * t1 + t35 * t3, t3, t6, t8, t10, t11; -m(1) * g(3) + t49, t49, t43, t42, m(5) * t16 - t47, t54;];
f_new = t2;
