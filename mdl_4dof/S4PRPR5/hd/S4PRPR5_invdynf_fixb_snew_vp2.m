% Calculate vector of cutting forces with Newton-Euler
% S4PRPR5
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x5]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4PRPR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR5_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR5_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR5_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR5_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR5_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR5_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR5_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:22:59
% EndTime: 2019-12-31 16:23:00
% DurationCPUTime: 0.24s
% Computational Cost: add. (1527->62), mult. (2520->89), div. (0->0), fcn. (1424->8), ass. (0->40)
t38 = sin(qJ(4));
t40 = cos(qJ(4));
t49 = qJD(2) * qJD(4);
t22 = t38 * qJDD(2) + t40 * t49;
t23 = t40 * qJDD(2) - t38 * t49;
t51 = qJD(2) * t38;
t27 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t51;
t50 = qJD(2) * t40;
t28 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t50;
t42 = qJD(2) ^ 2;
t35 = sin(pkin(6));
t37 = cos(pkin(6));
t26 = -t37 * g(1) - t35 * g(2);
t33 = -g(3) + qJDD(1);
t39 = sin(qJ(2));
t41 = cos(qJ(2));
t46 = -t39 * t26 + t41 * t33;
t17 = qJDD(2) * pkin(2) + t46;
t52 = t41 * t26 + t39 * t33;
t18 = -t42 * pkin(2) + t52;
t34 = sin(pkin(7));
t36 = cos(pkin(7));
t45 = t36 * t17 - t34 * t18;
t54 = (t38 * t27 - t40 * t28) * qJD(2) + m(5) * (-qJDD(2) * pkin(3) - t42 * pkin(5) - t45) - t23 * mrSges(5,1) + t22 * mrSges(5,2);
t53 = t34 * t17 + t36 * t18;
t14 = -t42 * pkin(3) + qJDD(2) * pkin(5) + t53;
t21 = (-mrSges(5,1) * t40 + mrSges(5,2) * t38) * qJD(2);
t25 = t35 * g(1) - t37 * g(2);
t24 = qJDD(3) - t25;
t11 = m(5) * (-t38 * t14 + t40 * t24) - t22 * mrSges(5,3) + qJDD(4) * mrSges(5,1) - t21 * t51 + qJD(4) * t28;
t12 = m(5) * (t40 * t14 + t38 * t24) + t23 * mrSges(5,3) - qJDD(4) * mrSges(5,2) + t21 * t50 - qJD(4) * t27;
t6 = m(4) * t53 - t42 * mrSges(4,1) - qJDD(2) * mrSges(4,2) - t38 * t11 + t40 * t12;
t8 = m(4) * t45 + qJDD(2) * mrSges(4,1) - t42 * mrSges(4,2) - t54;
t4 = m(3) * t46 + qJDD(2) * mrSges(3,1) - t42 * mrSges(3,2) + t34 * t6 + t36 * t8;
t5 = m(3) * t52 - t42 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t34 * t8 + t36 * t6;
t48 = m(2) * t33 + t39 * t5 + t41 * t4;
t47 = m(4) * t24 + t40 * t11 + t38 * t12;
t7 = (m(2) + m(3)) * t25 - t47;
t1 = m(2) * t26 - t39 * t4 + t41 * t5;
t2 = [-m(1) * g(1) + t37 * t1 - t35 * t7, t1, t5, t6, t12; -m(1) * g(2) + t35 * t1 + t37 * t7, t7, t4, t8, t11; -m(1) * g(3) + t48, t48, -m(3) * t25 + t47, t47, t54;];
f_new = t2;
