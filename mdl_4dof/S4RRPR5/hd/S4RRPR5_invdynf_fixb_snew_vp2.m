% Calculate vector of cutting forces with Newton-Euler
% S4RRPR5
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4RRPR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR5_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR5_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_invdynf_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR5_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR5_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR5_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:03:29
% EndTime: 2019-12-31 17:03:29
% DurationCPUTime: 0.20s
% Computational Cost: add. (1381->71), mult. (1663->89), div. (0->0), fcn. (722->6), ass. (0->43)
t54 = -m(3) - m(4);
t53 = -pkin(2) - pkin(6);
t29 = qJD(1) + qJD(2);
t30 = sin(qJ(4));
t52 = t29 * t30;
t33 = cos(qJ(4));
t51 = t29 * t33;
t50 = mrSges(3,1) - mrSges(4,2);
t49 = -mrSges(3,2) + mrSges(4,3);
t32 = sin(qJ(1));
t35 = cos(qJ(1));
t44 = t32 * g(1) - t35 * g(2);
t22 = qJDD(1) * pkin(1) + t44;
t36 = qJD(1) ^ 2;
t42 = -t35 * g(1) - t32 * g(2);
t23 = -t36 * pkin(1) + t42;
t31 = sin(qJ(2));
t34 = cos(qJ(2));
t48 = t31 * t22 + t34 * t23;
t47 = qJD(4) * t29;
t46 = -m(2) + t54;
t28 = qJDD(1) + qJDD(2);
t27 = t29 ^ 2;
t43 = t34 * t22 - t31 * t23;
t38 = -t27 * qJ(3) + qJDD(3) - t43;
t10 = t53 * t28 + t38;
t16 = (mrSges(5,1) * t30 + mrSges(5,2) * t33) * t29;
t18 = t33 * t28 - t30 * t47;
t24 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t52;
t6 = m(5) * (t30 * g(3) + t33 * t10) - t18 * mrSges(5,3) + qJDD(4) * mrSges(5,1) - t16 * t51 + qJD(4) * t24;
t17 = -t30 * t28 - t33 * t47;
t25 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t51;
t7 = m(5) * (-t33 * g(3) + t30 * t10) + t17 * mrSges(5,3) - qJDD(4) * mrSges(5,2) - t16 * t52 - qJD(4) * t25;
t45 = -t30 * t6 + t33 * t7;
t39 = t28 * qJ(3) + 0.2e1 * qJD(3) * t29 + t48;
t41 = -t17 * mrSges(5,1) + t24 * t52 + t25 * t51 + t18 * mrSges(5,2) + m(5) * (t53 * t27 + t39);
t40 = -m(4) * (-t28 * pkin(2) + t38) - t30 * t7 - t33 * t6;
t37 = -m(4) * (t27 * pkin(2) - t39) + t41;
t4 = m(3) * t48 - t50 * t27 + t49 * t28 + t37;
t3 = m(3) * t43 + t49 * t27 + t50 * t28 + t40;
t2 = m(2) * t42 - t36 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t31 * t3 + t34 * t4;
t1 = m(2) * t44 + qJDD(1) * mrSges(2,1) - t36 * mrSges(2,2) + t34 * t3 + t31 * t4;
t5 = [-m(1) * g(1) - t32 * t1 + t35 * t2, t2, t4, -m(4) * g(3) + t45, t7; -m(1) * g(2) + t35 * t1 + t32 * t2, t1, t3, -t27 * mrSges(4,2) - t28 * mrSges(4,3) - t37, t6; (-m(1) + t46) * g(3) + t45, t46 * g(3) + t45, t54 * g(3) + t45, t28 * mrSges(4,2) - t27 * mrSges(4,3) - t40, t41;];
f_new = t5;
