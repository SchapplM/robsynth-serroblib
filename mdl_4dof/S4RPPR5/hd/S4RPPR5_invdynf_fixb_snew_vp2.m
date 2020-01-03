% Calculate vector of cutting forces with Newton-Euler
% S4RPPR5
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
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4RPPR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR5_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR5_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_invdynf_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR5_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR5_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR5_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:39:43
% EndTime: 2019-12-31 16:39:44
% DurationCPUTime: 0.21s
% Computational Cost: add. (1197->71), mult. (2006->91), div. (0->0), fcn. (714->6), ass. (0->39)
t28 = sin(qJ(4));
t30 = cos(qJ(4));
t43 = qJD(1) * qJD(4);
t17 = -t28 * qJDD(1) - t30 * t43;
t18 = -t30 * qJDD(1) + t28 * t43;
t45 = qJD(1) * t28;
t19 = qJD(4) * mrSges(5,1) + mrSges(5,3) * t45;
t44 = qJD(1) * t30;
t20 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t44;
t32 = qJD(1) ^ 2;
t29 = sin(qJ(1));
t31 = cos(qJ(1));
t41 = -t31 * g(1) - t29 * g(2);
t37 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t41;
t48 = -pkin(1) - pkin(2);
t12 = t48 * t32 + t37;
t42 = t29 * g(1) - t31 * g(2);
t33 = -t32 * qJ(2) + qJDD(2) - t42;
t14 = t48 * qJDD(1) + t33;
t26 = sin(pkin(6));
t27 = cos(pkin(6));
t40 = -t26 * t12 + t27 * t14;
t50 = -(t28 * t19 - t30 * t20) * qJD(1) + m(5) * (qJDD(1) * pkin(3) - t32 * pkin(5) - t40) - t18 * mrSges(5,1) + t17 * mrSges(5,2);
t49 = -m(2) - m(3);
t47 = mrSges(2,1) + mrSges(3,1);
t46 = t27 * t12 + t26 * t14;
t16 = (mrSges(5,1) * t30 - mrSges(5,2) * t28) * qJD(1);
t25 = g(3) + qJDD(3);
t9 = -t32 * pkin(3) - qJDD(1) * pkin(5) + t46;
t6 = m(5) * (t30 * t25 - t28 * t9) - t17 * mrSges(5,3) + qJDD(4) * mrSges(5,1) + t16 * t45 + qJD(4) * t20;
t7 = m(5) * (t28 * t25 + t30 * t9) + t18 * mrSges(5,3) - qJDD(4) * mrSges(5,2) - t16 * t44 - qJD(4) * t19;
t4 = m(4) * t46 - t32 * mrSges(4,1) + qJDD(1) * mrSges(4,2) - t28 * t6 + t30 * t7;
t5 = m(4) * t40 - qJDD(1) * mrSges(4,1) - t32 * mrSges(4,2) - t50;
t38 = -t26 * t5 + t27 * t4 + m(3) * (-t32 * pkin(1) + t37) + qJDD(1) * mrSges(3,3);
t36 = -m(3) * (-qJDD(1) * pkin(1) + t33) - t26 * t4 - t27 * t5;
t35 = m(4) * t25 + t28 * t7 + t30 * t6;
t2 = m(2) * t42 + (-mrSges(2,2) + mrSges(3,3)) * t32 + t47 * qJDD(1) + t36;
t1 = m(2) * t41 - qJDD(1) * mrSges(2,2) - t47 * t32 + t38;
t3 = [-m(1) * g(1) + t31 * t1 - t29 * t2, t1, -t32 * mrSges(3,1) + t38, t4, t7; -m(1) * g(2) + t29 * t1 + t31 * t2, t2, -m(3) * g(3) - t35, t5, t6; (-m(1) + t49) * g(3) - t35, t49 * g(3) - t35, -qJDD(1) * mrSges(3,1) - t32 * mrSges(3,3) - t36, t35, t50;];
f_new = t3;
