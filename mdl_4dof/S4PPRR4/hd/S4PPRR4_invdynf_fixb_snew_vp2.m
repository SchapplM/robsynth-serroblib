% Calculate vector of cutting forces with Newton-Euler
% S4PPRR4
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
%   pkin=[a2,a3,a4,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-31 16:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4PPRR4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR4_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR4_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR4_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR4_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR4_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:18:28
% EndTime: 2019-12-31 16:18:28
% DurationCPUTime: 0.22s
% Computational Cost: add. (1351->55), mult. (2250->83), div. (0->0), fcn. (1424->8), ass. (0->39)
t33 = sin(qJ(4));
t35 = cos(qJ(4));
t45 = qJD(3) * qJD(4);
t21 = t33 * qJDD(3) + t35 * t45;
t22 = t35 * qJDD(3) - t33 * t45;
t47 = qJD(3) * t33;
t25 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t47;
t46 = qJD(3) * t35;
t26 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t46;
t37 = qJD(3) ^ 2;
t30 = sin(pkin(6));
t32 = cos(pkin(6));
t24 = -t32 * g(1) - t30 * g(2);
t28 = -g(3) + qJDD(1);
t29 = sin(pkin(7));
t31 = cos(pkin(7));
t17 = -t29 * t24 + t31 * t28;
t18 = t31 * t24 + t29 * t28;
t34 = sin(qJ(3));
t36 = cos(qJ(3));
t41 = t36 * t17 - t34 * t18;
t49 = (t33 * t25 - t35 * t26) * qJD(3) + m(5) * (-qJDD(3) * pkin(3) - t37 * pkin(5) - t41) - t22 * mrSges(5,1) + t21 * mrSges(5,2);
t48 = t34 * t17 + t36 * t18;
t14 = -t37 * pkin(3) + qJDD(3) * pkin(5) + t48;
t20 = (-mrSges(5,1) * t35 + mrSges(5,2) * t33) * qJD(3);
t42 = t30 * g(1) - t32 * g(2);
t23 = qJDD(2) - t42;
t11 = m(5) * (-t33 * t14 + t35 * t23) - t21 * mrSges(5,3) + qJDD(4) * mrSges(5,1) - t20 * t47 + qJD(4) * t26;
t12 = m(5) * (t35 * t14 + t33 * t23) + t22 * mrSges(5,3) - qJDD(4) * mrSges(5,2) + t20 * t46 - qJD(4) * t25;
t6 = m(4) * t48 - t37 * mrSges(4,1) - qJDD(3) * mrSges(4,2) - t33 * t11 + t35 * t12;
t8 = m(4) * t41 + qJDD(3) * mrSges(4,1) - t37 * mrSges(4,2) - t49;
t4 = m(3) * t17 + t34 * t6 + t36 * t8;
t5 = m(3) * t18 - t34 * t8 + t36 * t6;
t44 = m(2) * t28 + t29 * t5 + t31 * t4;
t43 = m(4) * t23 + t35 * t11 + t33 * t12;
t39 = m(3) * t23 + t43;
t7 = m(2) * t42 - t39;
t1 = m(2) * t24 - t29 * t4 + t31 * t5;
t2 = [-m(1) * g(1) + t32 * t1 - t30 * t7, t1, t5, t6, t12; -m(1) * g(2) + t30 * t1 + t32 * t7, t7, t4, t8, t11; -m(1) * g(3) + t44, t44, t39, t43, t49;];
f_new = t2;
