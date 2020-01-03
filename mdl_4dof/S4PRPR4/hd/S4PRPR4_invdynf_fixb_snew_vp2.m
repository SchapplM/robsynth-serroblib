% Calculate vector of cutting forces with Newton-Euler
% S4PRPR4
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
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% Datum: 2019-12-31 16:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4PRPR4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR4_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR4_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR4_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR4_invdynf_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR4_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR4_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR4_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:21:57
% EndTime: 2019-12-31 16:21:57
% DurationCPUTime: 0.18s
% Computational Cost: add. (832->61), mult. (1450->82), div. (0->0), fcn. (722->6), ass. (0->39)
t52 = -pkin(2) - pkin(5);
t51 = mrSges(3,1) - mrSges(4,2);
t50 = -mrSges(3,2) + mrSges(4,3);
t30 = sin(pkin(6));
t31 = cos(pkin(6));
t22 = t30 * g(1) - t31 * g(2);
t23 = -t31 * g(1) - t30 * g(2);
t33 = sin(qJ(2));
t35 = cos(qJ(2));
t49 = t33 * t22 + t35 * t23;
t32 = sin(qJ(4));
t48 = qJD(2) * t32;
t34 = cos(qJ(4));
t47 = qJD(2) * t34;
t46 = qJD(2) * qJD(4);
t45 = t35 * t22 - t33 * t23;
t29 = -g(3) + qJDD(1);
t36 = qJD(2) ^ 2;
t39 = -t36 * qJ(3) + qJDD(3) - t45;
t10 = t52 * qJDD(2) + t39;
t19 = (mrSges(5,1) * t32 + mrSges(5,2) * t34) * qJD(2);
t21 = t34 * qJDD(2) - t32 * t46;
t24 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t48;
t6 = m(5) * (t34 * t10 - t32 * t29) - t21 * mrSges(5,3) + qJDD(4) * mrSges(5,1) - t19 * t47 + qJD(4) * t24;
t20 = -t32 * qJDD(2) - t34 * t46;
t25 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t47;
t7 = m(5) * (t32 * t10 + t34 * t29) + t20 * mrSges(5,3) - qJDD(4) * mrSges(5,2) - t19 * t48 - qJD(4) * t25;
t44 = m(4) * t29 - t32 * t6 + t34 * t7;
t43 = m(3) * t29 + t44;
t42 = m(2) * t29 + t43;
t38 = qJDD(2) * qJ(3) + 0.2e1 * qJD(3) * qJD(2) + t49;
t41 = -t20 * mrSges(5,1) + t24 * t48 + t25 * t47 + t21 * mrSges(5,2) + m(5) * (t52 * t36 + t38);
t40 = -m(4) * (-qJDD(2) * pkin(2) + t39) - t32 * t7 - t34 * t6;
t37 = -m(4) * (t36 * pkin(2) - t38) + t41;
t4 = m(3) * t49 + t50 * qJDD(2) - t51 * t36 + t37;
t3 = m(3) * t45 + t51 * qJDD(2) + t50 * t36 + t40;
t2 = m(2) * t23 - t33 * t3 + t35 * t4;
t1 = m(2) * t22 + t35 * t3 + t33 * t4;
t5 = [-m(1) * g(1) - t30 * t1 + t31 * t2, t2, t4, t44, t7; -m(1) * g(2) + t31 * t1 + t30 * t2, t1, t3, -t36 * mrSges(4,2) - qJDD(2) * mrSges(4,3) - t37, t6; -m(1) * g(3) + t42, t42, t43, qJDD(2) * mrSges(4,2) - t36 * mrSges(4,3) - t40, t41;];
f_new = t5;
