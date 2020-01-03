% Calculate vector of cutting forces with Newton-Euler
% S4RPPR4
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
%   pkin=[a2,a3,a4,d1,d4,theta2]';
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

function f_new = S4RPPR4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR4_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR4_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR4_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR4_invdynf_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR4_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR4_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR4_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:38:49
% EndTime: 2019-12-31 16:38:50
% DurationCPUTime: 0.19s
% Computational Cost: add. (972->68), mult. (1663->88), div. (0->0), fcn. (722->6), ass. (0->40)
t53 = -pkin(2) - pkin(5);
t52 = mrSges(3,1) - mrSges(4,2);
t51 = -mrSges(3,2) + mrSges(4,3);
t33 = sin(qJ(1));
t35 = cos(qJ(1));
t46 = t33 * g(1) - t35 * g(2);
t19 = qJDD(1) * pkin(1) + t46;
t36 = qJD(1) ^ 2;
t43 = -t35 * g(1) - t33 * g(2);
t21 = -t36 * pkin(1) + t43;
t30 = sin(pkin(6));
t31 = cos(pkin(6));
t50 = t30 * t19 + t31 * t21;
t32 = sin(qJ(4));
t49 = qJD(1) * t32;
t34 = cos(qJ(4));
t48 = qJD(1) * t34;
t47 = qJD(1) * qJD(4);
t45 = t31 * t19 - t30 * t21;
t29 = -g(3) + qJDD(2);
t39 = -t36 * qJ(3) + qJDD(3) - t45;
t10 = t53 * qJDD(1) + t39;
t20 = (mrSges(5,1) * t32 + mrSges(5,2) * t34) * qJD(1);
t23 = t34 * qJDD(1) - t32 * t47;
t24 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t49;
t6 = m(5) * (t34 * t10 - t32 * t29) - t23 * mrSges(5,3) + qJDD(4) * mrSges(5,1) - t20 * t48 + qJD(4) * t24;
t22 = -t32 * qJDD(1) - t34 * t47;
t25 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t48;
t7 = m(5) * (t32 * t10 + t34 * t29) + t22 * mrSges(5,3) - qJDD(4) * mrSges(5,2) - t20 * t49 - qJD(4) * t25;
t44 = m(4) * t29 - t32 * t6 + t34 * t7;
t42 = m(3) * t29 + t44;
t38 = qJDD(1) * qJ(3) + 0.2e1 * qJD(3) * qJD(1) + t50;
t41 = -t22 * mrSges(5,1) + t24 * t49 + t25 * t48 + t23 * mrSges(5,2) + m(5) * (t53 * t36 + t38);
t40 = -m(4) * (-qJDD(1) * pkin(2) + t39) - t32 * t7 - t34 * t6;
t37 = -m(4) * (t36 * pkin(2) - t38) + t41;
t4 = m(3) * t50 + t51 * qJDD(1) - t52 * t36 + t37;
t3 = m(3) * t45 + t52 * qJDD(1) + t51 * t36 + t40;
t2 = m(2) * t43 - t36 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t30 * t3 + t31 * t4;
t1 = m(2) * t46 + qJDD(1) * mrSges(2,1) - t36 * mrSges(2,2) + t31 * t3 + t30 * t4;
t5 = [-m(1) * g(1) - t33 * t1 + t35 * t2, t2, t4, t44, t7; -m(1) * g(2) + t35 * t1 + t33 * t2, t1, t3, -t36 * mrSges(4,2) - qJDD(1) * mrSges(4,3) - t37, t6; (-m(1) - m(2)) * g(3) + t42, -m(2) * g(3) + t42, t42, qJDD(1) * mrSges(4,2) - t36 * mrSges(4,3) - t40, t41;];
f_new = t5;
