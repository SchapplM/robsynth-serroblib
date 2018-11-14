% Calculate vector of inverse dynamics joint torques for
% S4PPRR1
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
%   pkin=[a2,a3,a4,d3,d4,theta1]';
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:40
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau = S4PPRR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR1_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR1_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR1_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR1_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR1_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR1_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR1_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:40:08
% EndTime: 2018-11-14 13:40:09
% DurationCPUTime: 0.28s
% Computational Cost: add. (214->75), mult. (428->103), div. (0->0), fcn. (278->8), ass. (0->43)
t29 = qJD(3) + qJD(4);
t51 = m(5) * pkin(3);
t28 = qJDD(3) + qJDD(4);
t50 = pkin(3) * t28;
t36 = cos(qJ(3));
t34 = sin(qJ(3));
t44 = qJD(2) * qJD(3);
t41 = t34 * t44;
t21 = qJDD(2) * t36 - t41;
t16 = qJDD(3) * pkin(3) + t21;
t40 = t36 * t44;
t22 = qJDD(2) * t34 + t40;
t33 = sin(qJ(4));
t35 = cos(qJ(4));
t23 = qJD(3) * pkin(3) + qJD(2) * t36;
t45 = qJD(2) * t34;
t7 = t23 * t33 + t35 * t45;
t3 = -qJD(4) * t7 + t35 * t16 - t33 * t22;
t49 = mrSges(5,1) * t3 + Ifges(5,3) * t28;
t30 = qJ(3) + qJ(4);
t26 = sin(t30);
t27 = cos(t30);
t31 = sin(pkin(6));
t32 = cos(pkin(6));
t12 = -t26 * t31 - t27 * t32;
t13 = t26 * t32 - t27 * t31;
t48 = -mrSges(5,1) * t13 + mrSges(5,2) * t12;
t47 = mrSges(5,1) * t12 + mrSges(5,2) * t13;
t46 = t29 * mrSges(5,1);
t43 = pkin(3) * qJD(4) * t29;
t39 = -g(1) * t31 + g(2) * t32;
t38 = t31 * t36 - t32 * t34;
t17 = -t31 * t34 - t32 * t36;
t20 = t33 * t36 + t34 * t35;
t19 = -t33 * t34 + t35 * t36;
t6 = t23 * t35 - t33 * t45;
t37 = qJD(3) ^ 2;
t15 = t19 * qJD(2);
t14 = t20 * qJD(2);
t5 = t29 * t20;
t4 = t29 * t19;
t2 = qJD(4) * t6 + t33 * t16 + t35 * t22;
t1 = [(-g(3) + qJDD(1)) * (m(2) + m(3) + m(4) + m(5)); (-t20 * t28 - t29 * t4) * mrSges(5,2) + (-qJDD(3) * t34 - t36 * t37) * mrSges(4,2) + (t19 * t28 - t29 * t5) * mrSges(5,1) + (qJDD(3) * t36 - t34 * t37) * mrSges(4,1) + (t19 * t3 + t2 * t20 + t4 * t7 - t5 * t6 + t39) * m(5) + (t21 * t36 + t22 * t34 + t39) * m(4) + (qJDD(2) + t39) * m(3); -g(1) * (t38 * t51 + t48) - g(2) * (t17 * t51 + t47) + Ifges(4,3) * qJDD(3) + t14 * t46 - m(5) * (-t14 * t6 + t15 * t7) + (t2 * t33 + t3 * t35 + (-t33 * t6 + t35 * t7) * qJD(4)) * t51 + t49 + (-t33 * t43 + t35 * t50) * mrSges(5,1) + (-g(1) * t17 + g(2) * t38 - t22 + t40) * mrSges(4,2) + (-g(1) * t38 - g(2) * t17 + t21 + t41) * mrSges(4,1) + (t15 * t29 - t33 * t50 - t35 * t43 - t2) * mrSges(5,2); t7 * t46 - g(1) * t48 - g(2) * t47 + (t29 * t6 - t2) * mrSges(5,2) + t49;];
tau  = t1;
