% Calculate vector of inverse dynamics joint torques for
% S4PPRR2
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
%   pkin=[a2,a3,a4,d3,d4,theta2]';
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
% Datum: 2018-11-14 14:00
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau = S4PPRR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR2_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR2_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR2_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR2_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR2_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR2_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR2_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:59:32
% EndTime: 2018-11-14 13:59:33
% DurationCPUTime: 0.26s
% Computational Cost: add. (336->78), mult. (798->104), div. (0->0), fcn. (624->10), ass. (0->44)
t34 = sin(pkin(6));
t35 = cos(pkin(6));
t37 = sin(qJ(3));
t39 = cos(qJ(3));
t21 = -t37 * t34 + t39 * t35;
t32 = pkin(6) + qJ(3);
t30 = qJ(4) + t32;
t25 = sin(t30);
t50 = g(1) * t25;
t20 = t39 * t34 + t37 * t35;
t17 = t20 * qJD(1);
t36 = sin(qJ(4));
t49 = t36 * t17;
t38 = cos(qJ(4));
t47 = t38 * t17;
t45 = qJD(4) * t36;
t44 = qJD(4) * t38;
t43 = m(3) + m(4) + m(5);
t26 = cos(t30);
t19 = t20 * qJD(3);
t12 = -qJD(1) * t19 + t21 * qJDD(1);
t10 = qJDD(3) * pkin(3) + t12;
t18 = t21 * qJD(3);
t11 = qJD(1) * t18 + t20 * qJDD(1);
t16 = t21 * qJD(1);
t15 = qJD(3) * pkin(3) + t16;
t7 = t36 * t15 + t47;
t3 = -qJD(4) * t7 + t38 * t10 - t36 * t11;
t31 = qJDD(3) + qJDD(4);
t42 = Ifges(5,3) * t31 + (-g(2) * t26 + t3) * mrSges(5,1);
t28 = sin(t32);
t29 = cos(t32);
t41 = g(1) * t28 - g(2) * t29;
t6 = t38 * t15 - t49;
t14 = t38 * t20 + t36 * t21;
t13 = -t36 * t20 + t38 * t21;
t2 = qJD(4) * t6 + t36 * t10 + t38 * t11;
t40 = g(1) * t26 + g(2) * t25 - t2;
t33 = qJD(3) + qJD(4);
t9 = t38 * t16 - t49;
t8 = -t36 * t16 - t47;
t5 = -t14 * qJD(4) - t36 * t18 - t38 * t19;
t4 = t13 * qJD(4) + t38 * t18 - t36 * t19;
t1 = [(-t14 * t31 - t4 * t33) * mrSges(5,2) + (-t18 * qJD(3) - t20 * qJDD(3)) * mrSges(4,2) + (t13 * t31 + t5 * t33) * mrSges(5,1) + (-t19 * qJD(3) + t21 * qJDD(3)) * mrSges(4,1) + m(4) * (t11 * t20 + t12 * t21 - t16 * t19 + t17 * t18) + m(5) * (t3 * t13 + t2 * t14 + t7 * t4 + t6 * t5) + (-m(2) - t43) * g(2) + (m(2) + m(3) * (t34 ^ 2 + t35 ^ 2)) * qJDD(1); (-g(3) + qJDD(2)) * t43; Ifges(4,3) * qJDD(3) - m(5) * (t6 * t8 + t7 * t9) + (-t8 * t33 + t50) * mrSges(5,1) + (t9 * t33 + t40) * mrSges(5,2) + (g(1) * t29 + g(2) * t28 + t16 * qJD(3) - t11) * mrSges(4,2) + (t17 * qJD(3) + t12 + t41) * mrSges(4,1) + ((-t36 * t31 - t33 * t44) * mrSges(5,2) + (t38 * t31 - t33 * t45) * mrSges(5,1) + (t2 * t36 + t3 * t38 + t7 * t44 - t6 * t45 + t41) * m(5)) * pkin(3) + t42; (t7 * t33 + t50) * mrSges(5,1) + (t6 * t33 + t40) * mrSges(5,2) + t42;];
tau  = t1;
