% Calculate vector of inverse dynamics joint torques for
% S3PRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% qJDD [3x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';
% m_mdh [4x1]
%   mass of all robot links (including the base)
% mrSges [4x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [3x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:12
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau = S3PRR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRR1_invdynJ_fixb_slag_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRR1_invdynJ_fixb_slag_vp2: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3PRR1_invdynJ_fixb_slag_vp2: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3PRR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3PRR1_invdynJ_fixb_slag_vp2: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3PRR1_invdynJ_fixb_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3PRR1_invdynJ_fixb_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3PRR1_invdynJ_fixb_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:11:57
% EndTime: 2018-11-14 10:11:57
% DurationCPUTime: 0.23s
% Computational Cost: add. (182->64), mult. (381->89), div. (0->0), fcn. (230->6), ass. (0->34)
t22 = qJD(2) + qJD(3);
t23 = qJ(2) + qJ(3);
t19 = sin(t23);
t37 = g(1) * t19;
t27 = cos(qJ(2));
t36 = g(2) * t27;
t25 = sin(qJ(2));
t35 = qJD(1) * t25;
t24 = sin(qJ(3));
t34 = qJD(3) * t24;
t26 = cos(qJ(3));
t33 = qJD(3) * t26;
t32 = qJD(1) * qJD(2);
t31 = g(1) + t32;
t20 = cos(t23);
t21 = qJDD(2) + qJDD(3);
t13 = t27 * qJDD(1) - t25 * t32;
t10 = qJDD(2) * pkin(2) + t13;
t14 = t25 * qJDD(1) + t27 * t32;
t15 = qJD(2) * pkin(2) + t27 * qJD(1);
t7 = t24 * t15 + t26 * t35;
t3 = -qJD(3) * t7 + t26 * t10 - t24 * t14;
t30 = Ifges(4,3) * t21 + (-g(2) * t20 + t3) * mrSges(4,1);
t11 = t24 * t27 + t26 * t25;
t12 = -t24 * t25 + t26 * t27;
t6 = t26 * t15 - t24 * t35;
t2 = qJD(3) * t6 + t24 * t10 + t26 * t14;
t29 = g(1) * t20 + g(2) * t19 - t2;
t28 = qJD(2) ^ 2;
t9 = t12 * qJD(1);
t8 = t11 * qJD(1);
t5 = t22 * t11;
t4 = t22 * t12;
t1 = [m(2) * qJDD(1) + (-t11 * t21 - t4 * t22) * mrSges(4,2) + (-t25 * qJDD(2) - t28 * t27) * mrSges(3,2) + (t12 * t21 - t5 * t22) * mrSges(4,1) + (t27 * qJDD(2) - t28 * t25) * mrSges(3,1) + m(3) * (t13 * t27 + t14 * t25) + m(4) * (t2 * t11 + t3 * t12 + t7 * t4 - t6 * t5) + (-m(2) - m(3) - m(4)) * g(2); Ifges(3,3) * qJDD(2) - m(4) * (-t6 * t8 + t7 * t9) + (t8 * t22 + t37) * mrSges(4,1) + (t9 * t22 + t29) * mrSges(4,2) + (g(2) * t25 + t31 * t27 - t14) * mrSges(3,2) + (t31 * t25 + t13 - t36) * mrSges(3,1) + ((-t24 * t21 - t22 * t33) * mrSges(4,2) + (t26 * t21 - t22 * t34) * mrSges(4,1) + (g(1) * t25 + t2 * t24 + t26 * t3 + t7 * t33 - t6 * t34 - t36) * m(4)) * pkin(2) + t30; (t7 * t22 + t37) * mrSges(4,1) + (t6 * t22 + t29) * mrSges(4,2) + t30;];
tau  = t1;
