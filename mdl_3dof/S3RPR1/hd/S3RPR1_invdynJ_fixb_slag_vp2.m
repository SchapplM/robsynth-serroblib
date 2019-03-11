% Calculate vector of inverse dynamics joint torques for
% S3RPR1
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
%   pkin=[a2,a3,d1,d3]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S3RPR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPR1_invdynJ_fixb_slag_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPR1_invdynJ_fixb_slag_vp2: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3RPR1_invdynJ_fixb_slag_vp2: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RPR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RPR1_invdynJ_fixb_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3RPR1_invdynJ_fixb_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3RPR1_invdynJ_fixb_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3RPR1_invdynJ_fixb_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:05:51
% EndTime: 2019-03-08 18:05:52
% DurationCPUTime: 0.30s
% Computational Cost: add. (218->75), mult. (360->94), div. (0->0), fcn. (142->4), ass. (0->33)
t34 = qJD(1) - qJD(3);
t27 = -pkin(1) - pkin(2);
t21 = -qJDD(1) + qJDD(3);
t39 = Ifges(4,3) * t21;
t38 = -mrSges(2,1) - mrSges(3,1);
t24 = sin(qJ(1));
t26 = cos(qJ(1));
t37 = t26 * pkin(1) + t24 * qJ(2);
t36 = qJ(2) * qJD(1);
t35 = qJD(1) * qJD(2);
t13 = qJDD(1) * qJ(2) + t35;
t33 = t13 + t35;
t32 = -g(1) * t24 + g(2) * t26;
t23 = sin(qJ(3));
t25 = cos(qJ(3));
t10 = t25 * qJ(2) + t23 * t27;
t9 = -t23 * qJ(2) + t25 * t27;
t11 = t27 * qJDD(1) + qJDD(2);
t12 = t27 * qJD(1) + qJD(2);
t5 = t25 * t12 - t23 * t36;
t1 = qJD(3) * t5 + t23 * t11 + t25 * t13;
t7 = -t24 * t23 - t26 * t25;
t8 = t26 * t23 - t24 * t25;
t30 = -g(1) * t7 - g(2) * t8 - t1;
t6 = t23 * t12 + t25 * t36;
t2 = -qJD(3) * t6 + t25 * t11 - t23 * t13;
t29 = g(1) * t8 - g(2) * t7 + t2;
t28 = qJD(1) ^ 2;
t19 = t26 * qJ(2);
t17 = -qJDD(1) * pkin(1) + qJDD(2);
t4 = -t23 * qJD(2) - qJD(3) * t10;
t3 = t25 * qJD(2) + qJD(3) * t9;
t14 = [-t39 - t17 * mrSges(3,1) + m(3) * (-t17 * pkin(1) + qJ(2) * t33) + (-m(3) * t37 + t24 * mrSges(2,2) + t26 * t38) * g(2) + (-m(3) * t19 + t26 * mrSges(2,2) + (m(3) * pkin(1) - t38) * t24) * g(1) + (-g(1) * t26 - g(2) * t24 + t33) * mrSges(3,3) + (-g(1) * (t27 * t24 + t19) - g(2) * (t26 * pkin(2) + t37) + t1 * t10 + t2 * t9 + t6 * t3 + t5 * t4) * m(4) + (-t10 * t21 + t3 * t34 - t30) * mrSges(4,2) + (t9 * t21 - t34 * t4 - t29) * mrSges(4,1) + (pkin(1) * mrSges(3,1) + qJ(2) * mrSges(3,3) + Ifges(3,2) + Ifges(2,3)) * qJDD(1); -qJDD(1) * mrSges(3,1) - t28 * mrSges(3,3) + (t25 * mrSges(4,1) - t23 * mrSges(4,2)) * t21 + (t1 * t23 + t2 * t25 + t32 - t34 * (-t23 * t5 + t25 * t6)) * m(4) + (-t28 * qJ(2) + t17 + t32) * m(3) - (mrSges(4,1) * t23 + mrSges(4,2) * t25) * t34 ^ 2; t39 + (-t34 * t5 + t30) * mrSges(4,2) + (-t34 * t6 + t29) * mrSges(4,1);];
tau  = t14;
