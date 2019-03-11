% Calculate vector of inverse dynamics joint torques for
% S3RRP1
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
%   pkin=[a2,a3,d1,d2]';
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
% Datum: 2019-03-08 18:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S3RRP1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRP1_invdynJ_fixb_slag_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRP1_invdynJ_fixb_slag_vp2: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3RRP1_invdynJ_fixb_slag_vp2: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RRP1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RRP1_invdynJ_fixb_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3RRP1_invdynJ_fixb_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3RRP1_invdynJ_fixb_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3RRP1_invdynJ_fixb_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:06:51
% EndTime: 2019-03-08 18:06:51
% DurationCPUTime: 0.29s
% Computational Cost: add. (188->68), mult. (274->84), div. (0->0), fcn. (90->6), ass. (0->36)
t42 = mrSges(3,1) + mrSges(4,1);
t50 = m(4) * pkin(2) + t42;
t48 = mrSges(3,2) - mrSges(4,3);
t30 = cos(qJ(2));
t28 = sin(qJ(2));
t39 = qJD(2) * t28;
t5 = (qJD(1) * t39 - qJDD(1) * t30) * pkin(1);
t47 = t42 * t28;
t46 = mrSges(2,1) + (m(3) + m(4)) * pkin(1);
t25 = qJDD(1) + qJDD(2);
t44 = t25 * pkin(2);
t43 = t28 * pkin(1);
t41 = pkin(1) * qJD(1);
t38 = t30 * t41;
t6 = qJD(2) * t38 + qJDD(1) * t43;
t40 = pkin(1) * qJD(2);
t27 = qJ(1) + qJ(2);
t23 = sin(t27);
t3 = qJDD(3) - t44 + t5;
t37 = g(1) * t23 - t3;
t24 = cos(t27);
t36 = t48 * t24;
t26 = qJD(1) + qJD(2);
t35 = t25 * qJ(3) + t26 * qJD(3);
t2 = t35 + t6;
t34 = -t5 * mrSges(3,1) - t6 * mrSges(3,2) + t2 * mrSges(4,3) + (Ifges(4,2) + Ifges(3,3)) * t25;
t32 = -t50 * t24 + (-m(4) * qJ(3) + t48) * t23;
t31 = cos(qJ(1));
t29 = sin(qJ(1));
t19 = -t30 * pkin(1) - pkin(2);
t17 = qJ(3) + t43;
t12 = t24 * qJ(3);
t9 = t30 * t40 + qJD(3);
t8 = t26 * qJ(3) + t28 * t41;
t7 = -t26 * pkin(2) + qJD(3) - t38;
t1 = [-t3 * mrSges(4,1) + Ifges(2,3) * qJDD(1) + m(3) * (t28 * t6 - t30 * t5) * pkin(1) + m(4) * (t7 * pkin(1) * t39 + t2 * t17 + t3 * t19 + t8 * t9) + (t29 * mrSges(2,2) - t46 * t31 + t32) * g(2) + (-m(4) * t12 + t31 * mrSges(2,2) + t50 * t23 + t46 * t29 + t36) * g(1) + (t9 * mrSges(4,3) + (-t30 * mrSges(3,2) - t47) * t40) * t26 + (-t19 * mrSges(4,1) + t17 * mrSges(4,3) + (mrSges(3,1) * t30 - mrSges(3,2) * t28) * pkin(1)) * t25 + t34; m(4) * (-t3 * pkin(2) + t2 * qJ(3) + t8 * qJD(3)) + t32 * g(2) + (t23 * mrSges(3,1) - m(4) * (-t23 * pkin(2) + t12) + t36) * g(1) + t35 * mrSges(4,3) + (t37 + t44) * mrSges(4,1) + (-m(4) * (t28 * t7 + t30 * t8) + (t48 * t30 + t47) * t26) * t41 + t34; -t26 ^ 2 * mrSges(4,3) - t25 * mrSges(4,1) + (g(2) * t24 - t8 * t26 - t37) * m(4);];
tau  = t1;
