% Calculate vector of inverse dynamics joint torques for
% S3PPR1
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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d3]';
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
% Datum: 2019-03-08 18:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S3PPR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PPR1_invdynJ_fixb_slag_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PPR1_invdynJ_fixb_slag_vp2: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3PPR1_invdynJ_fixb_slag_vp2: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3PPR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PPR1_invdynJ_fixb_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3PPR1_invdynJ_fixb_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3PPR1_invdynJ_fixb_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3PPR1_invdynJ_fixb_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:01:50
% EndTime: 2019-03-08 18:01:50
% DurationCPUTime: 0.10s
% Computational Cost: add. (60->33), mult. (127->41), div. (0->0), fcn. (64->2), ass. (0->13)
t12 = m(2) + m(3);
t5 = sin(qJ(3));
t6 = cos(qJ(3));
t3 = -qJD(1) * t5 + qJD(2) * t6;
t11 = qJD(3) * t3;
t4 = qJD(1) * t6 + qJD(2) * t5;
t10 = qJD(3) * t4;
t7 = qJD(3) ^ 2;
t9 = t6 * qJDD(3) - t7 * t5;
t8 = -t5 * qJDD(3) - t7 * t6;
t2 = -qJDD(1) * t5 + qJDD(2) * t6 - t10;
t1 = qJDD(1) * t6 + qJDD(2) * t5 + t11;
t13 = [m(4) * (t1 * t6 - t2 * t5 + (-t3 * t6 - t4 * t5) * qJD(3)) + t12 * qJDD(1) - t9 * mrSges(4,2) + t8 * mrSges(4,1) + (-m(4) - t12) * g(2); t8 * mrSges(4,2) + t9 * mrSges(4,1) + (t1 * t5 + t2 * t6 + (-t3 * t5 + t4 * t6) * qJD(3) - g(1)) * m(4) + (qJDD(2) - g(1)) * m(3); Ifges(4,3) * qJDD(3) + (g(1) * t5 + g(2) * t6 - t1 + t11) * mrSges(4,2) + (-g(1) * t6 + g(2) * t5 + t10 + t2) * mrSges(4,1);];
tau  = t13;
