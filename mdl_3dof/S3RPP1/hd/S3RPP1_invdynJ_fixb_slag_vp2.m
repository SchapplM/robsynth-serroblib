% Calculate vector of inverse dynamics joint torques for
% S3RPP1
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
%   pkin=[a2,a3,d1]';
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
% Datum: 2018-11-14 10:14
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau = S3RPP1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPP1_invdynJ_fixb_slag_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPP1_invdynJ_fixb_slag_vp2: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3RPP1_invdynJ_fixb_slag_vp2: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RPP1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3RPP1_invdynJ_fixb_slag_vp2: pkin has to be [3x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3RPP1_invdynJ_fixb_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3RPP1_invdynJ_fixb_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3RPP1_invdynJ_fixb_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:13:38
% EndTime: 2018-11-14 10:13:38
% DurationCPUTime: 0.24s
% Computational Cost: add. (86->54), mult. (143->49), div. (0->0), fcn. (28->2), ass. (0->20)
t24 = -m(3) - m(4);
t23 = mrSges(3,2) - mrSges(4,3);
t17 = qJD(1) * qJD(2);
t4 = -qJDD(1) * qJ(2) - t17;
t22 = -mrSges(2,1) + t23;
t11 = sin(qJ(1));
t12 = cos(qJ(1));
t21 = -g(1) * t12 - g(2) * t11;
t20 = t21 - t4;
t10 = -pkin(1) - qJ(3);
t18 = qJDD(1) * pkin(1);
t16 = qJD(1) * qJD(3);
t14 = -g(1) * t11 + g(2) * t12;
t13 = qJD(1) ^ 2;
t6 = qJDD(2) - t18;
t5 = qJD(1) * qJ(2) + qJD(3);
t3 = t10 * qJD(1) + qJD(2);
t2 = qJDD(3) - t4;
t1 = t10 * qJDD(1) + qJDD(2) - t16;
t7 = [m(3) * (-pkin(1) * t6 + (-t4 + t17) * qJ(2)) + m(4) * (t2 * qJ(2) + t5 * qJD(2) - t3 * qJD(3) + t1 * t10) + (-t1 + t16) * mrSges(4,3) + (-t18 + t6) * mrSges(3,2) + (t11 * mrSges(2,2) + t24 * (t12 * pkin(1) + t11 * qJ(2)) + (-m(4) * qJ(3) + t22) * t12) * g(2) + ((m(3) * pkin(1) - m(4) * t10 - t22) * t11 + (t24 * qJ(2) + mrSges(2,2)) * t12) * g(1) + (t20 - t4) * mrSges(3,3) + (t2 + t20) * mrSges(4,2) + (-t10 * mrSges(4,3) + Ifges(3,1) + Ifges(4,1) + Ifges(2,3)) * qJDD(1); (-mrSges(3,3) - mrSges(4,2)) * t13 + t23 * qJDD(1) + (-t5 * qJD(1) + t1 + t14) * m(4) + (-t13 * qJ(2) + t14 + t6) * m(3); qJDD(1) * mrSges(4,2) - t13 * mrSges(4,3) + (t3 * qJD(1) + t2 + t21) * m(4);];
tau  = t7;
