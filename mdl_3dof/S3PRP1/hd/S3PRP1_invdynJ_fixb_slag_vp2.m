% Calculate vector of inverse dynamics joint torques for
% S3PRP1
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
%   pkin=[a2,a3,d2]';
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
% Datum: 2018-11-14 10:04
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau = S3PRP1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRP1_invdynJ_fixb_slag_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRP1_invdynJ_fixb_slag_vp2: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3PRP1_invdynJ_fixb_slag_vp2: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3PRP1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PRP1_invdynJ_fixb_slag_vp2: pkin has to be [3x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3PRP1_invdynJ_fixb_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3PRP1_invdynJ_fixb_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3PRP1_invdynJ_fixb_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:04:14
% EndTime: 2018-11-14 10:04:15
% DurationCPUTime: 0.16s
% Computational Cost: add. (84->50), mult. (170->60), div. (0->0), fcn. (62->2), ass. (0->17)
t10 = cos(qJ(2));
t17 = mrSges(3,2) - mrSges(4,3);
t18 = mrSges(3,1) + mrSges(4,1);
t9 = sin(qJ(2));
t19 = t17 * t10 + t18 * t9;
t16 = qJD(1) * qJD(2);
t4 = t9 * qJDD(1) + t10 * t16;
t15 = qJD(2) * qJD(3);
t14 = -m(4) * qJ(3) + mrSges(3,2);
t13 = m(4) * pkin(2) + t18;
t6 = qJD(2) * qJ(3) + t9 * qJD(1);
t12 = t10 * t6 + (-qJD(2) * pkin(2) - t10 * qJD(1) + qJD(3)) * t9;
t3 = t10 * qJDD(1) - t9 * t16;
t11 = qJD(2) ^ 2;
t2 = -qJDD(2) * pkin(2) + qJDD(3) - t3;
t1 = qJDD(2) * qJ(3) + t15 + t4;
t5 = [m(2) * qJDD(1) + m(3) * (t3 * t10 + t4 * t9) + m(4) * (t12 * qJD(2) + t1 * t9 - t2 * t10) + (-m(2) - m(3) - m(4)) * g(2) - t19 * t11 + (t18 * t10 - t17 * t9) * qJDD(2); -t4 * mrSges(3,2) + t3 * mrSges(3,1) - t2 * mrSges(4,1) + m(4) * (-t2 * pkin(2) + t1 * qJ(3) + t6 * qJD(3)) + (-t13 * t10 + t14 * t9) * g(2) + (t14 * t10 + t13 * t9) * g(1) + (-g(1) * t10 - g(2) * t9 + t1 + t15) * mrSges(4,3) + (pkin(2) * mrSges(4,1) + qJ(3) * mrSges(4,3) + Ifges(4,2) + Ifges(3,3)) * qJDD(2) + (-m(4) * t12 + t19 * qJD(2)) * qJD(1); -qJDD(2) * mrSges(4,1) - t11 * mrSges(4,3) + (-g(1) * t9 + g(2) * t10 - t6 * qJD(2) + t2) * m(4);];
tau  = t5;
