% Calculate vector of inverse dynamics joint torques for
% S4PPPR2
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d4,theta2]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PPPR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR2_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR2_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPPR2_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPPR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR2_invdynJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPPR2_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPPR2_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPPR2_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:10:07
% EndTime: 2019-03-08 18:10:07
% DurationCPUTime: 0.12s
% Computational Cost: add. (107->45), mult. (231->63), div. (0->0), fcn. (154->4), ass. (0->23)
t14 = cos(pkin(5));
t10 = -t14 * qJD(1) + qJD(3);
t15 = sin(qJ(4));
t16 = cos(qJ(4));
t13 = sin(pkin(5));
t21 = qJD(1) * t13;
t3 = t16 * t10 - t15 * t21;
t23 = t3 * qJD(4);
t4 = t15 * t10 + t16 * t21;
t22 = t4 * qJD(4);
t20 = qJDD(1) * t13;
t19 = m(3) + m(4) + m(5);
t18 = -g(1) * t13 + g(2) * t14;
t8 = -t13 * t16 + t14 * t15;
t7 = -t13 * t15 - t14 * t16;
t17 = qJD(4) ^ 2;
t12 = t13 ^ 2 * qJDD(1);
t9 = -t14 * qJDD(1) + qJDD(3);
t6 = t7 * qJD(4);
t5 = t8 * qJD(4);
t2 = -t15 * t20 + t16 * t9 - t22;
t1 = t15 * t9 + t16 * t20 + t23;
t11 = [m(2) * qJDD(1) + (-t6 * qJD(4) + t8 * qJDD(4)) * mrSges(5,2) + (t5 * qJD(4) + t7 * qJDD(4)) * mrSges(5,1) + m(3) * (t14 ^ 2 * qJDD(1) + t12) + m(4) * (-t9 * t14 + t12) + m(5) * (-t1 * t8 + t2 * t7 + t3 * t5 + t4 * t6) + (-m(2) - t19) * g(2); (-g(3) + qJDD(2)) * t19; (-t15 * qJDD(4) - t17 * t16) * mrSges(5,2) + (t16 * qJDD(4) - t17 * t15) * mrSges(5,1) + (t1 * t15 + t2 * t16 + (-t15 * t3 + t16 * t4) * qJD(4) + t18) * m(5) + (t18 + t9) * m(4); Ifges(5,3) * qJDD(4) + (-g(1) * t7 - g(2) * t8 - t1 + t23) * mrSges(5,2) + (g(1) * t8 - g(2) * t7 + t2 + t22) * mrSges(5,1);];
tau  = t11;
