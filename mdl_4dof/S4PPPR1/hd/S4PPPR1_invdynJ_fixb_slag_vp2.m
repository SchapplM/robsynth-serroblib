% Calculate vector of inverse dynamics joint torques for
% S4PPPR1
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
%   pkin=[a2,a3,a4,d4,theta1]';
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
% Datum: 2019-03-08 18:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PPPR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR1_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR1_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPPR1_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPPR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR1_invdynJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPPR1_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPPR1_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPPR1_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:08:53
% EndTime: 2019-03-08 18:08:54
% DurationCPUTime: 0.12s
% Computational Cost: add. (76->40), mult. (158->49), div. (0->0), fcn. (86->4), ass. (0->19)
t20 = m(3) + m(4);
t10 = cos(qJ(4));
t9 = sin(qJ(4));
t5 = -t9 * qJD(2) + t10 * qJD(3);
t19 = t5 * qJD(4);
t6 = t10 * qJD(2) + t9 * qJD(3);
t18 = t6 * qJD(4);
t7 = sin(pkin(5));
t8 = cos(pkin(5));
t16 = -g(1) * t8 - g(2) * t7;
t15 = -g(1) * t7 + g(2) * t8;
t11 = qJD(4) ^ 2;
t14 = -t9 * qJDD(4) - t11 * t10;
t13 = t10 * qJDD(4) - t11 * t9;
t4 = t7 * t10 + t8 * t9;
t3 = t8 * t10 - t7 * t9;
t2 = -t9 * qJDD(2) + t10 * qJDD(3) - t18;
t1 = t10 * qJDD(2) + t9 * qJDD(3) + t19;
t12 = [(-g(3) + qJDD(1)) * (m(2) + m(5) + t20); -t13 * mrSges(5,2) + t14 * mrSges(5,1) + (t1 * t10 - t2 * t9 + (-t10 * t5 - t6 * t9) * qJD(4) + t15) * m(5) + t20 * (qJDD(2) + t15); t14 * mrSges(5,2) + t13 * mrSges(5,1) + (t1 * t9 + t2 * t10 + (t10 * t6 - t5 * t9) * qJD(4) + t16) * m(5) + (qJDD(3) + t16) * m(4); Ifges(5,3) * qJDD(4) + (g(1) * t4 - g(2) * t3 - t1 + t19) * mrSges(5,2) + (-g(1) * t3 - g(2) * t4 + t18 + t2) * mrSges(5,1);];
tau  = t12;
