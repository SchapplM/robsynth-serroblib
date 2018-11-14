% Calculate vector of inverse dynamics joint torques for
% S4PPPR4
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
%   pkin=[a2,a3,a4,d4,theta3]';
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
% Datum: 2018-11-14 14:04
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau = S4PPPR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR4_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR4_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPPR4_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPPR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR4_invdynJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPPR4_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPPR4_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPPR4_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:04:11
% EndTime: 2018-11-14 14:04:11
% DurationCPUTime: 0.13s
% Computational Cost: add. (137->46), mult. (281->60), div. (0->0), fcn. (212->6), ass. (0->26)
t31 = m(4) + m(5);
t19 = sin(pkin(5));
t20 = cos(pkin(5));
t21 = sin(qJ(4));
t22 = cos(qJ(4));
t8 = -t21 * t19 + t22 * t20;
t11 = -t19 * qJD(1) + t20 * qJD(2);
t12 = t20 * qJD(1) + t19 * qJD(2);
t3 = t22 * t11 - t21 * t12;
t28 = t3 * qJD(4);
t4 = t21 * t11 + t22 * t12;
t27 = t4 * qJD(4);
t25 = m(3) + t31;
t7 = -t22 * t19 - t21 * t20;
t5 = t8 * qJD(4);
t24 = -t5 * qJD(4) + t7 * qJDD(4);
t6 = t7 * qJD(4);
t23 = t6 * qJD(4) + t8 * qJDD(4);
t18 = pkin(5) + qJ(4);
t17 = cos(t18);
t16 = sin(t18);
t10 = t20 * qJDD(1) + t19 * qJDD(2);
t9 = -t19 * qJDD(1) + t20 * qJDD(2);
t2 = -t21 * t10 + t22 * t9 - t27;
t1 = t22 * t10 + t21 * t9 + t28;
t13 = [(m(2) + m(3)) * qJDD(1) - t23 * mrSges(5,2) + t24 * mrSges(5,1) + m(4) * (t10 * t20 - t9 * t19) + m(5) * (t1 * t8 + t2 * t7 - t3 * t5 + t4 * t6) + (-m(2) - t25) * g(1); m(3) * qJDD(2) + t24 * mrSges(5,2) + t23 * mrSges(5,1) + m(4) * (t10 * t19 + t9 * t20) + m(5) * (-t1 * t7 + t2 * t8 + t3 * t6 + t4 * t5) + t25 * g(2); t31 * (g(3) + qJDD(3)); Ifges(5,3) * qJDD(4) + (g(1) * t17 - g(2) * t16 - t1 + t28) * mrSges(5,2) + (g(1) * t16 + g(2) * t17 + t2 + t27) * mrSges(5,1);];
tau  = t13;
