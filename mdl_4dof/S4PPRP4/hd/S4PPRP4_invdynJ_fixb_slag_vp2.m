% Calculate vector of inverse dynamics joint torques for
% S4PPRP4
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3]';
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
% Datum: 2018-11-14 14:06
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau = S4PPRP4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(4,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP4_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP4_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRP4_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRP4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PPRP4_invdynJ_fixb_slag_vp2: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP4_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRP4_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRP4_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:06:03
% EndTime: 2018-11-14 14:06:03
% DurationCPUTime: 0.20s
% Computational Cost: add. (132->52), mult. (271->58), div. (0->0), fcn. (134->2), ass. (0->25)
t14 = cos(qJ(3));
t13 = sin(qJ(3));
t23 = mrSges(4,2) + mrSges(5,2);
t18 = t23 * t13;
t24 = mrSges(4,1) + mrSges(5,1);
t31 = -t24 * t14 + t18;
t9 = -t13 * qJD(1) + t14 * qJD(2);
t30 = qJD(3) * t9;
t29 = t23 * t14;
t10 = t14 * qJD(1) + t13 * qJD(2);
t27 = qJD(3) * t10;
t4 = t14 * qJDD(1) + t13 * qJDD(2) + t30;
t26 = t4 * t13 + t14 * t27;
t25 = t10 * t13;
t22 = qJD(3) * t13;
t20 = m(3) + m(4) + m(5);
t19 = m(5) * pkin(3) + mrSges(5,1);
t17 = mrSges(4,1) + t19;
t16 = -t24 * t13 - t29;
t5 = -t13 * qJDD(1) + t14 * qJDD(2) - t27;
t15 = qJD(3) ^ 2;
t8 = qJD(3) * pkin(3) + t9;
t3 = qJDD(3) * pkin(3) + t5;
t2 = t4 * t14;
t1 = [(m(2) + m(3)) * qJDD(1) + m(4) * (-t5 * t13 + t2 + (-t14 * t9 - t25) * qJD(3)) + m(5) * (-t3 * t13 + t2 + (-t14 * t8 - t25) * qJD(3)) + t31 * t15 + t16 * qJDD(3) + (-m(2) - t20) * g(1); m(3) * qJDD(2) + m(4) * (t5 * t14 - t9 * t22 + t26) + m(5) * (t3 * t14 - t8 * t22 + t26) + t20 * g(2) + t16 * t15 - t31 * qJDD(3); t5 * mrSges(4,1) + t19 * t3 + (t17 * t14 - t18) * g(2) + (t17 * t13 + t29) * g(1) + (pkin(3) * mrSges(5,1) + Ifges(4,3) + Ifges(5,3)) * qJDD(3) + (-t4 + t30) * t23 + (-m(5) * (-t8 + t9) + t24 * qJD(3)) * t10; (g(3) + qJDD(4)) * m(5);];
tau  = t1;
