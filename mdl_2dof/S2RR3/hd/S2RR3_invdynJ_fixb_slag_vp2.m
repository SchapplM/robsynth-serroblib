% Calculate vector of inverse dynamics joint torques for
% S2RR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% qJDD [2x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,d1,d2]';
% m [3x1]
%   mass of all robot links (including the base)
% mrSges [3x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [3x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [2x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:14
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S2RR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR3_invdynJ_fixb_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR3_invdynJ_fixb_slag_vp2: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'S2RR3_invdynJ_fixb_slag_vp2: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S2RR3_invdynJ_fixb_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'S2RR3_invdynJ_fixb_slag_vp2: m has to be [3x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [3,3]), ...
  'S2RR3_invdynJ_fixb_slag_vp2: mrSges has to be [3x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [3 6]), ...
  'S2RR3_invdynJ_fixb_slag_vp2: Ifges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 09:14:25
% EndTime: 2020-06-19 09:14:25
% DurationCPUTime: 0.21s
% Computational Cost: add. (49->31), mult. (89->43), div. (0->0), fcn. (34->6), ass. (0->17)
t11 = sin(qJ(2));
t13 = cos(qJ(2));
t16 = qJD(1) * qJD(2);
t2 = (qJDD(1) * t13 - t11 * t16) * pkin(1);
t8 = qJDD(1) + qJDD(2);
t18 = t2 * mrSges(3,1) + Ifges(3,3) * t8;
t9 = qJD(1) + qJD(2);
t17 = qJD(2) * t9;
t15 = pkin(1) * qJD(1) * t9;
t14 = cos(qJ(1));
t12 = sin(qJ(1));
t10 = qJ(1) + qJ(2);
t7 = cos(t10);
t6 = sin(t10);
t4 = t7 * mrSges(3,1);
t3 = (qJDD(1) * t11 + t13 * t16) * pkin(1);
t1 = [-t3 * mrSges(3,2) + Ifges(2,3) * qJDD(1) + (-mrSges(2,1) * t14 + t12 * mrSges(2,2) + t6 * mrSges(3,2) - t4) * g(2) + (t12 * mrSges(2,1) + mrSges(3,1) * t6 + mrSges(2,2) * t14 + mrSges(3,2) * t7) * g(1) + ((-t11 * t8 - t13 * t17) * mrSges(3,2) + (-t11 * t17 + t13 * t8) * mrSges(3,1) + (g(1) * t12 - g(2) * t14 + t11 * t3 + t13 * t2) * m(3)) * pkin(1) + t18; -g(2) * t4 + (g(1) * t6 + t11 * t15) * mrSges(3,1) + (g(1) * t7 + g(2) * t6 + t13 * t15 - t3) * mrSges(3,2) + t18;];
tau = t1;
