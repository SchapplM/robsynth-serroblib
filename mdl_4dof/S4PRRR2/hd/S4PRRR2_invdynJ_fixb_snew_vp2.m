% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4PRRR2
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
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
% tauJ [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4PRRR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(2,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR2_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR2_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_invdynJ_fixb_snew_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR2_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR2_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR2_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:23
% EndTime: 2019-07-18 13:27:23
% DurationCPUTime: 0.07s
% Computational Cost: add. (338->41), mult. (455->49), div. (0->0), fcn. (246->6), ass. (0->28)
t143 = sin(qJ(2));
t146 = cos(qJ(2));
t150 = t146 * g(1) + t143 * g(3);
t127 = qJDD(2) * pkin(1) + t150;
t149 = t143 * g(1) - t146 * g(3);
t128 = -qJD(2) ^ 2 * pkin(1) + t149;
t142 = sin(qJ(3));
t145 = cos(qJ(3));
t122 = t145 * t127 - t142 * t128;
t138 = qJDD(2) + qJDD(3);
t119 = t138 * pkin(2) + t122;
t123 = t142 * t127 + t145 * t128;
t139 = qJD(2) + qJD(3);
t137 = t139 ^ 2;
t120 = -t137 * pkin(2) + t123;
t141 = sin(qJ(4));
t144 = cos(qJ(4));
t117 = t144 * t119 - t141 * t120;
t133 = qJD(4) + t139;
t131 = t133 ^ 2;
t132 = qJDD(4) + t138;
t114 = m(5) * t117 + t132 * mrSges(5,1) - t131 * mrSges(5,2);
t118 = t141 * t119 + t144 * t120;
t115 = m(5) * t118 - t131 * mrSges(5,1) - t132 * mrSges(5,2);
t151 = t144 * t114 + t141 * t115;
t148 = mrSges(5,1) * t117 - mrSges(5,2) * t118 + Ifges(5,3) * t132;
t147 = mrSges(4,1) * t122 - mrSges(4,2) * t123 + Ifges(4,3) * t138 + pkin(2) * t151 + t148;
t1 = [(m(2) + m(3) + m(4) + m(5)) * (g(2) + qJDD(1)); Ifges(3,3) * qJDD(2) + mrSges(3,1) * t150 - mrSges(3,2) * t149 + pkin(1) * (t142 * (m(4) * t123 - t137 * mrSges(4,1) - t138 * mrSges(4,2) - t141 * t114 + t144 * t115) + t145 * (m(4) * t122 + t138 * mrSges(4,1) - t137 * mrSges(4,2) + t151)) + t147; t147; t148;];
tauJ  = t1;
