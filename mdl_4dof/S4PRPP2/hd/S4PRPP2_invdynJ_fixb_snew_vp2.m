% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4PRPP2
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
%   pkin=[a2,a3,a4,d2,theta3]';
% m [5x1]
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
% Datum: 2021-01-15 16:10
% Revision: 9ee2a23c87189dabbba60bb1a629b0c16b06df77 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4PRPP2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP2_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP2_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPP2_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP2_invdynJ_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP2_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPP2_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPP2_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:10:18
% EndTime: 2021-01-15 16:10:18
% DurationCPUTime: 0.22s
% Computational Cost: add. (236->47), mult. (313->53), div. (0->0), fcn. (130->4), ass. (0->20)
t138 = -g(2) + qJDD(1);
t142 = sin(qJ(2));
t143 = cos(qJ(2));
t130 = t142 * g(1) + t143 * t138;
t128 = qJDD(2) * pkin(2) + t130;
t131 = -t143 * g(1) + t142 * t138;
t144 = qJD(2) ^ 2;
t129 = -t144 * pkin(2) + t131;
t140 = sin(pkin(5));
t141 = cos(pkin(5));
t125 = t140 * t128 + t141 * t129;
t122 = -t144 * pkin(3) + qJDD(2) * qJ(4) + (2 * qJD(4) * qJD(2)) + t125;
t145 = m(5) * t122 + qJDD(2) * mrSges(5,3);
t118 = m(4) * t125 - qJDD(2) * mrSges(4,2) + (-mrSges(4,1) - mrSges(5,1)) * t144 + t145;
t124 = t141 * t128 - t140 * t129;
t123 = -qJDD(2) * pkin(3) - t144 * qJ(4) + qJDD(4) - t124;
t120 = m(5) * t123 - qJDD(2) * mrSges(5,1) - t144 * mrSges(5,3);
t119 = m(4) * t124 + qJDD(2) * mrSges(4,1) - t144 * mrSges(4,2) - t120;
t146 = t140 * t118 + t141 * t119;
t1 = [m(2) * t138 + t142 * (m(3) * t131 - t144 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t141 * t118 - t140 * t119) + t143 * (m(3) * t130 + qJDD(2) * mrSges(3,1) - t144 * mrSges(3,2) + t146); mrSges(3,1) * t130 - mrSges(3,2) * t131 + mrSges(4,1) * t124 - mrSges(4,2) * t125 - mrSges(5,1) * t123 + mrSges(5,3) * t122 - pkin(3) * t120 + qJ(4) * (-t144 * mrSges(5,1) + t145) + pkin(2) * t146 + (Ifges(3,3) + Ifges(4,3) + Ifges(5,2)) * qJDD(2); (m(4) + m(5)) * (-g(3) + qJDD(3)); t120;];
tauJ = t1;
