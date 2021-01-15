% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4PRPP3
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
%   pkin=[a2,a3,a4,d2]';
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
% Datum: 2021-01-15 16:30
% Revision: 9ee2a23c87189dabbba60bb1a629b0c16b06df77 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4PRPP3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(4,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP3_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP3_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPP3_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PRPP3_invdynJ_fixb_snew_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP3_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPP3_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPP3_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:30:00
% EndTime: 2021-01-15 16:30:00
% DurationCPUTime: 0.14s
% Computational Cost: add. (139->43), mult. (174->44), div. (0->0), fcn. (40->2), ass. (0->18)
t149 = -pkin(2) - pkin(3);
t148 = -mrSges(4,1) - mrSges(5,1);
t138 = -g(2) + qJDD(1);
t141 = sin(qJ(2));
t142 = cos(qJ(2));
t130 = -t142 * g(1) + t141 * t138;
t143 = qJD(2) ^ 2;
t145 = qJDD(2) * qJ(3) + (2 * qJD(3) * qJD(2)) + t130;
t124 = t149 * t143 + t145;
t127 = -t143 * pkin(2) + t145;
t147 = m(4) * t127 + m(5) * t124 + (mrSges(5,2) + mrSges(4,3)) * qJDD(2);
t129 = t141 * g(1) + t142 * t138;
t144 = -t143 * qJ(3) + qJDD(3) - t129;
t126 = t149 * qJDD(2) + t144;
t146 = m(5) * t126 - qJDD(2) * mrSges(5,1) - t143 * mrSges(5,2);
t128 = -qJDD(2) * pkin(2) + t144;
t122 = m(4) * t128 - qJDD(2) * mrSges(4,1) - t143 * mrSges(4,3) + t146;
t1 = [m(2) * t138 + t141 * (m(3) * t130 - qJDD(2) * mrSges(3,2) + t147) + t142 * (m(3) * t129 + qJDD(2) * mrSges(3,1) - t122) + (t141 * (-mrSges(3,1) + t148) - t142 * mrSges(3,2)) * t143; mrSges(3,1) * t129 - mrSges(3,2) * t130 - mrSges(4,1) * t128 + mrSges(4,3) * t127 - mrSges(5,1) * t126 + mrSges(5,2) * t124 - pkin(3) * t146 - pkin(2) * t122 + qJ(3) * (t148 * t143 + t147) + (Ifges(3,3) + Ifges(4,2) + Ifges(5,3)) * qJDD(2); t122; m(5) * (g(3) + qJDD(4));];
tauJ = t1;
