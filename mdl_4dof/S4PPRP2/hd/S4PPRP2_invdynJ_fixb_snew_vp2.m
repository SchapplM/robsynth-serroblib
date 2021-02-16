% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4PPRP2
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
%   pkin=[a2,a3,a4,d3,theta2]';
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
% Datum: 2021-01-15 15:05
% Revision: 9ee2a23c87189dabbba60bb1a629b0c16b06df77 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4PPRP2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP2_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP2_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRP2_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRP2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP2_invdynJ_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP2_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRP2_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRP2_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:04:50
% EndTime: 2021-01-15 15:04:50
% DurationCPUTime: 0.11s
% Computational Cost: add. (151->38), mult. (193->44), div. (0->0), fcn. (100->4), ass. (0->17)
t123 = -g(2) + qJDD(1);
t125 = sin(pkin(5));
t126 = cos(pkin(5));
t118 = t125 * g(1) + t126 * t123;
t119 = -t126 * g(1) + t125 * t123;
t127 = sin(qJ(3));
t128 = cos(qJ(3));
t115 = t127 * t118 + t128 * t119;
t129 = qJD(3) ^ 2;
t112 = -t129 * pkin(3) + qJDD(3) * qJ(4) + (2 * qJD(4) * qJD(3)) + t115;
t130 = m(5) * t112 + qJDD(3) * mrSges(5,3);
t114 = t128 * t118 - t127 * t119;
t113 = -qJDD(3) * pkin(3) - t129 * qJ(4) + qJDD(4) - t114;
t110 = m(5) * t113 - qJDD(3) * mrSges(5,1) - t129 * mrSges(5,3);
t109 = m(4) * t114 + qJDD(3) * mrSges(4,1) - t129 * mrSges(4,2) - t110;
t108 = m(4) * t115 - qJDD(3) * mrSges(4,2) + (-mrSges(4,1) - mrSges(5,1)) * t129 + t130;
t1 = [m(2) * t123 + t125 * (m(3) * t119 + t128 * t108 - t127 * t109) + t126 * (m(3) * t118 + t127 * t108 + t128 * t109); (m(3) + m(4) + m(5)) * (-g(3) + qJDD(2)); mrSges(4,1) * t114 - mrSges(4,2) * t115 - mrSges(5,1) * t113 + mrSges(5,3) * t112 - pkin(3) * t110 + qJ(4) * (-t129 * mrSges(5,1) + t130) + (Ifges(4,3) + Ifges(5,2)) * qJDD(3); t110;];
tauJ = t1;
