% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S3RRR1
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% m [4x1]
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
% tauJ [3x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 13:58
% Revision: 9ee2a23c87189dabbba60bb1a629b0c16b06df77 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S3RRR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(5,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_invdynJ_fixb_snew_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRR1_invdynJ_fixb_snew_vp2: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3RRR1_invdynJ_fixb_snew_vp2: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RRR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_invdynJ_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3RRR1_invdynJ_fixb_snew_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3RRR1_invdynJ_fixb_snew_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3RRR1_invdynJ_fixb_snew_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 13:58:05
% EndTime: 2021-01-15 13:58:06
% DurationCPUTime: 0.20s
% Computational Cost: add. (331->37), mult. (451->48), div. (0->0), fcn. (246->6), ass. (0->28)
t121 = sin(qJ(3));
t124 = cos(qJ(3));
t123 = sin(qJ(1));
t126 = cos(qJ(1));
t130 = t123 * g(1) - t126 * g(2);
t110 = qJDD(1) * pkin(1) + t130;
t128 = -t126 * g(1) - t123 * g(2);
t111 = -qJD(1) ^ 2 * pkin(1) + t128;
t122 = sin(qJ(2));
t125 = cos(qJ(2));
t105 = t125 * t110 - t122 * t111;
t119 = qJDD(1) + qJDD(2);
t102 = t119 * pkin(2) + t105;
t106 = t122 * t110 + t125 * t111;
t120 = qJD(1) + qJD(2);
t118 = t120 ^ 2;
t103 = -t118 * pkin(2) + t106;
t100 = t124 * t102 - t121 * t103;
t116 = qJD(3) + t120;
t114 = t116 ^ 2;
t115 = qJDD(3) + t119;
t97 = m(4) * t100 + t115 * mrSges(4,1) - t114 * mrSges(4,2);
t101 = t121 * t102 + t124 * t103;
t98 = m(4) * t101 - t114 * mrSges(4,1) - t115 * mrSges(4,2);
t131 = t121 * t98 + t124 * t97;
t129 = mrSges(4,1) * t100 - mrSges(4,2) * t101 + Ifges(4,3) * t115;
t127 = mrSges(3,1) * t105 - mrSges(3,2) * t106 + Ifges(3,3) * t119 + pkin(2) * t131 + t129;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t130 - mrSges(2,2) * t128 + pkin(1) * (t122 * (m(3) * t106 - t118 * mrSges(3,1) - t119 * mrSges(3,2) - t121 * t97 + t124 * t98) + t125 * (m(3) * t105 + t119 * mrSges(3,1) - t118 * mrSges(3,2) + t131)) + t127; t127; t129;];
tauJ = t1;
