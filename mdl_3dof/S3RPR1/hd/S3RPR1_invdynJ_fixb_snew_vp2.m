% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S3RPR1
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
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
% Datum: 2021-01-15 13:33
% Revision: 9ee2a23c87189dabbba60bb1a629b0c16b06df77 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S3RPR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPR1_invdynJ_fixb_snew_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPR1_invdynJ_fixb_snew_vp2: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3RPR1_invdynJ_fixb_snew_vp2: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RPR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RPR1_invdynJ_fixb_snew_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3RPR1_invdynJ_fixb_snew_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3RPR1_invdynJ_fixb_snew_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3RPR1_invdynJ_fixb_snew_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 13:32:57
% EndTime: 2021-01-15 13:32:57
% DurationCPUTime: 0.13s
% Computational Cost: add. (214->38), mult. (304->45), div. (0->0), fcn. (94->4), ass. (0->25)
t116 = -pkin(1) - pkin(2);
t106 = sin(qJ(1));
t108 = cos(qJ(1));
t115 = t106 * g(1) - t108 * g(2);
t114 = -t108 * g(1) - t106 * g(2);
t105 = sin(qJ(3));
t107 = cos(qJ(3));
t103 = -qJD(1) + qJD(3);
t101 = t103 ^ 2;
t102 = -qJDD(1) + qJDD(3);
t109 = qJD(1) ^ 2;
t112 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t114;
t95 = t116 * t109 + t112;
t110 = -t109 * qJ(2) + qJDD(2) - t115;
t96 = t116 * qJDD(1) + t110;
t93 = -t105 * t95 + t107 * t96;
t91 = m(4) * t93 + t102 * mrSges(4,1) - t101 * mrSges(4,2);
t94 = t105 * t96 + t107 * t95;
t92 = m(4) * t94 - t101 * mrSges(4,1) - t102 * mrSges(4,2);
t113 = t105 * t92 + t107 * t91;
t111 = mrSges(4,1) * t93 - mrSges(4,2) * t94 + Ifges(4,3) * t102;
t98 = -qJDD(1) * pkin(1) + t110;
t97 = -t109 * pkin(1) + t112;
t90 = m(3) * t98 - qJDD(1) * mrSges(3,1) - t109 * mrSges(3,3) + t113;
t1 = [mrSges(2,1) * t115 - mrSges(2,2) * t114 - mrSges(3,1) * t98 + mrSges(3,3) * t97 - pkin(2) * t113 - pkin(1) * t90 + qJ(2) * (m(3) * t97 - t109 * mrSges(3,1) - t105 * t91 + t107 * t92) + (qJ(2) * mrSges(3,3) + Ifges(3,2) + Ifges(2,3)) * qJDD(1) - t111; t90; t111;];
tauJ = t1;
