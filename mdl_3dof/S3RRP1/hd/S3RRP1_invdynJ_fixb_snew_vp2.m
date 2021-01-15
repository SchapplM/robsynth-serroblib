% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S3RRP1
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
%   pkin=[a2,a3,d1,d2]';
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
% Datum: 2021-01-15 13:46
% Revision: 9ee2a23c87189dabbba60bb1a629b0c16b06df77 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S3RRP1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRP1_invdynJ_fixb_snew_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRP1_invdynJ_fixb_snew_vp2: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3RRP1_invdynJ_fixb_snew_vp2: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RRP1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RRP1_invdynJ_fixb_snew_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3RRP1_invdynJ_fixb_snew_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3RRP1_invdynJ_fixb_snew_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3RRP1_invdynJ_fixb_snew_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 13:45:53
% EndTime: 2021-01-15 13:45:53
% DurationCPUTime: 0.14s
% Computational Cost: add. (208->37), mult. (258->43), div. (0->0), fcn. (108->4), ass. (0->19)
t111 = sin(qJ(1));
t113 = cos(qJ(1));
t115 = -t113 * g(1) - t111 * g(2);
t100 = -qJD(1) ^ 2 * pkin(1) + t115;
t110 = sin(qJ(2));
t112 = cos(qJ(2));
t116 = t111 * g(1) - t113 * g(2);
t99 = qJDD(1) * pkin(1) + t116;
t96 = t112 * t100 + t110 * t99;
t108 = qJDD(1) + qJDD(2);
t109 = qJD(1) + qJD(2);
t107 = t109 ^ 2;
t92 = -t107 * pkin(2) + t108 * qJ(3) + 0.2e1 * qJD(3) * t109 + t96;
t117 = m(4) * t92 + t108 * mrSges(4,3);
t95 = -t110 * t100 + t112 * t99;
t93 = -t108 * pkin(2) - t107 * qJ(3) + qJDD(3) - t95;
t89 = m(4) * t93 - t108 * mrSges(4,1) - t107 * mrSges(4,3);
t114 = -mrSges(4,1) * t93 - mrSges(3,2) * t96 + qJ(3) * (-t107 * mrSges(4,1) + t117) - pkin(2) * t89 + mrSges(4,3) * t92 + mrSges(3,1) * t95 + (Ifges(3,3) + Ifges(4,2)) * t108;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t116 - mrSges(2,2) * t115 + pkin(1) * (t110 * (m(3) * t96 - t108 * mrSges(3,2) + (-mrSges(3,1) - mrSges(4,1)) * t107 + t117) + t112 * (m(3) * t95 + t108 * mrSges(3,1) - t107 * mrSges(3,2) - t89)) + t114; t114; t89;];
tauJ = t1;
