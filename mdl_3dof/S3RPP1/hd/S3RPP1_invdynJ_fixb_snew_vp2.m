% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S3RPP1
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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1]';
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
% Datum: 2021-01-15 13:21
% Revision: 9ee2a23c87189dabbba60bb1a629b0c16b06df77 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S3RPP1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPP1_invdynJ_fixb_snew_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPP1_invdynJ_fixb_snew_vp2: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3RPP1_invdynJ_fixb_snew_vp2: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RPP1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3RPP1_invdynJ_fixb_snew_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3RPP1_invdynJ_fixb_snew_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3RPP1_invdynJ_fixb_snew_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3RPP1_invdynJ_fixb_snew_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 13:21:15
% EndTime: 2021-01-15 13:21:15
% DurationCPUTime: 0.10s
% Computational Cost: add. (94->36), mult. (136->37), div. (0->0), fcn. (28->2), ass. (0->16)
t91 = qJD(1) ^ 2;
t89 = sin(qJ(1));
t90 = cos(qJ(1));
t96 = t89 * g(1) - t90 * g(2);
t93 = -t91 * qJ(2) + qJDD(2) - t96;
t97 = -pkin(1) - qJ(3);
t83 = -(2 * qJD(3) * qJD(1)) + t97 * qJDD(1) + t93;
t98 = m(4) * t83;
t94 = -t90 * g(1) - t89 * g(2);
t92 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t94;
t84 = t97 * t91 + qJDD(3) + t92;
t95 = m(4) * t84 + qJDD(1) * mrSges(4,2) - t91 * mrSges(4,3);
t86 = -qJDD(1) * pkin(1) + t93;
t85 = t91 * pkin(1) - t92;
t81 = m(3) * t86 + t98 + (-mrSges(4,2) - mrSges(3,3)) * t91 + (mrSges(3,2) - mrSges(4,3)) * qJDD(1);
t1 = [mrSges(2,1) * t96 - mrSges(2,2) * t94 + mrSges(3,2) * t86 - mrSges(3,3) * t85 + mrSges(4,2) * t84 - mrSges(4,3) * t83 - qJ(3) * (-t91 * mrSges(4,2) + t98) - pkin(1) * t81 + qJ(2) * (-m(3) * t85 + t91 * mrSges(3,2) + t95) + (qJ(2) * mrSges(3,3) + qJ(3) * mrSges(4,3) + Ifges(3,1) + Ifges(4,1) + Ifges(2,3)) * qJDD(1); t81; t95;];
tauJ = t1;
