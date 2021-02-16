% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S3PRP1
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
%   pkin=[a2,a3,d2]';
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
% Datum: 2021-01-15 12:54
% Revision: 9ee2a23c87189dabbba60bb1a629b0c16b06df77 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S3PRP1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRP1_invdynJ_fixb_snew_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRP1_invdynJ_fixb_snew_vp2: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3PRP1_invdynJ_fixb_snew_vp2: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3PRP1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PRP1_invdynJ_fixb_snew_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3PRP1_invdynJ_fixb_snew_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3PRP1_invdynJ_fixb_snew_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3PRP1_invdynJ_fixb_snew_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:53:57
% EndTime: 2021-01-15 12:53:57
% DurationCPUTime: 0.08s
% Computational Cost: add. (69->29), mult. (88->33), div. (0->0), fcn. (24->2), ass. (0->11)
t85 = -g(2) + qJDD(1);
t87 = sin(qJ(2));
t88 = cos(qJ(2));
t81 = -t88 * g(1) + t87 * t85;
t89 = qJD(2) ^ 2;
t78 = -t89 * pkin(2) + qJDD(2) * qJ(3) + (2 * qJD(3) * qJD(2)) + t81;
t90 = m(4) * t78 + qJDD(2) * mrSges(4,3);
t80 = t87 * g(1) + t88 * t85;
t79 = -qJDD(2) * pkin(2) - t89 * qJ(3) + qJDD(3) - t80;
t76 = m(4) * t79 - qJDD(2) * mrSges(4,1) - t89 * mrSges(4,3);
t1 = [m(2) * t85 + t87 * (m(3) * t81 - qJDD(2) * mrSges(3,2) + (-mrSges(3,1) - mrSges(4,1)) * t89 + t90) + t88 * (m(3) * t80 + qJDD(2) * mrSges(3,1) - t89 * mrSges(3,2) - t76); mrSges(3,1) * t80 - mrSges(3,2) * t81 - mrSges(4,1) * t79 + mrSges(4,3) * t78 - pkin(2) * t76 + qJ(3) * (-t89 * mrSges(4,1) + t90) + (Ifges(3,3) + Ifges(4,2)) * qJDD(2); t76;];
tauJ = t1;
