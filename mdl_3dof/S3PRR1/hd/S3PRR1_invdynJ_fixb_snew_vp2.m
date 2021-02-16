% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S3PRR1
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
%   pkin=[a2,a3,d2,d3]';
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
% Datum: 2021-01-15 13:07
% Revision: 9ee2a23c87189dabbba60bb1a629b0c16b06df77 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S3PRR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRR1_invdynJ_fixb_snew_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRR1_invdynJ_fixb_snew_vp2: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3PRR1_invdynJ_fixb_snew_vp2: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3PRR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3PRR1_invdynJ_fixb_snew_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3PRR1_invdynJ_fixb_snew_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3PRR1_invdynJ_fixb_snew_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3PRR1_invdynJ_fixb_snew_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 13:07:07
% EndTime: 2021-01-15 13:07:07
% DurationCPUTime: 0.09s
% Computational Cost: add. (124->29), mult. (149->38), div. (0->0), fcn. (76->4), ass. (0->20)
t102 = sin(qJ(3));
t104 = cos(qJ(3));
t101 = -g(2) + qJDD(1);
t103 = sin(qJ(2));
t105 = cos(qJ(2));
t92 = t103 * g(1) + t105 * t101;
t90 = qJDD(2) * pkin(2) + t92;
t106 = qJD(2) ^ 2;
t93 = -t105 * g(1) + t103 * t101;
t91 = -t106 * pkin(2) + t93;
t88 = -t102 * t91 + t104 * t90;
t100 = qJD(2) + qJD(3);
t98 = t100 ^ 2;
t99 = qJDD(2) + qJDD(3);
t85 = m(4) * t88 + t99 * mrSges(4,1) - t98 * mrSges(4,2);
t89 = t102 * t90 + t104 * t91;
t86 = m(4) * t89 - t98 * mrSges(4,1) - t99 * mrSges(4,2);
t108 = t102 * t86 + t104 * t85;
t107 = mrSges(4,1) * t88 - mrSges(4,2) * t89 + Ifges(4,3) * t99;
t1 = [m(2) * t101 + t103 * (m(3) * t93 - t106 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t102 * t85 + t104 * t86) + t105 * (m(3) * t92 + qJDD(2) * mrSges(3,1) - t106 * mrSges(3,2) + t108); mrSges(3,1) * t92 - mrSges(3,2) * t93 + Ifges(3,3) * qJDD(2) + pkin(2) * t108 + t107; t107;];
tauJ = t1;
