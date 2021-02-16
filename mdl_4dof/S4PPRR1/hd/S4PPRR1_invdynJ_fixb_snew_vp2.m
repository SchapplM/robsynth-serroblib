% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4PPRR1
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
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
% Datum: 2021-01-15 15:31
% Revision: 9ee2a23c87189dabbba60bb1a629b0c16b06df77 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4PPRR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR1_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR1_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR1_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR1_invdynJ_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR1_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR1_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR1_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:31:40
% EndTime: 2021-01-15 15:31:40
% DurationCPUTime: 0.20s
% Computational Cost: add. (180->35), mult. (251->43), div. (0->0), fcn. (174->6), ass. (0->23)
t124 = sin(pkin(6));
t125 = cos(pkin(6));
t117 = -t124 * g(1) + t125 * g(2) + qJDD(2);
t118 = -t125 * g(1) - t124 * g(2);
t127 = sin(qJ(3));
t129 = cos(qJ(3));
t112 = t129 * t117 - t127 * t118;
t110 = qJDD(3) * pkin(3) + t112;
t113 = t127 * t117 + t129 * t118;
t130 = qJD(3) ^ 2;
t111 = -t130 * pkin(3) + t113;
t126 = sin(qJ(4));
t128 = cos(qJ(4));
t108 = t128 * t110 - t126 * t111;
t122 = qJD(3) + qJD(4);
t120 = t122 ^ 2;
t121 = qJDD(3) + qJDD(4);
t105 = m(5) * t108 + t121 * mrSges(5,1) - t120 * mrSges(5,2);
t109 = t126 * t110 + t128 * t111;
t106 = m(5) * t109 - t120 * mrSges(5,1) - t121 * mrSges(5,2);
t132 = t128 * t105 + t126 * t106;
t131 = mrSges(5,1) * t108 - mrSges(5,2) * t109 + Ifges(5,3) * t121;
t1 = [(-m(2) - m(3) - m(4) - m(5)) * (g(3) - qJDD(1)); m(3) * t117 + t127 * (m(4) * t113 - t130 * mrSges(4,1) - qJDD(3) * mrSges(4,2) - t126 * t105 + t128 * t106) + t129 * (m(4) * t112 + qJDD(3) * mrSges(4,1) - t130 * mrSges(4,2) + t132); mrSges(4,1) * t112 - mrSges(4,2) * t113 + Ifges(4,3) * qJDD(3) + pkin(3) * t132 + t131; t131;];
tauJ = t1;
