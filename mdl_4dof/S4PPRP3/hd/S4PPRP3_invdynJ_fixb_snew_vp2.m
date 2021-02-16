% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4PPRP3
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
%   pkin=[a2,a3,a4,d3]';
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
% Datum: 2021-01-15 15:17
% Revision: 9ee2a23c87189dabbba60bb1a629b0c16b06df77 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4PPRP3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(4,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP3_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP3_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRP3_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRP3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PPRP3_invdynJ_fixb_snew_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP3_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRP3_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRP3_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:17:26
% EndTime: 2021-01-15 15:17:26
% DurationCPUTime: 0.08s
% Computational Cost: add. (83->28), mult. (87->30), div. (0->0), fcn. (30->2), ass. (0->14)
t119 = -mrSges(4,2) - mrSges(5,2);
t113 = -g(2) + qJDD(1);
t114 = -g(1) + qJDD(2);
t115 = sin(qJ(3));
t116 = cos(qJ(3));
t107 = -t115 * t113 + t116 * t114;
t105 = qJDD(3) * pkin(3) + t107;
t118 = m(5) * t105 + qJDD(3) * mrSges(5,1);
t108 = t116 * t113 + t115 * t114;
t117 = qJD(3) ^ 2;
t106 = -t117 * pkin(3) + t108;
t103 = m(4) * t108 + m(5) * t106 + (-mrSges(4,1) - mrSges(5,1)) * t117 + t119 * qJDD(3);
t102 = m(4) * t107 + qJDD(3) * mrSges(4,1) + t119 * t117 + t118;
t1 = [-t115 * t102 + t116 * t103 + (m(2) + m(3)) * t113; m(3) * t114 + t116 * t102 + t115 * t103; mrSges(4,1) * t107 - mrSges(4,2) * t108 + mrSges(5,1) * t105 - mrSges(5,2) * t106 + pkin(3) * (-t117 * mrSges(5,2) + t118) + (Ifges(4,3) + Ifges(5,3)) * qJDD(3); m(5) * (g(3) + qJDD(4));];
tauJ = t1;
