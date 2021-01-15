% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4PPRP1
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
%   pkin=[a2,a3,a4,d3,theta1]';
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
% Datum: 2021-01-15 14:52
% Revision: 9ee2a23c87189dabbba60bb1a629b0c16b06df77 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4PPRP1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP1_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP1_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRP1_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRP1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP1_invdynJ_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP1_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRP1_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRP1_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:52:12
% EndTime: 2021-01-15 14:52:12
% DurationCPUTime: 0.10s
% Computational Cost: add. (99->35), mult. (138->38), div. (0->0), fcn. (70->4), ass. (0->14)
t109 = sin(pkin(5));
t110 = cos(pkin(5));
t103 = -t109 * g(1) + t110 * g(2) + qJDD(2);
t104 = -t110 * g(1) - t109 * g(2);
t111 = sin(qJ(3));
t112 = cos(qJ(3));
t100 = t111 * t103 + t112 * t104;
t113 = qJD(3) ^ 2;
t97 = -t113 * pkin(3) + qJDD(3) * qJ(4) + (2 * qJD(4) * qJD(3)) + t100;
t114 = m(5) * t97 + qJDD(3) * mrSges(5,3);
t99 = t112 * t103 - t111 * t104;
t98 = -qJDD(3) * pkin(3) - t113 * qJ(4) + qJDD(4) - t99;
t95 = m(5) * t98 - qJDD(3) * mrSges(5,1) - t113 * mrSges(5,3);
t1 = [(-m(2) - m(3) - m(4) - m(5)) * (g(3) - qJDD(1)); m(3) * t103 + t111 * (m(4) * t100 - qJDD(3) * mrSges(4,2) + (-mrSges(4,1) - mrSges(5,1)) * t113 + t114) + t112 * (m(4) * t99 + qJDD(3) * mrSges(4,1) - t113 * mrSges(4,2) - t95); mrSges(4,1) * t99 - mrSges(4,2) * t100 - mrSges(5,1) * t98 + mrSges(5,3) * t97 - pkin(3) * t95 + qJ(4) * (-t113 * mrSges(5,1) + t114) + (Ifges(4,3) + Ifges(5,2)) * qJDD(3); t95;];
tauJ = t1;
