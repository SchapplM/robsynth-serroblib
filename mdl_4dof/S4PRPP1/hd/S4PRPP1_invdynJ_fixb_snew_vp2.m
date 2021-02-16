% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4PRPP1
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
%   pkin=[a2,a3,a4,d2,theta1]';
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
% Datum: 2021-01-15 15:57
% Revision: 9ee2a23c87189dabbba60bb1a629b0c16b06df77 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4PRPP1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP1_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP1_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPP1_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP1_invdynJ_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP1_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPP1_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPP1_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:57:41
% EndTime: 2021-01-15 15:57:42
% DurationCPUTime: 0.18s
% Computational Cost: add. (129->42), mult. (196->42), div. (0->0), fcn. (84->4), ass. (0->20)
t134 = qJD(2) ^ 2;
t130 = sin(pkin(5));
t131 = cos(pkin(5));
t126 = g(1) * t130 - g(2) * t131;
t127 = -g(1) * t131 - g(2) * t130;
t132 = sin(qJ(2));
t133 = cos(qJ(2));
t138 = t133 * t126 - t132 * t127;
t136 = -t134 * qJ(3) + qJDD(3) - t138;
t140 = -pkin(2) - qJ(4);
t119 = -(2 * qJD(4) * qJD(2)) + qJDD(2) * t140 + t136;
t141 = m(5) * t119;
t139 = t132 * t126 + t133 * t127;
t135 = qJDD(2) * qJ(3) + (2 * qJD(3) * qJD(2)) + t139;
t120 = t134 * t140 + qJDD(4) + t135;
t137 = m(5) * t120 + qJDD(2) * mrSges(5,2) - mrSges(5,3) * t134;
t122 = -qJDD(2) * pkin(2) + t136;
t121 = pkin(2) * t134 - t135;
t117 = m(4) * t122 + t141 + (-mrSges(5,2) - mrSges(4,3)) * t134 + (mrSges(4,2) - mrSges(5,3)) * qJDD(2);
t1 = [(m(2) + m(3) + m(4) + m(5)) * (-g(3) + qJDD(1)); mrSges(3,1) * t138 - mrSges(3,2) * t139 + mrSges(4,2) * t122 - mrSges(4,3) * t121 + mrSges(5,2) * t120 - mrSges(5,3) * t119 - qJ(4) * (-t134 * mrSges(5,2) + t141) - pkin(2) * t117 + qJ(3) * (-m(4) * t121 + mrSges(4,2) * t134 + t137) + (mrSges(4,3) * qJ(3) + mrSges(5,3) * qJ(4) + Ifges(4,1) + Ifges(5,1) + Ifges(3,3)) * qJDD(2); t117; t137;];
tauJ = t1;
