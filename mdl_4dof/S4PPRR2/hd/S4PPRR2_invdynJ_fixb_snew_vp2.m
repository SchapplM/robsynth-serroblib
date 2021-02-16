% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4PPRR2
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
%   pkin=[a2,a3,a4,d3,d4,theta2]';
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
% Datum: 2021-01-15 15:44
% Revision: 9ee2a23c87189dabbba60bb1a629b0c16b06df77 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4PPRR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR2_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR2_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR2_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR2_invdynJ_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR2_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR2_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR2_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:44:39
% EndTime: 2021-01-15 15:44:39
% DurationCPUTime: 0.13s
% Computational Cost: add. (295->38), mult. (366->49), div. (0->0), fcn. (252->6), ass. (0->26)
t137 = -g(2) + qJDD(1);
t138 = sin(pkin(6));
t139 = cos(pkin(6));
t130 = t138 * g(1) + t139 * t137;
t131 = -t139 * g(1) + t138 * t137;
t141 = sin(qJ(3));
t143 = cos(qJ(3));
t125 = t143 * t130 - t141 * t131;
t123 = qJDD(3) * pkin(3) + t125;
t126 = t141 * t130 + t143 * t131;
t144 = qJD(3) ^ 2;
t124 = -t144 * pkin(3) + t126;
t140 = sin(qJ(4));
t142 = cos(qJ(4));
t121 = t142 * t123 - t140 * t124;
t135 = qJD(3) + qJD(4);
t133 = t135 ^ 2;
t134 = qJDD(3) + qJDD(4);
t118 = m(5) * t121 + t134 * mrSges(5,1) - t133 * mrSges(5,2);
t122 = t140 * t123 + t142 * t124;
t119 = m(5) * t122 - t133 * mrSges(5,1) - t134 * mrSges(5,2);
t146 = t142 * t118 + t140 * t119;
t145 = mrSges(5,1) * t121 - mrSges(5,2) * t122 + Ifges(5,3) * t134;
t115 = m(4) * t126 - t144 * mrSges(4,1) - qJDD(3) * mrSges(4,2) - t140 * t118 + t142 * t119;
t114 = m(4) * t125 + qJDD(3) * mrSges(4,1) - t144 * mrSges(4,2) + t146;
t1 = [m(2) * t137 + t138 * (m(3) * t131 - t141 * t114 + t143 * t115) + t139 * (m(3) * t130 + t143 * t114 + t141 * t115); (m(3) + m(4) + m(5)) * (-g(3) + qJDD(2)); mrSges(4,1) * t125 - mrSges(4,2) * t126 + Ifges(4,3) * qJDD(3) + pkin(3) * t146 + t145; t145;];
tauJ = t1;
