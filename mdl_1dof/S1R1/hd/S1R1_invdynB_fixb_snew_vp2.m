% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S1R1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% qJDD [1x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1]';
% m [2x1]
%   mass of all robot links (including the base)
% mrSges [2x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [2x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:13
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S1R1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(1,1),zeros(3,1),zeros(1,1),zeros(2,1),zeros(2,3),zeros(2,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'S1R1_invdynB_fixb_snew_vp2: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'S1R1_invdynB_fixb_snew_vp2: qJD has to be [1x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [1 1]), ...
  'S1R1_invdynB_fixb_snew_vp2: qJDD has to be [1x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S1R1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S1R1_invdynB_fixb_snew_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [2 1]), ...
  'S1R1_invdynB_fixb_snew_vp2: m has to be [2x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [2,3]), ...
  'S1R1_invdynB_fixb_snew_vp2: mrSges has to be [2x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [2 6]), ...
  'S1R1_invdynB_fixb_snew_vp2: Ifges has to be [2x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 09:12:55
% EndTime: 2020-06-19 09:12:55
% DurationCPUTime: 0.07s
% Computational Cost: add. (61->29), mult. (107->41), div. (0->0), fcn. (40->2), ass. (0->12)
t29 = sin(qJ(1));
t30 = cos(qJ(1));
t27 = t29 * g(1) - t30 * g(2);
t31 = qJD(1) ^ 2;
t25 = m(2) * t27 + qJDD(1) * mrSges(2,1) - t31 * mrSges(2,2);
t28 = -t30 * g(1) - t29 * g(2);
t26 = m(2) * t28 - t31 * mrSges(2,1) - qJDD(1) * mrSges(2,2);
t33 = t30 * t25 + t29 * t26;
t32 = -t29 * t25 + t30 * t26;
t24 = -mrSges(2,2) * g(3) - mrSges(2,3) * t27 + Ifges(2,5) * qJDD(1) - t31 * Ifges(2,6);
t23 = mrSges(2,1) * g(3) + mrSges(2,3) * t28 + t31 * Ifges(2,5) + Ifges(2,6) * qJDD(1);
t1 = [-m(1) * g(1) + t32; -m(1) * g(2) + t33; (-m(1) - m(2)) * g(3); -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(1) * t33 - t29 * t23 + t30 * t24; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(1) * t32 + t30 * t23 + t29 * t24; -mrSges(1,1) * g(2) + mrSges(2,1) * t27 + mrSges(1,2) * g(1) - mrSges(2,2) * t28 + Ifges(2,3) * qJDD(1);];
tauB = t1;
