% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
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
% tauJB [(6+1)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:13
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S1R1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(1,1),zeros(3,1),zeros(1,1),zeros(2,1),zeros(2,3),zeros(2,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'S1R1_invdynJB_fixb_snew_vp2: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'S1R1_invdynJB_fixb_snew_vp2: qJD has to be [1x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [1 1]), ...
  'S1R1_invdynJB_fixb_snew_vp2: qJDD has to be [1x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S1R1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S1R1_invdynJB_fixb_snew_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [2 1]), ...
  'S1R1_invdynJB_fixb_snew_vp2: m has to be [2x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [2,3]), ...
  'S1R1_invdynJB_fixb_snew_vp2: mrSges has to be [2x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [2 6]), ...
  'S1R1_invdynJB_fixb_snew_vp2: Ifges has to be [2x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 09:12:55
% EndTime: 2020-06-19 09:12:55
% DurationCPUTime: 0.08s
% Computational Cost: add. (65->29), mult. (114->41), div. (0->0), fcn. (44->2), ass. (0->13)
t45 = sin(qJ(1));
t46 = cos(qJ(1));
t42 = t45 * g(1) - t46 * g(2);
t47 = qJD(1) ^ 2;
t39 = m(2) * t42 + qJDD(1) * mrSges(2,1) - t47 * mrSges(2,2);
t43 = -t46 * g(1) - t45 * g(2);
t40 = m(2) * t43 - t47 * mrSges(2,1) - qJDD(1) * mrSges(2,2);
t50 = t46 * t39 + t45 * t40;
t49 = -t45 * t39 + t46 * t40;
t48 = mrSges(2,1) * t42 - mrSges(2,2) * t43 + Ifges(2,3) * qJDD(1);
t38 = -mrSges(2,2) * g(3) - mrSges(2,3) * t42 + Ifges(2,5) * qJDD(1) - t47 * Ifges(2,6);
t37 = mrSges(2,1) * g(3) + mrSges(2,3) * t43 + t47 * Ifges(2,5) + Ifges(2,6) * qJDD(1);
t1 = [-m(1) * g(1) + t49; -m(1) * g(2) + t50; (-m(1) - m(2)) * g(3); -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(1) * t50 - t45 * t37 + t46 * t38; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(1) * t49 + t46 * t37 + t45 * t38; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t48; t48;];
tauJB = t1;
