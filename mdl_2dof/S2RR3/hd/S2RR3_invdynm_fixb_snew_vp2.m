% Calculate vector of cutting torques with Newton-Euler for
% S2RR3
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% qJDD [2x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,d1,d2]';
% m [3x1]
%   mass of all robot links (including the base)
% mrSges [3x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [3x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% m [3x3]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:14
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S2RR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR3_invdynm_fixb_snew_vp2: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR3_invdynm_fixb_snew_vp2: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'S2RR3_invdynm_fixb_snew_vp2: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S2RR3_invdynm_fixb_snew_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'S2RR3_invdynm_fixb_snew_vp2: m has to be [3x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [3,3]), ...
  'S2RR3_invdynm_fixb_snew_vp2: mrSges has to be [3x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [3 6]), ...
  'S2RR3_invdynm_fixb_snew_vp2: Ifges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 09:14:25
% EndTime: 2020-06-19 09:14:25
% DurationCPUTime: 0.23s
% Computational Cost: add. (492->57), mult. (717->74), div. (0->0), fcn. (336->4), ass. (0->27)
t47 = sin(qJ(1));
t49 = cos(qJ(1));
t38 = t47 * g(1) - t49 * g(2);
t35 = qJDD(1) * pkin(1) + t38;
t39 = -t49 * g(1) - t47 * g(2);
t50 = qJD(1) ^ 2;
t36 = -t50 * pkin(1) + t39;
t46 = sin(qJ(2));
t48 = cos(qJ(2));
t33 = t48 * t35 - t46 * t36;
t44 = qJD(1) + qJD(2);
t42 = t44 ^ 2;
t43 = qJDD(1) + qJDD(2);
t30 = m(3) * t33 + t43 * mrSges(3,1) - t42 * mrSges(3,2);
t34 = t46 * t35 + t48 * t36;
t31 = m(3) * t34 - t42 * mrSges(3,1) - t43 * mrSges(3,2);
t24 = t48 * t30 + t46 * t31;
t53 = -t46 * t30 + t48 * t31;
t52 = mrSges(3,1) * t33 - mrSges(3,2) * t34 + Ifges(3,3) * t43;
t51 = mrSges(2,1) * t38 - mrSges(2,2) * t39 + Ifges(2,3) * qJDD(1) + pkin(1) * t24 + t52;
t29 = -mrSges(3,2) * g(3) - mrSges(3,3) * t33 + Ifges(3,5) * t43 - t42 * Ifges(3,6);
t28 = mrSges(3,1) * g(3) + mrSges(3,3) * t34 + t42 * Ifges(3,5) + Ifges(3,6) * t43;
t22 = m(2) * t39 - t50 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t53;
t21 = m(2) * t38 + qJDD(1) * mrSges(2,1) - t50 * mrSges(2,2) + t24;
t20 = -mrSges(2,2) * g(3) - mrSges(2,3) * t38 + Ifges(2,5) * qJDD(1) - t50 * Ifges(2,6) - pkin(3) * t24 - t46 * t28 + t48 * t29;
t19 = Ifges(2,6) * qJDD(1) + t50 * Ifges(2,5) + mrSges(2,3) * t39 + t46 * t29 + t48 * t28 + pkin(3) * t53 + (pkin(1) * m(3) + mrSges(2,1)) * g(3);
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t49 * t20 - t47 * t19 - pkin(2) * (t49 * t21 + t47 * t22), t20, t29; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t47 * t20 + t49 * t19 + pkin(2) * (-t47 * t21 + t49 * t22), t19, t28; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t51, t51, t52;];
m_new = t1;
