% Calculate vector of cutting torques with Newton-Euler for
% S2RR1
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
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
% m_mdh [3x1]
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
% Datum: 2020-01-03 11:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S2RR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(1,1),zeros(3,1),zeros(3,3),zeros(3,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR1_invdynm_fixb_snew_vp2: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR1_invdynm_fixb_snew_vp2: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'S2RR1_invdynm_fixb_snew_vp2: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR1_invdynm_fixb_snew_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'S2RR1_invdynm_fixb_snew_vp2: m has to be [3x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [3,3]), ...
  'S2RR1_invdynm_fixb_snew_vp2: mrSges has to be [3x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [3 6]), ...
  'S2RR1_invdynm_fixb_snew_vp2: Ifges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:19:07
% EndTime: 2020-01-03 11:19:07
% DurationCPUTime: 0.14s
% Computational Cost: add. (439->68), mult. (837->99), div. (0->0), fcn. (396->4), ass. (0->29)
t54 = sin(qJ(2));
t63 = qJD(1) * t54;
t56 = cos(qJ(2));
t62 = qJD(1) * t56;
t61 = qJD(1) * qJD(2);
t55 = sin(qJ(1));
t57 = cos(qJ(1));
t50 = -t57 * g(1) + t55 * g(3);
t49 = -t55 * g(1) - t57 * g(3);
t43 = -qJDD(1) * pkin(1) + t49;
t35 = t56 * g(2) - t54 * t43;
t36 = t54 * g(2) + t56 * t43;
t39 = Ifges(3,6) * qJD(2) + (-Ifges(3,4) * t54 - Ifges(3,2) * t56) * qJD(1);
t40 = Ifges(3,5) * qJD(2) + (-Ifges(3,1) * t54 - Ifges(3,4) * t56) * qJD(1);
t46 = -t54 * qJDD(1) - t56 * t61;
t47 = -t56 * qJDD(1) + t54 * t61;
t60 = mrSges(3,1) * t35 - mrSges(3,2) * t36 + Ifges(3,5) * t46 + Ifges(3,6) * t47 + Ifges(3,3) * qJDD(2) - t39 * t63 + t40 * t62;
t38 = Ifges(3,3) * qJD(2) + (-Ifges(3,5) * t54 - Ifges(3,6) * t56) * qJD(1);
t58 = qJD(1) ^ 2;
t45 = -t58 * pkin(1) + t50;
t30 = -mrSges(3,1) * t45 + mrSges(3,3) * t36 + Ifges(3,4) * t46 + Ifges(3,2) * t47 + Ifges(3,6) * qJDD(2) + qJD(2) * t40 + t38 * t63;
t31 = mrSges(3,2) * t45 - mrSges(3,3) * t35 + Ifges(3,1) * t46 + Ifges(3,4) * t47 + Ifges(3,5) * qJDD(2) - qJD(2) * t39 - t38 * t62;
t44 = (mrSges(3,1) * t56 - mrSges(3,2) * t54) * qJD(1);
t32 = m(3) * t35 - t46 * mrSges(3,3) + qJDD(2) * mrSges(3,1) + t44 * t63 + qJD(2) * (-qJD(2) * mrSges(3,2) - mrSges(3,3) * t62);
t33 = m(3) * t36 + t47 * mrSges(3,3) - qJDD(2) * mrSges(3,2) - t44 * t62 - qJD(2) * (qJD(2) * mrSges(3,1) + mrSges(3,3) * t63);
t59 = -mrSges(2,2) * t49 - pkin(1) * (-t54 * t32 + t56 * t33) - t56 * t30 - t54 * t31 + mrSges(2,1) * t50 + Ifges(2,3) * qJDD(1);
t29 = mrSges(2,1) * g(2) + mrSges(2,3) * t49 + t58 * Ifges(2,5) + Ifges(2,6) * qJDD(1) + t60;
t27 = Ifges(2,5) * qJDD(1) - t58 * Ifges(2,6) - mrSges(2,2) * g(2) - mrSges(2,3) * t50 + t56 * t31 - t54 * t30 + pkin(1) * (-t56 * t32 - t54 * t33);
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t55 * t27 + t57 * t29, t27, t31; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t59, t29, t30; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t57 * t27 - t55 * t29, t59, t60;];
m_new = t1;
