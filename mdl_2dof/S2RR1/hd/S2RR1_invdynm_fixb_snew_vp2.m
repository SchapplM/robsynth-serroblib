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
% Datum: 2019-05-04 18:18
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
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
% StartTime: 2019-05-04 18:17:41
% EndTime: 2019-05-04 18:17:41
% DurationCPUTime: 0.14s
% Computational Cost: add. (439->68), mult. (837->99), div. (0->0), fcn. (396->4), ass. (0->29)
t59 = sin(qJ(1));
t61 = cos(qJ(1));
t51 = t59 * g(1) + t61 * g(3);
t58 = sin(qJ(2));
t67 = qJD(1) * t58;
t60 = cos(qJ(2));
t66 = qJD(1) * t60;
t65 = qJD(1) * qJD(2);
t52 = t61 * g(1) - t59 * g(3);
t45 = -qJDD(1) * pkin(1) + t51;
t37 = t60 * g(2) - t58 * t45;
t38 = t58 * g(2) + t60 * t45;
t41 = Ifges(3,6) * qJD(2) + (-Ifges(3,4) * t58 - Ifges(3,2) * t60) * qJD(1);
t42 = Ifges(3,5) * qJD(2) + (-Ifges(3,1) * t58 - Ifges(3,4) * t60) * qJD(1);
t48 = -t58 * qJDD(1) - t60 * t65;
t49 = -t60 * qJDD(1) + t58 * t65;
t64 = mrSges(3,1) * t37 - mrSges(3,2) * t38 + Ifges(3,5) * t48 + Ifges(3,6) * t49 + Ifges(3,3) * qJDD(2) - t41 * t67 + t42 * t66;
t40 = Ifges(3,3) * qJD(2) + (-Ifges(3,5) * t58 - Ifges(3,6) * t60) * qJD(1);
t62 = qJD(1) ^ 2;
t47 = -t62 * pkin(1) + t52;
t32 = -mrSges(3,1) * t47 + mrSges(3,3) * t38 + Ifges(3,4) * t48 + Ifges(3,2) * t49 + Ifges(3,6) * qJDD(2) + qJD(2) * t42 + t40 * t67;
t33 = mrSges(3,2) * t47 - mrSges(3,3) * t37 + Ifges(3,1) * t48 + Ifges(3,4) * t49 + Ifges(3,5) * qJDD(2) - qJD(2) * t41 - t40 * t66;
t46 = (mrSges(3,1) * t60 - mrSges(3,2) * t58) * qJD(1);
t34 = m(3) * t37 - t48 * mrSges(3,3) + qJDD(2) * mrSges(3,1) + t46 * t67 + qJD(2) * (-qJD(2) * mrSges(3,2) - mrSges(3,3) * t66);
t35 = m(3) * t38 + t49 * mrSges(3,3) - qJDD(2) * mrSges(3,2) - t46 * t66 - qJD(2) * (qJD(2) * mrSges(3,1) + mrSges(3,3) * t67);
t63 = -mrSges(2,2) * t51 - pkin(1) * (-t58 * t34 + t60 * t35) - t60 * t32 - t58 * t33 + mrSges(2,1) * t52 + Ifges(2,3) * qJDD(1);
t31 = mrSges(2,1) * g(2) + mrSges(2,3) * t51 + t62 * Ifges(2,5) + Ifges(2,6) * qJDD(1) + t64;
t29 = Ifges(2,5) * qJDD(1) - t62 * Ifges(2,6) - mrSges(2,2) * g(2) - mrSges(2,3) * t52 + t60 * t33 - t58 * t32 + pkin(1) * (-t60 * t34 - t58 * t35);
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - t59 * t29 - t61 * t31, t29, t33; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t63, t31, t32; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - t61 * t29 + t59 * t31, t63, t64;];
m_new  = t1;
