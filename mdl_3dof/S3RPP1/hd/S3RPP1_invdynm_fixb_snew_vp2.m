% Calculate vector of cutting torques with Newton-Euler for
% S3RPP1
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% qJDD [3x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1]';
% m_mdh [4x1]
%   mass of all robot links (including the base)
% mrSges [4x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% m [3x4]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S3RPP1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPP1_invdynm_fixb_snew_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPP1_invdynm_fixb_snew_vp2: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3RPP1_invdynm_fixb_snew_vp2: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RPP1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3RPP1_invdynm_fixb_snew_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3RPP1_invdynm_fixb_snew_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3RPP1_invdynm_fixb_snew_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3RPP1_invdynm_fixb_snew_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:27:22
% EndTime: 2019-05-04 18:27:22
% DurationCPUTime: 0.14s
% Computational Cost: add. (550->90), mult. (801->91), div. (0->0), fcn. (154->2), ass. (0->33)
t63 = qJD(1) ^ 2;
t80 = t63 * mrSges(4,3);
t79 = Ifges(3,4) - Ifges(4,5);
t78 = Ifges(3,5) - Ifges(2,6);
t77 = -pkin(1) - qJ(3);
t60 = sin(qJ(1));
t61 = cos(qJ(1));
t45 = -t61 * g(1) - t60 * g(2);
t66 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t45;
t39 = t63 * t77 + qJDD(3) + t66;
t76 = m(4) * t39 + qJDD(1) * mrSges(4,2);
t44 = t60 * g(1) - g(2) * t61;
t75 = mrSges(4,1) * t39 - Ifges(4,4) * qJDD(1);
t70 = -t63 * qJ(2) + qJDD(2) - t44;
t38 = -(2 * qJD(3) * qJD(1)) + qJDD(1) * t77 + t70;
t74 = mrSges(4,1) * t38 + mrSges(4,2) * g(3) + Ifges(4,4) * t63 + Ifges(4,5) * qJDD(1);
t33 = m(4) * t38 - mrSges(4,2) * t63 - qJDD(1) * mrSges(4,3);
t73 = -m(4) * qJ(3) + mrSges(3,2) - mrSges(4,3);
t72 = mrSges(4,2) * t39 - mrSges(4,3) * t38 + Ifges(4,1) * qJDD(1);
t40 = pkin(1) * t63 - t66;
t71 = -m(3) * t40 + mrSges(3,2) * t63 + qJDD(1) * mrSges(3,3) + t76;
t42 = -qJDD(1) * pkin(1) + t70;
t69 = -m(3) * t42 + mrSges(3,3) * t63 - t33;
t68 = mrSges(3,1) * t40 + pkin(2) * (-t76 + t80) - t75;
t67 = mrSges(3,1) * t42 + pkin(2) * t33 + t74;
t65 = mrSges(3,2) * t42 - mrSges(3,3) * t40 + Ifges(3,1) * qJDD(1) - qJ(3) * t33 + t72;
t64 = -mrSges(2,2) * t45 + qJ(2) * (t71 - t80) + pkin(1) * (-qJDD(1) * mrSges(3,2) + t69) + mrSges(2,1) * t44 + Ifges(2,3) * qJDD(1) + t65;
t46 = (-m(3) - m(4)) * g(3);
t30 = m(2) * t44 - t63 * mrSges(2,2) + (mrSges(2,1) - mrSges(3,2)) * qJDD(1) + t69;
t29 = m(2) * t45 - qJDD(1) * mrSges(2,2) + (-mrSges(2,1) - mrSges(4,3)) * t63 + t71;
t28 = -mrSges(2,3) * t44 - qJ(2) * t46 + t78 * t63 + (-Ifges(3,4) + Ifges(2,5)) * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) + t67;
t27 = mrSges(2,3) * t45 - pkin(1) * t46 - t78 * qJDD(1) + (Ifges(2,5) - t79) * t63 + (mrSges(2,1) - t73) * g(3) - t68;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t61 * t28 - t60 * t27 - pkin(3) * (t29 * t60 + t30 * t61), t28, t65, t72; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t60 * t28 + t61 * t27 + pkin(3) * (t29 * t61 - t30 * t60), t27, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - t63 * Ifges(3,5) - t67, -mrSges(4,3) * g(3) - t63 * Ifges(4,5) - t75; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t64, t64, Ifges(3,5) * qJDD(1) + g(3) * t73 + t63 * t79 + t68, t74;];
m_new  = t1;
