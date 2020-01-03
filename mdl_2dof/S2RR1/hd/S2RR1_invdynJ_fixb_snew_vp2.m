% Calculate vector of inverse dynamics joint torques for with Newton-Euler
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
% tauJ [2x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S2RR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(1,1),zeros(3,1),zeros(3,3),zeros(3,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR1_invdynJ_fixb_snew_vp2: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR1_invdynJ_fixb_snew_vp2: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'S2RR1_invdynJ_fixb_snew_vp2: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR1_invdynJ_fixb_snew_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'S2RR1_invdynJ_fixb_snew_vp2: m has to be [3x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [3,3]), ...
  'S2RR1_invdynJ_fixb_snew_vp2: mrSges has to be [3x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [3 6]), ...
  'S2RR1_invdynJ_fixb_snew_vp2: Ifges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:19:07
% EndTime: 2020-01-03 11:19:07
% DurationCPUTime: 0.10s
% Computational Cost: add. (80->46), mult. (155->73), div. (0->0), fcn. (76->4), ass. (0->21)
t75 = sin(qJ(2));
t84 = qJD(1) * t75;
t77 = cos(qJ(2));
t83 = qJD(1) * t77;
t82 = qJD(1) * qJD(2);
t76 = sin(qJ(1));
t78 = cos(qJ(1));
t81 = -t78 * g(1) + t76 * g(3);
t80 = -t76 * g(1) - t78 * g(3);
t79 = t77 * mrSges(3,1) - t75 * mrSges(3,2);
t73 = -t77 * qJDD(1) + t75 * t82;
t72 = -t75 * qJDD(1) - t77 * t82;
t71 = -qJD(1) ^ 2 * pkin(1) + t81;
t70 = t79 * qJD(1);
t69 = -qJDD(1) * pkin(1) + t80;
t68 = Ifges(3,5) * qJD(2) + (-t75 * Ifges(3,1) - t77 * Ifges(3,4)) * qJD(1);
t67 = Ifges(3,6) * qJD(2) + (-t75 * Ifges(3,4) - t77 * Ifges(3,2)) * qJD(1);
t66 = Ifges(3,3) * qJD(2) + (-t75 * Ifges(3,5) - t77 * Ifges(3,6)) * qJD(1);
t65 = t75 * g(2) + t77 * t69;
t64 = t77 * g(2) - t75 * t69;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t81 - mrSges(2,2) * t80 - t75 * (mrSges(3,2) * t71 - mrSges(3,3) * t64 + Ifges(3,1) * t72 + Ifges(3,4) * t73 + Ifges(3,5) * qJDD(2) - qJD(2) * t67 - t66 * t83) - t77 * (-mrSges(3,1) * t71 + mrSges(3,3) * t65 + Ifges(3,4) * t72 + Ifges(3,2) * t73 + Ifges(3,6) * qJDD(2) + qJD(2) * t68 + t66 * t84) - pkin(1) * (t77 * (m(3) * t65 - qJDD(2) * mrSges(3,2) + t73 * mrSges(3,3) - t70 * t83) - t75 * (m(3) * t64 + qJDD(2) * mrSges(3,1) - t72 * mrSges(3,3) + t70 * t84) - t79 * qJD(2) ^ 2); mrSges(3,1) * t64 - mrSges(3,2) * t65 + Ifges(3,5) * t72 + Ifges(3,6) * t73 + Ifges(3,3) * qJDD(2) + (-t75 * t67 + t77 * t68) * qJD(1);];
tauJ = t1;
