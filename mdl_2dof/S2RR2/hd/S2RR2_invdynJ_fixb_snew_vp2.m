% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S2RR2
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
% tauJ [2x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:16
% Revision: 9ee2a23c87189dabbba60bb1a629b0c16b06df77 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S2RR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(1,1),zeros(3,1),zeros(3,3),zeros(3,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR2_invdynJ_fixb_snew_vp2: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR2_invdynJ_fixb_snew_vp2: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'S2RR2_invdynJ_fixb_snew_vp2: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR2_invdynJ_fixb_snew_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'S2RR2_invdynJ_fixb_snew_vp2: m has to be [3x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [3,3]), ...
  'S2RR2_invdynJ_fixb_snew_vp2: mrSges has to be [3x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [3 6]), ...
  'S2RR2_invdynJ_fixb_snew_vp2: Ifges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:15:48
% EndTime: 2021-01-15 12:15:48
% DurationCPUTime: 0.20s
% Computational Cost: add. (80->46), mult. (155->73), div. (0->0), fcn. (76->4), ass. (0->21)
t76 = sin(qJ(2));
t85 = qJD(1) * t76;
t78 = cos(qJ(2));
t84 = qJD(1) * t78;
t83 = qJD(1) * qJD(2);
t77 = sin(qJ(1));
t79 = cos(qJ(1));
t82 = -t79 * g(1) + t77 * g(3);
t81 = t77 * g(1) + t79 * g(3);
t80 = -t78 * mrSges(3,1) + t76 * mrSges(3,2);
t74 = t78 * qJDD(1) - t76 * t83;
t73 = t76 * qJDD(1) + t78 * t83;
t72 = -qJD(1) ^ 2 * pkin(1) - t81;
t71 = t80 * qJD(1);
t70 = qJDD(1) * pkin(1) + t82;
t69 = Ifges(3,5) * qJD(2) + (t76 * Ifges(3,1) + t78 * Ifges(3,4)) * qJD(1);
t68 = Ifges(3,6) * qJD(2) + (t76 * Ifges(3,4) + t78 * Ifges(3,2)) * qJD(1);
t67 = Ifges(3,3) * qJD(2) + (t76 * Ifges(3,5) + t78 * Ifges(3,6)) * qJD(1);
t66 = -t76 * g(2) + t78 * t70;
t65 = -t78 * g(2) - t76 * t70;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t81 - mrSges(2,2) * t82 + t76 * (mrSges(3,2) * t72 - mrSges(3,3) * t65 + Ifges(3,1) * t73 + Ifges(3,4) * t74 + Ifges(3,5) * qJDD(2) - qJD(2) * t68 + t67 * t84) + t78 * (-mrSges(3,1) * t72 + mrSges(3,3) * t66 + Ifges(3,4) * t73 + Ifges(3,2) * t74 + Ifges(3,6) * qJDD(2) + qJD(2) * t69 - t67 * t85) + pkin(1) * (t78 * (m(3) * t66 - qJDD(2) * mrSges(3,2) + t74 * mrSges(3,3) + t71 * t84) - t76 * (m(3) * t65 + qJDD(2) * mrSges(3,1) - t73 * mrSges(3,3) - t71 * t85) + t80 * qJD(2) ^ 2); mrSges(3,1) * t65 - mrSges(3,2) * t66 + Ifges(3,5) * t73 + Ifges(3,6) * t74 + Ifges(3,3) * qJDD(2) + (t76 * t68 - t78 * t69) * qJD(1);];
tauJ = t1;
