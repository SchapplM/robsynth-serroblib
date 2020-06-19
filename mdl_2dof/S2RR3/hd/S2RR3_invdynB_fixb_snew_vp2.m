% Calculate vector of inverse dynamics base forces with Newton-Euler for
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:14
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S2RR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR3_invdynB_fixb_snew_vp2: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR3_invdynB_fixb_snew_vp2: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'S2RR3_invdynB_fixb_snew_vp2: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S2RR3_invdynB_fixb_snew_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'S2RR3_invdynB_fixb_snew_vp2: m has to be [3x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [3,3]), ...
  'S2RR3_invdynB_fixb_snew_vp2: mrSges has to be [3x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [3 6]), ...
  'S2RR3_invdynB_fixb_snew_vp2: Ifges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 09:14:25
% EndTime: 2020-06-19 09:14:25
% DurationCPUTime: 0.19s
% Computational Cost: add. (426->61), mult. (622->77), div. (0->0), fcn. (290->4), ass. (0->27)
t100 = qJD(1) ^ 2;
t97 = sin(qJ(1));
t99 = cos(qJ(1));
t90 = t97 * g(1) - t99 * g(2);
t88 = qJDD(1) * pkin(1) + t90;
t91 = -t99 * g(1) - t97 * g(2);
t89 = -t100 * pkin(1) + t91;
t96 = sin(qJ(2));
t98 = cos(qJ(2));
t86 = t98 * t88 - t96 * t89;
t95 = qJD(1) + qJD(2);
t93 = t95 ^ 2;
t94 = qJDD(1) + qJDD(2);
t84 = m(3) * t86 + t94 * mrSges(3,1) - t93 * mrSges(3,2);
t87 = t96 * t88 + t98 * t89;
t85 = m(3) * t87 - t93 * mrSges(3,1) - t94 * mrSges(3,2);
t78 = t98 * t84 + t96 * t85;
t76 = m(2) * t90 + qJDD(1) * mrSges(2,1) - t100 * mrSges(2,2) + t78;
t101 = -t96 * t84 + t98 * t85;
t77 = m(2) * t91 - t100 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t101;
t103 = t99 * t76 + t97 * t77;
t102 = -t97 * t76 + t99 * t77;
t83 = -mrSges(3,2) * g(3) - mrSges(3,3) * t86 + Ifges(3,5) * t94 - t93 * Ifges(3,6);
t82 = mrSges(3,1) * g(3) + mrSges(3,3) * t87 + t93 * Ifges(3,5) + Ifges(3,6) * t94;
t72 = -mrSges(2,2) * g(3) - mrSges(2,3) * t90 + Ifges(2,5) * qJDD(1) - t100 * Ifges(2,6) - pkin(3) * t78 - t96 * t82 + t98 * t83;
t71 = Ifges(2,6) * qJDD(1) + t100 * Ifges(2,5) + mrSges(2,3) * t91 + t96 * t83 + t98 * t82 + pkin(3) * t101 + (pkin(1) * m(3) + mrSges(2,1)) * g(3);
t1 = [-m(1) * g(1) + t102; -m(1) * g(2) + t103; (-m(1) - m(2) - m(3)) * g(3); -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(2) * t103 - t97 * t71 + t99 * t72; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(2) * t102 + t99 * t71 + t97 * t72; -mrSges(1,1) * g(2) + mrSges(2,1) * t90 + mrSges(3,1) * t86 + mrSges(1,2) * g(1) - mrSges(2,2) * t91 - mrSges(3,2) * t87 + Ifges(2,3) * qJDD(1) + Ifges(3,3) * t94 + pkin(1) * t78;];
tauB = t1;
