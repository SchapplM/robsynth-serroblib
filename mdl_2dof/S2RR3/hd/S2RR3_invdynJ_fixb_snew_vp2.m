% Calculate vector of inverse dynamics joint torques for with Newton-Euler
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
% tauJ [2x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:14
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S2RR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR3_invdynJ_fixb_snew_vp2: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR3_invdynJ_fixb_snew_vp2: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'S2RR3_invdynJ_fixb_snew_vp2: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S2RR3_invdynJ_fixb_snew_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'S2RR3_invdynJ_fixb_snew_vp2: m has to be [3x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [3,3]), ...
  'S2RR3_invdynJ_fixb_snew_vp2: mrSges has to be [3x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [3 6]), ...
  'S2RR3_invdynJ_fixb_snew_vp2: Ifges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 09:14:25
% EndTime: 2020-06-19 09:14:25
% DurationCPUTime: 0.07s
% Computational Cost: add. (51->19), mult. (78->27), div. (0->0), fcn. (42->4), ass. (0->15)
t65 = sin(qJ(1));
t67 = cos(qJ(1));
t70 = t65 * g(1) - t67 * g(2);
t57 = qJDD(1) * pkin(1) + t70;
t68 = -t67 * g(1) - t65 * g(2);
t58 = -qJD(1) ^ 2 * pkin(1) + t68;
t64 = sin(qJ(2));
t66 = cos(qJ(2));
t55 = t66 * t57 - t64 * t58;
t56 = t64 * t57 + t66 * t58;
t62 = qJDD(1) + qJDD(2);
t69 = mrSges(3,1) * t55 - mrSges(3,2) * t56 + Ifges(3,3) * t62;
t63 = qJD(1) + qJD(2);
t61 = t63 ^ 2;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t70 - mrSges(2,2) * t68 + pkin(1) * (t64 * (m(3) * t56 - t61 * mrSges(3,1) - t62 * mrSges(3,2)) + t66 * (m(3) * t55 + t62 * mrSges(3,1) - t61 * mrSges(3,2))) + t69; t69;];
tauJ = t1;
