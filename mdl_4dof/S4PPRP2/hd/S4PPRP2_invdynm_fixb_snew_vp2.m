% Calculate vector of cutting torques with Newton-Euler for
% S4PPRP2
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% m [3x5]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4PPRP2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP2_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP2_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRP2_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRP2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP2_invdynm_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP2_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRP2_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRP2_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:44:19
% EndTime: 2019-05-04 18:44:19
% DurationCPUTime: 0.22s
% Computational Cost: add. (1994->91), mult. (2500->102), div. (0->0), fcn. (1236->4), ass. (0->40)
t82 = -g(2) + qJDD(1);
t85 = sin(pkin(5));
t86 = cos(pkin(5));
t71 = t85 * g(1) + t86 * t82;
t72 = -t86 * g(1) + t85 * t82;
t87 = sin(qJ(3));
t88 = cos(qJ(3));
t66 = t87 * t71 + t88 * t72;
t91 = qJD(3) ^ 2;
t61 = -t91 * pkin(3) + qJDD(3) * qJ(4) + (2 * qJD(4) * qJD(3)) + t66;
t107 = mrSges(5,2) * t61;
t106 = -mrSges(4,1) - mrSges(5,1);
t105 = m(5) * t61 + qJDD(3) * mrSges(5,3);
t53 = m(4) * t66 - qJDD(3) * mrSges(4,2) + t106 * t91 + t105;
t65 = t88 * t71 - t87 * t72;
t63 = -qJDD(3) * pkin(3) - t91 * qJ(4) + qJDD(4) - t65;
t99 = -m(5) * t63 + qJDD(3) * mrSges(5,1) + t91 * mrSges(5,3);
t54 = m(4) * t65 + qJDD(3) * mrSges(4,1) - t91 * mrSges(4,2) + t99;
t47 = t87 * t53 + t88 * t54;
t81 = -g(3) + qJDD(2);
t104 = (m(4) + m(5)) * t81;
t103 = mrSges(5,2) * t63 + Ifges(5,4) * qJDD(3) + t91 * Ifges(5,6);
t44 = m(3) * t71 + t47;
t101 = t88 * t53 - t87 * t54;
t45 = m(3) * t72 + t101;
t102 = -t85 * t44 + t86 * t45;
t100 = m(3) * t81 + t104;
t98 = qJ(1) * m(2) + mrSges(1,2) + mrSges(2,3);
t97 = -mrSges(5,1) * t63 + mrSges(5,3) * t61 + Ifges(5,2) * qJDD(3);
t48 = t107 + mrSges(4,3) * t66 + (Ifges(5,4) + Ifges(4,5)) * t91 + (Ifges(4,6) - Ifges(5,6)) * qJDD(3) + (-m(5) * pkin(3) + t106) * t81;
t49 = -mrSges(4,3) * t65 + Ifges(4,5) * qJDD(3) - t91 * Ifges(4,6) + (-m(5) * qJ(4) + mrSges(4,2) - mrSges(5,3)) * t81 + t103;
t37 = -mrSges(3,1) * t81 + mrSges(3,3) * t72 - pkin(2) * t104 + pkin(4) * t101 + t88 * t48 + t87 * t49;
t40 = mrSges(3,2) * t81 - mrSges(3,3) * t71 - pkin(4) * t47 - t87 * t48 + t88 * t49;
t96 = mrSges(2,1) * g(3) + mrSges(2,2) * g(1) - pkin(1) * t100 + qJ(2) * t102 + t86 * t37 + t85 * t40;
t42 = t86 * t44 + t85 * t45;
t95 = mrSges(2,2) * t82 - qJ(2) * t42 - t85 * t37 + t86 * t40;
t94 = -mrSges(4,2) * t66 + qJ(4) * (-t91 * mrSges(5,1) + t105) + pkin(3) * t99 + mrSges(4,1) * t65 + Ifges(4,3) * qJDD(3) + t97;
t93 = mrSges(3,1) * t71 - mrSges(3,2) * t72 + pkin(2) * t47 + t94;
t92 = mrSges(2,1) * t82 + pkin(1) * t42 + t93;
t1 = [mrSges(1,3) * g(2) - t98 * g(3) + qJ(1) * t100 + t95, -mrSges(2,3) * g(3) + t95, t40, t49, -mrSges(5,3) * t81 + t103; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t96, -mrSges(2,3) * g(1) - t92, t37, t48, t97; -mrSges(1,1) * g(2) + t98 * g(1) - qJ(1) * t102 + t92, t96, t93, t94, mrSges(5,1) * t81 - t91 * Ifges(5,4) + Ifges(5,6) * qJDD(3) - t107;];
m_new  = t1;
