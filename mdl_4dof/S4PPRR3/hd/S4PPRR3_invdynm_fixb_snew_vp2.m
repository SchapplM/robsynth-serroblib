% Calculate vector of cutting torques with Newton-Euler for
% S4PPRR3
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4PPRR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR3_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR3_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR3_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR3_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR3_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:17:23
% EndTime: 2019-12-31 16:17:24
% DurationCPUTime: 0.45s
% Computational Cost: add. (3629->124), mult. (6125->159), div. (0->0), fcn. (3300->6), ass. (0->55)
t101 = qJD(3) ^ 2;
t100 = cos(qJ(3));
t95 = sin(pkin(6));
t96 = cos(pkin(6));
t90 = g(1) * t95 - g(2) * t96;
t88 = qJDD(2) - t90;
t91 = -g(1) * t96 - g(2) * t95;
t98 = sin(qJ(3));
t74 = t100 * t91 + t98 * t88;
t71 = -pkin(3) * t101 + qJDD(3) * pkin(5) + t74;
t94 = g(3) - qJDD(1);
t97 = sin(qJ(4));
t99 = cos(qJ(4));
t68 = -t71 * t97 + t94 * t99;
t69 = t71 * t99 + t94 * t97;
t76 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t97 + Ifges(5,2) * t99) * qJD(3);
t77 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t97 + Ifges(5,4) * t99) * qJD(3);
t111 = qJD(3) * qJD(4);
t84 = qJDD(3) * t97 + t111 * t99;
t85 = qJDD(3) * t99 - t111 * t97;
t114 = mrSges(5,1) * t68 - mrSges(5,2) * t69 + Ifges(5,5) * t84 + Ifges(5,6) * t85 + Ifges(5,3) * qJDD(4) + (t76 * t97 - t77 * t99) * qJD(3);
t113 = qJD(3) * t97;
t112 = qJD(3) * t99;
t83 = (-mrSges(5,1) * t99 + mrSges(5,2) * t97) * qJD(3);
t93 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t112;
t65 = m(5) * t68 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t84 + qJD(4) * t93 - t113 * t83;
t92 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t113;
t66 = m(5) * t69 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t85 - qJD(4) * t92 + t112 * t83;
t60 = -t65 * t97 + t99 * t66;
t56 = m(4) * t74 - mrSges(4,1) * t101 - qJDD(3) * mrSges(4,2) + t60;
t73 = t100 * t88 - t91 * t98;
t70 = -qJDD(3) * pkin(3) - pkin(5) * t101 - t73;
t67 = -m(5) * t70 + t85 * mrSges(5,1) - mrSges(5,2) * t84 + t93 * t112 - t113 * t92;
t63 = m(4) * t73 + qJDD(3) * mrSges(4,1) - mrSges(4,2) * t101 + t67;
t54 = t100 * t56 - t63 * t98;
t110 = m(3) * t91 + t54;
t59 = t99 * t65 + t97 * t66;
t53 = t100 * t63 + t56 * t98;
t75 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t97 + Ifges(5,6) * t99) * qJD(3);
t61 = -mrSges(5,1) * t70 + mrSges(5,3) * t69 + Ifges(5,4) * t84 + Ifges(5,2) * t85 + Ifges(5,6) * qJDD(4) + qJD(4) * t77 - t113 * t75;
t62 = mrSges(5,2) * t70 - mrSges(5,3) * t68 + Ifges(5,1) * t84 + Ifges(5,4) * t85 + Ifges(5,5) * qJDD(4) - qJD(4) * t76 + t112 * t75;
t47 = mrSges(4,2) * t94 - mrSges(4,3) * t73 + Ifges(4,5) * qJDD(3) - Ifges(4,6) * t101 - pkin(5) * t59 - t61 * t97 + t62 * t99;
t48 = -mrSges(4,1) * t94 + mrSges(4,3) * t74 + t101 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t59 - t114;
t108 = mrSges(3,2) * t88 - pkin(4) * t53 + t100 * t47 - t98 * t48;
t107 = -m(3) * t88 - t53;
t106 = -pkin(2) * (-m(4) * t94 - t59) - pkin(4) * t54 - t100 * t48 - t98 * t47;
t104 = mrSges(4,1) * t73 - mrSges(4,2) * t74 + Ifges(4,3) * qJDD(3) + pkin(3) * t67 + pkin(5) * t60 + t61 * t99 + t62 * t97;
t103 = -mrSges(3,1) * t88 + mrSges(3,3) * t91 - pkin(2) * t53 - t104;
t102 = mrSges(2,1) * t90 - mrSges(2,2) * t91 + pkin(1) * t107 + qJ(2) * t110 + t103;
t57 = (-m(3) - m(4)) * t94 - t59;
t50 = m(2) * t91 + t110;
t49 = m(2) * t90 + t107;
t45 = -mrSges(2,3) * t90 - qJ(2) * t57 + (-mrSges(2,2) + mrSges(3,3)) * t94 + t108;
t44 = -pkin(1) * t57 + (mrSges(2,1) + mrSges(3,1)) * t94 + (mrSges(3,2) + mrSges(2,3)) * t91 + t106;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t96 * t45 - t95 * t44 - qJ(1) * (t49 * t96 + t50 * t95), t45, mrSges(3,3) * t94 + t108, t47, t62; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t95 * t45 + t96 * t44 + qJ(1) * (-t49 * t95 + t50 * t96), t44, t103, t48, t61; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t102, t102, -mrSges(3,1) * t94 - mrSges(3,2) * t91 - t106, t104, t114;];
m_new = t1;
