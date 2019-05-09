% Calculate vector of cutting torques with Newton-Euler for
% S4PPRR1
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
% Datum: 2019-05-04 18:48
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4PPRR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR1_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR1_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR1_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR1_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR1_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR1_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR1_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:47:36
% EndTime: 2019-05-04 18:47:36
% DurationCPUTime: 0.29s
% Computational Cost: add. (3378->94), mult. (4618->112), div. (0->0), fcn. (3026->6), ass. (0->47)
t102 = -m(4) - m(5);
t86 = sin(pkin(6));
t87 = cos(pkin(6));
t79 = t86 * g(1) - t87 * g(2);
t76 = qJDD(2) - t79;
t80 = -t87 * g(1) - t86 * g(2);
t89 = sin(qJ(3));
t91 = cos(qJ(3));
t67 = t91 * t76 - t89 * t80;
t64 = qJDD(3) * pkin(3) + t67;
t68 = t89 * t76 + t91 * t80;
t92 = qJD(3) ^ 2;
t65 = -t92 * pkin(3) + t68;
t88 = sin(qJ(4));
t90 = cos(qJ(4));
t62 = t90 * t64 - t88 * t65;
t84 = qJD(3) + qJD(4);
t82 = t84 ^ 2;
t83 = qJDD(3) + qJDD(4);
t58 = m(5) * t62 + t83 * mrSges(5,1) - t82 * mrSges(5,2);
t63 = t88 * t64 + t90 * t65;
t59 = m(5) * t63 - t82 * mrSges(5,1) - t83 * mrSges(5,2);
t52 = t90 * t58 + t88 * t59;
t101 = mrSges(5,1) * t62 - mrSges(5,2) * t63 + Ifges(5,3) * t83;
t50 = m(4) * t67 + qJDD(3) * mrSges(4,1) - t92 * mrSges(4,2) + t52;
t100 = -t88 * t58 + t90 * t59;
t51 = m(4) * t68 - t92 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t100;
t45 = -t89 * t50 + t91 * t51;
t99 = m(3) * t80 + t45;
t44 = t91 * t50 + t89 * t51;
t85 = g(3) - qJDD(1);
t53 = -mrSges(5,1) * t85 + mrSges(5,3) * t63 + t82 * Ifges(5,5) + Ifges(5,6) * t83;
t54 = mrSges(5,2) * t85 - mrSges(5,3) * t62 + Ifges(5,5) * t83 - t82 * Ifges(5,6);
t46 = Ifges(4,6) * qJDD(3) + t92 * Ifges(4,5) + mrSges(4,3) * t68 + t88 * t54 + t90 * t53 + pkin(5) * t100 + (-pkin(3) * m(5) - mrSges(4,1)) * t85;
t48 = mrSges(4,2) * t85 - mrSges(4,3) * t67 + Ifges(4,5) * qJDD(3) - t92 * Ifges(4,6) - pkin(5) * t52 - t88 * t53 + t90 * t54;
t98 = mrSges(3,2) * t76 - pkin(4) * t44 - t89 * t46 + t91 * t48;
t97 = -m(3) * t76 - t44;
t96 = -pkin(2) * t102 * t85 - pkin(4) * t45 - t91 * t46 - t89 * t48;
t95 = mrSges(4,1) * t67 - mrSges(4,2) * t68 + Ifges(4,3) * qJDD(3) + pkin(3) * t52 + t101;
t94 = -mrSges(3,1) * t76 + mrSges(3,3) * t80 - pkin(2) * t44 - t95;
t93 = mrSges(2,1) * t79 - mrSges(2,2) * t80 + pkin(1) * t97 + qJ(2) * t99 + t94;
t69 = (-m(3) + t102) * t85;
t41 = m(2) * t80 + t99;
t40 = m(2) * t79 + t97;
t39 = -mrSges(2,3) * t79 - qJ(2) * t69 + (-mrSges(2,2) + mrSges(3,3)) * t85 + t98;
t38 = -pkin(1) * t69 + (mrSges(2,1) + mrSges(3,1)) * t85 + (mrSges(3,2) + mrSges(2,3)) * t80 + t96;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t87 * t39 - t86 * t38 - qJ(1) * (t87 * t40 + t86 * t41), t39, mrSges(3,3) * t85 + t98, t48, t54; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t86 * t39 + t87 * t38 + qJ(1) * (-t86 * t40 + t87 * t41), t38, t94, t46, t53; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t93, t93, -mrSges(3,1) * t85 - mrSges(3,2) * t80 - t96, t95, t101;];
m_new  = t1;
