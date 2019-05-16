% Calculate vector of cutting torques with Newton-Euler for
% S4RPPR2
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
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2019-05-04 19:11
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RPPR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR2_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR2_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR2_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR2_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR2_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:10:44
% EndTime: 2019-05-04 19:10:44
% DurationCPUTime: 0.45s
% Computational Cost: add. (6270->120), mult. (9189->138), div. (0->0), fcn. (3026->6), ass. (0->54)
t127 = -pkin(1) - pkin(2);
t126 = mrSges(2,1) + mrSges(3,1);
t109 = sin(qJ(4));
t111 = cos(qJ(4));
t107 = sin(pkin(6));
t108 = cos(pkin(6));
t114 = qJD(1) ^ 2;
t110 = sin(qJ(1));
t112 = cos(qJ(1));
t90 = -t112 * g(1) - t110 * g(2);
t122 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t90;
t79 = t127 * t114 + t122;
t89 = t110 * g(1) - t112 * g(2);
t121 = -t114 * qJ(2) + qJDD(2) - t89;
t82 = t127 * qJDD(1) + t121;
t74 = -t107 * t79 + t108 * t82;
t71 = -qJDD(1) * pkin(3) + t74;
t75 = t107 * t82 + t108 * t79;
t72 = -t114 * pkin(3) + t75;
t69 = -t109 * t72 + t111 * t71;
t96 = -qJD(1) + qJD(4);
t94 = t96 ^ 2;
t95 = -qJDD(1) + qJDD(4);
t65 = m(5) * t69 + t95 * mrSges(5,1) - (t94 * mrSges(5,2));
t70 = t109 * t71 + t111 * t72;
t66 = m(5) * t70 - t94 * mrSges(5,1) - t95 * mrSges(5,2);
t59 = t109 * t66 + t111 * t65;
t125 = mrSges(5,1) * t69 - mrSges(5,2) * t70 + Ifges(5,3) * t95;
t104 = g(3) + qJDD(3);
t88 = (-m(4) - m(5)) * t104;
t57 = m(4) * t74 - qJDD(1) * mrSges(4,1) - t114 * mrSges(4,2) + t59;
t124 = -t109 * t65 + t111 * t66;
t58 = m(4) * t75 - t114 * mrSges(4,1) + qJDD(1) * mrSges(4,2) + t124;
t52 = -t107 * t57 + t108 * t58;
t51 = t107 * t58 + t108 * t57;
t83 = -t114 * pkin(1) + t122;
t123 = m(3) * t83 + qJDD(1) * mrSges(3,3) + t52;
t85 = -qJDD(1) * pkin(1) + t121;
t120 = -m(3) * t85 + qJDD(1) * mrSges(3,1) + t114 * mrSges(3,3) - t51;
t119 = mrSges(4,1) * t74 - mrSges(4,2) * t75 - Ifges(4,3) * qJDD(1) + pkin(3) * t59 + t125;
t60 = -mrSges(5,1) * t104 + mrSges(5,3) * t70 + (t94 * Ifges(5,5)) + Ifges(5,6) * t95;
t61 = mrSges(5,2) * t104 - mrSges(5,3) * t69 + Ifges(5,5) * t95 - t94 * Ifges(5,6);
t53 = -Ifges(4,6) * qJDD(1) + t114 * Ifges(4,5) + mrSges(4,3) * t75 + t109 * t61 + t111 * t60 + pkin(5) * t124 + (-pkin(3) * m(5) - mrSges(4,1)) * t104;
t55 = mrSges(4,2) * t104 - mrSges(4,3) * t74 - Ifges(4,5) * qJDD(1) - t114 * Ifges(4,6) - pkin(5) * t59 - t109 * t60 + t111 * t61;
t118 = mrSges(3,2) * t85 + mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) + t114 * Ifges(3,6) - qJ(3) * t51 - t107 * t53 + t108 * t55;
t117 = mrSges(3,2) * t83 - pkin(2) * t88 - qJ(3) * t52 - t107 * t55 - t108 * t53;
t116 = -mrSges(3,1) * t85 + mrSges(3,3) * t83 + Ifges(3,2) * qJDD(1) - pkin(2) * t51 - t119;
t115 = -mrSges(2,2) * t90 + Ifges(2,3) * qJDD(1) + t116 + qJ(2) * (-t114 * mrSges(3,1) + t123) + pkin(1) * t120 + mrSges(2,1) * t89;
t86 = -m(3) * g(3) + t88;
t48 = m(2) * t89 + qJDD(1) * mrSges(2,1) - t114 * mrSges(2,2) + t120;
t47 = m(2) * t90 - qJDD(1) * mrSges(2,2) - t126 * t114 + t123;
t46 = -mrSges(2,2) * g(3) - mrSges(2,3) * t89 + Ifges(2,5) * qJDD(1) - t114 * Ifges(2,6) - qJ(2) * t86 + t118;
t45 = mrSges(2,3) * t90 - pkin(1) * t86 + (Ifges(3,4) + Ifges(2,5)) * t114 + (Ifges(2,6) - Ifges(3,6)) * qJDD(1) + t126 * g(3) + t117;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t112 * t46 - t110 * t45 - pkin(4) * (t110 * t47 + t112 * t48), t46, t118, t55, t61; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t110 * t46 + t112 * t45 + pkin(4) * (-t110 * t48 + t112 * t47), t45, t116, t53, t60; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t115, t115, -mrSges(3,1) * g(3) - t114 * Ifges(3,4) + Ifges(3,6) * qJDD(1) - t117, t119, t125;];
m_new  = t1;
