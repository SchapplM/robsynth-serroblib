% Calculate vector of cutting torques with Newton-Euler for
% S5PPRPR5
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% m [3x6]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PPRPR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR5_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR5_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR5_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR5_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR5_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:26
% EndTime: 2019-12-31 17:33:27
% DurationCPUTime: 0.60s
% Computational Cost: add. (6096->156), mult. (9620->184), div. (0->0), fcn. (4768->6), ass. (0->70)
t163 = -m(4) - m(5);
t162 = -pkin(3) - pkin(6);
t161 = mrSges(4,1) - mrSges(5,2);
t160 = -Ifges(5,4) + Ifges(4,5);
t159 = Ifges(5,5) - Ifges(4,6);
t137 = sin(qJ(5));
t139 = cos(qJ(5));
t118 = (mrSges(6,1) * t137 + mrSges(6,2) * t139) * qJD(3);
t156 = qJD(3) * qJD(5);
t120 = t139 * qJDD(3) - t137 * t156;
t158 = qJD(3) * t137;
t127 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t158;
t157 = qJD(3) * t139;
t132 = g(3) - qJDD(1);
t135 = sin(pkin(7));
t136 = cos(pkin(7));
t125 = t135 * g(1) - t136 * g(2);
t123 = qJDD(2) - t125;
t126 = -t136 * g(1) - t135 * g(2);
t138 = sin(qJ(3));
t140 = cos(qJ(3));
t102 = t140 * t123 - t138 * t126;
t141 = qJD(3) ^ 2;
t151 = -t141 * qJ(4) + qJDD(4) - t102;
t97 = t162 * qJDD(3) + t151;
t92 = -t137 * t132 + t139 * t97;
t89 = m(6) * t92 + qJDD(5) * mrSges(6,1) - t120 * mrSges(6,3) + qJD(5) * t127 - t118 * t157;
t119 = -t137 * qJDD(3) - t139 * t156;
t128 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t157;
t93 = t139 * t132 + t137 * t97;
t90 = m(6) * t93 - qJDD(5) * mrSges(6,2) + t119 * mrSges(6,3) - qJD(5) * t128 - t118 * t158;
t78 = -t137 * t89 + t139 * t90;
t103 = t138 * t123 + t140 * t126;
t100 = -qJDD(3) * pkin(3) + t151;
t77 = t137 * t90 + t139 * t89;
t149 = -m(5) * t100 + t141 * mrSges(5,3) - t77;
t71 = m(4) * t102 - t141 * mrSges(4,2) + t161 * qJDD(3) + t149;
t150 = qJDD(3) * qJ(4) + 0.2e1 * qJD(4) * qJD(3) + t103;
t95 = t162 * t141 + t150;
t88 = -m(6) * t95 + t119 * mrSges(6,1) - t120 * mrSges(6,2) - t127 * t158 - t128 * t157;
t98 = t141 * pkin(3) - t150;
t85 = -m(5) * t98 + t141 * mrSges(5,2) + qJDD(3) * mrSges(5,3) - t88;
t80 = m(4) * t103 - t141 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t85;
t70 = -t138 * t71 + t140 * t80;
t155 = m(3) * t126 + t70;
t69 = t138 * t80 + t140 * t71;
t106 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t139 - Ifges(6,6) * t137) * qJD(3);
t108 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t139 - Ifges(6,4) * t137) * qJD(3);
t83 = -mrSges(6,1) * t95 + mrSges(6,3) * t93 + Ifges(6,4) * t120 + Ifges(6,2) * t119 + Ifges(6,6) * qJDD(5) + qJD(5) * t108 - t106 * t157;
t107 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t139 - Ifges(6,2) * t137) * qJD(3);
t84 = mrSges(6,2) * t95 - mrSges(6,3) * t92 + Ifges(6,1) * t120 + Ifges(6,4) * t119 + Ifges(6,5) * qJDD(5) - qJD(5) * t107 - t106 * t158;
t154 = mrSges(5,2) * t100 - mrSges(5,3) * t98 + Ifges(5,1) * qJDD(3) - pkin(6) * t77 - t137 * t83 + t139 * t84;
t146 = -mrSges(5,1) * t98 - pkin(4) * t88 - pkin(6) * t78 - t137 * t84 - t139 * t83;
t75 = m(5) * t132 + t78;
t62 = mrSges(4,3) * t103 - pkin(3) * t75 - t159 * qJDD(3) - t161 * t132 + t160 * t141 + t146;
t148 = mrSges(6,1) * t92 - mrSges(6,2) * t93 + Ifges(6,5) * t120 + Ifges(6,6) * t119 + Ifges(6,3) * qJDD(5) + t107 * t157 + t108 * t158;
t145 = mrSges(5,1) * t100 + pkin(4) * t77 + t148;
t64 = -mrSges(4,3) * t102 - qJ(4) * t75 + t159 * t141 + (mrSges(4,2) - mrSges(5,3)) * t132 + t160 * qJDD(3) + t145;
t153 = mrSges(3,2) * t123 - pkin(5) * t69 - t138 * t62 + t140 * t64;
t152 = -m(3) * t123 - t69;
t147 = -pkin(2) * (t163 * t132 - t78) - pkin(5) * t70 - t138 * t64 - t140 * t62;
t144 = mrSges(4,1) * t102 + Ifges(4,3) * qJDD(3) + pkin(3) * (-qJDD(3) * mrSges(5,2) + t149) + qJ(4) * t85 - mrSges(4,2) * t103 + t154;
t143 = -mrSges(3,1) * t123 + mrSges(3,3) * t126 - pkin(2) * t69 - t144;
t142 = mrSges(2,1) * t125 - mrSges(2,2) * t126 + pkin(1) * t152 + qJ(2) * t155 + t143;
t73 = (-m(3) + t163) * t132 - t78;
t66 = m(2) * t126 + t155;
t65 = m(2) * t125 + t152;
t61 = -mrSges(2,3) * t125 - qJ(2) * t73 + (-mrSges(2,2) + mrSges(3,3)) * t132 + t153;
t60 = -pkin(1) * t73 + (mrSges(2,1) + mrSges(3,1)) * t132 + (mrSges(3,2) + mrSges(2,3)) * t126 + t147;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t136 * t61 - t135 * t60 - qJ(1) * (t135 * t66 + t136 * t65), t61, mrSges(3,3) * t132 + t153, t64, t154, t84; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t135 * t61 + t136 * t60 + qJ(1) * (-t135 * t65 + t136 * t66), t60, t143, t62, mrSges(5,3) * t132 + Ifges(5,4) * qJDD(3) - t141 * Ifges(5,5) - t145, t83; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t142, t142, -mrSges(3,1) * t132 - mrSges(3,2) * t126 - t147, t144, -mrSges(5,2) * t132 + t141 * Ifges(5,4) + Ifges(5,5) * qJDD(3) - t146, t148;];
m_new = t1;
