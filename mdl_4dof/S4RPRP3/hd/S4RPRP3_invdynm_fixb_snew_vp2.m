% Calculate vector of cutting torques with Newton-Euler for
% S4RPRP3
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
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2019-12-31 16:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RPRP3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP3_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP3_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP3_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP3_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP3_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:42:33
% EndTime: 2019-12-31 16:42:35
% DurationCPUTime: 0.88s
% Computational Cost: add. (6743->189), mult. (12764->235), div. (0->0), fcn. (5708->6), ass. (0->74)
t167 = sin(qJ(1));
t169 = cos(qJ(1));
t151 = t167 * g(1) - t169 * g(2);
t138 = qJDD(1) * pkin(1) + t151;
t152 = -t169 * g(1) - t167 * g(2);
t170 = qJD(1) ^ 2;
t141 = -t170 * pkin(1) + t152;
t164 = sin(pkin(6));
t165 = cos(pkin(6));
t113 = t164 * t138 + t165 * t141;
t110 = -t170 * pkin(2) + qJDD(1) * pkin(5) + t113;
t163 = -g(3) + qJDD(2);
t168 = cos(qJ(3));
t154 = t168 * t163;
t166 = sin(qJ(3));
t106 = -t166 * t110 + t154;
t107 = t168 * t110 + t166 * t163;
t124 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t166 + Ifges(4,2) * t168) * qJD(1);
t125 = Ifges(5,5) * qJD(3) + (Ifges(5,1) * t166 + Ifges(5,4) * t168) * qJD(1);
t126 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t166 + Ifges(4,4) * t168) * qJD(1);
t186 = qJD(1) * qJD(3);
t182 = t168 * t186;
t142 = t166 * qJDD(1) + t182;
t143 = t168 * qJDD(1) - t166 * t186;
t185 = qJD(1) * qJD(4);
t193 = pkin(3) * t170;
t102 = qJDD(3) * pkin(3) + t154 + (-t142 + t182) * qJ(4) + (t168 * t193 - t110 - 0.2e1 * t185) * t166;
t188 = qJD(1) * t166;
t146 = qJD(3) * pkin(3) - qJ(4) * t188;
t162 = t168 ^ 2;
t103 = t143 * qJ(4) - qJD(3) * t146 - t162 * t193 + 0.2e1 * t168 * t185 + t107;
t123 = Ifges(5,6) * qJD(3) + (Ifges(5,4) * t166 + Ifges(5,2) * t168) * qJD(1);
t177 = -mrSges(5,1) * t102 + mrSges(5,2) * t103 - Ifges(5,5) * t142 - Ifges(5,6) * t143 - Ifges(5,3) * qJDD(3) - t123 * t188;
t139 = (-mrSges(5,1) * t168 + mrSges(5,2) * t166) * qJD(1);
t187 = qJD(1) * t168;
t149 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t187;
t184 = m(5) * t102 + qJDD(3) * mrSges(5,1) + qJD(3) * t149;
t97 = -t142 * mrSges(5,3) - t139 * t188 + t184;
t195 = mrSges(4,1) * t106 - mrSges(4,2) * t107 + Ifges(4,5) * t142 + Ifges(4,6) * t143 + Ifges(4,3) * qJDD(3) + pkin(3) * t97 - (-t166 * t124 + (t125 + t126) * t168) * qJD(1) - t177;
t192 = -mrSges(4,2) - mrSges(5,2);
t140 = (-mrSges(4,1) * t168 + mrSges(4,2) * t166) * qJD(1);
t150 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t187;
t94 = m(4) * t106 + qJDD(3) * mrSges(4,1) + qJD(3) * t150 + (-mrSges(4,3) - mrSges(5,3)) * t142 + (-t139 - t140) * t188 + t184;
t183 = m(5) * t103 + t143 * mrSges(5,3) + t139 * t187;
t147 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t188;
t189 = -qJD(3) * mrSges(4,1) + mrSges(4,3) * t188 - t147;
t95 = m(4) * t107 + t143 * mrSges(4,3) + t189 * qJD(3) + t192 * qJDD(3) + t140 * t187 + t183;
t180 = -t166 * t94 + t168 * t95;
t86 = m(3) * t113 - t170 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t180;
t112 = t165 * t138 - t164 * t141;
t178 = -qJDD(1) * pkin(2) - t112;
t109 = -t170 * pkin(5) + t178;
t105 = t146 * t188 - t143 * pkin(3) + qJDD(4) + (-qJ(4) * t162 - pkin(5)) * t170 + t178;
t179 = -m(5) * t105 + t143 * mrSges(5,1) + t149 * t187;
t172 = -m(4) * t109 + t143 * mrSges(4,1) + t192 * t142 + t150 * t187 + t189 * t188 + t179;
t90 = m(3) * t112 + qJDD(1) * mrSges(3,1) - t170 * mrSges(3,2) + t172;
t79 = t164 * t86 + t165 * t90;
t88 = t166 * t95 + t168 * t94;
t181 = -t164 * t90 + t165 * t86;
t176 = -mrSges(5,1) * t105 + mrSges(5,3) * t103 + Ifges(5,4) * t142 + Ifges(5,2) * t143 + Ifges(5,6) * qJDD(3) + qJD(3) * t125;
t121 = Ifges(5,3) * qJD(3) + (Ifges(5,5) * t166 + Ifges(5,6) * t168) * qJD(1);
t175 = mrSges(5,2) * t105 - mrSges(5,3) * t102 + Ifges(5,1) * t142 + Ifges(5,4) * t143 + Ifges(5,5) * qJDD(3) + t121 * t187;
t122 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t166 + Ifges(4,6) * t168) * qJD(1);
t81 = Ifges(4,4) * t142 + Ifges(4,2) * t143 + Ifges(4,6) * qJDD(3) + qJD(3) * t126 - mrSges(4,1) * t109 + mrSges(4,3) * t107 - pkin(3) * (t142 * mrSges(5,2) - t179) + qJ(4) * (-qJDD(3) * mrSges(5,2) - qJD(3) * t147 + t183) + (-pkin(3) * t147 - t121 - t122) * t188 + t176;
t83 = t122 * t187 + mrSges(4,2) * t109 - mrSges(4,3) * t106 + Ifges(4,1) * t142 + Ifges(4,4) * t143 + Ifges(4,5) * qJDD(3) - qJ(4) * t97 + (-t123 - t124) * qJD(3) + t175;
t174 = mrSges(3,1) * t112 - mrSges(3,2) * t113 + Ifges(3,3) * qJDD(1) + pkin(2) * t172 + pkin(5) * t180 + t166 * t83 + t168 * t81;
t173 = mrSges(2,1) * t151 - mrSges(2,2) * t152 + Ifges(2,3) * qJDD(1) + pkin(1) * t79 + t174;
t77 = m(2) * t152 - t170 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t181;
t76 = m(2) * t151 + qJDD(1) * mrSges(2,1) - t170 * mrSges(2,2) + t79;
t75 = -mrSges(3,1) * t163 + mrSges(3,3) * t113 + t170 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t88 - t195;
t74 = mrSges(3,2) * t163 - mrSges(3,3) * t112 + Ifges(3,5) * qJDD(1) - t170 * Ifges(3,6) - pkin(5) * t88 - t166 * t81 + t168 * t83;
t73 = -mrSges(2,2) * g(3) - mrSges(2,3) * t151 + Ifges(2,5) * qJDD(1) - t170 * Ifges(2,6) - qJ(2) * t79 - t164 * t75 + t165 * t74;
t72 = Ifges(2,6) * qJDD(1) + t170 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t152 + t164 * t74 + t165 * t75 - pkin(1) * (m(3) * t163 + t88) + qJ(2) * t181;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t169 * t73 - t167 * t72 - pkin(4) * (t167 * t77 + t169 * t76), t73, t74, t83, -qJD(3) * t123 + t175; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t167 * t73 + t169 * t72 + pkin(4) * (-t167 * t76 + t169 * t77), t72, t75, t81, -t121 * t188 + t176; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t173, t173, t174, t195, -t125 * t187 - t177;];
m_new = t1;
