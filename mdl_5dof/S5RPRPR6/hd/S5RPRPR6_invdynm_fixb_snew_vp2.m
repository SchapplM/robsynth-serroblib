% Calculate vector of cutting torques with Newton-Euler for
% S5RPRPR6
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRPR6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR6_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR6_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR6_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR6_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR6_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:17:44
% EndTime: 2019-12-31 18:17:46
% DurationCPUTime: 1.35s
% Computational Cost: add. (20954->180), mult. (27659->219), div. (0->0), fcn. (13400->8), ass. (0->81)
t191 = -pkin(3) - pkin(7);
t190 = mrSges(4,1) - mrSges(5,2);
t189 = -Ifges(5,4) + Ifges(4,5);
t188 = Ifges(5,5) - Ifges(4,6);
t162 = sin(pkin(8));
t163 = cos(pkin(8));
t166 = sin(qJ(1));
t169 = cos(qJ(1));
t145 = t166 * g(1) - t169 * g(2);
t140 = qJDD(1) * pkin(1) + t145;
t146 = -t169 * g(1) - t166 * g(2);
t170 = qJD(1) ^ 2;
t141 = -t170 * pkin(1) + t146;
t122 = t163 * t140 - t162 * t141;
t165 = sin(qJ(3));
t168 = cos(qJ(3));
t119 = qJDD(1) * pkin(2) + t122;
t123 = t162 * t140 + t163 * t141;
t120 = -t170 * pkin(2) + t123;
t114 = t168 * t119 - t165 * t120;
t157 = qJD(1) + qJD(3);
t155 = t157 ^ 2;
t156 = qJDD(1) + qJDD(3);
t180 = -t155 * qJ(4) + qJDD(4) - t114;
t112 = -t156 * pkin(3) + t180;
t109 = t191 * t156 + t180;
t161 = -g(3) + qJDD(2);
t164 = sin(qJ(5));
t167 = cos(qJ(5));
t105 = t167 * t109 - t164 * t161;
t134 = (mrSges(6,1) * t164 + mrSges(6,2) * t167) * t157;
t185 = qJD(5) * t157;
t136 = t167 * t156 - t164 * t185;
t187 = t157 * t164;
t142 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t187;
t186 = t157 * t167;
t102 = m(6) * t105 + qJDD(5) * mrSges(6,1) - t136 * mrSges(6,3) + qJD(5) * t142 - t134 * t186;
t106 = t164 * t109 + t167 * t161;
t135 = -t164 * t156 - t167 * t185;
t143 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t186;
t103 = m(6) * t106 - qJDD(5) * mrSges(6,2) + t135 * mrSges(6,3) - qJD(5) * t143 - t134 * t187;
t94 = t167 * t102 + t164 * t103;
t179 = -m(5) * t112 + t155 * mrSges(5,3) - t94;
t88 = m(4) * t114 - t155 * mrSges(4,2) + t190 * t156 + t179;
t115 = t165 * t119 + t168 * t120;
t181 = t156 * qJ(4) + 0.2e1 * qJD(4) * t157 + t115;
t108 = t191 * t155 + t181;
t100 = -m(6) * t108 + t135 * mrSges(6,1) - t136 * mrSges(6,2) - t142 * t187 - t143 * t186;
t110 = t155 * pkin(3) - t181;
t176 = -m(5) * t110 + t155 * mrSges(5,2) + t156 * mrSges(5,3) - t100;
t92 = m(4) * t115 - t155 * mrSges(4,1) - t156 * mrSges(4,2) + t176;
t86 = t165 * t92 + t168 * t88;
t83 = m(3) * t122 + qJDD(1) * mrSges(3,1) - t170 * mrSges(3,2) + t86;
t183 = -t165 * t88 + t168 * t92;
t84 = m(3) * t123 - t170 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t183;
t77 = t162 * t84 + t163 * t83;
t184 = -t162 * t83 + t163 * t84;
t95 = -t164 * t102 + t167 * t103;
t93 = m(5) * t161 + t95;
t182 = m(4) * t161 + t93;
t127 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t167 - Ifges(6,2) * t164) * t157;
t128 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t167 - Ifges(6,4) * t164) * t157;
t178 = mrSges(6,1) * t105 - mrSges(6,2) * t106 + Ifges(6,5) * t136 + Ifges(6,6) * t135 + Ifges(6,3) * qJDD(5) + t127 * t186 + t128 * t187;
t126 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t167 - Ifges(6,6) * t164) * t157;
t97 = -mrSges(6,1) * t108 + mrSges(6,3) * t106 + Ifges(6,4) * t136 + Ifges(6,2) * t135 + Ifges(6,6) * qJDD(5) + qJD(5) * t128 - t126 * t186;
t98 = mrSges(6,2) * t108 - mrSges(6,3) * t105 + Ifges(6,1) * t136 + Ifges(6,4) * t135 + Ifges(6,5) * qJDD(5) - qJD(5) * t127 - t126 * t187;
t177 = mrSges(5,2) * t112 - mrSges(5,3) * t110 + Ifges(5,1) * t156 - pkin(7) * t94 - t164 * t97 + t167 * t98;
t175 = -mrSges(5,1) * t110 - pkin(4) * t100 - pkin(7) * t95 - t164 * t98 - t167 * t97;
t174 = mrSges(5,1) * t112 + pkin(4) * t94 + t178;
t173 = -mrSges(4,2) * t115 + mrSges(4,1) * t114 + Ifges(4,3) * t156 + t177 + pkin(3) * (-t156 * mrSges(5,2) + t179) + qJ(4) * t176;
t172 = mrSges(3,1) * t122 - mrSges(3,2) * t123 + Ifges(3,3) * qJDD(1) + pkin(2) * t86 + t173;
t171 = mrSges(2,1) * t145 - mrSges(2,2) * t146 + Ifges(2,3) * qJDD(1) + pkin(1) * t77 + t172;
t79 = -mrSges(4,3) * t114 - qJ(4) * t93 + (mrSges(4,2) - mrSges(5,3)) * t161 + t189 * t156 + t188 * t155 + t174;
t78 = mrSges(4,3) * t115 - pkin(3) * t93 + t189 * t155 - t188 * t156 - t190 * t161 + t175;
t75 = m(2) * t146 - t170 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t184;
t74 = m(2) * t145 + qJDD(1) * mrSges(2,1) - t170 * mrSges(2,2) + t77;
t73 = mrSges(3,2) * t161 - mrSges(3,3) * t122 + Ifges(3,5) * qJDD(1) - t170 * Ifges(3,6) - pkin(6) * t86 - t165 * t78 + t168 * t79;
t72 = -mrSges(3,1) * t161 + mrSges(3,3) * t123 + t170 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t182 + pkin(6) * t183 + t165 * t79 + t168 * t78;
t71 = -mrSges(2,2) * g(3) - mrSges(2,3) * t145 + Ifges(2,5) * qJDD(1) - t170 * Ifges(2,6) - qJ(2) * t77 - t162 * t72 + t163 * t73;
t70 = Ifges(2,6) * qJDD(1) + t170 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t146 + t162 * t73 + t163 * t72 - pkin(1) * (m(3) * t161 + t182) + qJ(2) * t184;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t169 * t71 - t166 * t70 - pkin(5) * (t166 * t75 + t169 * t74), t71, t73, t79, t177, t98; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t166 * t71 + t169 * t70 + pkin(5) * (-t166 * t74 + t169 * t75), t70, t72, t78, mrSges(5,3) * t161 + Ifges(5,4) * t156 - t155 * Ifges(5,5) - t174, t97; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t171, t171, t172, t173, -mrSges(5,2) * t161 + t155 * Ifges(5,4) + Ifges(5,5) * t156 - t175, t178;];
m_new = t1;
