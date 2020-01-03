% Calculate vector of cutting torques with Newton-Euler for
% S5PPRPR4
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PPRPR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR4_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR4_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR4_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR4_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR4_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:32:19
% EndTime: 2019-12-31 17:32:20
% DurationCPUTime: 1.65s
% Computational Cost: add. (18389->179), mult. (37304->224), div. (0->0), fcn. (23645->8), ass. (0->85)
t167 = qJD(3) ^ 2;
t160 = sin(pkin(7));
t162 = cos(pkin(7));
t146 = t160 * g(1) - t162 * g(2);
t144 = qJDD(2) - t146;
t147 = -t162 * g(1) - t160 * g(2);
t164 = sin(qJ(3));
t166 = cos(qJ(3));
t129 = t164 * t144 + t166 * t147;
t124 = -t167 * pkin(3) + qJDD(3) * qJ(4) + t129;
t159 = sin(pkin(8));
t158 = g(3) - qJDD(1);
t161 = cos(pkin(8));
t189 = qJD(3) * qJD(4);
t191 = t161 * t158 - 0.2e1 * t159 * t189;
t112 = -t159 * t124 + t191;
t113 = t159 * t158 + (t124 + 0.2e1 * t189) * t161;
t194 = pkin(4) * t161;
t108 = (-pkin(6) * qJDD(3) + t167 * t194 - t124) * t159 + t191;
t188 = qJDD(3) * t161;
t156 = t161 ^ 2;
t192 = t156 * t167;
t109 = -pkin(4) * t192 + pkin(6) * t188 + t113;
t163 = sin(qJ(5));
t165 = cos(qJ(5));
t106 = t165 * t108 - t163 * t109;
t107 = t163 * t108 + t165 * t109;
t179 = -t159 * t163 + t161 * t165;
t132 = t179 * qJD(3);
t180 = t159 * t165 + t161 * t163;
t133 = t180 * qJD(3);
t115 = Ifges(6,4) * t133 + Ifges(6,2) * t132 + Ifges(6,6) * qJD(5);
t116 = Ifges(6,1) * t133 + Ifges(6,4) * t132 + Ifges(6,5) * qJD(5);
t125 = -t133 * qJD(5) + t179 * qJDD(3);
t126 = t132 * qJD(5) + t180 * qJDD(3);
t173 = -mrSges(6,1) * t106 + mrSges(6,2) * t107 - Ifges(6,5) * t126 - Ifges(6,6) * t125 - Ifges(6,3) * qJDD(5) - t133 * t115 + t132 * t116;
t184 = Ifges(5,4) * t159 + Ifges(5,2) * t161;
t185 = Ifges(5,1) * t159 + Ifges(5,4) * t161;
t120 = -t132 * mrSges(6,1) + t133 * mrSges(6,2);
t130 = -qJD(5) * mrSges(6,2) + t132 * mrSges(6,3);
t103 = m(6) * t106 + qJDD(5) * mrSges(6,1) - t126 * mrSges(6,3) + qJD(5) * t130 - t133 * t120;
t131 = qJD(5) * mrSges(6,1) - t133 * mrSges(6,3);
t104 = m(6) * t107 - qJDD(5) * mrSges(6,2) + t125 * mrSges(6,3) - qJD(5) * t131 + t132 * t120;
t95 = t165 * t103 + t163 * t104;
t195 = -mrSges(5,1) * t112 + mrSges(5,2) * t113 - pkin(4) * t95 - (t159 * t184 - t161 * t185) * t167 + t173;
t193 = mrSges(5,2) * t159;
t183 = Ifges(5,5) * t159 + Ifges(5,6) * t161;
t190 = t167 * t183;
t178 = mrSges(5,3) * qJDD(3) + t167 * (-mrSges(5,1) * t161 + t193);
t93 = m(5) * t112 - t178 * t159 + t95;
t187 = -t163 * t103 + t165 * t104;
t94 = m(5) * t113 + t178 * t161 + t187;
t91 = -t159 * t93 + t161 * t94;
t87 = m(4) * t129 - t167 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t91;
t128 = t166 * t144 - t164 * t147;
t182 = qJDD(4) - t128;
t119 = -qJDD(3) * pkin(3) - t167 * qJ(4) + t182;
t155 = t159 ^ 2;
t111 = (-pkin(3) - t194) * qJDD(3) + (-qJ(4) + (-t155 - t156) * pkin(6)) * t167 + t182;
t174 = m(6) * t111 - t125 * mrSges(6,1) + t126 * mrSges(6,2) - t132 * t130 + t133 * t131;
t172 = -m(5) * t119 + mrSges(5,1) * t188 - t174 + (t155 * t167 + t192) * mrSges(5,3);
t98 = m(4) * t128 - t167 * mrSges(4,2) + (mrSges(4,1) - t193) * qJDD(3) + t172;
t84 = -t164 * t98 + t166 * t87;
t186 = m(3) * t147 + t84;
t90 = t159 * t94 + t161 * t93;
t83 = t164 * t87 + t166 * t98;
t114 = Ifges(6,5) * t133 + Ifges(6,6) * t132 + Ifges(6,3) * qJD(5);
t96 = -mrSges(6,1) * t111 + mrSges(6,3) * t107 + Ifges(6,4) * t126 + Ifges(6,2) * t125 + Ifges(6,6) * qJDD(5) + qJD(5) * t116 - t133 * t114;
t97 = mrSges(6,2) * t111 - mrSges(6,3) * t106 + Ifges(6,1) * t126 + Ifges(6,4) * t125 + Ifges(6,5) * qJDD(5) - qJD(5) * t115 + t132 * t114;
t78 = -mrSges(5,1) * t119 + mrSges(5,3) * t113 - pkin(4) * t174 + pkin(6) * t187 + t184 * qJDD(3) - t159 * t190 + t163 * t97 + t165 * t96;
t85 = mrSges(5,2) * t119 - mrSges(5,3) * t112 - pkin(6) * t95 + t185 * qJDD(3) + t161 * t190 - t163 * t96 + t165 * t97;
t76 = mrSges(4,2) * t158 - mrSges(4,3) * t128 + Ifges(4,5) * qJDD(3) - t167 * Ifges(4,6) - qJ(4) * t90 - t159 * t78 + t161 * t85;
t77 = (Ifges(4,6) - t183) * qJDD(3) + t167 * Ifges(4,5) - mrSges(4,1) * t158 + mrSges(4,3) * t129 - pkin(3) * t90 + t195;
t177 = mrSges(3,2) * t144 - pkin(5) * t83 - t164 * t77 + t166 * t76;
t176 = -m(3) * t144 - t83;
t175 = -pkin(2) * (-m(4) * t158 - t90) - pkin(5) * t84 - t164 * t76 - t166 * t77;
t171 = mrSges(4,1) * t128 + Ifges(4,3) * qJDD(3) + pkin(3) * (-qJDD(3) * t193 + t172) + qJ(4) * t91 + t159 * t85 + t161 * t78 - mrSges(4,2) * t129;
t169 = -mrSges(3,1) * t144 + mrSges(3,3) * t147 - pkin(2) * t83 - t171;
t168 = mrSges(2,1) * t146 - mrSges(2,2) * t147 + pkin(1) * t176 + qJ(2) * t186 + t169;
t88 = (-m(3) - m(4)) * t158 - t90;
t80 = m(2) * t147 + t186;
t79 = m(2) * t146 + t176;
t74 = -mrSges(2,3) * t146 - qJ(2) * t88 + (-mrSges(2,2) + mrSges(3,3)) * t158 + t177;
t73 = -pkin(1) * t88 + (mrSges(2,1) + mrSges(3,1)) * t158 + (mrSges(3,2) + mrSges(2,3)) * t147 + t175;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t162 * t74 - t160 * t73 - qJ(1) * (t160 * t80 + t162 * t79), t74, mrSges(3,3) * t158 + t177, t76, t85, t97; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t160 * t74 + t162 * t73 + qJ(1) * (-t160 * t79 + t162 * t80), t73, t169, t77, t78, t96; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t168, t168, -mrSges(3,1) * t158 - mrSges(3,2) * t147 - t175, t171, t183 * qJDD(3) - t195, -t173;];
m_new = t1;
