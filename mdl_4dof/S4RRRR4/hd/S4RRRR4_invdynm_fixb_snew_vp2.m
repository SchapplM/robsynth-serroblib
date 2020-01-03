% Calculate vector of cutting torques with Newton-Euler for
% S4RRRR4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RRRR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR4_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR4_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR4_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR4_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR4_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR4_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR4_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:25:44
% EndTime: 2019-12-31 17:25:48
% DurationCPUTime: 2.60s
% Computational Cost: add. (30978->237), mult. (63270->306), div. (0->0), fcn. (40681->8), ass. (0->97)
t191 = sin(qJ(3));
t192 = sin(qJ(2));
t195 = cos(qJ(3));
t196 = cos(qJ(2));
t167 = (t191 * t192 - t195 * t196) * qJD(1);
t193 = sin(qJ(1));
t197 = cos(qJ(1));
t182 = -t197 * g(1) - t193 * g(2);
t198 = qJD(1) ^ 2;
t170 = -t198 * pkin(1) + qJDD(1) * pkin(5) + t182;
t215 = t192 * t170;
t158 = -t196 * g(3) - t215;
t159 = -t192 * g(3) + t196 * t170;
t165 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t192 + Ifges(3,2) * t196) * qJD(1);
t166 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t192 + Ifges(3,4) * t196) * qJD(1);
t212 = qJD(1) * qJD(2);
t175 = t192 * qJDD(1) + t196 * t212;
t176 = t196 * qJDD(1) - t192 * t212;
t168 = (t191 * t196 + t192 * t195) * qJD(1);
t145 = -t168 * qJD(3) - t191 * t175 + t195 * t176;
t146 = -t167 * qJD(3) + t195 * t175 + t191 * t176;
t214 = qJD(1) * t192;
t180 = qJD(2) * pkin(2) - pkin(6) * t214;
t189 = t196 ^ 2;
t181 = t193 * g(1) - t197 * g(2);
t206 = -qJDD(1) * pkin(1) - t181;
t147 = -t176 * pkin(2) + t180 * t214 + (-pkin(6) * t189 - pkin(5)) * t198 + t206;
t187 = qJD(2) + qJD(3);
t122 = (t167 * t187 - t146) * pkin(7) + (t168 * t187 - t145) * pkin(3) + t147;
t216 = pkin(2) * t198;
t139 = qJDD(2) * pkin(2) - t175 * pkin(6) - t215 + (pkin(6) * t212 + t192 * t216 - g(3)) * t196;
t140 = t176 * pkin(6) - qJD(2) * t180 - t189 * t216 + t159;
t127 = t191 * t139 + t195 * t140;
t155 = t167 * pkin(3) - t168 * pkin(7);
t185 = t187 ^ 2;
t186 = qJDD(2) + qJDD(3);
t124 = -t185 * pkin(3) + t186 * pkin(7) - t167 * t155 + t127;
t190 = sin(qJ(4));
t194 = cos(qJ(4));
t121 = t190 * t122 + t194 * t124;
t126 = t195 * t139 - t191 * t140;
t123 = -t186 * pkin(3) - t185 * pkin(7) + t168 * t155 - t126;
t157 = t194 * t168 + t190 * t187;
t129 = -t157 * qJD(4) - t190 * t146 + t194 * t186;
t156 = -t190 * t168 + t194 * t187;
t130 = t156 * qJD(4) + t194 * t146 + t190 * t186;
t163 = qJD(4) + t167;
t131 = Ifges(5,5) * t157 + Ifges(5,6) * t156 + Ifges(5,3) * t163;
t133 = Ifges(5,1) * t157 + Ifges(5,4) * t156 + Ifges(5,5) * t163;
t144 = qJDD(4) - t145;
t109 = -mrSges(5,1) * t123 + mrSges(5,3) * t121 + Ifges(5,4) * t130 + Ifges(5,2) * t129 + Ifges(5,6) * t144 - t157 * t131 + t163 * t133;
t120 = t194 * t122 - t190 * t124;
t132 = Ifges(5,4) * t157 + Ifges(5,2) * t156 + Ifges(5,6) * t163;
t110 = mrSges(5,2) * t123 - mrSges(5,3) * t120 + Ifges(5,1) * t130 + Ifges(5,4) * t129 + Ifges(5,5) * t144 + t156 * t131 - t163 * t132;
t151 = Ifges(4,4) * t168 - Ifges(4,2) * t167 + Ifges(4,6) * t187;
t152 = Ifges(4,1) * t168 - Ifges(4,4) * t167 + Ifges(4,5) * t187;
t148 = -t163 * mrSges(5,2) + t156 * mrSges(5,3);
t149 = t163 * mrSges(5,1) - t157 * mrSges(5,3);
t204 = -m(5) * t123 + t129 * mrSges(5,1) - t130 * mrSges(5,2) + t156 * t148 - t157 * t149;
t137 = -t156 * mrSges(5,1) + t157 * mrSges(5,2);
t116 = m(5) * t120 + t144 * mrSges(5,1) - t130 * mrSges(5,3) - t157 * t137 + t163 * t148;
t117 = m(5) * t121 - t144 * mrSges(5,2) + t129 * mrSges(5,3) + t156 * t137 - t163 * t149;
t209 = -t190 * t116 + t194 * t117;
t202 = -mrSges(4,1) * t126 + mrSges(4,2) * t127 - Ifges(4,5) * t146 - Ifges(4,6) * t145 - Ifges(4,3) * t186 - pkin(3) * t204 - pkin(7) * t209 - t194 * t109 - t190 * t110 - t168 * t151 - t167 * t152;
t154 = t167 * mrSges(4,1) + t168 * mrSges(4,2);
t161 = t187 * mrSges(4,1) - t168 * mrSges(4,3);
t103 = m(4) * t127 - t186 * mrSges(4,2) + t145 * mrSges(4,3) - t167 * t154 - t187 * t161 + t209;
t160 = -t187 * mrSges(4,2) - t167 * mrSges(4,3);
t112 = m(4) * t126 + t186 * mrSges(4,1) - t146 * mrSges(4,3) - t168 * t154 + t187 * t160 + t204;
t98 = t191 * t103 + t195 * t112;
t217 = mrSges(3,1) * t158 - mrSges(3,2) * t159 + Ifges(3,5) * t175 + Ifges(3,6) * t176 + Ifges(3,3) * qJDD(2) + pkin(2) * t98 + (t192 * t165 - t196 * t166) * qJD(1) - t202;
t105 = t194 * t116 + t190 * t117;
t213 = qJD(1) * t196;
t174 = (-mrSges(3,1) * t196 + mrSges(3,2) * t192) * qJD(1);
t179 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t213;
t96 = m(3) * t158 + qJDD(2) * mrSges(3,1) - t175 * mrSges(3,3) + qJD(2) * t179 - t174 * t214 + t98;
t178 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t214;
t210 = t195 * t103 - t191 * t112;
t97 = m(3) * t159 - qJDD(2) * mrSges(3,2) + t176 * mrSges(3,3) - qJD(2) * t178 + t174 * t213 + t210;
t211 = -t192 * t96 + t196 * t97;
t169 = -t198 * pkin(5) + t206;
t203 = m(4) * t147 - t145 * mrSges(4,1) + t146 * mrSges(4,2) + t167 * t160 + t168 * t161 + t105;
t200 = -m(3) * t169 + t176 * mrSges(3,1) - t175 * mrSges(3,2) - t178 * t214 + t179 * t213 - t203;
t164 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t192 + Ifges(3,6) * t196) * qJD(1);
t150 = Ifges(4,5) * t168 - Ifges(4,6) * t167 + Ifges(4,3) * t187;
t93 = mrSges(4,2) * t147 - mrSges(4,3) * t126 + Ifges(4,1) * t146 + Ifges(4,4) * t145 + Ifges(4,5) * t186 - pkin(7) * t105 - t190 * t109 + t194 * t110 - t167 * t150 - t187 * t151;
t201 = mrSges(5,1) * t120 - mrSges(5,2) * t121 + Ifges(5,5) * t130 + Ifges(5,6) * t129 + Ifges(5,3) * t144 + t157 * t132 - t156 * t133;
t94 = -mrSges(4,1) * t147 + mrSges(4,3) * t127 + Ifges(4,4) * t146 + Ifges(4,2) * t145 + Ifges(4,6) * t186 - pkin(3) * t105 - t168 * t150 + t187 * t152 - t201;
t87 = -mrSges(3,1) * t169 + mrSges(3,3) * t159 + Ifges(3,4) * t175 + Ifges(3,2) * t176 + Ifges(3,6) * qJDD(2) - pkin(2) * t203 + pkin(6) * t210 + qJD(2) * t166 - t164 * t214 + t191 * t93 + t195 * t94;
t89 = mrSges(3,2) * t169 - mrSges(3,3) * t158 + Ifges(3,1) * t175 + Ifges(3,4) * t176 + Ifges(3,5) * qJDD(2) - pkin(6) * t98 - qJD(2) * t165 + t164 * t213 - t191 * t94 + t195 * t93;
t205 = mrSges(2,1) * t181 - mrSges(2,2) * t182 + Ifges(2,3) * qJDD(1) + pkin(1) * t200 + pkin(5) * t211 + t192 * t89 + t196 * t87;
t99 = m(2) * t181 + qJDD(1) * mrSges(2,1) - t198 * mrSges(2,2) + t200;
t92 = t192 * t97 + t196 * t96;
t90 = m(2) * t182 - t198 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t211;
t85 = mrSges(2,1) * g(3) + mrSges(2,3) * t182 + t198 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t92 - t217;
t84 = -mrSges(2,2) * g(3) - mrSges(2,3) * t181 + Ifges(2,5) * qJDD(1) - t198 * Ifges(2,6) - pkin(5) * t92 - t192 * t87 + t196 * t89;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t197 * t84 - t193 * t85 - pkin(4) * (t193 * t90 + t197 * t99), t84, t89, t93, t110; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t193 * t84 + t197 * t85 + pkin(4) * (-t193 * t99 + t197 * t90), t85, t87, t94, t109; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t205, t205, t217, -t202, t201;];
m_new = t1;
