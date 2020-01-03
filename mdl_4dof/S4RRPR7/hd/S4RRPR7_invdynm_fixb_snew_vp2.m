% Calculate vector of cutting torques with Newton-Euler for
% S4RRPR7
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RRPR7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR7_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR7_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR7_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR7_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR7_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:06:00
% EndTime: 2019-12-31 17:06:04
% DurationCPUTime: 2.46s
% Computational Cost: add. (26641->237), mult. (60012->308), div. (0->0), fcn. (37749->8), ass. (0->97)
t188 = sin(pkin(7));
t189 = cos(pkin(7));
t191 = sin(qJ(2));
t194 = cos(qJ(2));
t165 = (t188 * t191 - t189 * t194) * qJD(1);
t192 = sin(qJ(1));
t195 = cos(qJ(1));
t183 = -t195 * g(1) - t192 * g(2);
t197 = qJD(1) ^ 2;
t171 = -t197 * pkin(1) + qJDD(1) * pkin(5) + t183;
t215 = t191 * t171;
t158 = -t194 * g(3) - t215;
t159 = -t191 * g(3) + t194 * t171;
t168 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t191 + Ifges(3,2) * t194) * qJD(1);
t169 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t191 + Ifges(3,4) * t194) * qJD(1);
t212 = qJD(1) * qJD(2);
t176 = t191 * qJDD(1) + t194 * t212;
t177 = t194 * qJDD(1) - t191 * t212;
t216 = pkin(2) * t197;
t138 = qJDD(2) * pkin(2) - t176 * qJ(3) - t215 + (qJ(3) * t212 + t191 * t216 - g(3)) * t194;
t214 = qJD(1) * t191;
t179 = qJD(2) * pkin(2) - qJ(3) * t214;
t187 = t194 ^ 2;
t139 = t177 * qJ(3) - qJD(2) * t179 - t187 * t216 + t159;
t217 = 2 * qJD(3);
t127 = t188 * t138 + t189 * t139 - t165 * t217;
t166 = (t188 * t194 + t189 * t191) * qJD(1);
t149 = t165 * pkin(3) - t166 * pkin(6);
t196 = qJD(2) ^ 2;
t123 = -t196 * pkin(3) + qJDD(2) * pkin(6) - t165 * t149 + t127;
t182 = t192 * g(1) - t195 * g(2);
t205 = -qJDD(1) * pkin(1) - t182;
t141 = -t177 * pkin(2) + qJDD(3) + t179 * t214 + (-qJ(3) * t187 - pkin(5)) * t197 + t205;
t154 = -t188 * t176 + t189 * t177;
t155 = t189 * t176 + t188 * t177;
t124 = (qJD(2) * t165 - t155) * pkin(6) + (qJD(2) * t166 - t154) * pkin(3) + t141;
t190 = sin(qJ(4));
t193 = cos(qJ(4));
t121 = t193 * t123 + t190 * t124;
t208 = -t189 * t138 + t188 * t139;
t122 = -qJDD(2) * pkin(3) - t196 * pkin(6) + (t217 + t149) * t166 + t208;
t156 = t193 * qJD(2) - t190 * t166;
t157 = t190 * qJD(2) + t193 * t166;
t164 = qJD(4) + t165;
t128 = Ifges(5,5) * t157 + Ifges(5,6) * t156 + Ifges(5,3) * t164;
t130 = Ifges(5,1) * t157 + Ifges(5,4) * t156 + Ifges(5,5) * t164;
t133 = -t157 * qJD(4) + t193 * qJDD(2) - t190 * t155;
t134 = t156 * qJD(4) + t190 * qJDD(2) + t193 * t155;
t153 = qJDD(4) - t154;
t109 = -mrSges(5,1) * t122 + mrSges(5,3) * t121 + Ifges(5,4) * t134 + Ifges(5,2) * t133 + Ifges(5,6) * t153 - t157 * t128 + t164 * t130;
t120 = -t190 * t123 + t193 * t124;
t129 = Ifges(5,4) * t157 + Ifges(5,2) * t156 + Ifges(5,6) * t164;
t110 = mrSges(5,2) * t122 - mrSges(5,3) * t120 + Ifges(5,1) * t134 + Ifges(5,4) * t133 + Ifges(5,5) * t153 + t156 * t128 - t164 * t129;
t126 = -0.2e1 * qJD(3) * t166 - t208;
t145 = Ifges(4,4) * t166 - Ifges(4,2) * t165 + Ifges(4,6) * qJD(2);
t146 = Ifges(4,1) * t166 - Ifges(4,4) * t165 + Ifges(4,5) * qJD(2);
t142 = -t164 * mrSges(5,2) + t156 * mrSges(5,3);
t143 = t164 * mrSges(5,1) - t157 * mrSges(5,3);
t203 = -m(5) * t122 + t133 * mrSges(5,1) - t134 * mrSges(5,2) + t156 * t142 - t157 * t143;
t140 = -t156 * mrSges(5,1) + t157 * mrSges(5,2);
t116 = m(5) * t120 + t153 * mrSges(5,1) - t134 * mrSges(5,3) - t157 * t140 + t164 * t142;
t117 = m(5) * t121 - t153 * mrSges(5,2) + t133 * mrSges(5,3) + t156 * t140 - t164 * t143;
t209 = -t190 * t116 + t193 * t117;
t201 = -mrSges(4,1) * t126 + mrSges(4,2) * t127 - Ifges(4,5) * t155 - Ifges(4,6) * t154 - Ifges(4,3) * qJDD(2) - pkin(3) * t203 - pkin(6) * t209 - t193 * t109 - t190 * t110 - t166 * t145 - t165 * t146;
t148 = t165 * mrSges(4,1) + t166 * mrSges(4,2);
t161 = qJD(2) * mrSges(4,1) - t166 * mrSges(4,3);
t103 = m(4) * t127 - qJDD(2) * mrSges(4,2) + t154 * mrSges(4,3) - qJD(2) * t161 - t165 * t148 + t209;
t160 = -qJD(2) * mrSges(4,2) - t165 * mrSges(4,3);
t112 = m(4) * t126 + qJDD(2) * mrSges(4,1) - t155 * mrSges(4,3) + qJD(2) * t160 - t166 * t148 + t203;
t98 = t188 * t103 + t189 * t112;
t218 = mrSges(3,1) * t158 - mrSges(3,2) * t159 + Ifges(3,5) * t176 + Ifges(3,6) * t177 + Ifges(3,3) * qJDD(2) + pkin(2) * t98 + (t191 * t168 - t194 * t169) * qJD(1) - t201;
t105 = t193 * t116 + t190 * t117;
t213 = qJD(1) * t194;
t175 = (-mrSges(3,1) * t194 + mrSges(3,2) * t191) * qJD(1);
t181 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t213;
t96 = m(3) * t158 + qJDD(2) * mrSges(3,1) - t176 * mrSges(3,3) + qJD(2) * t181 - t175 * t214 + t98;
t180 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t214;
t210 = t189 * t103 - t188 * t112;
t97 = m(3) * t159 - qJDD(2) * mrSges(3,2) + t177 * mrSges(3,3) - qJD(2) * t180 + t175 * t213 + t210;
t211 = -t191 * t96 + t194 * t97;
t170 = -t197 * pkin(5) + t205;
t202 = m(4) * t141 - t154 * mrSges(4,1) + t155 * mrSges(4,2) + t165 * t160 + t166 * t161 + t105;
t199 = -m(3) * t170 + t177 * mrSges(3,1) - t176 * mrSges(3,2) - t180 * t214 + t181 * t213 - t202;
t167 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t191 + Ifges(3,6) * t194) * qJD(1);
t144 = Ifges(4,5) * t166 - Ifges(4,6) * t165 + Ifges(4,3) * qJD(2);
t93 = mrSges(4,2) * t141 - mrSges(4,3) * t126 + Ifges(4,1) * t155 + Ifges(4,4) * t154 + Ifges(4,5) * qJDD(2) - pkin(6) * t105 - qJD(2) * t145 - t190 * t109 + t193 * t110 - t165 * t144;
t200 = mrSges(5,1) * t120 - mrSges(5,2) * t121 + Ifges(5,5) * t134 + Ifges(5,6) * t133 + Ifges(5,3) * t153 + t157 * t129 - t156 * t130;
t94 = -mrSges(4,1) * t141 + mrSges(4,3) * t127 + Ifges(4,4) * t155 + Ifges(4,2) * t154 + Ifges(4,6) * qJDD(2) - pkin(3) * t105 + qJD(2) * t146 - t166 * t144 - t200;
t87 = -mrSges(3,1) * t170 + mrSges(3,3) * t159 + Ifges(3,4) * t176 + Ifges(3,2) * t177 + Ifges(3,6) * qJDD(2) - pkin(2) * t202 + qJ(3) * t210 + qJD(2) * t169 - t167 * t214 + t188 * t93 + t189 * t94;
t89 = mrSges(3,2) * t170 - mrSges(3,3) * t158 + Ifges(3,1) * t176 + Ifges(3,4) * t177 + Ifges(3,5) * qJDD(2) - qJ(3) * t98 - qJD(2) * t168 + t167 * t213 - t188 * t94 + t189 * t93;
t204 = mrSges(2,1) * t182 - mrSges(2,2) * t183 + Ifges(2,3) * qJDD(1) + pkin(1) * t199 + pkin(5) * t211 + t191 * t89 + t194 * t87;
t99 = m(2) * t182 + qJDD(1) * mrSges(2,1) - t197 * mrSges(2,2) + t199;
t92 = t191 * t97 + t194 * t96;
t90 = m(2) * t183 - t197 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t211;
t85 = mrSges(2,1) * g(3) + mrSges(2,3) * t183 + t197 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t92 - t218;
t84 = -mrSges(2,2) * g(3) - mrSges(2,3) * t182 + Ifges(2,5) * qJDD(1) - t197 * Ifges(2,6) - pkin(5) * t92 - t191 * t87 + t194 * t89;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t195 * t84 - t192 * t85 - pkin(4) * (t192 * t90 + t195 * t99), t84, t89, t93, t110; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t192 * t84 + t195 * t85 + pkin(4) * (-t192 * t99 + t195 * t90), t85, t87, t94, t109; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t204, t204, t218, -t201, t200;];
m_new = t1;
