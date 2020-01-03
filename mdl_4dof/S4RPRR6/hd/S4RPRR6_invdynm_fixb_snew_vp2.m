% Calculate vector of cutting torques with Newton-Euler for
% S4RPRR6
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RPRR6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR6_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR6_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR6_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR6_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR6_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:52:31
% EndTime: 2019-12-31 16:52:35
% DurationCPUTime: 2.82s
% Computational Cost: add. (30080->216), mult. (72642->275), div. (0->0), fcn. (50185->8), ass. (0->96)
t198 = qJD(1) ^ 2;
t194 = sin(qJ(1));
t197 = cos(qJ(1));
t176 = -t197 * g(1) - t194 * g(2);
t169 = -t198 * pkin(1) + qJDD(1) * qJ(2) + t176;
t190 = sin(pkin(7));
t191 = cos(pkin(7));
t219 = qJD(1) * qJD(2);
t217 = -t191 * g(3) - 0.2e1 * t190 * t219;
t159 = -t190 * t169 + t217;
t160 = -t190 * g(3) + (t169 + 0.2e1 * t219) * t191;
t224 = pkin(2) * t191;
t144 = (-pkin(5) * qJDD(1) + t198 * t224 - t169) * t190 + t217;
t218 = qJDD(1) * t191;
t186 = t191 ^ 2;
t222 = t186 * t198;
t145 = -pkin(2) * t222 + pkin(5) * t218 + t160;
t193 = sin(qJ(3));
t196 = cos(qJ(3));
t130 = t196 * t144 - t193 * t145;
t207 = t190 * t196 + t191 * t193;
t206 = -t190 * t193 + t191 * t196;
t167 = t206 * qJD(1);
t220 = t167 * qJD(3);
t158 = t207 * qJDD(1) + t220;
t168 = t207 * qJD(1);
t118 = (-t158 + t220) * pkin(6) + (t167 * t168 + qJDD(3)) * pkin(3) + t130;
t131 = t193 * t144 + t196 * t145;
t157 = -t168 * qJD(3) + t206 * qJDD(1);
t163 = qJD(3) * pkin(3) - t168 * pkin(6);
t166 = t167 ^ 2;
t119 = -t166 * pkin(3) + t157 * pkin(6) - qJD(3) * t163 + t131;
t192 = sin(qJ(4));
t195 = cos(qJ(4));
t116 = t195 * t118 - t192 * t119;
t149 = t195 * t167 - t192 * t168;
t128 = t149 * qJD(4) + t192 * t157 + t195 * t158;
t150 = t192 * t167 + t195 * t168;
t137 = -t149 * mrSges(5,1) + t150 * mrSges(5,2);
t187 = qJD(3) + qJD(4);
t141 = -t187 * mrSges(5,2) + t149 * mrSges(5,3);
t184 = qJDD(3) + qJDD(4);
t113 = m(5) * t116 + t184 * mrSges(5,1) - t128 * mrSges(5,3) - t150 * t137 + t187 * t141;
t117 = t192 * t118 + t195 * t119;
t127 = -t150 * qJD(4) + t195 * t157 - t192 * t158;
t142 = t187 * mrSges(5,1) - t150 * mrSges(5,3);
t114 = m(5) * t117 - t184 * mrSges(5,2) + t127 * mrSges(5,3) + t149 * t137 - t187 * t142;
t105 = t195 * t113 + t192 * t114;
t147 = Ifges(4,4) * t168 + Ifges(4,2) * t167 + Ifges(4,6) * qJD(3);
t148 = Ifges(4,1) * t168 + Ifges(4,4) * t167 + Ifges(4,5) * qJD(3);
t133 = Ifges(5,4) * t150 + Ifges(5,2) * t149 + Ifges(5,6) * t187;
t134 = Ifges(5,1) * t150 + Ifges(5,4) * t149 + Ifges(5,5) * t187;
t203 = -mrSges(5,1) * t116 + mrSges(5,2) * t117 - Ifges(5,5) * t128 - Ifges(5,6) * t127 - Ifges(5,3) * t184 - t150 * t133 + t149 * t134;
t200 = -mrSges(4,1) * t130 + mrSges(4,2) * t131 - Ifges(4,5) * t158 - Ifges(4,6) * t157 - Ifges(4,3) * qJDD(3) - pkin(3) * t105 - t168 * t147 + t167 * t148 + t203;
t211 = Ifges(3,4) * t190 + Ifges(3,2) * t191;
t212 = Ifges(3,1) * t190 + Ifges(3,4) * t191;
t152 = -t167 * mrSges(4,1) + t168 * mrSges(4,2);
t161 = -qJD(3) * mrSges(4,2) + t167 * mrSges(4,3);
t102 = m(4) * t130 + qJDD(3) * mrSges(4,1) - t158 * mrSges(4,3) + qJD(3) * t161 - t168 * t152 + t105;
t162 = qJD(3) * mrSges(4,1) - t168 * mrSges(4,3);
t214 = -t192 * t113 + t195 * t114;
t103 = m(4) * t131 - qJDD(3) * mrSges(4,2) + t157 * mrSges(4,3) - qJD(3) * t162 + t167 * t152 + t214;
t98 = t196 * t102 + t193 * t103;
t225 = -mrSges(3,1) * t159 + mrSges(3,2) * t160 - pkin(2) * t98 - (t190 * t211 - t191 * t212) * t198 + t200;
t223 = mrSges(3,2) * t190;
t210 = Ifges(3,5) * t190 + Ifges(3,6) * t191;
t221 = t198 * t210;
t175 = t194 * g(1) - t197 * g(2);
t205 = mrSges(3,3) * qJDD(1) + t198 * (-mrSges(3,1) * t191 + t223);
t96 = m(3) * t159 - t205 * t190 + t98;
t215 = -t193 * t102 + t196 * t103;
t97 = m(3) * t160 + t205 * t191 + t215;
t216 = -t190 * t96 + t191 * t97;
t213 = qJDD(2) - t175;
t185 = t190 ^ 2;
t156 = (-pkin(1) - t224) * qJDD(1) + (-qJ(2) + (-t185 - t186) * pkin(5)) * t198 + t213;
t122 = -t157 * pkin(3) - t166 * pkin(6) + t168 * t163 + t156;
t209 = m(5) * t122 - t127 * mrSges(5,1) + t128 * mrSges(5,2) - t149 * t141 + t150 * t142;
t165 = -qJDD(1) * pkin(1) - t198 * qJ(2) + t213;
t202 = m(4) * t156 - t157 * mrSges(4,1) + t158 * mrSges(4,2) - t167 * t161 + t168 * t162 + t209;
t201 = -m(3) * t165 + mrSges(3,1) * t218 - t202 + (t185 * t198 + t222) * mrSges(3,3);
t132 = Ifges(5,5) * t150 + Ifges(5,6) * t149 + Ifges(5,3) * t187;
t106 = -mrSges(5,1) * t122 + mrSges(5,3) * t117 + Ifges(5,4) * t128 + Ifges(5,2) * t127 + Ifges(5,6) * t184 - t150 * t132 + t187 * t134;
t107 = mrSges(5,2) * t122 - mrSges(5,3) * t116 + Ifges(5,1) * t128 + Ifges(5,4) * t127 + Ifges(5,5) * t184 + t149 * t132 - t187 * t133;
t146 = Ifges(4,5) * t168 + Ifges(4,6) * t167 + Ifges(4,3) * qJD(3);
t93 = -mrSges(4,1) * t156 + mrSges(4,3) * t131 + Ifges(4,4) * t158 + Ifges(4,2) * t157 + Ifges(4,6) * qJDD(3) - pkin(3) * t209 + pkin(6) * t214 + qJD(3) * t148 + t195 * t106 + t192 * t107 - t168 * t146;
t94 = mrSges(4,2) * t156 - mrSges(4,3) * t130 + Ifges(4,1) * t158 + Ifges(4,4) * t157 + Ifges(4,5) * qJDD(3) - pkin(6) * t105 - qJD(3) * t147 - t192 * t106 + t195 * t107 + t167 * t146;
t87 = -mrSges(3,1) * t165 + mrSges(3,3) * t160 - pkin(2) * t202 + pkin(5) * t215 + t211 * qJDD(1) - t190 * t221 + t193 * t94 + t196 * t93;
t89 = mrSges(3,2) * t165 - mrSges(3,3) * t159 - pkin(5) * t98 + t212 * qJDD(1) + t191 * t221 - t193 * t93 + t196 * t94;
t204 = -mrSges(2,2) * t176 + qJ(2) * t216 + t190 * t89 + t191 * t87 + pkin(1) * (-qJDD(1) * t223 + t201) + mrSges(2,1) * t175 + Ifges(2,3) * qJDD(1);
t108 = (mrSges(2,1) - t223) * qJDD(1) + t201 - t198 * mrSges(2,2) + m(2) * t175;
t92 = t190 * t97 + t191 * t96;
t90 = m(2) * t176 - t198 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t216;
t85 = mrSges(2,1) * g(3) + (Ifges(2,6) - t210) * qJDD(1) + t198 * Ifges(2,5) + mrSges(2,3) * t176 - pkin(1) * t92 + t225;
t84 = -mrSges(2,2) * g(3) - mrSges(2,3) * t175 + Ifges(2,5) * qJDD(1) - t198 * Ifges(2,6) - qJ(2) * t92 - t190 * t87 + t191 * t89;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t197 * t84 - t194 * t85 - pkin(4) * (t197 * t108 + t194 * t90), t84, t89, t94, t107; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t194 * t84 + t197 * t85 + pkin(4) * (-t194 * t108 + t197 * t90), t85, t87, t93, t106; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t204, t204, t210 * qJDD(1) - t225, -t200, -t203;];
m_new = t1;
