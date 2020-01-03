% Calculate vector of cutting torques with Newton-Euler for
% S5RRPPR1
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2020-01-03 11:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPPR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR1_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR1_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR1_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR1_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:55:49
% EndTime: 2020-01-03 11:55:52
% DurationCPUTime: 3.56s
% Computational Cost: add. (65636->204), mult. (89300->259), div. (0->0), fcn. (51705->10), ass. (0->97)
t191 = qJD(1) + qJD(2);
t187 = t191 ^ 2;
t201 = sin(qJ(1));
t204 = cos(qJ(1));
t177 = -t204 * g(2) - t201 * g(3);
t171 = qJDD(1) * pkin(1) + t177;
t176 = -t201 * g(2) + t204 * g(3);
t205 = qJD(1) ^ 2;
t172 = -t205 * pkin(1) + t176;
t200 = sin(qJ(2));
t203 = cos(qJ(2));
t158 = t203 * t171 - t200 * t172;
t188 = qJDD(1) + qJDD(2);
t155 = t188 * pkin(2) + t158;
t159 = t200 * t171 + t203 * t172;
t156 = -t187 * pkin(2) + t159;
t196 = sin(pkin(8));
t198 = cos(pkin(8));
t140 = t196 * t155 + t198 * t156;
t137 = -t187 * pkin(3) + t188 * qJ(4) + t140;
t195 = sin(pkin(9));
t194 = -g(1) + qJDD(3);
t197 = cos(pkin(9));
t226 = qJD(4) * t191;
t227 = t197 * t194 - 0.2e1 * t195 * t226;
t232 = pkin(4) * t197;
t130 = (-pkin(7) * t188 + t187 * t232 - t137) * t195 + t227;
t134 = t195 * t194 + (t137 + 0.2e1 * t226) * t197;
t229 = t188 * t197;
t190 = t197 ^ 2;
t230 = t187 * t190;
t131 = -pkin(4) * t230 + pkin(7) * t229 + t134;
t199 = sin(qJ(5));
t202 = cos(qJ(5));
t128 = t202 * t130 - t199 * t131;
t213 = -t195 * t199 + t197 * t202;
t162 = t213 * t191;
t214 = t195 * t202 + t197 * t199;
t163 = t214 * t191;
t149 = -t162 * mrSges(6,1) + t163 * mrSges(6,2);
t151 = t162 * qJD(5) + t188 * t214;
t160 = -qJD(5) * mrSges(6,2) + t162 * mrSges(6,3);
t125 = m(6) * t128 + qJDD(5) * mrSges(6,1) - t151 * mrSges(6,3) + qJD(5) * t160 - t163 * t149;
t129 = t199 * t130 + t202 * t131;
t150 = -t163 * qJD(5) + t188 * t213;
t161 = qJD(5) * mrSges(6,1) - t163 * mrSges(6,3);
t126 = m(6) * t129 - qJDD(5) * mrSges(6,2) + t150 * mrSges(6,3) - qJD(5) * t161 + t162 * t149;
t116 = t202 * t125 + t199 * t126;
t133 = -t195 * t137 + t227;
t143 = Ifges(6,4) * t163 + Ifges(6,2) * t162 + Ifges(6,6) * qJD(5);
t144 = Ifges(6,1) * t163 + Ifges(6,4) * t162 + Ifges(6,5) * qJD(5);
t210 = -mrSges(6,1) * t128 + mrSges(6,2) * t129 - Ifges(6,5) * t151 - Ifges(6,6) * t150 - Ifges(6,3) * qJDD(5) - t163 * t143 + t162 * t144;
t219 = Ifges(5,4) * t195 + Ifges(5,2) * t197;
t220 = Ifges(5,1) * t195 + Ifges(5,4) * t197;
t233 = -mrSges(5,1) * t133 + mrSges(5,2) * t134 - pkin(4) * t116 - (t195 * t219 - t197 * t220) * t187 + t210;
t231 = mrSges(5,2) * t195;
t216 = mrSges(5,3) * t188 + (-mrSges(5,1) * t197 + t231) * t187;
t114 = m(5) * t133 - t195 * t216 + t116;
t221 = -t199 * t125 + t202 * t126;
t115 = m(5) * t134 + t197 * t216 + t221;
t222 = -t195 * t114 + t197 * t115;
t108 = m(4) * t140 - t187 * mrSges(4,1) - t188 * mrSges(4,2) + t222;
t139 = t198 * t155 - t196 * t156;
t217 = qJDD(4) - t139;
t136 = -t188 * pkin(3) - t187 * qJ(4) + t217;
t189 = t195 ^ 2;
t132 = (-pkin(3) - t232) * t188 + (-qJ(4) + (-t189 - t190) * pkin(7)) * t187 + t217;
t211 = m(6) * t132 - t150 * mrSges(6,1) + t151 * mrSges(6,2) - t162 * t160 + t163 * t161;
t208 = -m(5) * t136 + mrSges(5,1) * t229 - t211 + (t187 * t189 + t230) * mrSges(5,3);
t120 = m(4) * t139 - t187 * mrSges(4,2) + (mrSges(4,1) - t231) * t188 + t208;
t103 = t196 * t108 + t198 * t120;
t98 = m(3) * t158 + t188 * mrSges(3,1) - t187 * mrSges(3,2) + t103;
t223 = t198 * t108 - t196 * t120;
t99 = m(3) * t159 - t187 * mrSges(3,1) - t188 * mrSges(3,2) + t223;
t93 = t200 * t99 + t203 * t98;
t218 = Ifges(5,5) * t195 + Ifges(5,6) * t197;
t228 = t187 * t218;
t110 = t197 * t114 + t195 * t115;
t225 = m(4) * t194 + t110;
t224 = -t200 * t98 + t203 * t99;
t142 = Ifges(6,5) * t163 + Ifges(6,6) * t162 + Ifges(6,3) * qJD(5);
t117 = -mrSges(6,1) * t132 + mrSges(6,3) * t129 + Ifges(6,4) * t151 + Ifges(6,2) * t150 + Ifges(6,6) * qJDD(5) + qJD(5) * t144 - t163 * t142;
t118 = mrSges(6,2) * t132 - mrSges(6,3) * t128 + Ifges(6,1) * t151 + Ifges(6,4) * t150 + Ifges(6,5) * qJDD(5) - qJD(5) * t143 + t162 * t142;
t101 = -mrSges(5,1) * t136 + mrSges(5,3) * t134 - pkin(4) * t211 + pkin(7) * t221 + t202 * t117 + t199 * t118 + t219 * t188 - t195 * t228;
t105 = mrSges(5,2) * t136 - mrSges(5,3) * t133 - pkin(7) * t116 - t199 * t117 + t202 * t118 + t188 * t220 + t197 * t228;
t212 = -mrSges(4,2) * t140 + qJ(4) * t222 + t197 * t101 + t195 * t105 + pkin(3) * (-t188 * t231 + t208) + mrSges(4,1) * t139 + Ifges(4,3) * t188;
t209 = mrSges(3,1) * t158 - mrSges(3,2) * t159 + Ifges(3,3) * t188 + pkin(2) * t103 + t212;
t206 = mrSges(2,1) * t177 - mrSges(2,2) * t176 + Ifges(2,3) * qJDD(1) + pkin(1) * t93 + t209;
t94 = (Ifges(4,6) - t218) * t188 - mrSges(4,1) * t194 + t187 * Ifges(4,5) + mrSges(4,3) * t140 - pkin(3) * t110 + t233;
t91 = m(2) * t177 + qJDD(1) * mrSges(2,1) - t205 * mrSges(2,2) + t93;
t90 = m(2) * t176 - t205 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t224;
t89 = mrSges(4,2) * t194 - mrSges(4,3) * t139 + Ifges(4,5) * t188 - t187 * Ifges(4,6) - qJ(4) * t110 - t195 * t101 + t197 * t105;
t88 = -mrSges(3,2) * g(1) - mrSges(3,3) * t158 + Ifges(3,5) * t188 - t187 * Ifges(3,6) - qJ(3) * t103 - t196 * t94 + t198 * t89;
t87 = mrSges(3,1) * g(1) + mrSges(3,3) * t159 + t187 * Ifges(3,5) + Ifges(3,6) * t188 - pkin(2) * t225 + qJ(3) * t223 + t196 * t89 + t198 * t94;
t86 = -mrSges(2,2) * g(1) - mrSges(2,3) * t177 + Ifges(2,5) * qJDD(1) - t205 * Ifges(2,6) - pkin(6) * t93 - t200 * t87 + t203 * t88;
t85 = Ifges(2,6) * qJDD(1) + t205 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t176 + t200 * t88 + t203 * t87 - pkin(1) * (-m(3) * g(1) + t225) + pkin(6) * t224;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t206, t86, t88, t89, t105, t118; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t201 * t86 + t204 * t85 - pkin(5) * (t201 * t91 - t204 * t90), t85, t87, t94, t101, t117; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - t204 * t86 + t201 * t85 + pkin(5) * (t201 * t90 + t204 * t91), t206, t209, t212, t188 * t218 - t233, -t210;];
m_new = t1;
