% Calculate vector of cutting torques with Newton-Euler for
% S5RPRPR2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% m [6x1]
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
% Datum: 2022-01-23 09:19
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRPR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR2_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR2_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR2_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR2_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:18:53
% EndTime: 2022-01-23 09:18:56
% DurationCPUTime: 3.64s
% Computational Cost: add. (62114->204), mult. (89300->259), div. (0->0), fcn. (51705->10), ass. (0->97)
t191 = qJD(1) + qJD(3);
t187 = t191 ^ 2;
t202 = sin(qJ(1));
t205 = cos(qJ(1));
t177 = t202 * g(1) - t205 * g(2);
t172 = qJDD(1) * pkin(1) + t177;
t178 = -t205 * g(1) - t202 * g(2);
t206 = qJD(1) ^ 2;
t173 = -t206 * pkin(1) + t178;
t197 = sin(pkin(8));
t199 = cos(pkin(8));
t159 = t199 * t172 - t197 * t173;
t156 = qJDD(1) * pkin(2) + t159;
t160 = t197 * t172 + t199 * t173;
t157 = -t206 * pkin(2) + t160;
t201 = sin(qJ(3));
t204 = cos(qJ(3));
t141 = t201 * t156 + t204 * t157;
t188 = qJDD(1) + qJDD(3);
t138 = -t187 * pkin(3) + t188 * qJ(4) + t141;
t196 = sin(pkin(9));
t195 = -g(3) + qJDD(2);
t198 = cos(pkin(9));
t227 = qJD(4) * t191;
t228 = t198 * t195 - 0.2e1 * t196 * t227;
t233 = pkin(4) * t198;
t131 = (-pkin(7) * t188 + t187 * t233 - t138) * t196 + t228;
t135 = t196 * t195 + (t138 + 0.2e1 * t227) * t198;
t230 = t188 * t198;
t190 = t198 ^ 2;
t231 = t187 * t190;
t132 = -pkin(4) * t231 + pkin(7) * t230 + t135;
t200 = sin(qJ(5));
t203 = cos(qJ(5));
t129 = t203 * t131 - t200 * t132;
t214 = -t196 * t200 + t198 * t203;
t163 = t214 * t191;
t215 = t196 * t203 + t198 * t200;
t164 = t215 * t191;
t150 = -t163 * mrSges(6,1) + t164 * mrSges(6,2);
t152 = t163 * qJD(5) + t215 * t188;
t161 = -qJD(5) * mrSges(6,2) + t163 * mrSges(6,3);
t126 = m(6) * t129 + qJDD(5) * mrSges(6,1) - t152 * mrSges(6,3) + qJD(5) * t161 - t164 * t150;
t130 = t200 * t131 + t203 * t132;
t151 = -t164 * qJD(5) + t214 * t188;
t162 = qJD(5) * mrSges(6,1) - t164 * mrSges(6,3);
t127 = m(6) * t130 - qJDD(5) * mrSges(6,2) + t151 * mrSges(6,3) - qJD(5) * t162 + t163 * t150;
t117 = t203 * t126 + t200 * t127;
t134 = -t196 * t138 + t228;
t144 = Ifges(6,4) * t164 + Ifges(6,2) * t163 + Ifges(6,6) * qJD(5);
t145 = Ifges(6,1) * t164 + Ifges(6,4) * t163 + Ifges(6,5) * qJD(5);
t211 = -mrSges(6,1) * t129 + mrSges(6,2) * t130 - Ifges(6,5) * t152 - Ifges(6,6) * t151 - Ifges(6,3) * qJDD(5) - t164 * t144 + t163 * t145;
t220 = Ifges(5,4) * t196 + Ifges(5,2) * t198;
t221 = Ifges(5,1) * t196 + Ifges(5,4) * t198;
t234 = -mrSges(5,1) * t134 + mrSges(5,2) * t135 - pkin(4) * t117 - (t196 * t220 - t198 * t221) * t187 + t211;
t232 = mrSges(5,2) * t196;
t217 = mrSges(5,3) * t188 + (-mrSges(5,1) * t198 + t232) * t187;
t115 = m(5) * t134 - t217 * t196 + t117;
t222 = -t200 * t126 + t203 * t127;
t116 = m(5) * t135 + t217 * t198 + t222;
t223 = -t196 * t115 + t198 * t116;
t109 = m(4) * t141 - t187 * mrSges(4,1) - t188 * mrSges(4,2) + t223;
t140 = t204 * t156 - t201 * t157;
t218 = qJDD(4) - t140;
t137 = -t188 * pkin(3) - t187 * qJ(4) + t218;
t189 = t196 ^ 2;
t133 = (-pkin(3) - t233) * t188 + (-qJ(4) + (-t189 - t190) * pkin(7)) * t187 + t218;
t212 = m(6) * t133 - t151 * mrSges(6,1) + t152 * mrSges(6,2) - t163 * t161 + t164 * t162;
t209 = -m(5) * t137 + mrSges(5,1) * t230 - t212 + (t187 * t189 + t231) * mrSges(5,3);
t121 = m(4) * t140 - t187 * mrSges(4,2) + (mrSges(4,1) - t232) * t188 + t209;
t224 = t204 * t109 - t201 * t121;
t100 = m(3) * t160 - t206 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t224;
t104 = t201 * t109 + t204 * t121;
t99 = m(3) * t159 + qJDD(1) * mrSges(3,1) - t206 * mrSges(3,2) + t104;
t94 = t197 * t100 + t199 * t99;
t219 = Ifges(5,5) * t196 + Ifges(5,6) * t198;
t229 = t187 * t219;
t111 = t198 * t115 + t196 * t116;
t226 = m(4) * t195 + t111;
t225 = t199 * t100 - t197 * t99;
t143 = Ifges(6,5) * t164 + Ifges(6,6) * t163 + Ifges(6,3) * qJD(5);
t118 = -mrSges(6,1) * t133 + mrSges(6,3) * t130 + Ifges(6,4) * t152 + Ifges(6,2) * t151 + Ifges(6,6) * qJDD(5) + qJD(5) * t145 - t164 * t143;
t119 = mrSges(6,2) * t133 - mrSges(6,3) * t129 + Ifges(6,1) * t152 + Ifges(6,4) * t151 + Ifges(6,5) * qJDD(5) - qJD(5) * t144 + t163 * t143;
t102 = -mrSges(5,1) * t137 + mrSges(5,3) * t135 - pkin(4) * t212 + pkin(7) * t222 + t203 * t118 + t200 * t119 + t220 * t188 - t196 * t229;
t106 = mrSges(5,2) * t137 - mrSges(5,3) * t134 - pkin(7) * t117 - t200 * t118 + t203 * t119 + t221 * t188 + t198 * t229;
t213 = -mrSges(4,2) * t141 + qJ(4) * t223 + t198 * t102 + t196 * t106 + pkin(3) * (-t188 * t232 + t209) + mrSges(4,1) * t140 + Ifges(4,3) * t188;
t210 = mrSges(3,1) * t159 - mrSges(3,2) * t160 + Ifges(3,3) * qJDD(1) + pkin(2) * t104 + t213;
t207 = mrSges(2,1) * t177 - mrSges(2,2) * t178 + Ifges(2,3) * qJDD(1) + pkin(1) * t94 + t210;
t95 = (Ifges(4,6) - t219) * t188 - mrSges(4,1) * t195 + t187 * Ifges(4,5) + mrSges(4,3) * t141 - pkin(3) * t111 + t234;
t92 = m(2) * t178 - t206 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t225;
t91 = m(2) * t177 + qJDD(1) * mrSges(2,1) - t206 * mrSges(2,2) + t94;
t90 = mrSges(4,2) * t195 - mrSges(4,3) * t140 + Ifges(4,5) * t188 - t187 * Ifges(4,6) - qJ(4) * t111 - t196 * t102 + t198 * t106;
t89 = mrSges(3,2) * t195 - mrSges(3,3) * t159 + Ifges(3,5) * qJDD(1) - t206 * Ifges(3,6) - pkin(6) * t104 - t201 * t95 + t204 * t90;
t88 = -mrSges(3,1) * t195 + mrSges(3,3) * t160 + t206 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t226 + pkin(6) * t224 + t201 * t90 + t204 * t95;
t87 = -mrSges(2,2) * g(3) - mrSges(2,3) * t177 + Ifges(2,5) * qJDD(1) - t206 * Ifges(2,6) - qJ(2) * t94 - t197 * t88 + t199 * t89;
t86 = Ifges(2,6) * qJDD(1) + t206 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t178 + t197 * t89 + t199 * t88 - pkin(1) * (m(3) * t195 + t226) + qJ(2) * t225;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t205 * t87 - t202 * t86 - pkin(5) * (t202 * t92 + t205 * t91), t87, t89, t90, t106, t119; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t202 * t87 + t205 * t86 + pkin(5) * (-t202 * t91 + t205 * t92), t86, t88, t95, t102, t118; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t207, t207, t210, t213, t219 * t188 - t234, -t211;];
m_new = t1;
