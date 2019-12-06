% Calculate vector of cutting torques with Newton-Euler for
% S5RRRRR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-05 18:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRRR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR5_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR5_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR5_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR5_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:58:16
% EndTime: 2019-12-05 18:58:22
% DurationCPUTime: 4.25s
% Computational Cost: add. (94882->228), mult. (97595->290), div. (0->0), fcn. (55165->10), ass. (0->100)
t197 = qJDD(1) + qJDD(2);
t190 = qJDD(3) + t197;
t203 = sin(qJ(4));
t208 = cos(qJ(4));
t199 = qJD(1) + qJD(2);
t191 = qJD(3) + t199;
t226 = qJD(4) * t191;
t173 = t203 * t190 + t208 * t226;
t206 = sin(qJ(1));
t211 = cos(qJ(1));
t185 = t211 * g(2) + t206 * g(3);
t181 = qJDD(1) * pkin(1) + t185;
t184 = t206 * g(2) - t211 * g(3);
t212 = qJD(1) ^ 2;
t182 = -t212 * pkin(1) + t184;
t205 = sin(qJ(2));
t210 = cos(qJ(2));
t161 = t210 * t181 - t205 * t182;
t158 = t197 * pkin(2) + t161;
t162 = t205 * t181 + t210 * t182;
t195 = t199 ^ 2;
t159 = -t195 * pkin(2) + t162;
t204 = sin(qJ(3));
t209 = cos(qJ(3));
t143 = t204 * t158 + t209 * t159;
t189 = t191 ^ 2;
t140 = -t189 * pkin(3) + t190 * pkin(8) + t143;
t227 = t203 * t140;
t230 = pkin(4) * t189;
t133 = qJDD(4) * pkin(4) - t173 * pkin(9) - t227 + (pkin(9) * t226 + t203 * t230 - g(1)) * t208;
t137 = -t203 * g(1) + t208 * t140;
t174 = t208 * t190 - t203 * t226;
t229 = t191 * t203;
t180 = qJD(4) * pkin(4) - pkin(9) * t229;
t201 = t208 ^ 2;
t134 = t174 * pkin(9) - qJD(4) * t180 - t201 * t230 + t137;
t202 = sin(qJ(5));
t207 = cos(qJ(5));
t131 = t207 * t133 - t202 * t134;
t168 = (-t202 * t203 + t207 * t208) * t191;
t149 = t168 * qJD(5) + t207 * t173 + t202 * t174;
t169 = (t202 * t208 + t203 * t207) * t191;
t154 = -t168 * mrSges(6,1) + t169 * mrSges(6,2);
t198 = qJD(4) + qJD(5);
t163 = -t198 * mrSges(6,2) + t168 * mrSges(6,3);
t196 = qJDD(4) + qJDD(5);
t128 = m(6) * t131 + t196 * mrSges(6,1) - t149 * mrSges(6,3) - t169 * t154 + t198 * t163;
t132 = t202 * t133 + t207 * t134;
t148 = -t169 * qJD(5) - t202 * t173 + t207 * t174;
t164 = t198 * mrSges(6,1) - t169 * mrSges(6,3);
t129 = m(6) * t132 - t196 * mrSges(6,2) + t148 * mrSges(6,3) + t168 * t154 - t198 * t164;
t119 = t207 * t128 + t202 * t129;
t136 = -t208 * g(1) - t227;
t166 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t203 + Ifges(5,2) * t208) * t191;
t167 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t203 + Ifges(5,4) * t208) * t191;
t151 = Ifges(6,4) * t169 + Ifges(6,2) * t168 + Ifges(6,6) * t198;
t152 = Ifges(6,1) * t169 + Ifges(6,4) * t168 + Ifges(6,5) * t198;
t217 = -mrSges(6,1) * t131 + mrSges(6,2) * t132 - Ifges(6,5) * t149 - Ifges(6,6) * t148 - Ifges(6,3) * t196 - t169 * t151 + t168 * t152;
t231 = mrSges(5,1) * t136 - mrSges(5,2) * t137 + Ifges(5,5) * t173 + Ifges(5,6) * t174 + Ifges(5,3) * qJDD(4) + pkin(4) * t119 + (t203 * t166 - t208 * t167) * t191 - t217;
t228 = t191 * t208;
t172 = (-mrSges(5,1) * t208 + mrSges(5,2) * t203) * t191;
t179 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t228;
t117 = m(5) * t136 + qJDD(4) * mrSges(5,1) - t173 * mrSges(5,3) + qJD(4) * t179 - t172 * t229 + t119;
t178 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t229;
t222 = -t202 * t128 + t207 * t129;
t118 = m(5) * t137 - qJDD(4) * mrSges(5,2) + t174 * mrSges(5,3) - qJD(4) * t178 + t172 * t228 + t222;
t223 = -t203 * t117 + t208 * t118;
t111 = m(4) * t143 - t189 * mrSges(4,1) - t190 * mrSges(4,2) + t223;
t142 = t209 * t158 - t204 * t159;
t220 = -t190 * pkin(3) - t142;
t139 = -t189 * pkin(8) + t220;
t135 = t180 * t229 - t174 * pkin(4) + (-pkin(9) * t201 - pkin(8)) * t189 + t220;
t218 = m(6) * t135 - t148 * mrSges(6,1) + t149 * mrSges(6,2) - t168 * t163 + t169 * t164;
t214 = -m(5) * t139 + t174 * mrSges(5,1) - t173 * mrSges(5,2) - t178 * t229 + t179 * t228 - t218;
t123 = m(4) * t142 + t190 * mrSges(4,1) - t189 * mrSges(4,2) + t214;
t106 = t204 * t111 + t209 * t123;
t103 = m(3) * t161 + t197 * mrSges(3,1) - t195 * mrSges(3,2) + t106;
t224 = t209 * t111 - t204 * t123;
t104 = m(3) * t162 - t195 * mrSges(3,1) - t197 * mrSges(3,2) + t224;
t96 = t210 * t103 + t205 * t104;
t113 = t208 * t117 + t203 * t118;
t225 = -t205 * t103 + t210 * t104;
t150 = Ifges(6,5) * t169 + Ifges(6,6) * t168 + Ifges(6,3) * t198;
t120 = -mrSges(6,1) * t135 + mrSges(6,3) * t132 + Ifges(6,4) * t149 + Ifges(6,2) * t148 + Ifges(6,6) * t196 - t169 * t150 + t198 * t152;
t121 = mrSges(6,2) * t135 - mrSges(6,3) * t131 + Ifges(6,1) * t149 + Ifges(6,4) * t148 + Ifges(6,5) * t196 + t168 * t150 - t198 * t151;
t165 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t203 + Ifges(5,6) * t208) * t191;
t108 = mrSges(5,2) * t139 - mrSges(5,3) * t136 + Ifges(5,1) * t173 + Ifges(5,4) * t174 + Ifges(5,5) * qJDD(4) - pkin(9) * t119 - qJD(4) * t166 - t202 * t120 + t207 * t121 + t165 * t228;
t99 = -mrSges(5,1) * t139 + mrSges(5,3) * t137 + Ifges(5,4) * t173 + Ifges(5,2) * t174 + Ifges(5,6) * qJDD(4) - pkin(4) * t218 + pkin(9) * t222 + qJD(4) * t167 + t207 * t120 + t202 * t121 - t165 * t229;
t219 = mrSges(4,1) * t142 - mrSges(4,2) * t143 + Ifges(4,3) * t190 + pkin(3) * t214 + pkin(8) * t223 + t203 * t108 + t208 * t99;
t216 = mrSges(3,1) * t161 - mrSges(3,2) * t162 + Ifges(3,3) * t197 + pkin(2) * t106 + t219;
t215 = mrSges(2,1) * t185 - mrSges(2,2) * t184 + Ifges(2,3) * qJDD(1) + pkin(1) * t96 + t216;
t97 = mrSges(4,1) * g(1) + mrSges(4,3) * t143 + t189 * Ifges(4,5) + Ifges(4,6) * t190 - pkin(3) * t113 - t231;
t94 = m(2) * t185 + qJDD(1) * mrSges(2,1) - t212 * mrSges(2,2) + t96;
t93 = m(2) * t184 - t212 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t225;
t92 = -mrSges(4,2) * g(1) - mrSges(4,3) * t142 + Ifges(4,5) * t190 - t189 * Ifges(4,6) - pkin(8) * t113 + t208 * t108 - t203 * t99;
t91 = -mrSges(3,2) * g(1) - mrSges(3,3) * t161 + Ifges(3,5) * t197 - t195 * Ifges(3,6) - pkin(7) * t106 - t204 * t97 + t209 * t92;
t90 = Ifges(3,6) * t197 + t195 * Ifges(3,5) + mrSges(3,1) * g(1) + mrSges(3,3) * t162 + t204 * t92 + t209 * t97 - pkin(2) * (-m(4) * g(1) + t113) + pkin(7) * t224;
t89 = -mrSges(2,2) * g(1) - mrSges(2,3) * t185 + Ifges(2,5) * qJDD(1) - t212 * Ifges(2,6) - pkin(6) * t96 - t205 * t90 + t210 * t91;
t88 = Ifges(2,6) * qJDD(1) + t212 * Ifges(2,5) + mrSges(2,3) * t184 + t205 * t91 + t210 * t90 - pkin(1) * t113 + pkin(6) * t225 + (mrSges(2,1) - pkin(1) * (-m(3) - m(4))) * g(1);
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t215, t89, t91, t92, t108, t121; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - t206 * t89 - t211 * t88 - pkin(5) * (-t206 * t94 + t211 * t93), t88, t90, t97, t99, t120; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t211 * t89 - t206 * t88 + pkin(5) * (-t206 * t93 - t211 * t94), t215, t216, t219, t231, -t217;];
m_new = t1;
