% Calculate vector of cutting torques with Newton-Euler for
% S5RRPRR4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPRR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR4_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR4_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR4_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR4_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:31:56
% EndTime: 2019-12-05 18:32:00
% DurationCPUTime: 4.12s
% Computational Cost: add. (75833->225), mult. (97595->288), div. (0->0), fcn. (55165->10), ass. (0->98)
t208 = sin(qJ(1));
t212 = cos(qJ(1));
t185 = t212 * g(2) + t208 * g(3);
t178 = qJDD(1) * pkin(1) + t185;
t184 = t208 * g(2) - t212 * g(3);
t213 = qJD(1) ^ 2;
t179 = -t213 * pkin(1) + t184;
t207 = sin(qJ(2));
t211 = cos(qJ(2));
t161 = t211 * t178 - t207 * t179;
t197 = qJDD(1) + qJDD(2);
t158 = t197 * pkin(2) + t161;
t162 = t207 * t178 + t211 * t179;
t199 = qJD(1) + qJD(2);
t195 = t199 ^ 2;
t159 = -t195 * pkin(2) + t162;
t203 = sin(pkin(9));
t204 = cos(pkin(9));
t143 = t203 * t158 + t204 * t159;
t140 = -t195 * pkin(3) + t197 * pkin(7) + t143;
t202 = -g(1) + qJDD(3);
t206 = sin(qJ(4));
t210 = cos(qJ(4));
t136 = -t206 * t140 + t210 * t202;
t229 = qJD(4) * t199;
t227 = t210 * t229;
t173 = t206 * t197 + t227;
t133 = (-t173 + t227) * pkin(8) + (t195 * t206 * t210 + qJDD(4)) * pkin(4) + t136;
t137 = t210 * t140 + t206 * t202;
t174 = t210 * t197 - t206 * t229;
t231 = t199 * t206;
t182 = qJD(4) * pkin(4) - pkin(8) * t231;
t201 = t210 ^ 2;
t134 = -t201 * t195 * pkin(4) + t174 * pkin(8) - qJD(4) * t182 + t137;
t205 = sin(qJ(5));
t209 = cos(qJ(5));
t131 = t209 * t133 - t205 * t134;
t168 = (-t205 * t206 + t209 * t210) * t199;
t149 = t168 * qJD(5) + t209 * t173 + t205 * t174;
t169 = (t205 * t210 + t206 * t209) * t199;
t154 = -t168 * mrSges(6,1) + t169 * mrSges(6,2);
t198 = qJD(4) + qJD(5);
t163 = -t198 * mrSges(6,2) + t168 * mrSges(6,3);
t196 = qJDD(4) + qJDD(5);
t128 = m(6) * t131 + t196 * mrSges(6,1) - t149 * mrSges(6,3) - t169 * t154 + t198 * t163;
t132 = t205 * t133 + t209 * t134;
t148 = -t169 * qJD(5) - t205 * t173 + t209 * t174;
t164 = t198 * mrSges(6,1) - t169 * mrSges(6,3);
t129 = m(6) * t132 - t196 * mrSges(6,2) + t148 * mrSges(6,3) + t168 * t154 - t198 * t164;
t119 = t209 * t128 + t205 * t129;
t166 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t206 + Ifges(5,2) * t210) * t199;
t167 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t206 + Ifges(5,4) * t210) * t199;
t151 = Ifges(6,4) * t169 + Ifges(6,2) * t168 + Ifges(6,6) * t198;
t152 = Ifges(6,1) * t169 + Ifges(6,4) * t168 + Ifges(6,5) * t198;
t218 = -mrSges(6,1) * t131 + mrSges(6,2) * t132 - Ifges(6,5) * t149 - Ifges(6,6) * t148 - Ifges(6,3) * t196 - t169 * t151 + t168 * t152;
t232 = mrSges(5,1) * t136 - mrSges(5,2) * t137 + Ifges(5,5) * t173 + Ifges(5,6) * t174 + Ifges(5,3) * qJDD(4) + pkin(4) * t119 + (t206 * t166 - t210 * t167) * t199 - t218;
t230 = t199 * t210;
t172 = (-mrSges(5,1) * t210 + mrSges(5,2) * t206) * t199;
t181 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t230;
t117 = m(5) * t136 + qJDD(4) * mrSges(5,1) - t173 * mrSges(5,3) + qJD(4) * t181 - t172 * t231 + t119;
t180 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t231;
t223 = -t205 * t128 + t209 * t129;
t118 = m(5) * t137 - qJDD(4) * mrSges(5,2) + t174 * mrSges(5,3) - qJD(4) * t180 + t172 * t230 + t223;
t224 = -t206 * t117 + t210 * t118;
t111 = m(4) * t143 - t195 * mrSges(4,1) - t197 * mrSges(4,2) + t224;
t142 = t204 * t158 - t203 * t159;
t221 = -t197 * pkin(3) - t142;
t139 = -t195 * pkin(7) + t221;
t135 = t182 * t231 - t174 * pkin(4) + (-pkin(8) * t201 - pkin(7)) * t195 + t221;
t219 = m(6) * t135 - t148 * mrSges(6,1) + t149 * mrSges(6,2) - t168 * t163 + t169 * t164;
t215 = -m(5) * t139 + t174 * mrSges(5,1) - t173 * mrSges(5,2) - t180 * t231 + t181 * t230 - t219;
t123 = m(4) * t142 + t197 * mrSges(4,1) - t195 * mrSges(4,2) + t215;
t106 = t203 * t111 + t204 * t123;
t103 = m(3) * t161 + t197 * mrSges(3,1) - t195 * mrSges(3,2) + t106;
t225 = t204 * t111 - t203 * t123;
t104 = m(3) * t162 - t195 * mrSges(3,1) - t197 * mrSges(3,2) + t225;
t96 = t211 * t103 + t207 * t104;
t113 = t210 * t117 + t206 * t118;
t228 = m(4) * t202 + t113;
t226 = -t207 * t103 + t211 * t104;
t150 = Ifges(6,5) * t169 + Ifges(6,6) * t168 + Ifges(6,3) * t198;
t120 = -mrSges(6,1) * t135 + mrSges(6,3) * t132 + Ifges(6,4) * t149 + Ifges(6,2) * t148 + Ifges(6,6) * t196 - t169 * t150 + t198 * t152;
t121 = mrSges(6,2) * t135 - mrSges(6,3) * t131 + Ifges(6,1) * t149 + Ifges(6,4) * t148 + Ifges(6,5) * t196 + t168 * t150 - t198 * t151;
t165 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t206 + Ifges(5,6) * t210) * t199;
t108 = mrSges(5,2) * t139 - mrSges(5,3) * t136 + Ifges(5,1) * t173 + Ifges(5,4) * t174 + Ifges(5,5) * qJDD(4) - pkin(8) * t119 - qJD(4) * t166 - t205 * t120 + t209 * t121 + t165 * t230;
t99 = -mrSges(5,1) * t139 + mrSges(5,3) * t137 + Ifges(5,4) * t173 + Ifges(5,2) * t174 + Ifges(5,6) * qJDD(4) - pkin(4) * t219 + pkin(8) * t223 + qJD(4) * t167 + t209 * t120 + t205 * t121 - t165 * t231;
t220 = mrSges(4,1) * t142 - mrSges(4,2) * t143 + Ifges(4,3) * t197 + pkin(3) * t215 + pkin(7) * t224 + t206 * t108 + t210 * t99;
t217 = mrSges(3,1) * t161 - mrSges(3,2) * t162 + Ifges(3,3) * t197 + pkin(2) * t106 + t220;
t216 = mrSges(2,1) * t185 - mrSges(2,2) * t184 + Ifges(2,3) * qJDD(1) + pkin(1) * t96 + t217;
t97 = -mrSges(4,1) * t202 + mrSges(4,3) * t143 + t195 * Ifges(4,5) + Ifges(4,6) * t197 - pkin(3) * t113 - t232;
t94 = m(2) * t185 + qJDD(1) * mrSges(2,1) - t213 * mrSges(2,2) + t96;
t93 = m(2) * t184 - t213 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t226;
t92 = mrSges(4,2) * t202 - mrSges(4,3) * t142 + Ifges(4,5) * t197 - t195 * Ifges(4,6) - pkin(7) * t113 + t210 * t108 - t206 * t99;
t91 = -mrSges(3,2) * g(1) - mrSges(3,3) * t161 + Ifges(3,5) * t197 - t195 * Ifges(3,6) - qJ(3) * t106 - t203 * t97 + t204 * t92;
t90 = mrSges(3,1) * g(1) + mrSges(3,3) * t162 + t195 * Ifges(3,5) + Ifges(3,6) * t197 - pkin(2) * t228 + qJ(3) * t225 + t203 * t92 + t204 * t97;
t89 = -mrSges(2,2) * g(1) - mrSges(2,3) * t185 + Ifges(2,5) * qJDD(1) - t213 * Ifges(2,6) - pkin(6) * t96 - t207 * t90 + t211 * t91;
t88 = Ifges(2,6) * qJDD(1) + t213 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t184 + t207 * t91 + t211 * t90 - pkin(1) * (-m(3) * g(1) + t228) + pkin(6) * t226;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t216, t89, t91, t92, t108, t121; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - t208 * t89 - t212 * t88 - pkin(5) * (-t208 * t94 + t212 * t93), t88, t90, t97, t99, t120; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t212 * t89 - t208 * t88 + pkin(5) * (-t208 * t93 - t212 * t94), t216, t217, t220, t232, -t218;];
m_new = t1;
