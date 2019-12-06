% Calculate vector of cutting torques with Newton-Euler for
% S5RPRRR5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-05 18:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRRR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR5_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR5_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR5_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR5_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:16:15
% EndTime: 2019-12-05 18:16:20
% DurationCPUTime: 4.09s
% Computational Cost: add. (72311->225), mult. (97595->288), div. (0->0), fcn. (55165->10), ass. (0->98)
t209 = sin(qJ(1));
t213 = cos(qJ(1));
t186 = t213 * g(2) + t209 * g(3);
t179 = qJDD(1) * pkin(1) + t186;
t185 = t209 * g(2) - g(3) * t213;
t214 = qJD(1) ^ 2;
t180 = -pkin(1) * t214 + t185;
t204 = sin(pkin(9));
t205 = cos(pkin(9));
t162 = t205 * t179 - t180 * t204;
t159 = qJDD(1) * pkin(2) + t162;
t163 = t204 * t179 + t205 * t180;
t160 = -pkin(2) * t214 + t163;
t208 = sin(qJ(3));
t212 = cos(qJ(3));
t144 = t208 * t159 + t212 * t160;
t199 = qJD(1) + qJD(3);
t195 = t199 ^ 2;
t197 = qJDD(1) + qJDD(3);
t141 = -pkin(3) * t195 + pkin(7) * t197 + t144;
t203 = -g(1) + qJDD(2);
t207 = sin(qJ(4));
t211 = cos(qJ(4));
t137 = -t207 * t141 + t211 * t203;
t230 = qJD(4) * t199;
t228 = t211 * t230;
t174 = t197 * t207 + t228;
t134 = (-t174 + t228) * pkin(8) + (t195 * t207 * t211 + qJDD(4)) * pkin(4) + t137;
t138 = t211 * t141 + t207 * t203;
t175 = t197 * t211 - t207 * t230;
t232 = t199 * t207;
t183 = qJD(4) * pkin(4) - pkin(8) * t232;
t202 = t211 ^ 2;
t135 = -pkin(4) * t195 * t202 + pkin(8) * t175 - qJD(4) * t183 + t138;
t206 = sin(qJ(5));
t210 = cos(qJ(5));
t132 = t134 * t210 - t135 * t206;
t169 = (-t206 * t207 + t210 * t211) * t199;
t150 = qJD(5) * t169 + t174 * t210 + t175 * t206;
t170 = (t206 * t211 + t207 * t210) * t199;
t155 = -mrSges(6,1) * t169 + mrSges(6,2) * t170;
t198 = qJD(4) + qJD(5);
t164 = -mrSges(6,2) * t198 + mrSges(6,3) * t169;
t196 = qJDD(4) + qJDD(5);
t129 = m(6) * t132 + mrSges(6,1) * t196 - mrSges(6,3) * t150 - t155 * t170 + t164 * t198;
t133 = t134 * t206 + t135 * t210;
t149 = -qJD(5) * t170 - t174 * t206 + t175 * t210;
t165 = mrSges(6,1) * t198 - mrSges(6,3) * t170;
t130 = m(6) * t133 - mrSges(6,2) * t196 + mrSges(6,3) * t149 + t155 * t169 - t165 * t198;
t120 = t129 * t210 + t130 * t206;
t167 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t207 + Ifges(5,2) * t211) * t199;
t168 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t207 + Ifges(5,4) * t211) * t199;
t152 = Ifges(6,4) * t170 + Ifges(6,2) * t169 + Ifges(6,6) * t198;
t153 = Ifges(6,1) * t170 + Ifges(6,4) * t169 + Ifges(6,5) * t198;
t219 = -mrSges(6,1) * t132 + mrSges(6,2) * t133 - Ifges(6,5) * t150 - Ifges(6,6) * t149 - Ifges(6,3) * t196 - t170 * t152 + t169 * t153;
t233 = mrSges(5,1) * t137 - mrSges(5,2) * t138 + Ifges(5,5) * t174 + Ifges(5,6) * t175 + Ifges(5,3) * qJDD(4) + pkin(4) * t120 + (t167 * t207 - t168 * t211) * t199 - t219;
t231 = t199 * t211;
t173 = (-mrSges(5,1) * t211 + mrSges(5,2) * t207) * t199;
t182 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t231;
t118 = m(5) * t137 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t174 + qJD(4) * t182 - t173 * t232 + t120;
t181 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t232;
t224 = -t129 * t206 + t130 * t210;
t119 = m(5) * t138 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t175 - qJD(4) * t181 + t173 * t231 + t224;
t225 = -t118 * t207 + t119 * t211;
t112 = m(4) * t144 - mrSges(4,1) * t195 - mrSges(4,2) * t197 + t225;
t143 = t212 * t159 - t208 * t160;
t222 = -t197 * pkin(3) - t143;
t140 = -pkin(7) * t195 + t222;
t136 = t183 * t232 - t175 * pkin(4) + (-pkin(8) * t202 - pkin(7)) * t195 + t222;
t220 = m(6) * t136 - t149 * mrSges(6,1) + mrSges(6,2) * t150 - t169 * t164 + t165 * t170;
t216 = -m(5) * t140 + t175 * mrSges(5,1) - mrSges(5,2) * t174 - t181 * t232 + t182 * t231 - t220;
t124 = m(4) * t143 + mrSges(4,1) * t197 - mrSges(4,2) * t195 + t216;
t107 = t112 * t208 + t124 * t212;
t104 = m(3) * t162 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t214 + t107;
t226 = t112 * t212 - t124 * t208;
t105 = m(3) * t163 - mrSges(3,1) * t214 - qJDD(1) * mrSges(3,2) + t226;
t97 = t104 * t205 + t105 * t204;
t114 = t118 * t211 + t119 * t207;
t229 = m(4) * t203 + t114;
t227 = -t104 * t204 + t105 * t205;
t151 = Ifges(6,5) * t170 + Ifges(6,6) * t169 + Ifges(6,3) * t198;
t121 = -mrSges(6,1) * t136 + mrSges(6,3) * t133 + Ifges(6,4) * t150 + Ifges(6,2) * t149 + Ifges(6,6) * t196 - t151 * t170 + t153 * t198;
t122 = mrSges(6,2) * t136 - mrSges(6,3) * t132 + Ifges(6,1) * t150 + Ifges(6,4) * t149 + Ifges(6,5) * t196 + t151 * t169 - t152 * t198;
t166 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t207 + Ifges(5,6) * t211) * t199;
t100 = -mrSges(5,1) * t140 + mrSges(5,3) * t138 + Ifges(5,4) * t174 + Ifges(5,2) * t175 + Ifges(5,6) * qJDD(4) - pkin(4) * t220 + pkin(8) * t224 + qJD(4) * t168 + t210 * t121 + t206 * t122 - t166 * t232;
t109 = mrSges(5,2) * t140 - mrSges(5,3) * t137 + Ifges(5,1) * t174 + Ifges(5,4) * t175 + Ifges(5,5) * qJDD(4) - pkin(8) * t120 - qJD(4) * t167 - t121 * t206 + t122 * t210 + t166 * t231;
t221 = mrSges(4,1) * t143 - mrSges(4,2) * t144 + Ifges(4,3) * t197 + pkin(3) * t216 + pkin(7) * t225 + t100 * t211 + t109 * t207;
t218 = mrSges(3,1) * t162 - mrSges(3,2) * t163 + Ifges(3,3) * qJDD(1) + pkin(2) * t107 + t221;
t217 = mrSges(2,1) * t186 - mrSges(2,2) * t185 + Ifges(2,3) * qJDD(1) + pkin(1) * t97 + t218;
t98 = -mrSges(4,1) * t203 + mrSges(4,3) * t144 + t195 * Ifges(4,5) + Ifges(4,6) * t197 - pkin(3) * t114 - t233;
t95 = m(2) * t186 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t214 + t97;
t94 = m(2) * t185 - mrSges(2,1) * t214 - qJDD(1) * mrSges(2,2) + t227;
t93 = mrSges(4,2) * t203 - mrSges(4,3) * t143 + Ifges(4,5) * t197 - Ifges(4,6) * t195 - pkin(7) * t114 - t100 * t207 + t109 * t211;
t92 = mrSges(3,2) * t203 - mrSges(3,3) * t162 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t214 - pkin(6) * t107 - t208 * t98 + t212 * t93;
t91 = -mrSges(3,1) * t203 + mrSges(3,3) * t163 + t214 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t229 + pkin(6) * t226 + t208 * t93 + t212 * t98;
t90 = -mrSges(2,2) * g(1) - mrSges(2,3) * t186 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t214 - qJ(2) * t97 - t204 * t91 + t205 * t92;
t89 = Ifges(2,6) * qJDD(1) + t214 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t185 + t204 * t92 + t205 * t91 - pkin(1) * (m(3) * t203 + t229) + qJ(2) * t227;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t217, t90, t92, t93, t109, t122; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - t209 * t90 - t213 * t89 - pkin(5) * (-t209 * t95 + t213 * t94), t89, t91, t98, t100, t121; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t213 * t90 - t209 * t89 + pkin(5) * (-t209 * t94 - t213 * t95), t217, t218, t221, t233, -t219;];
m_new = t1;
