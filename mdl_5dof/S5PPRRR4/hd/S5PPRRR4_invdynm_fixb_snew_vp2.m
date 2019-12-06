% Calculate vector of cutting torques with Newton-Euler for
% S5PPRRR4
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PPRRR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR4_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR4_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR4_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR4_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR4_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:18:57
% EndTime: 2019-12-05 15:19:08
% DurationCPUTime: 7.88s
% Computational Cost: add. (142135->217), mult. (251938->293), div. (0->0), fcn. (196231->14), ass. (0->108)
t192 = sin(pkin(10));
t196 = cos(pkin(10));
t186 = -t196 * g(1) - t192 * g(2);
t191 = sin(pkin(11));
t195 = cos(pkin(11));
t185 = t192 * g(1) - t196 * g(2);
t190 = -g(3) + qJDD(1);
t194 = sin(pkin(5));
t198 = cos(pkin(5));
t212 = t185 * t198 + t190 * t194;
t158 = -t191 * t186 + t212 * t195;
t159 = t195 * t186 + t212 * t191;
t171 = -t194 * t185 + t198 * t190 + qJDD(2);
t204 = cos(qJ(3));
t197 = cos(pkin(6));
t201 = sin(qJ(3));
t222 = t197 * t201;
t193 = sin(pkin(6));
t223 = t193 * t201;
t152 = t158 * t222 + t204 * t159 + t171 * t223;
t206 = qJD(3) ^ 2;
t150 = -t206 * pkin(3) + qJDD(3) * pkin(8) + t152;
t154 = -t193 * t158 + t197 * t171;
t200 = sin(qJ(4));
t203 = cos(qJ(4));
t146 = t203 * t150 + t200 * t154;
t181 = (-pkin(4) * t203 - pkin(9) * t200) * qJD(3);
t205 = qJD(4) ^ 2;
t219 = t203 * qJD(3);
t144 = -t205 * pkin(4) + qJDD(4) * pkin(9) + t181 * t219 + t146;
t151 = -t201 * t159 + (t158 * t197 + t171 * t193) * t204;
t149 = -qJDD(3) * pkin(3) - t206 * pkin(8) - t151;
t218 = qJD(3) * qJD(4);
t216 = t203 * t218;
t182 = t200 * qJDD(3) + t216;
t217 = t200 * t218;
t183 = t203 * qJDD(3) - t217;
t147 = (-t182 - t216) * pkin(9) + (-t183 + t217) * pkin(4) + t149;
t199 = sin(qJ(5));
t202 = cos(qJ(5));
t140 = -t199 * t144 + t202 * t147;
t220 = qJD(3) * t200;
t178 = t202 * qJD(4) - t199 * t220;
t166 = t178 * qJD(5) + t199 * qJDD(4) + t202 * t182;
t179 = t199 * qJD(4) + t202 * t220;
t167 = -t178 * mrSges(6,1) + t179 * mrSges(6,2);
t189 = qJD(5) - t219;
t169 = -t189 * mrSges(6,2) + t178 * mrSges(6,3);
t177 = qJDD(5) - t183;
t138 = m(6) * t140 + t177 * mrSges(6,1) - t166 * mrSges(6,3) - t179 * t167 + t189 * t169;
t141 = t202 * t144 + t199 * t147;
t165 = -t179 * qJD(5) + t202 * qJDD(4) - t199 * t182;
t170 = t189 * mrSges(6,1) - t179 * mrSges(6,3);
t139 = m(6) * t141 - t177 * mrSges(6,2) + t165 * mrSges(6,3) + t178 * t167 - t189 * t170;
t132 = -t199 * t138 + t202 * t139;
t180 = (-mrSges(5,1) * t203 + mrSges(5,2) * t200) * qJD(3);
t187 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t220;
t130 = m(5) * t146 - qJDD(4) * mrSges(5,2) + t183 * mrSges(5,3) - qJD(4) * t187 + t180 * t219 + t132;
t221 = t203 * t154;
t143 = -qJDD(4) * pkin(4) - t205 * pkin(9) - t221 + (qJD(3) * t181 + t150) * t200;
t142 = -m(6) * t143 + t165 * mrSges(6,1) - t166 * mrSges(6,2) + t178 * t169 - t179 * t170;
t145 = -t200 * t150 + t221;
t188 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t219;
t136 = m(5) * t145 + qJDD(4) * mrSges(5,1) - t182 * mrSges(5,3) + qJD(4) * t188 - t180 * t220 + t142;
t215 = t203 * t130 - t200 * t136;
t121 = m(4) * t152 - t206 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t215;
t124 = t200 * t130 + t203 * t136;
t123 = m(4) * t154 + t124;
t131 = t202 * t138 + t199 * t139;
t209 = -m(5) * t149 + t183 * mrSges(5,1) - t182 * mrSges(5,2) - t187 * t220 + t188 * t219 - t131;
t127 = m(4) * t151 + qJDD(3) * mrSges(4,1) - t206 * mrSges(4,2) + t209;
t224 = t127 * t204;
t111 = t121 * t222 - t193 * t123 + t197 * t224;
t108 = m(3) * t158 + t111;
t115 = t204 * t121 - t201 * t127;
t114 = m(3) * t159 + t115;
t232 = t108 * t195 + t114 * t191;
t103 = -t191 * t108 + t195 * t114;
t160 = Ifges(6,5) * t179 + Ifges(6,6) * t178 + Ifges(6,3) * t189;
t162 = Ifges(6,1) * t179 + Ifges(6,4) * t178 + Ifges(6,5) * t189;
t133 = -mrSges(6,1) * t143 + mrSges(6,3) * t141 + Ifges(6,4) * t166 + Ifges(6,2) * t165 + Ifges(6,6) * t177 - t179 * t160 + t189 * t162;
t161 = Ifges(6,4) * t179 + Ifges(6,2) * t178 + Ifges(6,6) * t189;
t134 = mrSges(6,2) * t143 - mrSges(6,3) * t140 + Ifges(6,1) * t166 + Ifges(6,4) * t165 + Ifges(6,5) * t177 + t178 * t160 - t189 * t161;
t172 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t200 + Ifges(5,6) * t203) * qJD(3);
t173 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t200 + Ifges(5,2) * t203) * qJD(3);
t116 = mrSges(5,2) * t149 - mrSges(5,3) * t145 + Ifges(5,1) * t182 + Ifges(5,4) * t183 + Ifges(5,5) * qJDD(4) - pkin(9) * t131 - qJD(4) * t173 - t199 * t133 + t202 * t134 + t172 * t219;
t174 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t200 + Ifges(5,4) * t203) * qJD(3);
t208 = mrSges(6,1) * t140 - mrSges(6,2) * t141 + Ifges(6,5) * t166 + Ifges(6,6) * t165 + Ifges(6,3) * t177 + t179 * t161 - t178 * t162;
t117 = -mrSges(5,1) * t149 + mrSges(5,3) * t146 + Ifges(5,4) * t182 + Ifges(5,2) * t183 + Ifges(5,6) * qJDD(4) - pkin(4) * t131 + qJD(4) * t174 - t172 * t220 - t208;
t104 = mrSges(4,1) * t151 - mrSges(4,2) * t152 + Ifges(4,3) * qJDD(3) + pkin(3) * t209 + pkin(8) * t215 + t200 * t116 + t203 * t117;
t110 = t121 * t223 + t197 * t123 + t193 * t224;
t105 = mrSges(4,2) * t154 - mrSges(4,3) * t151 + Ifges(4,5) * qJDD(3) - t206 * Ifges(4,6) - pkin(8) * t124 + t203 * t116 - t200 * t117;
t230 = mrSges(5,1) * t145 - mrSges(5,2) * t146 + Ifges(5,5) * t182 + Ifges(5,6) * t183 + Ifges(5,3) * qJDD(4) + pkin(4) * t142 + pkin(9) * t132 + t202 * t133 + t199 * t134 + (t200 * t173 - t203 * t174) * qJD(3);
t106 = -mrSges(4,1) * t154 + mrSges(4,3) * t152 + t206 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t124 - t230;
t211 = pkin(7) * t115 + t105 * t201 + t106 * t204;
t94 = -mrSges(3,1) * t171 + mrSges(3,3) * t159 - pkin(2) * t110 - t193 * t104 + t211 * t197;
t96 = mrSges(3,2) * t171 - mrSges(3,3) * t158 + t204 * t105 - t201 * t106 + (-t110 * t193 - t111 * t197) * pkin(7);
t231 = qJ(2) * t103 + t191 * t96 + t195 * t94;
t109 = m(3) * t171 + t110;
t100 = -t194 * t109 + t232 * t198;
t92 = mrSges(3,1) * t158 - mrSges(3,2) * t159 + pkin(2) * t111 + t197 * t104 + t211 * t193;
t210 = mrSges(2,1) * t185 - mrSges(2,2) * t186 + pkin(1) * t100 + t231 * t194 + t198 * t92;
t101 = m(2) * t186 + t103;
t99 = t198 * t109 + t232 * t194;
t97 = m(2) * t185 + t100;
t90 = mrSges(2,2) * t190 - mrSges(2,3) * t185 - t191 * t94 + t195 * t96 + (-t100 * t198 - t194 * t99) * qJ(2);
t89 = -mrSges(2,1) * t190 + mrSges(2,3) * t186 - pkin(1) * t99 - t194 * t92 + t231 * t198;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t196 * t90 - t192 * t89 - qJ(1) * (t192 * t101 + t196 * t97), t90, t96, t105, t116, t134; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t192 * t90 + t196 * t89 + qJ(1) * (t196 * t101 - t192 * t97), t89, t94, t106, t117, t133; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t210, t210, t92, t104, t230, t208;];
m_new = t1;
