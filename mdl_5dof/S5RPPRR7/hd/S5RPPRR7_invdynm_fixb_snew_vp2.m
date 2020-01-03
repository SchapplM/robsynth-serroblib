% Calculate vector of cutting torques with Newton-Euler for
% S5RPPRR7
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPPRR7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR7_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR7_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR7_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR7_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR7_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:34
% EndTime: 2019-12-31 17:59:36
% DurationCPUTime: 1.82s
% Computational Cost: add. (22552->223), mult. (40467->272), div. (0->0), fcn. (21025->8), ass. (0->96)
t198 = sin(qJ(1));
t201 = cos(qJ(1));
t178 = t198 * g(1) - t201 * g(2);
t169 = qJDD(1) * pkin(1) + t178;
t179 = -t201 * g(1) - t198 * g(2);
t203 = qJD(1) ^ 2;
t171 = -t203 * pkin(1) + t179;
t194 = sin(pkin(8));
t195 = cos(pkin(8));
t149 = t194 * t169 + t195 * t171;
t229 = -qJDD(1) * qJ(3) - (2 * qJD(3) * qJD(1)) - t149;
t228 = -pkin(2) - pkin(6);
t227 = mrSges(3,1) - mrSges(4,2);
t226 = -Ifges(4,4) + Ifges(3,5);
t225 = Ifges(4,5) - Ifges(3,6);
t135 = t228 * t203 - t229;
t197 = sin(qJ(4));
t200 = cos(qJ(4));
t221 = qJD(1) * qJD(4);
t217 = t200 * t221;
t173 = -t197 * qJDD(1) - t217;
t218 = t197 * t221;
t174 = t200 * qJDD(1) - t218;
t127 = (-t174 + t218) * pkin(7) + (-t173 + t217) * pkin(4) + t135;
t148 = t195 * t169 - t194 * t171;
t214 = -t203 * qJ(3) + qJDD(3) - t148;
t136 = t228 * qJDD(1) + t214;
t191 = -g(3) + qJDD(2);
t132 = t197 * t136 + t200 * t191;
t172 = (pkin(4) * t197 - pkin(7) * t200) * qJD(1);
t202 = qJD(4) ^ 2;
t222 = t197 * qJD(1);
t129 = -t202 * pkin(4) + qJDD(4) * pkin(7) - t172 * t222 + t132;
t196 = sin(qJ(5));
t199 = cos(qJ(5));
t125 = t199 * t127 - t196 * t129;
t223 = qJD(1) * t200;
t167 = t199 * qJD(4) - t196 * t223;
t146 = t167 * qJD(5) + t196 * qJDD(4) + t199 * t174;
t168 = t196 * qJD(4) + t199 * t223;
t150 = -t167 * mrSges(6,1) + t168 * mrSges(6,2);
t180 = qJD(5) + t222;
t151 = -t180 * mrSges(6,2) + t167 * mrSges(6,3);
t166 = qJDD(5) - t173;
t121 = m(6) * t125 + t166 * mrSges(6,1) - t146 * mrSges(6,3) - t168 * t150 + t180 * t151;
t126 = t196 * t127 + t199 * t129;
t145 = -t168 * qJD(5) + t199 * qJDD(4) - t196 * t174;
t152 = t180 * mrSges(6,1) - t168 * mrSges(6,3);
t122 = m(6) * t126 - t166 * mrSges(6,2) + t145 * mrSges(6,3) + t167 * t150 - t180 * t152;
t111 = t199 * t121 + t196 * t122;
t176 = -(qJD(4) * mrSges(5,2)) - mrSges(5,3) * t222;
t177 = (qJD(4) * mrSges(5,1)) - mrSges(5,3) * t223;
t108 = -m(5) * t135 + t173 * mrSges(5,1) - t174 * mrSges(5,2) - t176 * t222 - t177 * t223 - t111;
t137 = t203 * pkin(2) + t229;
t208 = -m(4) * t137 + t203 * mrSges(4,2) + qJDD(1) * mrSges(4,3) - t108;
t105 = m(3) * t149 - t203 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t208;
t170 = (mrSges(5,1) * t197 + mrSges(5,2) * t200) * qJD(1);
t215 = -t196 * t121 + t199 * t122;
t109 = m(5) * t132 - qJDD(4) * mrSges(5,2) + t173 * mrSges(5,3) - qJD(4) * t177 - t170 * t222 + t215;
t224 = t197 * t191;
t131 = t200 * t136 - t224;
t128 = -qJDD(4) * pkin(4) - t202 * pkin(7) + t224 + (qJD(1) * t172 - t136) * t200;
t212 = -m(6) * t128 + t145 * mrSges(6,1) - t146 * mrSges(6,2) + t167 * t151 - t168 * t152;
t117 = m(5) * t131 + qJDD(4) * mrSges(5,1) - t174 * mrSges(5,3) + qJD(4) * t176 - t170 * t223 + t212;
t101 = t197 * t109 + t200 * t117;
t139 = -qJDD(1) * pkin(2) + t214;
t213 = -m(4) * t139 + t203 * mrSges(4,3) - t101;
t98 = m(3) * t148 - t203 * mrSges(3,2) + t227 * qJDD(1) + t213;
t93 = t194 * t105 + t195 * t98;
t216 = t195 * t105 - t194 * t98;
t102 = t200 * t109 - t197 * t117;
t100 = m(4) * t191 + t102;
t140 = Ifges(6,5) * t168 + Ifges(6,6) * t167 + Ifges(6,3) * t180;
t142 = Ifges(6,1) * t168 + Ifges(6,4) * t167 + Ifges(6,5) * t180;
t115 = -mrSges(6,1) * t128 + mrSges(6,3) * t126 + Ifges(6,4) * t146 + Ifges(6,2) * t145 + Ifges(6,6) * t166 - t168 * t140 + t180 * t142;
t141 = Ifges(6,4) * t168 + Ifges(6,2) * t167 + Ifges(6,6) * t180;
t116 = mrSges(6,2) * t128 - mrSges(6,3) * t125 + Ifges(6,1) * t146 + Ifges(6,4) * t145 + Ifges(6,5) * t166 + t167 * t140 - t180 * t141;
t158 = (Ifges(5,3) * qJD(4)) + (Ifges(5,5) * t200 - Ifges(5,6) * t197) * qJD(1);
t159 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t200 - Ifges(5,2) * t197) * qJD(1);
t95 = mrSges(5,2) * t135 - mrSges(5,3) * t131 + Ifges(5,1) * t174 + Ifges(5,4) * t173 + Ifges(5,5) * qJDD(4) - pkin(7) * t111 - qJD(4) * t159 - t196 * t115 + t199 * t116 - t158 * t222;
t160 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t200 - Ifges(5,4) * t197) * qJD(1);
t205 = mrSges(6,1) * t125 - mrSges(6,2) * t126 + Ifges(6,5) * t146 + Ifges(6,6) * t145 + Ifges(6,3) * t166 + t168 * t141 - t167 * t142;
t96 = -mrSges(5,1) * t135 + mrSges(5,3) * t132 + Ifges(5,4) * t174 + Ifges(5,2) * t173 + Ifges(5,6) * qJDD(4) - pkin(4) * t111 + qJD(4) * t160 - t158 * t223 - t205;
t211 = mrSges(4,2) * t139 - mrSges(4,3) * t137 + Ifges(4,1) * qJDD(1) - pkin(6) * t101 - t197 * t96 + t200 * t95;
t210 = -mrSges(4,1) * t137 - pkin(3) * t108 - pkin(6) * t102 - t197 * t95 - t200 * t96;
t209 = mrSges(5,1) * t131 - mrSges(5,2) * t132 + Ifges(5,5) * t174 + Ifges(5,6) * t173 + Ifges(5,3) * qJDD(4) + pkin(4) * t212 + pkin(7) * t215 + t199 * t115 + t196 * t116 + t159 * t223 + t160 * t222;
t207 = -mrSges(3,2) * t149 + qJ(3) * t208 + mrSges(3,1) * t148 + Ifges(3,3) * qJDD(1) + t211 + pkin(2) * (-qJDD(1) * mrSges(4,2) + t213);
t206 = mrSges(4,1) * t139 + pkin(3) * t101 + t209;
t204 = mrSges(2,1) * t178 - mrSges(2,2) * t179 + Ifges(2,3) * qJDD(1) + pkin(1) * t93 + t207;
t91 = m(2) * t179 - t203 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t216;
t90 = m(2) * t178 + qJDD(1) * mrSges(2,1) - t203 * mrSges(2,2) + t93;
t89 = t225 * t203 + (mrSges(3,2) - mrSges(4,3)) * t191 + t226 * qJDD(1) - mrSges(3,3) * t148 - qJ(3) * t100 + t206;
t88 = mrSges(3,3) * t149 - pkin(2) * t100 - t225 * qJDD(1) - t227 * t191 + t226 * t203 + t210;
t87 = -mrSges(2,2) * g(3) - mrSges(2,3) * t178 + Ifges(2,5) * qJDD(1) - t203 * Ifges(2,6) - qJ(2) * t93 - t194 * t88 + t195 * t89;
t86 = Ifges(2,6) * qJDD(1) + t203 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t179 + t194 * t89 + t195 * t88 - pkin(1) * (m(3) * t191 + t100) + qJ(2) * t216;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t201 * t87 - t198 * t86 - pkin(5) * (t198 * t91 + t201 * t90), t87, t89, t211, t95, t116; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t198 * t87 + t201 * t86 + pkin(5) * (-t198 * t90 + t201 * t91), t86, t88, mrSges(4,3) * t191 + Ifges(4,4) * qJDD(1) - t203 * Ifges(4,5) - t206, t96, t115; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t204, t204, t207, -mrSges(4,2) * t191 + t203 * Ifges(4,4) + Ifges(4,5) * qJDD(1) - t210, t209, t205;];
m_new = t1;
