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
% Datum: 2020-01-03 12:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:13:09
% EndTime: 2020-01-03 12:13:15
% DurationCPUTime: 4.22s
% Computational Cost: add. (94882->228), mult. (97595->290), div. (0->0), fcn. (55165->10), ass. (0->100)
t193 = qJDD(1) + qJDD(2);
t188 = qJDD(3) + t193;
t199 = sin(qJ(4));
t204 = cos(qJ(4));
t195 = qJD(1) + qJD(2);
t189 = qJD(3) + t195;
t222 = qJD(4) * t189;
t171 = t199 * t188 + t204 * t222;
t202 = sin(qJ(1));
t207 = cos(qJ(1));
t183 = -t207 * g(2) - t202 * g(3);
t179 = qJDD(1) * pkin(1) + t183;
t182 = -t202 * g(2) + t207 * g(3);
t208 = qJD(1) ^ 2;
t180 = -t208 * pkin(1) + t182;
t201 = sin(qJ(2));
t206 = cos(qJ(2));
t159 = t206 * t179 - t201 * t180;
t156 = t193 * pkin(2) + t159;
t160 = t201 * t179 + t206 * t180;
t191 = t195 ^ 2;
t157 = -t191 * pkin(2) + t160;
t200 = sin(qJ(3));
t205 = cos(qJ(3));
t141 = t200 * t156 + t205 * t157;
t187 = t189 ^ 2;
t138 = -t187 * pkin(3) + t188 * pkin(8) + t141;
t223 = t199 * t138;
t226 = pkin(4) * t187;
t131 = qJDD(4) * pkin(4) - t171 * pkin(9) - t223 + (pkin(9) * t222 + t199 * t226 - g(1)) * t204;
t135 = -t199 * g(1) + t204 * t138;
t172 = t204 * t188 - t199 * t222;
t225 = t189 * t199;
t178 = qJD(4) * pkin(4) - pkin(9) * t225;
t197 = t204 ^ 2;
t132 = t172 * pkin(9) - qJD(4) * t178 - t197 * t226 + t135;
t198 = sin(qJ(5));
t203 = cos(qJ(5));
t129 = t203 * t131 - t198 * t132;
t166 = (-t198 * t199 + t203 * t204) * t189;
t147 = t166 * qJD(5) + t203 * t171 + t198 * t172;
t167 = (t198 * t204 + t199 * t203) * t189;
t152 = -t166 * mrSges(6,1) + t167 * mrSges(6,2);
t194 = qJD(4) + qJD(5);
t161 = -t194 * mrSges(6,2) + t166 * mrSges(6,3);
t192 = qJDD(4) + qJDD(5);
t126 = m(6) * t129 + t192 * mrSges(6,1) - t147 * mrSges(6,3) - t167 * t152 + t194 * t161;
t130 = t198 * t131 + t203 * t132;
t146 = -t167 * qJD(5) - t198 * t171 + t203 * t172;
t162 = t194 * mrSges(6,1) - t167 * mrSges(6,3);
t127 = m(6) * t130 - t192 * mrSges(6,2) + t146 * mrSges(6,3) + t166 * t152 - t194 * t162;
t117 = t203 * t126 + t198 * t127;
t134 = -t204 * g(1) - t223;
t164 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t199 + Ifges(5,2) * t204) * t189;
t165 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t199 + Ifges(5,4) * t204) * t189;
t149 = Ifges(6,4) * t167 + Ifges(6,2) * t166 + Ifges(6,6) * t194;
t150 = Ifges(6,1) * t167 + Ifges(6,4) * t166 + Ifges(6,5) * t194;
t213 = -mrSges(6,1) * t129 + mrSges(6,2) * t130 - Ifges(6,5) * t147 - Ifges(6,6) * t146 - Ifges(6,3) * t192 - t167 * t149 + t166 * t150;
t227 = mrSges(5,1) * t134 - mrSges(5,2) * t135 + Ifges(5,5) * t171 + Ifges(5,6) * t172 + Ifges(5,3) * qJDD(4) + pkin(4) * t117 + (t199 * t164 - t204 * t165) * t189 - t213;
t170 = (-mrSges(5,1) * t204 + mrSges(5,2) * t199) * t189;
t224 = t189 * t204;
t177 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t224;
t115 = m(5) * t134 + qJDD(4) * mrSges(5,1) - t171 * mrSges(5,3) + qJD(4) * t177 - t170 * t225 + t117;
t176 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t225;
t218 = -t198 * t126 + t203 * t127;
t116 = m(5) * t135 - qJDD(4) * mrSges(5,2) + t172 * mrSges(5,3) - qJD(4) * t176 + t170 * t224 + t218;
t219 = -t199 * t115 + t204 * t116;
t109 = m(4) * t141 - t187 * mrSges(4,1) - t188 * mrSges(4,2) + t219;
t140 = t205 * t156 - t200 * t157;
t216 = -t188 * pkin(3) - t140;
t137 = -t187 * pkin(8) + t216;
t133 = t178 * t225 - t172 * pkin(4) + (-pkin(9) * t197 - pkin(8)) * t187 + t216;
t214 = m(6) * t133 - t146 * mrSges(6,1) + t147 * mrSges(6,2) - t166 * t161 + t167 * t162;
t210 = -m(5) * t137 + t172 * mrSges(5,1) - t171 * mrSges(5,2) - t176 * t225 + t177 * t224 - t214;
t121 = m(4) * t140 + t188 * mrSges(4,1) - t187 * mrSges(4,2) + t210;
t104 = t200 * t109 + t205 * t121;
t101 = m(3) * t159 + t193 * mrSges(3,1) - t191 * mrSges(3,2) + t104;
t220 = t205 * t109 - t200 * t121;
t102 = m(3) * t160 - t191 * mrSges(3,1) - t193 * mrSges(3,2) + t220;
t94 = t206 * t101 + t201 * t102;
t111 = t204 * t115 + t199 * t116;
t221 = -t201 * t101 + t206 * t102;
t148 = Ifges(6,5) * t167 + Ifges(6,6) * t166 + Ifges(6,3) * t194;
t118 = -mrSges(6,1) * t133 + mrSges(6,3) * t130 + Ifges(6,4) * t147 + Ifges(6,2) * t146 + Ifges(6,6) * t192 - t167 * t148 + t194 * t150;
t119 = mrSges(6,2) * t133 - mrSges(6,3) * t129 + Ifges(6,1) * t147 + Ifges(6,4) * t146 + Ifges(6,5) * t192 + t166 * t148 - t194 * t149;
t163 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t199 + Ifges(5,6) * t204) * t189;
t106 = mrSges(5,2) * t137 - mrSges(5,3) * t134 + Ifges(5,1) * t171 + Ifges(5,4) * t172 + Ifges(5,5) * qJDD(4) - pkin(9) * t117 - qJD(4) * t164 - t198 * t118 + t203 * t119 + t163 * t224;
t97 = -mrSges(5,1) * t137 + mrSges(5,3) * t135 + Ifges(5,4) * t171 + Ifges(5,2) * t172 + Ifges(5,6) * qJDD(4) - pkin(4) * t214 + pkin(9) * t218 + qJD(4) * t165 + t203 * t118 + t198 * t119 - t163 * t225;
t215 = mrSges(4,1) * t140 - mrSges(4,2) * t141 + Ifges(4,3) * t188 + pkin(3) * t210 + pkin(8) * t219 + t199 * t106 + t204 * t97;
t212 = mrSges(3,1) * t159 - mrSges(3,2) * t160 + Ifges(3,3) * t193 + pkin(2) * t104 + t215;
t211 = mrSges(2,1) * t183 - mrSges(2,2) * t182 + Ifges(2,3) * qJDD(1) + pkin(1) * t94 + t212;
t95 = mrSges(4,1) * g(1) + mrSges(4,3) * t141 + t187 * Ifges(4,5) + Ifges(4,6) * t188 - pkin(3) * t111 - t227;
t92 = m(2) * t183 + qJDD(1) * mrSges(2,1) - t208 * mrSges(2,2) + t94;
t91 = m(2) * t182 - t208 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t221;
t90 = -mrSges(4,2) * g(1) - mrSges(4,3) * t140 + Ifges(4,5) * t188 - t187 * Ifges(4,6) - pkin(8) * t111 + t204 * t106 - t199 * t97;
t89 = -mrSges(3,2) * g(1) - mrSges(3,3) * t159 + Ifges(3,5) * t193 - t191 * Ifges(3,6) - pkin(7) * t104 - t200 * t95 + t205 * t90;
t88 = Ifges(3,6) * t193 + t191 * Ifges(3,5) + mrSges(3,1) * g(1) + mrSges(3,3) * t160 + t200 * t90 + t205 * t95 - pkin(2) * (-m(4) * g(1) + t111) + pkin(7) * t220;
t87 = -mrSges(2,2) * g(1) - mrSges(2,3) * t183 + Ifges(2,5) * qJDD(1) - t208 * Ifges(2,6) - pkin(6) * t94 - t201 * t88 + t206 * t89;
t86 = Ifges(2,6) * qJDD(1) + t208 * Ifges(2,5) + mrSges(2,3) * t182 + t201 * t89 + t206 * t88 - pkin(1) * t111 + pkin(6) * t221 + (mrSges(2,1) - pkin(1) * (-m(3) - m(4))) * g(1);
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t211, t87, t89, t90, t106, t119; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t202 * t87 + t207 * t86 - pkin(5) * (t202 * t92 - t207 * t91), t86, t88, t95, t97, t118; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - t207 * t87 + t202 * t86 + pkin(5) * (t202 * t91 + t207 * t92), t211, t212, t215, t227, -t213;];
m_new = t1;
