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
% Datum: 2022-01-20 10:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 10:48:13
% EndTime: 2022-01-20 10:48:19
% DurationCPUTime: 4.62s
% Computational Cost: add. (75833->225), mult. (97595->288), div. (0->0), fcn. (55165->10), ass. (0->98)
t204 = sin(qJ(1));
t208 = cos(qJ(1));
t182 = t204 * g(1) - t208 * g(2);
t176 = qJDD(1) * pkin(1) + t182;
t183 = -t208 * g(1) - t204 * g(2);
t209 = qJD(1) ^ 2;
t177 = -t209 * pkin(1) + t183;
t203 = sin(qJ(2));
t207 = cos(qJ(2));
t159 = t207 * t176 - t203 * t177;
t193 = qJDD(1) + qJDD(2);
t156 = t193 * pkin(2) + t159;
t160 = t203 * t176 + t207 * t177;
t195 = qJD(1) + qJD(2);
t191 = t195 ^ 2;
t157 = -t191 * pkin(2) + t160;
t199 = sin(pkin(9));
t200 = cos(pkin(9));
t141 = t199 * t156 + t200 * t157;
t138 = -t191 * pkin(3) + t193 * pkin(7) + t141;
t198 = -g(3) + qJDD(3);
t202 = sin(qJ(4));
t206 = cos(qJ(4));
t134 = -t202 * t138 + t206 * t198;
t225 = qJD(4) * t195;
t223 = t206 * t225;
t171 = t202 * t193 + t223;
t131 = (-t171 + t223) * pkin(8) + (t191 * t202 * t206 + qJDD(4)) * pkin(4) + t134;
t135 = t206 * t138 + t202 * t198;
t172 = t206 * t193 - t202 * t225;
t227 = t195 * t202;
t180 = qJD(4) * pkin(4) - pkin(8) * t227;
t197 = t206 ^ 2;
t132 = -t197 * t191 * pkin(4) + t172 * pkin(8) - qJD(4) * t180 + t135;
t201 = sin(qJ(5));
t205 = cos(qJ(5));
t129 = t205 * t131 - t201 * t132;
t166 = (-t201 * t202 + t205 * t206) * t195;
t147 = t166 * qJD(5) + t205 * t171 + t201 * t172;
t167 = (t201 * t206 + t202 * t205) * t195;
t152 = -t166 * mrSges(6,1) + t167 * mrSges(6,2);
t194 = qJD(4) + qJD(5);
t161 = -t194 * mrSges(6,2) + t166 * mrSges(6,3);
t192 = qJDD(4) + qJDD(5);
t126 = m(6) * t129 + t192 * mrSges(6,1) - t147 * mrSges(6,3) - t167 * t152 + t194 * t161;
t130 = t201 * t131 + t205 * t132;
t146 = -t167 * qJD(5) - t201 * t171 + t205 * t172;
t162 = t194 * mrSges(6,1) - t167 * mrSges(6,3);
t127 = m(6) * t130 - t192 * mrSges(6,2) + t146 * mrSges(6,3) + t166 * t152 - t194 * t162;
t117 = t205 * t126 + t201 * t127;
t164 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t202 + Ifges(5,2) * t206) * t195;
t165 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t202 + Ifges(5,4) * t206) * t195;
t149 = Ifges(6,4) * t167 + Ifges(6,2) * t166 + Ifges(6,6) * t194;
t150 = Ifges(6,1) * t167 + Ifges(6,4) * t166 + Ifges(6,5) * t194;
t214 = -mrSges(6,1) * t129 + mrSges(6,2) * t130 - Ifges(6,5) * t147 - Ifges(6,6) * t146 - Ifges(6,3) * t192 - t167 * t149 + t166 * t150;
t228 = mrSges(5,1) * t134 - mrSges(5,2) * t135 + Ifges(5,5) * t171 + Ifges(5,6) * t172 + Ifges(5,3) * qJDD(4) + pkin(4) * t117 + (t202 * t164 - t206 * t165) * t195 - t214;
t170 = (-mrSges(5,1) * t206 + mrSges(5,2) * t202) * t195;
t226 = t195 * t206;
t179 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t226;
t115 = m(5) * t134 + qJDD(4) * mrSges(5,1) - t171 * mrSges(5,3) + qJD(4) * t179 - t170 * t227 + t117;
t178 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t227;
t219 = -t201 * t126 + t205 * t127;
t116 = m(5) * t135 - qJDD(4) * mrSges(5,2) + t172 * mrSges(5,3) - qJD(4) * t178 + t170 * t226 + t219;
t220 = -t202 * t115 + t206 * t116;
t109 = m(4) * t141 - t191 * mrSges(4,1) - t193 * mrSges(4,2) + t220;
t140 = t200 * t156 - t199 * t157;
t217 = -t193 * pkin(3) - t140;
t137 = -t191 * pkin(7) + t217;
t133 = t180 * t227 - t172 * pkin(4) + (-pkin(8) * t197 - pkin(7)) * t191 + t217;
t215 = m(6) * t133 - t146 * mrSges(6,1) + t147 * mrSges(6,2) - t166 * t161 + t167 * t162;
t211 = -m(5) * t137 + t172 * mrSges(5,1) - t171 * mrSges(5,2) - t178 * t227 + t179 * t226 - t215;
t121 = m(4) * t140 + t193 * mrSges(4,1) - t191 * mrSges(4,2) + t211;
t104 = t199 * t109 + t200 * t121;
t101 = m(3) * t159 + t193 * mrSges(3,1) - t191 * mrSges(3,2) + t104;
t221 = t200 * t109 - t199 * t121;
t102 = m(3) * t160 - t191 * mrSges(3,1) - t193 * mrSges(3,2) + t221;
t94 = t207 * t101 + t203 * t102;
t111 = t206 * t115 + t202 * t116;
t224 = m(4) * t198 + t111;
t222 = -t203 * t101 + t207 * t102;
t148 = Ifges(6,5) * t167 + Ifges(6,6) * t166 + Ifges(6,3) * t194;
t118 = -mrSges(6,1) * t133 + mrSges(6,3) * t130 + Ifges(6,4) * t147 + Ifges(6,2) * t146 + Ifges(6,6) * t192 - t167 * t148 + t194 * t150;
t119 = mrSges(6,2) * t133 - mrSges(6,3) * t129 + Ifges(6,1) * t147 + Ifges(6,4) * t146 + Ifges(6,5) * t192 + t166 * t148 - t194 * t149;
t163 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t202 + Ifges(5,6) * t206) * t195;
t106 = mrSges(5,2) * t137 - mrSges(5,3) * t134 + Ifges(5,1) * t171 + Ifges(5,4) * t172 + Ifges(5,5) * qJDD(4) - pkin(8) * t117 - qJD(4) * t164 - t201 * t118 + t205 * t119 + t163 * t226;
t97 = -mrSges(5,1) * t137 + mrSges(5,3) * t135 + Ifges(5,4) * t171 + Ifges(5,2) * t172 + Ifges(5,6) * qJDD(4) - pkin(4) * t215 + pkin(8) * t219 + qJD(4) * t165 + t205 * t118 + t201 * t119 - t163 * t227;
t216 = mrSges(4,1) * t140 - mrSges(4,2) * t141 + Ifges(4,3) * t193 + pkin(3) * t211 + pkin(7) * t220 + t202 * t106 + t206 * t97;
t213 = mrSges(3,1) * t159 - mrSges(3,2) * t160 + Ifges(3,3) * t193 + pkin(2) * t104 + t216;
t212 = mrSges(2,1) * t182 - mrSges(2,2) * t183 + Ifges(2,3) * qJDD(1) + pkin(1) * t94 + t213;
t95 = -mrSges(4,1) * t198 + mrSges(4,3) * t141 + t191 * Ifges(4,5) + Ifges(4,6) * t193 - pkin(3) * t111 - t228;
t92 = m(2) * t183 - t209 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t222;
t91 = m(2) * t182 + qJDD(1) * mrSges(2,1) - t209 * mrSges(2,2) + t94;
t90 = mrSges(4,2) * t198 - mrSges(4,3) * t140 + Ifges(4,5) * t193 - t191 * Ifges(4,6) - pkin(7) * t111 + t206 * t106 - t202 * t97;
t89 = -mrSges(3,2) * g(3) - mrSges(3,3) * t159 + Ifges(3,5) * t193 - t191 * Ifges(3,6) - qJ(3) * t104 - t199 * t95 + t200 * t90;
t88 = mrSges(3,1) * g(3) + mrSges(3,3) * t160 + t191 * Ifges(3,5) + Ifges(3,6) * t193 - pkin(2) * t224 + qJ(3) * t221 + t199 * t90 + t200 * t95;
t87 = -mrSges(2,2) * g(3) - mrSges(2,3) * t182 + Ifges(2,5) * qJDD(1) - t209 * Ifges(2,6) - pkin(6) * t94 - t203 * t88 + t207 * t89;
t86 = Ifges(2,6) * qJDD(1) + t209 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t183 + t203 * t89 + t207 * t88 - pkin(1) * (-m(3) * g(3) + t224) + pkin(6) * t222;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t208 * t87 - t204 * t86 - pkin(5) * (t204 * t92 + t208 * t91), t87, t89, t90, t106, t119; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t204 * t87 + t208 * t86 + pkin(5) * (-t204 * t91 + t208 * t92), t86, t88, t95, t97, t118; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t212, t212, t213, t216, t228, -t214;];
m_new = t1;
