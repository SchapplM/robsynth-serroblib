% Calculate vector of cutting torques with Newton-Euler for
% S5RRPPR2
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
% Datum: 2022-01-20 10:06
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPPR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR2_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR2_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR2_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:05:35
% EndTime: 2022-01-20 10:05:39
% DurationCPUTime: 3.67s
% Computational Cost: add. (55820->203), mult. (73380->265), div. (0->0), fcn. (40823->10), ass. (0->101)
t227 = 2 * qJD(4);
t191 = sin(qJ(1));
t194 = cos(qJ(1));
t170 = t191 * g(1) - t194 * g(2);
t164 = qJDD(1) * pkin(1) + t170;
t171 = -t194 * g(1) - t191 * g(2);
t195 = qJD(1) ^ 2;
t165 = -t195 * pkin(1) + t171;
t190 = sin(qJ(2));
t193 = cos(qJ(2));
t148 = t193 * t164 - t190 * t165;
t181 = qJDD(1) + qJDD(2);
t142 = t181 * pkin(2) + t148;
t149 = t190 * t164 + t193 * t165;
t182 = (qJD(1) + qJD(2));
t180 = t182 ^ 2;
t143 = -t180 * pkin(2) + t149;
t186 = sin(pkin(8));
t188 = cos(pkin(8));
t138 = t186 * t142 + t188 * t143;
t135 = -t180 * pkin(3) + t181 * qJ(4) + t138;
t226 = (t182 * t227) + t135;
t185 = sin(pkin(9));
t221 = t182 * t185;
t187 = cos(pkin(9));
t219 = t187 * t182;
t184 = -g(3) + qJDD(3);
t131 = t185 * t184 + t226 * t187;
t209 = -pkin(4) * t187 - pkin(7) * t185;
t160 = t209 * t182;
t129 = t160 * t219 + t131;
t137 = t188 * t142 - t186 * t143;
t202 = -t180 * qJ(4) + qJDD(4) - t137;
t132 = (-pkin(3) + t209) * t181 + t202;
t189 = sin(qJ(5));
t192 = cos(qJ(5));
t125 = -t189 * t129 + t192 * t132;
t168 = qJD(5) - t219;
t215 = t189 * t221;
t151 = -t168 * mrSges(6,2) - mrSges(6,3) * t215;
t153 = (mrSges(6,1) * t189 + mrSges(6,2) * t192) * t221;
t216 = qJD(5) * t182;
t155 = (t181 * t192 - t189 * t216) * t185;
t220 = t187 * t181;
t167 = qJDD(5) - t220;
t214 = t192 * t221;
t123 = m(6) * t125 + t167 * mrSges(6,1) - t155 * mrSges(6,3) + t168 * t151 - t153 * t214;
t126 = t192 * t129 + t189 * t132;
t152 = t168 * mrSges(6,1) - mrSges(6,3) * t214;
t154 = (-t181 * t189 - t192 * t216) * t185;
t124 = m(6) * t126 - t167 * mrSges(6,2) + t154 * mrSges(6,3) - t168 * t152 - t153 * t215;
t117 = -t189 * t123 + t192 * t124;
t218 = t187 * t184;
t128 = -t218 + (t135 + (t227 + t160) * t182) * t185;
t144 = Ifges(6,3) * t168 + (Ifges(6,5) * t192 - Ifges(6,6) * t189) * t221;
t146 = Ifges(6,5) * t168 + (Ifges(6,1) * t192 - Ifges(6,4) * t189) * t221;
t118 = -mrSges(6,1) * t128 + mrSges(6,3) * t126 + Ifges(6,4) * t155 + Ifges(6,2) * t154 + Ifges(6,6) * t167 - t144 * t214 + t168 * t146;
t145 = Ifges(6,6) * t168 + (Ifges(6,4) * t192 - Ifges(6,2) * t189) * t221;
t119 = mrSges(6,2) * t128 - mrSges(6,3) * t125 + Ifges(6,1) * t155 + Ifges(6,4) * t154 + Ifges(6,5) * t167 - t144 * t215 - t168 * t145;
t130 = -t226 * t185 + t218;
t203 = -m(6) * t128 + t154 * mrSges(6,1) - t155 * mrSges(6,2);
t205 = -t151 * t189 - t152 * t192;
t208 = Ifges(5,1) * t185 + Ifges(5,4) * t187;
t225 = -((Ifges(5,4) * t185 + Ifges(5,2) * t187) * t221 - t208 * t219) * t182 - mrSges(5,1) * t130 + mrSges(5,2) * t131 - pkin(4) * (t205 * t221 + t203) - pkin(7) * t117 - t192 * t118 - t189 * t119;
t224 = mrSges(5,2) * t185;
t156 = (-mrSges(5,1) * t187 + t224) * t182;
t222 = mrSges(5,3) * t181;
t114 = m(5) * t131 + (t156 * t182 + t222) * t187 + t117;
t121 = m(5) * t130 + (-t222 + (-t156 + t205) * t182) * t185 + t203;
t210 = t187 * t114 - t185 * t121;
t107 = m(4) * t138 - t180 * mrSges(4,1) - t181 * mrSges(4,2) + t210;
t116 = t192 * t123 + t189 * t124;
t134 = -t181 * pkin(3) + t202;
t200 = -m(5) * t134 + mrSges(5,1) * t220 - t116 + (t185 ^ 2 + t187 ^ 2) * mrSges(5,3) * t180;
t111 = m(4) * t137 - t180 * mrSges(4,2) + (mrSges(4,1) - t224) * t181 + t200;
t100 = t186 * t107 + t188 * t111;
t97 = m(3) * t148 + t181 * mrSges(3,1) - t180 * mrSges(3,2) + t100;
t211 = t188 * t107 - t186 * t111;
t98 = m(3) * t149 - t180 * mrSges(3,1) - t181 * mrSges(3,2) + t211;
t92 = t190 * t98 + t193 * t97;
t109 = t185 * t114 + t187 * t121;
t213 = m(4) * t184 + t109;
t212 = -t190 * t97 + t193 * t98;
t207 = Ifges(5,5) * t185 + Ifges(5,6) * t187;
t206 = t145 * t192 + t146 * t189;
t157 = t207 * t182;
t102 = mrSges(5,2) * t134 - mrSges(5,3) * t130 - pkin(7) * t116 - t189 * t118 + t192 * t119 + t157 * t219 + t208 * t181;
t199 = mrSges(6,1) * t125 - mrSges(6,2) * t126 + Ifges(6,5) * t155 + Ifges(6,6) * t154 + Ifges(6,3) * t167;
t104 = Ifges(5,2) * t220 - mrSges(5,1) * t134 + mrSges(5,3) * t131 - pkin(4) * t116 + (Ifges(5,4) * t181 + (-t157 - t206) * t182) * t185 - t199;
t201 = -mrSges(4,2) * t138 + qJ(4) * t210 + t185 * t102 + t187 * t104 + pkin(3) * (-t181 * t224 + t200) + mrSges(4,1) * t137 + Ifges(4,3) * t181;
t198 = mrSges(3,1) * t148 - mrSges(3,2) * t149 + Ifges(3,3) * t181 + pkin(2) * t100 + t201;
t196 = mrSges(2,1) * t170 - mrSges(2,2) * t171 + Ifges(2,3) * qJDD(1) + pkin(1) * t92 + t198;
t93 = -mrSges(4,1) * t184 + mrSges(4,3) * t138 + t180 * Ifges(4,5) - pkin(3) * t109 + (Ifges(4,6) - t207) * t181 + t225;
t90 = m(2) * t171 - t195 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t212;
t89 = m(2) * t170 + qJDD(1) * mrSges(2,1) - t195 * mrSges(2,2) + t92;
t88 = mrSges(4,2) * t184 - mrSges(4,3) * t137 + Ifges(4,5) * t181 - t180 * Ifges(4,6) - qJ(4) * t109 + t187 * t102 - t185 * t104;
t87 = -mrSges(3,2) * g(3) - mrSges(3,3) * t148 + Ifges(3,5) * t181 - t180 * Ifges(3,6) - qJ(3) * t100 - t186 * t93 + t188 * t88;
t86 = mrSges(3,1) * g(3) + mrSges(3,3) * t149 + t180 * Ifges(3,5) + Ifges(3,6) * t181 - pkin(2) * t213 + qJ(3) * t211 + t186 * t88 + t188 * t93;
t85 = -mrSges(2,2) * g(3) - mrSges(2,3) * t170 + Ifges(2,5) * qJDD(1) - t195 * Ifges(2,6) - pkin(6) * t92 - t190 * t86 + t193 * t87;
t84 = Ifges(2,6) * qJDD(1) + t195 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t171 + t190 * t87 + t193 * t86 - pkin(1) * (-m(3) * g(3) + t213) + pkin(6) * t212;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t194 * t85 - t191 * t84 - pkin(5) * (t191 * t90 + t194 * t89), t85, t87, t88, t102, t119; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t191 * t85 + t194 * t84 + pkin(5) * (-t191 * t89 + t194 * t90), t84, t86, t93, t104, t118; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t196, t196, t198, t201, t207 * t181 - t225, t206 * t221 + t199;];
m_new = t1;
