% Calculate vector of cutting torques with Newton-Euler for
% S5PRRPR2
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRRPR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR2_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR2_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR2_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR2_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR2_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:17:24
% EndTime: 2019-12-05 16:17:28
% DurationCPUTime: 3.19s
% Computational Cost: add. (46066->192), mult. (63628->254), div. (0->0), fcn. (40823->10), ass. (0->99)
t223 = 2 * qJD(4);
t182 = sin(pkin(8));
t184 = cos(pkin(8));
t165 = t182 * g(1) - t184 * g(2);
t167 = -t184 * g(1) - t182 * g(2);
t187 = sin(qJ(2));
t190 = cos(qJ(2));
t148 = t190 * t165 - t187 * t167;
t145 = qJDD(2) * pkin(2) + t148;
t149 = t187 * t165 + t190 * t167;
t191 = qJD(2) ^ 2;
t146 = -t191 * pkin(2) + t149;
t186 = sin(qJ(3));
t189 = cos(qJ(3));
t138 = t186 * t145 + t189 * t146;
t178 = (qJD(2) + qJD(3));
t176 = t178 ^ 2;
t177 = qJDD(2) + qJDD(3);
t135 = -t176 * pkin(3) + t177 * qJ(4) + t138;
t222 = (t178 * t223) + t135;
t181 = sin(pkin(9));
t217 = t178 * t181;
t183 = cos(pkin(9));
t215 = t183 * t178;
t180 = -g(3) + qJDD(1);
t131 = t181 * t180 + t222 * t183;
t205 = -pkin(4) * t183 - pkin(7) * t181;
t160 = t205 * t178;
t129 = t160 * t215 + t131;
t137 = t189 * t145 - t186 * t146;
t198 = -t176 * qJ(4) + qJDD(4) - t137;
t132 = (-pkin(3) + t205) * t177 + t198;
t185 = sin(qJ(5));
t188 = cos(qJ(5));
t126 = -t185 * t129 + t188 * t132;
t168 = qJD(5) - t215;
t211 = t185 * t217;
t151 = -t168 * mrSges(6,2) - mrSges(6,3) * t211;
t153 = (mrSges(6,1) * t185 + mrSges(6,2) * t188) * t217;
t212 = qJD(5) * t178;
t155 = (t177 * t188 - t185 * t212) * t181;
t216 = t183 * t177;
t166 = qJDD(5) - t216;
t210 = t188 * t217;
t123 = m(6) * t126 + t166 * mrSges(6,1) - t155 * mrSges(6,3) + t168 * t151 - t153 * t210;
t127 = t188 * t129 + t185 * t132;
t152 = t168 * mrSges(6,1) - mrSges(6,3) * t210;
t154 = (-t177 * t185 - t188 * t212) * t181;
t124 = m(6) * t127 - t166 * mrSges(6,2) + t154 * mrSges(6,3) - t168 * t152 - t153 * t211;
t117 = -t185 * t123 + t188 * t124;
t214 = t183 * t180;
t128 = -t214 + (t135 + (t223 + t160) * t178) * t181;
t139 = Ifges(6,3) * t168 + (Ifges(6,5) * t188 - Ifges(6,6) * t185) * t217;
t141 = Ifges(6,5) * t168 + (Ifges(6,1) * t188 - Ifges(6,4) * t185) * t217;
t118 = -mrSges(6,1) * t128 + mrSges(6,3) * t127 + Ifges(6,4) * t155 + Ifges(6,2) * t154 + Ifges(6,6) * t166 - t139 * t210 + t168 * t141;
t140 = Ifges(6,6) * t168 + (Ifges(6,4) * t188 - Ifges(6,2) * t185) * t217;
t119 = mrSges(6,2) * t128 - mrSges(6,3) * t126 + Ifges(6,1) * t155 + Ifges(6,4) * t154 + Ifges(6,5) * t166 - t139 * t211 - t168 * t140;
t130 = -t222 * t181 + t214;
t199 = -m(6) * t128 + t154 * mrSges(6,1) - t155 * mrSges(6,2);
t201 = -t151 * t185 - t152 * t188;
t204 = Ifges(5,1) * t181 + Ifges(5,4) * t183;
t221 = -((Ifges(5,4) * t181 + Ifges(5,2) * t183) * t217 - t204 * t215) * t178 - mrSges(5,1) * t130 + mrSges(5,2) * t131 - pkin(4) * (t201 * t217 + t199) - pkin(7) * t117 - t188 * t118 - t185 * t119;
t220 = mrSges(5,2) * t181;
t156 = (-mrSges(5,1) * t183 + t220) * t178;
t218 = mrSges(5,3) * t177;
t114 = m(5) * t131 + (t156 * t178 + t218) * t183 + t117;
t121 = m(5) * t130 + (-t218 + (-t156 + t201) * t178) * t181 + t199;
t206 = t183 * t114 - t181 * t121;
t107 = m(4) * t138 - t176 * mrSges(4,1) - t177 * mrSges(4,2) + t206;
t116 = t188 * t123 + t185 * t124;
t134 = -t177 * pkin(3) + t198;
t196 = -m(5) * t134 + mrSges(5,1) * t216 - t116 + (t181 ^ 2 + t183 ^ 2) * mrSges(5,3) * t176;
t111 = m(4) * t137 - t176 * mrSges(4,2) + (mrSges(4,1) - t220) * t177 + t196;
t100 = t186 * t107 + t189 * t111;
t97 = m(3) * t148 + qJDD(2) * mrSges(3,1) - t191 * mrSges(3,2) + t100;
t207 = t189 * t107 - t186 * t111;
t98 = m(3) * t149 - t191 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t207;
t92 = t187 * t98 + t190 * t97;
t109 = t181 * t114 + t183 * t121;
t209 = m(4) * t180 + t109;
t208 = -t187 * t97 + t190 * t98;
t203 = Ifges(5,5) * t181 + Ifges(5,6) * t183;
t202 = t140 * t188 + t141 * t185;
t157 = t203 * t178;
t102 = mrSges(5,2) * t134 - mrSges(5,3) * t130 - pkin(7) * t116 - t185 * t118 + t188 * t119 + t157 * t215 + t204 * t177;
t195 = mrSges(6,1) * t126 - mrSges(6,2) * t127 + Ifges(6,5) * t155 + Ifges(6,6) * t154 + Ifges(6,3) * t166;
t104 = Ifges(5,2) * t216 - mrSges(5,1) * t134 + mrSges(5,3) * t131 - pkin(4) * t116 + (Ifges(5,4) * t177 + (-t157 - t202) * t178) * t181 - t195;
t197 = -mrSges(4,2) * t138 + qJ(4) * t206 + t181 * t102 + t183 * t104 + pkin(3) * (-t177 * t220 + t196) + mrSges(4,1) * t137 + Ifges(4,3) * t177;
t194 = mrSges(3,1) * t148 - mrSges(3,2) * t149 + Ifges(3,3) * qJDD(2) + pkin(2) * t100 + t197;
t192 = mrSges(2,1) * t165 - mrSges(2,2) * t167 + pkin(1) * t92 + t194;
t93 = -mrSges(4,1) * t180 + mrSges(4,3) * t138 + t176 * Ifges(4,5) - pkin(3) * t109 + (Ifges(4,6) - t203) * t177 + t221;
t90 = m(2) * t167 + t208;
t89 = m(2) * t165 + t92;
t88 = mrSges(4,2) * t180 - mrSges(4,3) * t137 + Ifges(4,5) * t177 - t176 * Ifges(4,6) - qJ(4) * t109 + t183 * t102 - t181 * t104;
t87 = mrSges(3,2) * t180 - mrSges(3,3) * t148 + Ifges(3,5) * qJDD(2) - t191 * Ifges(3,6) - pkin(6) * t100 - t186 * t93 + t189 * t88;
t86 = -mrSges(3,1) * t180 + mrSges(3,3) * t149 + t191 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t209 + pkin(6) * t207 + t186 * t88 + t189 * t93;
t85 = mrSges(2,2) * t180 - mrSges(2,3) * t165 - pkin(5) * t92 - t187 * t86 + t190 * t87;
t84 = -mrSges(2,1) * t180 + mrSges(2,3) * t167 + t187 * t87 + t190 * t86 - pkin(1) * (m(3) * t180 + t209) + pkin(5) * t208;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t184 * t85 - t182 * t84 - qJ(1) * (t182 * t90 + t184 * t89), t85, t87, t88, t102, t119; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t182 * t85 + t184 * t84 + qJ(1) * (-t182 * t89 + t184 * t90), t84, t86, t93, t104, t118; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t192, t192, t194, t197, t203 * t177 - t221, t202 * t217 + t195;];
m_new = t1;
