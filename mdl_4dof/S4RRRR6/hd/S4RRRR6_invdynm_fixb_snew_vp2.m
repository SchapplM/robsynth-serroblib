% Calculate vector of cutting torques with Newton-Euler for
% S4RRRR6
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% m [3x5]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RRRR6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(8,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR6_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR6_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR6_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR6_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR6_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:29:28
% EndTime: 2019-12-31 17:29:38
% DurationCPUTime: 5.18s
% Computational Cost: add. (64747->249), mult. (140137->336), div. (0->0), fcn. (103859->10), ass. (0->109)
t195 = sin(pkin(4));
t199 = sin(qJ(2));
t203 = cos(qJ(2));
t214 = qJD(1) * qJD(2);
t183 = (-qJDD(1) * t203 + t199 * t214) * t195;
t223 = pkin(6) * t195;
t196 = cos(pkin(4));
t222 = t196 * g(3);
t221 = t195 * t199;
t220 = t195 * t203;
t219 = t196 * t199;
t218 = t196 * t203;
t216 = qJD(1) * t195;
t181 = (-pkin(2) * t203 - pkin(7) * t199) * t216;
t191 = t196 * qJD(1) + qJD(2);
t189 = t191 ^ 2;
t190 = t196 * qJDD(1) + qJDD(2);
t215 = qJD(1) * t203;
t200 = sin(qJ(1));
t204 = cos(qJ(1));
t187 = t200 * g(1) - t204 * g(2);
t205 = qJD(1) ^ 2;
t178 = qJDD(1) * pkin(1) + t205 * t223 + t187;
t188 = -t204 * g(1) - t200 * g(2);
t179 = -t205 * pkin(1) + qJDD(1) * t223 + t188;
t217 = t178 * t219 + t203 * t179;
t143 = -t189 * pkin(2) + t190 * pkin(7) + (-g(3) * t199 + t181 * t215) * t195 + t217;
t182 = (qJDD(1) * t199 + t203 * t214) * t195;
t144 = t183 * pkin(2) - t182 * pkin(7) - t222 + (-t178 + (pkin(2) * t199 - pkin(7) * t203) * t191 * qJD(1)) * t195;
t198 = sin(qJ(3));
t202 = cos(qJ(3));
t132 = t202 * t143 + t198 * t144;
t213 = t199 * t216;
t171 = t202 * t191 - t198 * t213;
t172 = t198 * t191 + t202 * t213;
t159 = -t171 * pkin(3) - t172 * pkin(8);
t175 = qJDD(3) + t183;
t212 = t195 * t215;
t186 = qJD(3) - t212;
t185 = t186 ^ 2;
t129 = -t185 * pkin(3) + t175 * pkin(8) + t171 * t159 + t132;
t156 = -g(3) * t220 + t178 * t218 - t199 * t179;
t142 = -t190 * pkin(2) - t189 * pkin(7) + t181 * t213 - t156;
t154 = -t172 * qJD(3) - t198 * t182 + t202 * t190;
t155 = t171 * qJD(3) + t202 * t182 + t198 * t190;
t130 = (-t171 * t186 - t155) * pkin(8) + (t172 * t186 - t154) * pkin(3) + t142;
t197 = sin(qJ(4));
t201 = cos(qJ(4));
t126 = -t197 * t129 + t201 * t130;
t160 = -t197 * t172 + t201 * t186;
t135 = t160 * qJD(4) + t201 * t155 + t197 * t175;
t161 = t201 * t172 + t197 * t186;
t145 = -t160 * mrSges(5,1) + t161 * mrSges(5,2);
t170 = qJD(4) - t171;
t146 = -t170 * mrSges(5,2) + t160 * mrSges(5,3);
t152 = qJDD(4) - t154;
t123 = m(5) * t126 + t152 * mrSges(5,1) - t135 * mrSges(5,3) - t161 * t145 + t170 * t146;
t127 = t201 * t129 + t197 * t130;
t134 = -t161 * qJD(4) - t197 * t155 + t201 * t175;
t147 = t170 * mrSges(5,1) - t161 * mrSges(5,3);
t124 = m(5) * t127 - t152 * mrSges(5,2) + t134 * mrSges(5,3) + t160 * t145 - t170 * t147;
t117 = -t197 * t123 + t201 * t124;
t158 = -t171 * mrSges(4,1) + t172 * mrSges(4,2);
t163 = t186 * mrSges(4,1) - t172 * mrSges(4,3);
t115 = m(4) * t132 - t175 * mrSges(4,2) + t154 * mrSges(4,3) + t171 * t158 - t186 * t163 + t117;
t131 = -t198 * t143 + t202 * t144;
t128 = -t175 * pkin(3) - t185 * pkin(8) + t172 * t159 - t131;
t125 = -m(5) * t128 + t134 * mrSges(5,1) - t135 * mrSges(5,2) + t160 * t146 - t161 * t147;
t162 = -t186 * mrSges(4,2) + t171 * mrSges(4,3);
t121 = m(4) * t131 + t175 * mrSges(4,1) - t155 * mrSges(4,3) - t172 * t158 + t186 * t162 + t125;
t110 = t198 * t115 + t202 * t121;
t157 = -g(3) * t221 + t217;
t176 = t191 * mrSges(3,1) - mrSges(3,3) * t213;
t180 = (-mrSges(3,1) * t203 + mrSges(3,2) * t199) * t216;
t211 = t202 * t115 - t198 * t121;
t108 = m(3) * t157 - t190 * mrSges(3,2) - t183 * mrSges(3,3) - t191 * t176 + t180 * t212 + t211;
t177 = -t191 * mrSges(3,2) + mrSges(3,3) * t212;
t116 = t201 * t123 + t197 * t124;
t208 = -m(4) * t142 + t154 * mrSges(4,1) - t155 * mrSges(4,2) + t171 * t162 - t172 * t163 - t116;
t112 = m(3) * t156 + t190 * mrSges(3,1) - t182 * mrSges(3,3) + t191 * t177 - t180 * t213 + t208;
t103 = t203 * t108 - t199 * t112;
t167 = -t195 * t178 - t222;
t109 = m(3) * t167 + t183 * mrSges(3,1) + t182 * mrSges(3,2) + (t176 * t199 - t177 * t203) * t216 + t110;
t100 = t108 * t219 - t195 * t109 + t112 * t218;
t136 = Ifges(5,5) * t161 + Ifges(5,6) * t160 + Ifges(5,3) * t170;
t138 = Ifges(5,1) * t161 + Ifges(5,4) * t160 + Ifges(5,5) * t170;
t118 = -mrSges(5,1) * t128 + mrSges(5,3) * t127 + Ifges(5,4) * t135 + Ifges(5,2) * t134 + Ifges(5,6) * t152 - t161 * t136 + t170 * t138;
t137 = Ifges(5,4) * t161 + Ifges(5,2) * t160 + Ifges(5,6) * t170;
t119 = mrSges(5,2) * t128 - mrSges(5,3) * t126 + Ifges(5,1) * t135 + Ifges(5,4) * t134 + Ifges(5,5) * t152 + t160 * t136 - t170 * t137;
t148 = Ifges(4,5) * t172 + Ifges(4,6) * t171 + Ifges(4,3) * t186;
t149 = Ifges(4,4) * t172 + Ifges(4,2) * t171 + Ifges(4,6) * t186;
t104 = mrSges(4,2) * t142 - mrSges(4,3) * t131 + Ifges(4,1) * t155 + Ifges(4,4) * t154 + Ifges(4,5) * t175 - pkin(8) * t116 - t197 * t118 + t201 * t119 + t171 * t148 - t186 * t149;
t150 = Ifges(4,1) * t172 + Ifges(4,4) * t171 + Ifges(4,5) * t186;
t207 = mrSges(5,1) * t126 - mrSges(5,2) * t127 + Ifges(5,5) * t135 + Ifges(5,6) * t134 + Ifges(5,3) * t152 + t161 * t137 - t160 * t138;
t105 = -mrSges(4,1) * t142 + mrSges(4,3) * t132 + Ifges(4,4) * t155 + Ifges(4,2) * t154 + Ifges(4,6) * t175 - pkin(3) * t116 - t172 * t148 + t186 * t150 - t207;
t165 = Ifges(3,6) * t191 + (Ifges(3,4) * t199 + Ifges(3,2) * t203) * t216;
t166 = Ifges(3,5) * t191 + (Ifges(3,1) * t199 + Ifges(3,4) * t203) * t216;
t92 = Ifges(3,5) * t182 - Ifges(3,6) * t183 + Ifges(3,3) * t190 + mrSges(3,1) * t156 - mrSges(3,2) * t157 + t198 * t104 + t202 * t105 + pkin(2) * t208 + pkin(7) * t211 + (t165 * t199 - t166 * t203) * t216;
t164 = Ifges(3,3) * t191 + (Ifges(3,5) * t199 + Ifges(3,6) * t203) * t216;
t94 = mrSges(3,2) * t167 - mrSges(3,3) * t156 + Ifges(3,1) * t182 - Ifges(3,4) * t183 + Ifges(3,5) * t190 - pkin(7) * t110 + t202 * t104 - t198 * t105 + t164 * t212 - t191 * t165;
t206 = mrSges(4,1) * t131 - mrSges(4,2) * t132 + Ifges(4,5) * t155 + Ifges(4,6) * t154 + Ifges(4,3) * t175 + pkin(3) * t125 + pkin(8) * t117 + t201 * t118 + t197 * t119 + t172 * t149 - t171 * t150;
t96 = -mrSges(3,1) * t167 + mrSges(3,3) * t157 + Ifges(3,4) * t182 - Ifges(3,2) * t183 + Ifges(3,6) * t190 - pkin(2) * t110 - t164 * t213 + t191 * t166 - t206;
t209 = mrSges(2,1) * t187 - mrSges(2,2) * t188 + Ifges(2,3) * qJDD(1) + pkin(1) * t100 + t103 * t223 + t196 * t92 + t96 * t220 + t94 * t221;
t101 = m(2) * t188 - t205 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t103;
t99 = t196 * t109 + (t108 * t199 + t112 * t203) * t195;
t97 = m(2) * t187 + qJDD(1) * mrSges(2,1) - t205 * mrSges(2,2) + t100;
t90 = -mrSges(2,2) * g(3) - mrSges(2,3) * t187 + Ifges(2,5) * qJDD(1) - t205 * Ifges(2,6) - t199 * t96 + t203 * t94 + (-t100 * t196 - t195 * t99) * pkin(6);
t89 = mrSges(2,1) * g(3) + mrSges(2,3) * t188 + t205 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t99 - t195 * t92 + (pkin(6) * t103 + t199 * t94 + t203 * t96) * t196;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t204 * t90 - t200 * t89 - pkin(5) * (t200 * t101 + t204 * t97), t90, t94, t104, t119; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t200 * t90 + t204 * t89 + pkin(5) * (t204 * t101 - t200 * t97), t89, t96, t105, t118; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t209, t209, t92, t206, t207;];
m_new = t1;
