% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRRR14V3
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% m [3x7]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 03:56
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRRR14V3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(1,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14V3_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR14V3_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR14V3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_invdynm_fixb_snew_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14V3_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR14V3_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR14V3_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 03:48:26
% EndTime: 2019-05-07 03:48:33
% DurationCPUTime: 2.67s
% Computational Cost: add. (30805->278), mult. (58333->349), div. (0->0), fcn. (42155->10), ass. (0->95)
t209 = sin(qJ(2));
t214 = cos(qJ(2));
t226 = qJD(1) * qJD(2);
t225 = t214 * t226;
t194 = t209 * qJDD(1) + t225;
t210 = sin(qJ(1));
t215 = cos(qJ(1));
t198 = t210 * g(1) - t215 * g(2);
t228 = qJD(1) * t209;
t162 = -0.2e1 * qJD(3) * t228 + qJ(3) * (-t194 - t225) - t198;
t199 = -t215 * g(1) - t210 * g(2);
t177 = -t209 * g(3) + t214 * t199;
t216 = qJD(1) ^ 2;
t232 = t214 * t216;
t166 = 0.2e1 * qJD(3) * qJD(2) + (-t209 * t232 + qJDD(2)) * qJ(3) + t177;
t208 = sin(qJ(4));
t213 = cos(qJ(4));
t145 = t208 * t162 + t213 * t166;
t176 = -t214 * g(3) - t209 * t199;
t171 = qJDD(3) + (-t209 ^ 2 * t216 - qJD(2) ^ 2) * qJ(3) - t176;
t207 = sin(qJ(5));
t212 = cos(qJ(5));
t136 = t212 * t145 + t207 * t171;
t144 = -t213 * t162 + t208 * t166;
t206 = sin(qJ(6));
t211 = cos(qJ(6));
t128 = -t206 * t136 + t211 * t144;
t191 = t213 * qJD(2) - t208 * t228;
t164 = t191 * qJD(4) + t208 * qJDD(2) + t213 * t194;
t192 = t208 * qJD(2) + t213 * t228;
t227 = t214 * qJD(1);
t200 = qJD(4) - t227;
t172 = -t207 * t192 + t212 * t200;
t195 = t214 * qJDD(1) - t209 * t226;
t190 = qJDD(4) - t195;
t140 = t172 * qJD(5) + t212 * t164 + t207 * t190;
t173 = t212 * t192 + t207 * t200;
t185 = qJD(5) - t191;
t151 = -t206 * t173 + t211 * t185;
t163 = -t192 * qJD(4) + t213 * qJDD(2) - t208 * t194;
t161 = qJDD(5) - t163;
t131 = t151 * qJD(6) + t211 * t140 + t206 * t161;
t152 = t211 * t173 + t206 * t185;
t137 = -t151 * mrSges(7,1) + t152 * mrSges(7,2);
t139 = -t173 * qJD(5) - t207 * t164 + t212 * t190;
t138 = qJDD(6) - t139;
t170 = qJD(6) - t172;
t141 = -t170 * mrSges(7,2) + t151 * mrSges(7,3);
t126 = m(7) * t128 + t138 * mrSges(7,1) - t131 * mrSges(7,3) - t152 * t137 + t170 * t141;
t129 = t211 * t136 + t206 * t144;
t130 = -t152 * qJD(6) - t206 * t140 + t211 * t161;
t142 = t170 * mrSges(7,1) - t152 * mrSges(7,3);
t127 = m(7) * t129 - t138 * mrSges(7,2) + t130 * mrSges(7,3) + t151 * t137 - t170 * t142;
t135 = t207 * t145 - t212 * t171;
t149 = -t172 * mrSges(6,1) + t173 * mrSges(6,2);
t153 = -t185 * mrSges(6,2) + t172 * mrSges(6,3);
t154 = t185 * mrSges(6,1) - t173 * mrSges(6,3);
t169 = -t191 * mrSges(5,1) + t192 * mrSges(5,2);
t116 = m(5) * t145 + t163 * mrSges(5,3) - t190 * mrSges(5,2) + t191 * t169 - t200 * (t200 * mrSges(5,1) - t192 * mrSges(5,3)) + t212 * (m(6) * t136 - t161 * mrSges(6,2) + t139 * mrSges(6,3) - t206 * t126 + t211 * t127 + t172 * t149 - t185 * t154) - t207 * (t161 * mrSges(6,1) + t130 * mrSges(7,1) - t131 * mrSges(7,2) - t140 * mrSges(6,3) + t151 * t141 - t152 * t142 - t173 * t149 + t185 * t153 + (-m(6) - m(7)) * t135);
t121 = -t164 * mrSges(5,3) + t190 * mrSges(5,1) - t192 * t169 + t200 * (-t200 * mrSges(5,2) + t191 * mrSges(5,3)) - t140 * mrSges(6,2) + t139 * mrSges(6,1) - t173 * t154 + t172 * t153 - t206 * t127 - t211 * t126 + (-m(5) - m(6)) * t144;
t178 = Ifges(4,6) * qJD(2) + (Ifges(4,5) * t209 - Ifges(4,3) * t214) * qJD(1);
t181 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t209 + Ifges(3,2) * t214) * qJD(1);
t197 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t228;
t132 = Ifges(7,5) * t152 + Ifges(7,6) * t151 + Ifges(7,3) * t170;
t134 = Ifges(7,1) * t152 + Ifges(7,4) * t151 + Ifges(7,5) * t170;
t124 = -mrSges(7,1) * t135 + mrSges(7,3) * t129 + Ifges(7,4) * t131 + Ifges(7,2) * t130 + Ifges(7,6) * t138 - t152 * t132 + t170 * t134;
t133 = Ifges(7,4) * t152 + Ifges(7,2) * t151 + Ifges(7,6) * t170;
t125 = mrSges(7,2) * t135 - mrSges(7,3) * t128 + Ifges(7,1) * t131 + Ifges(7,4) * t130 + Ifges(7,5) * t138 + t151 * t132 - t170 * t133;
t146 = Ifges(6,5) * t173 + Ifges(6,6) * t172 + Ifges(6,3) * t185;
t147 = Ifges(6,4) * t173 + Ifges(6,2) * t172 + Ifges(6,6) * t185;
t120 = mrSges(6,2) * t144 + mrSges(6,3) * t135 + Ifges(6,1) * t140 + Ifges(6,4) * t139 + Ifges(6,5) * t161 - t206 * t124 + t211 * t125 + t172 * t146 - t185 * t147;
t148 = Ifges(6,1) * t173 + Ifges(6,4) * t172 + Ifges(6,5) * t185;
t219 = mrSges(7,1) * t128 - mrSges(7,2) * t129 + Ifges(7,5) * t131 + Ifges(7,6) * t130 + Ifges(7,3) * t138 + t152 * t133 - t151 * t134;
t123 = -mrSges(6,1) * t144 + mrSges(6,3) * t136 + Ifges(6,4) * t140 + Ifges(6,2) * t139 + Ifges(6,6) * t161 - t173 * t146 + t185 * t148 - t219;
t156 = Ifges(5,5) * t192 + Ifges(5,6) * t191 + Ifges(5,3) * t200;
t157 = Ifges(5,4) * t192 + Ifges(5,2) * t191 + Ifges(5,6) * t200;
t115 = mrSges(5,2) * t171 + mrSges(5,3) * t144 + Ifges(5,1) * t164 + Ifges(5,4) * t163 + Ifges(5,5) * t190 + t212 * t120 - t207 * t123 + t191 * t156 - t200 * t157;
t158 = Ifges(5,1) * t192 + Ifges(5,4) * t191 + Ifges(5,5) * t200;
t218 = mrSges(6,1) * t135 + mrSges(6,2) * t136 - Ifges(6,5) * t140 - Ifges(6,6) * t139 - Ifges(6,3) * t161 - t211 * t124 - t206 * t125 - t173 * t147 + t172 * t148;
t118 = -mrSges(5,1) * t171 + mrSges(5,3) * t145 + Ifges(5,4) * t164 + Ifges(5,2) * t163 + Ifges(5,6) * t190 - t192 * t156 + t200 * t158 + t218;
t223 = -mrSges(4,1) * t171 + mrSges(4,3) * t166 + Ifges(4,4) * t194 + Ifges(4,2) * qJDD(2) - Ifges(4,6) * t195 + t208 * t115 + t213 * t118;
t182 = Ifges(4,4) * qJD(2) + (Ifges(4,1) * t209 - Ifges(4,5) * t214) * qJD(1);
t229 = t182 + Ifges(3,5) * qJD(2) + (Ifges(3,1) * t209 + Ifges(3,4) * t214) * qJD(1);
t233 = -((t178 - t181) * t209 + t229 * t214) * qJD(1) + mrSges(3,1) * t176 - mrSges(3,2) * t177 + Ifges(3,5) * t194 + Ifges(3,6) * t195 + Ifges(3,3) * qJDD(2) + qJ(3) * (m(4) * t166 + t195 * mrSges(4,2) + qJDD(2) * mrSges(4,3) + qJD(2) * t197 + t213 * t116 - t208 * t121 + (-mrSges(4,1) * t214 - mrSges(4,3) * t209) * t232) + t223;
t179 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t209 + Ifges(3,6) * t214) * qJD(1);
t180 = Ifges(4,2) * qJD(2) + (Ifges(4,4) * t209 - Ifges(4,6) * t214) * qJD(1);
t221 = mrSges(4,2) * t171 - mrSges(4,3) * t162 + Ifges(4,1) * t194 + Ifges(4,4) * qJDD(2) - Ifges(4,5) * t195 + qJD(2) * t178 + t213 * t115 - t208 * t118 + t180 * t227;
t109 = Ifges(3,5) * qJDD(2) + t221 + t179 * t227 - qJ(3) * (m(4) * t162 - t194 * mrSges(4,3) - t195 * mrSges(4,1) - t197 * t228 - (mrSges(4,2) * t227 + qJD(2) * mrSges(4,3)) * t227 + t208 * t116 + t213 * t121) + Ifges(3,1) * t194 + Ifges(3,4) * t195 - mrSges(3,2) * t198 - qJD(2) * t181 - mrSges(3,3) * t176;
t222 = mrSges(5,1) * t144 + mrSges(5,2) * t145 - Ifges(5,5) * t164 - Ifges(5,6) * t163 - Ifges(5,3) * t190 - t207 * t120 - t212 * t123 - t192 * t157 + t191 * t158;
t220 = -mrSges(4,1) * t162 + mrSges(4,2) * t166 + t222;
t112 = (-t179 - t180) * t228 + t220 + (Ifges(3,2) + Ifges(4,3)) * t195 + (Ifges(3,4) - Ifges(4,5)) * t194 + (Ifges(3,6) - Ifges(4,6)) * qJDD(2) + t229 * qJD(2) + mrSges(3,1) * t198 + mrSges(3,3) * t177;
t224 = mrSges(2,1) * t198 - mrSges(2,2) * t199 + Ifges(2,3) * qJDD(1) + t209 * t109 + t214 * t112;
t107 = mrSges(2,1) * g(3) + mrSges(2,3) * t199 + t216 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - t233;
t106 = -mrSges(2,2) * g(3) - mrSges(2,3) * t198 + Ifges(2,5) * qJDD(1) - t216 * Ifges(2,6) + t214 * t109 - t209 * t112;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t215 * t106 - t210 * t107, t106, t109, t221, t115, t120, t125; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t210 * t106 + t215 * t107, t107, t112, (-t209 * t178 - t214 * t182) * qJD(1) + t223, t118, t123, t124; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t224, t224, t233, Ifges(4,5) * t194 + Ifges(4,6) * qJDD(2) - Ifges(4,3) * t195 - qJD(2) * t182 + t180 * t228 - t220, -t222, -t218, t219;];
m_new  = t1;
