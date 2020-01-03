% Calculate vector of cutting torques with Newton-Euler for
% S4RRPP5
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RRPP5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP5_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP5_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_invdynm_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP5_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP5_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP5_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:00:11
% EndTime: 2019-12-31 17:00:14
% DurationCPUTime: 1.03s
% Computational Cost: add. (4660->241), mult. (9538->281), div. (0->0), fcn. (3784->4), ass. (0->86)
t196 = sin(qJ(1));
t198 = cos(qJ(1));
t180 = -t198 * g(1) - t196 * g(2);
t200 = qJD(1) ^ 2;
t144 = -t200 * pkin(1) + qJDD(1) * pkin(5) + t180;
t195 = sin(qJ(2));
t197 = cos(qJ(2));
t123 = -t197 * g(3) - t195 * t144;
t124 = -t195 * g(3) + t197 * t144;
t135 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t195 + Ifges(3,2) * t197) * qJD(1);
t136 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t195 + Ifges(3,4) * t197) * qJD(1);
t137 = Ifges(5,5) * qJD(2) + (-Ifges(5,6) * t197 + Ifges(5,3) * t195) * qJD(1);
t221 = qJD(1) * qJD(2);
t218 = t197 * t221;
t166 = t195 * qJDD(1) + t218;
t219 = t195 * t221;
t167 = t197 * qJDD(1) - t219;
t222 = qJD(1) * t197;
t176 = -mrSges(4,1) * t222 - qJD(2) * mrSges(4,3);
t165 = (-mrSges(5,2) * t195 - mrSges(5,3) * t197) * qJD(1);
t162 = (-pkin(2) * t197 - qJ(3) * t195) * qJD(1);
t199 = qJD(2) ^ 2;
t223 = qJD(1) * t195;
t122 = -qJDD(2) * pkin(2) - t199 * qJ(3) + t162 * t223 + qJDD(3) - t123;
t233 = -2 * qJD(4);
t117 = qJD(2) * t233 + (-t195 * t197 * t200 - qJDD(2)) * qJ(4) + (t166 - t218) * pkin(3) + t122;
t177 = mrSges(5,1) * t222 + qJD(2) * mrSges(5,2);
t215 = -m(5) * t117 + qJDD(2) * mrSges(5,3) + qJD(2) * t177;
t230 = t166 * mrSges(5,1);
t109 = t165 * t223 - t215 + t230;
t206 = -t199 * pkin(2) + qJDD(2) * qJ(3) + t162 * t222 + t124;
t234 = -2 * qJD(3);
t120 = qJD(2) * t234 - t206;
t140 = Ifges(4,4) * qJD(2) + (-Ifges(4,2) * t195 - Ifges(4,6) * t197) * qJD(1);
t174 = pkin(3) * t223 - qJD(2) * qJ(4);
t194 = t197 ^ 2;
t118 = -t194 * t200 * qJ(4) + t167 * pkin(3) + qJDD(4) + ((2 * qJD(3)) + t174) * qJD(2) + t206;
t213 = mrSges(5,2) * t118 - mrSges(5,3) * t117 + Ifges(5,1) * qJDD(2) - Ifges(5,4) * t167 + Ifges(5,5) * t166;
t204 = mrSges(4,2) * t122 - mrSges(4,3) * t120 + Ifges(4,1) * qJDD(2) - Ifges(4,4) * t166 - Ifges(4,5) * t167 - qJ(4) * t109 + t140 * t222 + t213;
t163 = (mrSges(4,2) * t197 - mrSges(4,3) * t195) * qJD(1);
t178 = mrSges(4,1) * t223 + qJD(2) * mrSges(4,2);
t175 = mrSges(5,1) * t223 - qJD(2) * mrSges(5,3);
t216 = m(5) * t118 + qJDD(2) * mrSges(5,2) + qJD(2) * t175 + t165 * t222;
t207 = -m(4) * t120 + qJDD(2) * mrSges(4,3) + qJD(2) * t178 + t163 * t222 + t216;
t211 = -m(4) * t122 - t166 * mrSges(4,1) + t215;
t224 = -t163 - t165;
t138 = Ifges(4,5) * qJD(2) + (-Ifges(4,6) * t195 - Ifges(4,3) * t197) * qJD(1);
t139 = Ifges(5,4) * qJD(2) + (-Ifges(5,2) * t197 + Ifges(5,6) * t195) * qJD(1);
t225 = -t138 - t139;
t237 = ((t135 + t225) * t195 - (t136 + t137) * t197) * qJD(1) + mrSges(3,1) * t123 - mrSges(3,2) * t124 + Ifges(3,5) * t166 + Ifges(3,6) * t167 + Ifges(3,3) * qJDD(2) + pkin(2) * (-qJDD(2) * mrSges(4,2) - qJD(2) * t176 + t224 * t223 + t211 - t230) + qJ(3) * ((mrSges(4,1) + mrSges(5,1)) * t167 + t207) + t204;
t236 = pkin(2) * t219 + t223 * t234;
t232 = -mrSges(5,1) - mrSges(3,3);
t231 = Ifges(3,4) + Ifges(4,6);
t229 = t197 * t137;
t142 = Ifges(4,1) * qJD(2) + (-Ifges(4,4) * t195 - Ifges(4,5) * t197) * qJD(1);
t228 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t195 + Ifges(3,6) * t197) * qJD(1) + t142;
t226 = -t137 + t140;
t179 = t196 * g(1) - t198 * g(2);
t164 = (-mrSges(3,1) * t197 + mrSges(3,2) * t195) * qJD(1);
t172 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t223;
t103 = t164 * t222 + m(3) * t124 - qJDD(2) * mrSges(3,2) - qJD(2) * t172 + (mrSges(4,1) - t232) * t167 + t207;
t173 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t222;
t104 = m(3) * t123 + t232 * t166 + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + (t173 - t176) * qJD(2) + (-t164 + t224) * t223 + t211;
t217 = t197 * t103 - t195 * t104;
t214 = -qJDD(1) * pkin(1) - t179;
t113 = -t166 * qJ(3) + (-pkin(3) * t194 - pkin(5)) * t200 + (-pkin(2) - qJ(4)) * t167 + (-t174 * t195 + (-qJ(3) * qJD(2) + t233) * t197) * qJD(1) + t214 + t236;
t108 = m(5) * t113 - t166 * mrSges(5,2) - t167 * mrSges(5,3) - t175 * t223 - t177 * t222;
t141 = Ifges(5,1) * qJD(2) + (-Ifges(5,4) * t197 + Ifges(5,5) * t195) * qJD(1);
t212 = mrSges(5,1) * t118 - mrSges(5,3) * t113 - Ifges(5,4) * qJDD(2) + Ifges(5,2) * t167 - Ifges(5,6) * t166 - t141 * t223;
t143 = -t200 * pkin(5) + t214;
t119 = -t167 * pkin(2) + (-t166 - t218) * qJ(3) + t143 + t236;
t208 = -m(4) * t119 - t167 * mrSges(4,2) + t178 * t223 - t108;
t202 = -m(3) * t143 + t173 * t222 + t167 * mrSges(3,1) + (-mrSges(3,2) + mrSges(4,3)) * t166 + (-t172 * t195 - t176 * t197) * qJD(1) + t208;
t105 = -t166 * mrSges(4,3) + t176 * t222 - t208;
t203 = mrSges(4,1) * t120 - mrSges(4,2) * t119 + pkin(3) * (-t167 * mrSges(5,1) - t216) + qJ(4) * t108 - t212;
t94 = -t203 - t228 * t223 - pkin(2) * t105 + mrSges(3,3) * t124 - mrSges(3,1) * t143 + (t136 - t226) * qJD(2) + (Ifges(3,6) - Ifges(4,5)) * qJDD(2) + t231 * t166 + (Ifges(3,2) + Ifges(4,3)) * t167;
t209 = mrSges(5,1) * t117 - mrSges(5,2) * t113 + Ifges(5,5) * qJDD(2) - Ifges(5,6) * t167 + Ifges(5,3) * t166 + qJD(2) * t139 + t141 * t222;
t205 = mrSges(4,1) * t122 - mrSges(4,3) * t119 + pkin(3) * t109 + t209;
t96 = t205 + t228 * t222 - qJ(3) * t105 - mrSges(3,3) * t123 + mrSges(3,2) * t143 + (-t135 + t138) * qJD(2) + (Ifges(3,5) - Ifges(4,4)) * qJDD(2) + (Ifges(3,1) + Ifges(4,2)) * t166 + t231 * t167;
t210 = mrSges(2,1) * t179 - mrSges(2,2) * t180 + Ifges(2,3) * qJDD(1) + pkin(1) * t202 + pkin(5) * t217 + t195 * t96 + t197 * t94;
t100 = m(2) * t179 + qJDD(1) * mrSges(2,1) - t200 * mrSges(2,2) + t202;
t99 = t195 * t103 + t197 * t104;
t97 = m(2) * t180 - t200 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t217;
t92 = mrSges(2,1) * g(3) + mrSges(2,3) * t180 + t200 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t99 - t237;
t91 = -mrSges(2,2) * g(3) - mrSges(2,3) * t179 + Ifges(2,5) * qJDD(1) - t200 * Ifges(2,6) - pkin(5) * t99 - t195 * t94 + t197 * t96;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t198 * t91 - t196 * t92 - pkin(4) * (t198 * t100 + t196 * t97), t91, t96, t204 + (t225 * t195 - t229) * qJD(1), (-t195 * t139 - t229) * qJD(1) + t213; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t196 * t91 + t198 * t92 + pkin(4) * (-t196 * t100 + t198 * t97), t92, t94, Ifges(4,4) * qJDD(2) - Ifges(4,2) * t166 - Ifges(4,6) * t167 - qJD(2) * t138 - t142 * t222 - t205, -qJD(2) * t137 - t212; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t210, t210, t237, Ifges(4,5) * qJDD(2) - Ifges(4,6) * t166 - Ifges(4,3) * t167 + t226 * qJD(2) + t142 * t223 + t203, t209;];
m_new = t1;
