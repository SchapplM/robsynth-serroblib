% Calculate vector of cutting torques with Newton-Euler for
% S4RRRP4
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RRRP4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP4_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP4_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP4_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP4_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP4_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:12
% EndTime: 2019-12-31 17:15:16
% DurationCPUTime: 1.73s
% Computational Cost: add. (15091->232), mult. (31856->288), div. (0->0), fcn. (18978->6), ass. (0->86)
t198 = sin(qJ(3));
t199 = sin(qJ(2));
t201 = cos(qJ(3));
t202 = cos(qJ(2));
t172 = (-t198 * t199 + t201 * t202) * qJD(1);
t221 = qJD(1) * qJD(2);
t180 = t199 * qJDD(1) + t202 * t221;
t181 = t202 * qJDD(1) - t199 * t221;
t145 = t172 * qJD(3) + t201 * t180 + t198 * t181;
t173 = (t198 * t202 + t199 * t201) * qJD(1);
t157 = -t172 * mrSges(5,1) + t173 * mrSges(5,2);
t200 = sin(qJ(1));
t203 = cos(qJ(1));
t187 = -t203 * g(1) - t200 * g(2);
t204 = qJD(1) ^ 2;
t175 = -t204 * pkin(1) + qJDD(1) * pkin(5) + t187;
t225 = t199 * t175;
t226 = pkin(2) * t204;
t129 = qJDD(2) * pkin(2) - t180 * pkin(6) - t225 + (pkin(6) * t221 + t199 * t226 - g(3)) * t202;
t161 = -t199 * g(3) + t202 * t175;
t223 = qJD(1) * t199;
t185 = qJD(2) * pkin(2) - pkin(6) * t223;
t197 = t202 ^ 2;
t130 = t181 * pkin(6) - qJD(2) * t185 - t197 * t226 + t161;
t121 = t201 * t129 - t198 * t130;
t194 = qJDD(2) + qJDD(3);
t195 = qJD(2) + qJD(3);
t113 = -0.2e1 * qJD(4) * t173 + (t172 * t195 - t145) * qJ(4) + (t172 * t173 + t194) * pkin(3) + t121;
t162 = -t195 * mrSges(5,2) + t172 * mrSges(5,3);
t220 = m(5) * t113 + t194 * mrSges(5,1) + t195 * t162;
t110 = -t145 * mrSges(5,3) - t173 * t157 + t220;
t122 = t198 * t129 + t201 * t130;
t144 = -t173 * qJD(3) - t198 * t180 + t201 * t181;
t151 = Ifges(4,4) * t173 + Ifges(4,2) * t172 + Ifges(4,6) * t195;
t152 = Ifges(5,1) * t173 + Ifges(5,4) * t172 + Ifges(5,5) * t195;
t153 = Ifges(4,1) * t173 + Ifges(4,4) * t172 + Ifges(4,5) * t195;
t164 = t195 * pkin(3) - t173 * qJ(4);
t168 = t172 ^ 2;
t116 = -t168 * pkin(3) + t144 * qJ(4) + 0.2e1 * qJD(4) * t172 - t195 * t164 + t122;
t150 = Ifges(5,4) * t173 + Ifges(5,2) * t172 + Ifges(5,6) * t195;
t211 = -mrSges(5,1) * t113 + mrSges(5,2) * t116 - Ifges(5,5) * t145 - Ifges(5,6) * t144 - Ifges(5,3) * t194 - t173 * t150;
t229 = mrSges(4,1) * t121 - mrSges(4,2) * t122 + Ifges(4,5) * t145 + Ifges(4,6) * t144 + Ifges(4,3) * t194 + pkin(3) * t110 + t173 * t151 - t211 + (-t153 - t152) * t172;
t158 = -t172 * mrSges(4,1) + t173 * mrSges(4,2);
t163 = -t195 * mrSges(4,2) + t172 * mrSges(4,3);
t105 = m(4) * t121 + t194 * mrSges(4,1) + t195 * t163 + (-t157 - t158) * t173 + (-mrSges(4,3) - mrSges(5,3)) * t145 + t220;
t165 = t195 * mrSges(5,1) - t173 * mrSges(5,3);
t166 = t195 * mrSges(4,1) - t173 * mrSges(4,3);
t219 = m(5) * t116 + t144 * mrSges(5,3) + t172 * t157;
t108 = m(4) * t122 + t144 * mrSges(4,3) + t172 * t158 + (-t165 - t166) * t195 + (-mrSges(4,2) - mrSges(5,2)) * t194 + t219;
t101 = t201 * t105 + t198 * t108;
t160 = -t202 * g(3) - t225;
t170 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t199 + Ifges(3,2) * t202) * qJD(1);
t171 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t199 + Ifges(3,4) * t202) * qJD(1);
t228 = mrSges(3,1) * t160 - mrSges(3,2) * t161 + Ifges(3,5) * t180 + Ifges(3,6) * t181 + Ifges(3,3) * qJDD(2) + pkin(2) * t101 + (t199 * t170 - t202 * t171) * qJD(1) + t229;
t222 = qJD(1) * t202;
t186 = t200 * g(1) - t203 * g(2);
t179 = (-mrSges(3,1) * t202 + mrSges(3,2) * t199) * qJD(1);
t183 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t223;
t216 = -t198 * t105 + t201 * t108;
t100 = m(3) * t161 - qJDD(2) * mrSges(3,2) + t181 * mrSges(3,3) - qJD(2) * t183 + t179 * t222 + t216;
t184 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t222;
t99 = m(3) * t160 + qJDD(2) * mrSges(3,1) - t180 * mrSges(3,3) + qJD(2) * t184 - t179 * t223 + t101;
t218 = t202 * t100 - t199 * t99;
t213 = -qJDD(1) * pkin(1) - t186;
t146 = -t181 * pkin(2) + t185 * t223 + (-pkin(6) * t197 - pkin(5)) * t204 + t213;
t119 = -t144 * pkin(3) - t168 * qJ(4) + t173 * t164 + qJDD(4) + t146;
t215 = m(5) * t119 - t144 * mrSges(5,1) + t145 * mrSges(5,2) - t172 * t162 + t173 * t165;
t212 = -mrSges(5,1) * t119 + mrSges(5,3) * t116 + Ifges(5,4) * t145 + Ifges(5,2) * t144 + Ifges(5,6) * t194 + t195 * t152;
t148 = Ifges(5,5) * t173 + Ifges(5,6) * t172 + Ifges(5,3) * t195;
t210 = mrSges(5,2) * t119 - mrSges(5,3) * t113 + Ifges(5,1) * t145 + Ifges(5,4) * t144 + Ifges(5,5) * t194 + t172 * t148;
t174 = -t204 * pkin(5) + t213;
t208 = m(4) * t146 - t144 * mrSges(4,1) + t145 * mrSges(4,2) - t172 * t163 + t173 * t166 + t215;
t206 = -m(3) * t174 + t181 * mrSges(3,1) - t180 * mrSges(3,2) - t183 * t223 + t184 * t222 - t208;
t169 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t199 + Ifges(3,6) * t202) * qJD(1);
t149 = Ifges(4,5) * t173 + Ifges(4,6) * t172 + Ifges(4,3) * t195;
t96 = Ifges(4,4) * t145 + Ifges(4,2) * t144 + Ifges(4,6) * t194 + t195 * t153 - mrSges(4,1) * t146 + mrSges(4,3) * t122 - pkin(3) * t215 + qJ(4) * (-t194 * mrSges(5,2) - t195 * t165 + t219) + (-t149 - t148) * t173 + t212;
t97 = mrSges(4,2) * t146 - mrSges(4,3) * t121 + Ifges(4,1) * t145 + Ifges(4,4) * t144 + Ifges(4,5) * t194 - qJ(4) * t110 + t172 * t149 + (-t150 - t151) * t195 + t210;
t90 = -mrSges(3,1) * t174 + mrSges(3,3) * t161 + Ifges(3,4) * t180 + Ifges(3,2) * t181 + Ifges(3,6) * qJDD(2) - pkin(2) * t208 + pkin(6) * t216 + qJD(2) * t171 - t169 * t223 + t198 * t97 + t201 * t96;
t92 = mrSges(3,2) * t174 - mrSges(3,3) * t160 + Ifges(3,1) * t180 + Ifges(3,4) * t181 + Ifges(3,5) * qJDD(2) - pkin(6) * t101 - qJD(2) * t170 + t169 * t222 - t198 * t96 + t201 * t97;
t209 = mrSges(2,1) * t186 - mrSges(2,2) * t187 + Ifges(2,3) * qJDD(1) + pkin(1) * t206 + pkin(5) * t218 + t199 * t92 + t202 * t90;
t102 = m(2) * t186 + qJDD(1) * mrSges(2,1) - t204 * mrSges(2,2) + t206;
t95 = t199 * t100 + t202 * t99;
t93 = m(2) * t187 - t204 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t218;
t88 = mrSges(2,1) * g(3) + mrSges(2,3) * t187 + t204 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t95 - t228;
t87 = -mrSges(2,2) * g(3) - mrSges(2,3) * t186 + Ifges(2,5) * qJDD(1) - t204 * Ifges(2,6) - pkin(5) * t95 - t199 * t90 + t202 * t92;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t203 * t87 - t200 * t88 - pkin(4) * (t203 * t102 + t200 * t93), t87, t92, t97, -t195 * t150 + t210; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t200 * t87 + t203 * t88 + pkin(4) * (-t200 * t102 + t203 * t93), t88, t90, t96, -t173 * t148 + t212; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t209, t209, t228, t229, -t172 * t152 - t211;];
m_new = t1;
