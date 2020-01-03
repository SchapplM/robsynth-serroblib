% Calculate vector of cutting torques with Newton-Euler for
% S5RPPRR11
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPPRR11_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR11_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR11_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR11_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR11_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR11_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR11_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:05:34
% EndTime: 2019-12-31 18:05:36
% DurationCPUTime: 1.23s
% Computational Cost: add. (13641->225), mult. (24784->262), div. (0->0), fcn. (11616->6), ass. (0->88)
t207 = sin(qJ(1));
t210 = cos(qJ(1));
t181 = t207 * g(1) - t210 * g(2);
t213 = qJD(1) ^ 2;
t161 = -qJDD(1) * pkin(1) - t213 * qJ(2) + qJDD(2) - t181;
t153 = -qJDD(1) * qJ(3) - (2 * qJD(3) * qJD(1)) + t161;
t182 = -t210 * g(1) - t207 * g(2);
t244 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t182;
t206 = sin(qJ(4));
t243 = t206 * g(3);
t242 = mrSges(3,2) - mrSges(4,3);
t241 = Ifges(3,4) - Ifges(4,5);
t240 = Ifges(2,6) - Ifges(3,5);
t239 = mrSges(4,3) * t213;
t154 = qJDD(3) + (-pkin(1) - qJ(3)) * t213 + t244;
t149 = -qJDD(1) * pkin(6) + t154;
t209 = cos(qJ(4));
t144 = -g(3) * t209 + t206 * t149;
t174 = (mrSges(5,1) * t206 + mrSges(5,2) * t209) * qJD(1);
t236 = qJD(1) * qJD(4);
t232 = t209 * t236;
t176 = -qJDD(1) * t206 - t232;
t237 = qJD(1) * t209;
t180 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t237;
t148 = -t213 * pkin(6) - t153;
t233 = t206 * t236;
t177 = qJDD(1) * t209 - t233;
t131 = (-t177 + t233) * pkin(7) + (-t176 + t232) * pkin(4) + t148;
t175 = (pkin(4) * t206 - pkin(7) * t209) * qJD(1);
t212 = qJD(4) ^ 2;
t238 = qJD(1) * t206;
t133 = -pkin(4) * t212 + qJDD(4) * pkin(7) - t175 * t238 + t144;
t205 = sin(qJ(5));
t208 = cos(qJ(5));
t129 = t131 * t208 - t133 * t205;
t172 = qJD(4) * t208 - t205 * t237;
t142 = qJD(5) * t172 + qJDD(4) * t205 + t177 * t208;
t173 = qJD(4) * t205 + t208 * t237;
t146 = -mrSges(6,1) * t172 + mrSges(6,2) * t173;
t183 = qJD(5) + t238;
t155 = -mrSges(6,2) * t183 + mrSges(6,3) * t172;
t171 = qJDD(5) - t176;
t125 = m(6) * t129 + mrSges(6,1) * t171 - mrSges(6,3) * t142 - t146 * t173 + t155 * t183;
t130 = t131 * t205 + t133 * t208;
t141 = -qJD(5) * t173 + qJDD(4) * t208 - t177 * t205;
t156 = mrSges(6,1) * t183 - mrSges(6,3) * t173;
t126 = m(6) * t130 - mrSges(6,2) * t171 + mrSges(6,3) * t141 + t146 * t172 - t156 * t183;
t230 = -t125 * t205 + t208 * t126;
t112 = m(5) * t144 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t176 - qJD(4) * t180 - t174 * t238 + t230;
t143 = t149 * t209 + t243;
t179 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t238;
t132 = -qJDD(4) * pkin(4) - t212 * pkin(7) - t243 + (qJD(1) * t175 - t149) * t209;
t223 = -m(6) * t132 + t141 * mrSges(6,1) - mrSges(6,2) * t142 + t172 * t155 - t156 * t173;
t121 = m(5) * t143 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t177 + qJD(4) * t179 - t174 * t237 + t223;
t103 = t206 * t112 + t209 * t121;
t114 = t208 * t125 + t205 * t126;
t231 = t209 * t112 - t206 * t121;
t229 = m(4) * t154 + qJDD(1) * mrSges(4,2) + t103;
t226 = -m(5) * t148 + t176 * mrSges(5,1) - t177 * mrSges(5,2) - t179 * t238 - t180 * t237 - t114;
t134 = Ifges(6,5) * t173 + Ifges(6,6) * t172 + Ifges(6,3) * t183;
t136 = Ifges(6,1) * t173 + Ifges(6,4) * t172 + Ifges(6,5) * t183;
t118 = -mrSges(6,1) * t132 + mrSges(6,3) * t130 + Ifges(6,4) * t142 + Ifges(6,2) * t141 + Ifges(6,6) * t171 - t134 * t173 + t136 * t183;
t135 = Ifges(6,4) * t173 + Ifges(6,2) * t172 + Ifges(6,6) * t183;
t119 = mrSges(6,2) * t132 - mrSges(6,3) * t129 + Ifges(6,1) * t142 + Ifges(6,4) * t141 + Ifges(6,5) * t171 + t134 * t172 - t135 * t183;
t162 = (Ifges(5,3) * qJD(4)) + (Ifges(5,5) * t209 - Ifges(5,6) * t206) * qJD(1);
t163 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t209 - Ifges(5,2) * t206) * qJD(1);
t94 = mrSges(5,2) * t148 - mrSges(5,3) * t143 + Ifges(5,1) * t177 + Ifges(5,4) * t176 + Ifges(5,5) * qJDD(4) - pkin(7) * t114 - qJD(4) * t163 - t118 * t205 + t119 * t208 - t162 * t238;
t164 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t209 - Ifges(5,4) * t206) * qJD(1);
t216 = mrSges(6,1) * t129 - mrSges(6,2) * t130 + Ifges(6,5) * t142 + Ifges(6,6) * t141 + Ifges(6,3) * t171 + t135 * t173 - t136 * t172;
t96 = -mrSges(5,1) * t148 + mrSges(5,3) * t144 + Ifges(5,4) * t177 + Ifges(5,2) * t176 + Ifges(5,6) * qJDD(4) - pkin(4) * t114 + qJD(4) * t164 - t162 * t237 - t216;
t225 = mrSges(4,1) * t153 + mrSges(4,2) * g(3) + t213 * Ifges(4,4) + Ifges(4,5) * qJDD(1) + pkin(3) * t226 + pkin(6) * t231 + t206 * t94 + t209 * t96;
t159 = pkin(1) * t213 - t244;
t224 = -m(3) * t159 + t213 * mrSges(3,2) + qJDD(1) * mrSges(3,3) + t229;
t222 = mrSges(4,2) * t154 - mrSges(4,3) * t153 + Ifges(4,1) * qJDD(1) - pkin(6) * t103 - t206 * t96 + t209 * t94;
t221 = mrSges(5,1) * t143 - mrSges(5,2) * t144 + Ifges(5,5) * t177 + Ifges(5,6) * t176 + Ifges(5,3) * qJDD(4) + pkin(4) * t223 + pkin(7) * t230 + t208 * t118 + t205 * t119 + t163 * t237 + t164 * t238;
t107 = m(4) * t153 - t213 * mrSges(4,2) - qJDD(1) * mrSges(4,3) + t226;
t220 = mrSges(3,1) * t161 + pkin(2) * t107 + t225;
t219 = -m(3) * t161 + t213 * mrSges(3,3) - t107;
t218 = mrSges(4,1) * t154 - Ifges(4,4) * qJDD(1) + pkin(3) * t103 + t221;
t217 = mrSges(3,2) * t161 - mrSges(3,3) * t159 + Ifges(3,1) * qJDD(1) - qJ(3) * t107 + t222;
t215 = -mrSges(2,2) * t182 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t219) + mrSges(2,1) * t181 + Ifges(2,3) * qJDD(1) + t217 + qJ(2) * (t224 - t239);
t214 = mrSges(3,1) * t159 + pkin(2) * (-t229 + t239) + qJ(3) * (-m(4) * g(3) + t231) - t218;
t105 = t219 + (mrSges(2,1) - mrSges(3,2)) * qJDD(1) - t213 * mrSges(2,2) + m(2) * t181;
t100 = (-m(3) - m(4)) * g(3) + t231;
t97 = m(2) * t182 - qJDD(1) * mrSges(2,2) + (-mrSges(2,1) - mrSges(4,3)) * t213 + t224;
t91 = -t214 + (Ifges(2,5) - t241) * t213 + t240 * qJDD(1) + (mrSges(2,1) - t242) * g(3) + mrSges(2,3) * t182 - pkin(1) * t100;
t90 = t220 - t240 * t213 + (-mrSges(2,2) + mrSges(3,3)) * g(3) + (-Ifges(3,4) + Ifges(2,5)) * qJDD(1) - mrSges(2,3) * t181 - qJ(2) * t100;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t210 * t90 - t207 * t91 - pkin(5) * (t105 * t210 + t207 * t97), t90, t217, t222, t94, t119; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t207 * t90 + t210 * t91 + pkin(5) * (-t105 * t207 + t210 * t97), t91, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - t213 * Ifges(3,5) - t220, -mrSges(4,3) * g(3) - t213 * Ifges(4,5) - t218, t96, t118; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t215, t215, Ifges(3,5) * qJDD(1) + t242 * g(3) + t241 * t213 + t214, t225, t221, t216;];
m_new = t1;
