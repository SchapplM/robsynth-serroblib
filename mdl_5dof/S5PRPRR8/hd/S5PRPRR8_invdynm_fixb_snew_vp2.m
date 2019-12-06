% Calculate vector of cutting torques with Newton-Euler for
% S5PRPRR8
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRPRR8_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR8_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR8_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR8_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR8_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR8_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR8_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:02:32
% EndTime: 2019-12-05 16:02:39
% DurationCPUTime: 2.91s
% Computational Cost: add. (35211->222), mult. (63585->280), div. (0->0), fcn. (39227->10), ass. (0->103)
t205 = sin(pkin(9));
t207 = cos(pkin(9));
t191 = t205 * g(1) - t207 * g(2);
t192 = -t207 * g(1) - t205 * g(2);
t202 = -g(3) + qJDD(1);
t206 = sin(pkin(5));
t208 = cos(pkin(5));
t211 = sin(qJ(2));
t214 = cos(qJ(2));
t155 = -t211 * t192 + (t191 * t208 + t202 * t206) * t214;
t216 = qJD(2) ^ 2;
t222 = -t216 * qJ(3) + qJDD(3) - t155;
t245 = -pkin(2) - pkin(7);
t151 = t245 * qJDD(2) + t222;
t169 = -t206 * t191 + t208 * t202;
t210 = sin(qJ(4));
t213 = cos(qJ(4));
t147 = t210 * t151 + t213 * t169;
t186 = (mrSges(5,1) * t210 + mrSges(5,2) * t213) * qJD(2);
t233 = qJD(2) * qJD(4);
t230 = t213 * t233;
t188 = -t210 * qJDD(2) - t230;
t235 = qJD(2) * t213;
t194 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t235;
t187 = (pkin(4) * t210 - pkin(8) * t213) * qJD(2);
t215 = qJD(4) ^ 2;
t234 = t210 * qJD(2);
t143 = -t215 * pkin(4) + qJDD(4) * pkin(8) - t187 * t234 + t147;
t237 = t208 * t211;
t238 = t206 * t211;
t156 = t191 * t237 + t214 * t192 + t202 * t238;
t246 = -qJDD(2) * qJ(3) - 0.2e1 * qJD(3) * qJD(2) - t156;
t150 = t245 * t216 - t246;
t231 = t210 * t233;
t189 = t213 * qJDD(2) - t231;
t144 = (-t189 + t231) * pkin(8) + (-t188 + t230) * pkin(4) + t150;
t209 = sin(qJ(5));
t212 = cos(qJ(5));
t140 = -t209 * t143 + t212 * t144;
t184 = t212 * qJD(4) - t209 * t235;
t163 = t184 * qJD(5) + t209 * qJDD(4) + t212 * t189;
t185 = t209 * qJD(4) + t212 * t235;
t164 = -t184 * mrSges(6,1) + t185 * mrSges(6,2);
t196 = qJD(5) + t234;
t166 = -t196 * mrSges(6,2) + t184 * mrSges(6,3);
t181 = qJDD(5) - t188;
t136 = m(6) * t140 + t181 * mrSges(6,1) - t163 * mrSges(6,3) - t185 * t164 + t196 * t166;
t141 = t212 * t143 + t209 * t144;
t162 = -t185 * qJD(5) + t212 * qJDD(4) - t209 * t189;
t167 = t196 * mrSges(6,1) - t185 * mrSges(6,3);
t137 = m(6) * t141 - t181 * mrSges(6,2) + t162 * mrSges(6,3) + t184 * t164 - t196 * t167;
t229 = -t209 * t136 + t212 * t137;
t123 = m(5) * t147 - qJDD(4) * mrSges(5,2) + t188 * mrSges(5,3) - qJD(4) * t194 - t186 * t234 + t229;
t236 = t210 * t169;
t146 = t213 * t151 - t236;
t193 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t234;
t142 = -qJDD(4) * pkin(4) - t215 * pkin(8) + t236 + (qJD(2) * t187 - t151) * t213;
t224 = -m(6) * t142 + t162 * mrSges(6,1) - t163 * mrSges(6,2) + t184 * t166 - t185 * t167;
t132 = m(5) * t146 + qJDD(4) * mrSges(5,1) - t189 * mrSges(5,3) + qJD(4) * t193 - t186 * t235 + t224;
t118 = t213 * t123 - t210 * t132;
t116 = m(4) * t169 + t118;
t126 = t212 * t136 + t209 * t137;
t157 = Ifges(6,5) * t185 + Ifges(6,6) * t184 + Ifges(6,3) * t196;
t159 = Ifges(6,1) * t185 + Ifges(6,4) * t184 + Ifges(6,5) * t196;
t130 = -mrSges(6,1) * t142 + mrSges(6,3) * t141 + Ifges(6,4) * t163 + Ifges(6,2) * t162 + Ifges(6,6) * t181 - t185 * t157 + t196 * t159;
t158 = Ifges(6,4) * t185 + Ifges(6,2) * t184 + Ifges(6,6) * t196;
t131 = mrSges(6,2) * t142 - mrSges(6,3) * t140 + Ifges(6,1) * t163 + Ifges(6,4) * t162 + Ifges(6,5) * t181 + t184 * t157 - t196 * t158;
t173 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t213 - Ifges(5,6) * t210) * qJD(2);
t174 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t213 - Ifges(5,2) * t210) * qJD(2);
t111 = mrSges(5,2) * t150 - mrSges(5,3) * t146 + Ifges(5,1) * t189 + Ifges(5,4) * t188 + Ifges(5,5) * qJDD(4) - pkin(8) * t126 - qJD(4) * t174 - t209 * t130 + t212 * t131 - t173 * t234;
t175 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t213 - Ifges(5,4) * t210) * qJD(2);
t217 = mrSges(6,1) * t140 - mrSges(6,2) * t141 + Ifges(6,5) * t163 + Ifges(6,6) * t162 + Ifges(6,3) * t181 + t185 * t158 - t184 * t159;
t112 = -mrSges(5,1) * t150 + mrSges(5,3) * t147 + Ifges(5,4) * t189 + Ifges(5,2) * t188 + Ifges(5,6) * qJDD(4) - pkin(4) * t126 + qJD(4) * t175 - t173 * t235 - t217;
t124 = -m(5) * t150 + t188 * mrSges(5,1) - t189 * mrSges(5,2) - t193 * t234 - t194 * t235 - t126;
t152 = t216 * pkin(2) + t246;
t221 = -mrSges(4,1) * t152 - pkin(3) * t124 - pkin(7) * t118 - t210 * t111 - t213 * t112;
t241 = Ifges(4,5) - Ifges(3,6);
t242 = -Ifges(4,4) + Ifges(3,5);
t243 = mrSges(3,1) - mrSges(4,2);
t100 = mrSges(3,3) * t156 - pkin(2) * t116 - t241 * qJDD(2) - t243 * t169 + t242 * t216 + t221;
t117 = t210 * t123 + t213 * t132;
t154 = -qJDD(2) * pkin(2) + t222;
t225 = -m(4) * t154 + t216 * mrSges(4,3) - t117;
t114 = m(3) * t155 - t216 * mrSges(3,2) + t243 * qJDD(2) + t225;
t219 = -m(4) * t152 + t216 * mrSges(4,2) + qJDD(2) * mrSges(4,3) - t124;
t121 = m(3) * t156 - t216 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t219;
t109 = -t211 * t114 + t214 * t121;
t247 = pkin(6) * t109 + t100 * t214;
t239 = t114 * t214;
t115 = m(3) * t169 + t116;
t106 = -t206 * t115 + t121 * t237 + t208 * t239;
t220 = mrSges(5,1) * t146 - mrSges(5,2) * t147 + Ifges(5,5) * t189 + Ifges(5,6) * t188 + Ifges(5,3) * qJDD(4) + pkin(4) * t224 + pkin(8) * t229 + t212 * t130 + t209 * t131 + t174 * t235 + t175 * t234;
t218 = mrSges(4,1) * t154 + pkin(3) * t117 + t220;
t102 = t218 - mrSges(3,3) * t155 - qJ(3) * t116 + t242 * qJDD(2) + (mrSges(3,2) - mrSges(4,3)) * t169 + t241 * t216;
t223 = mrSges(4,2) * t154 - mrSges(4,3) * t152 + Ifges(4,1) * qJDD(2) - pkin(7) * t117 + t213 * t111 - t210 * t112;
t98 = Ifges(3,3) * qJDD(2) + mrSges(3,1) * t155 - mrSges(3,2) * t156 + pkin(2) * (-qJDD(2) * mrSges(4,2) + t225) + qJ(3) * t219 + t223;
t226 = mrSges(2,1) * t191 - mrSges(2,2) * t192 + pkin(1) * t106 + t102 * t238 + t247 * t206 + t208 * t98;
t107 = m(2) * t192 + t109;
t105 = t208 * t115 + (t121 * t211 + t239) * t206;
t103 = m(2) * t191 + t106;
t96 = mrSges(2,2) * t202 - mrSges(2,3) * t191 - t211 * t100 + t214 * t102 + (-t105 * t206 - t106 * t208) * pkin(6);
t95 = -mrSges(2,1) * t202 + mrSges(2,3) * t192 - pkin(1) * t105 - t206 * t98 + (t102 * t211 + t247) * t208;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t207 * t96 - t205 * t95 - qJ(1) * (t207 * t103 + t205 * t107), t96, t102, t223, t111, t131; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t205 * t96 + t207 * t95 + qJ(1) * (-t205 * t103 + t207 * t107), t95, t100, mrSges(4,3) * t169 + Ifges(4,4) * qJDD(2) - t216 * Ifges(4,5) - t218, t112, t130; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t226, t226, t98, -mrSges(4,2) * t169 + t216 * Ifges(4,4) + Ifges(4,5) * qJDD(2) - t221, t220, t217;];
m_new = t1;
