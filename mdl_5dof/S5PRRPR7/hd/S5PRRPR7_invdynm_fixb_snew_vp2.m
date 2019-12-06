% Calculate vector of cutting torques with Newton-Euler for
% S5PRRPR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRRPR7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR7_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR7_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR7_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR7_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR7_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:35:18
% EndTime: 2019-12-05 16:35:34
% DurationCPUTime: 6.50s
% Computational Cost: add. (101577->267), mult. (207362->348), div. (0->0), fcn. (141119->12), ass. (0->120)
t220 = sin(pkin(9));
t223 = cos(pkin(9));
t211 = t220 * g(1) - t223 * g(2);
t212 = -t223 * g(1) - t220 * g(2);
t218 = -g(3) + qJDD(1);
t230 = cos(qJ(2));
t224 = cos(pkin(5));
t227 = sin(qJ(2));
t249 = t224 * t227;
t221 = sin(pkin(5));
t250 = t221 * t227;
t176 = t211 * t249 + t230 * t212 + t218 * t250;
t232 = qJD(2) ^ 2;
t172 = -t232 * pkin(2) + qJDD(2) * pkin(7) + t176;
t192 = -t221 * t211 + t224 * t218;
t226 = sin(qJ(3));
t229 = cos(qJ(3));
t160 = t229 * t172 + t226 * t192;
t206 = (-pkin(3) * t229 - qJ(4) * t226) * qJD(2);
t231 = qJD(3) ^ 2;
t246 = qJD(2) * t229;
t156 = -t231 * pkin(3) + qJDD(3) * qJ(4) + t206 * t246 + t160;
t175 = -t227 * t212 + (t211 * t224 + t218 * t221) * t230;
t171 = -qJDD(2) * pkin(2) - t232 * pkin(7) - t175;
t245 = qJD(2) * qJD(3);
t243 = t229 * t245;
t208 = t226 * qJDD(2) + t243;
t244 = t226 * t245;
t209 = t229 * qJDD(2) - t244;
t158 = (-t208 - t243) * qJ(4) + (-t209 + t244) * pkin(3) + t171;
t219 = sin(pkin(10));
t222 = cos(pkin(10));
t247 = qJD(2) * t226;
t201 = t222 * qJD(3) - t219 * t247;
t255 = 2 * qJD(4);
t152 = t222 * t156 + t219 * t158 + t201 * t255;
t202 = t219 * qJD(3) + t222 * t247;
t182 = -t201 * pkin(4) - t202 * pkin(8);
t251 = t229 ^ 2 * t232;
t150 = -pkin(4) * t251 - t209 * pkin(8) + t201 * t182 + t152;
t190 = t222 * qJDD(3) - t219 * t208;
t191 = t219 * qJDD(3) + t222 * t208;
t169 = t226 * t172;
t238 = -qJDD(3) * pkin(3) - t231 * qJ(4) + t206 * t247 + qJDD(4) + t169;
t153 = -t190 * pkin(4) - t191 * pkin(8) + (-t192 + (-pkin(4) * t202 + pkin(8) * t201) * qJD(2)) * t229 + t238;
t225 = sin(qJ(5));
t228 = cos(qJ(5));
t147 = -t225 * t150 + t228 * t153;
t183 = -t225 * t202 - t228 * t246;
t167 = t183 * qJD(5) + t228 * t191 - t225 * t209;
t184 = t228 * t202 - t225 * t246;
t168 = -t183 * mrSges(6,1) + t184 * mrSges(6,2);
t199 = qJD(5) - t201;
t173 = -t199 * mrSges(6,2) + t183 * mrSges(6,3);
t187 = qJDD(5) - t190;
t144 = m(6) * t147 + t187 * mrSges(6,1) - t167 * mrSges(6,3) - t184 * t168 + t199 * t173;
t148 = t228 * t150 + t225 * t153;
t166 = -t184 * qJD(5) - t225 * t191 - t228 * t209;
t174 = t199 * mrSges(6,1) - t184 * mrSges(6,3);
t145 = m(6) * t148 - t187 * mrSges(6,2) + t166 * mrSges(6,3) + t183 * t168 - t199 * t174;
t139 = -t225 * t144 + t228 * t145;
t181 = -t201 * mrSges(5,1) + t202 * mrSges(5,2);
t189 = -mrSges(5,1) * t246 - t202 * mrSges(5,3);
t136 = m(5) * t152 + t209 * mrSges(5,2) + t190 * mrSges(5,3) + t201 * t181 + t189 * t246 + t139;
t241 = t219 * t156 - t222 * t158;
t149 = -pkin(8) * t251 + t209 * pkin(4) + (t255 + t182) * t202 + t241;
t146 = -m(6) * t149 + t166 * mrSges(6,1) - t167 * mrSges(6,2) + t183 * t173 - t184 * t174;
t151 = -0.2e1 * qJD(4) * t202 - t241;
t188 = mrSges(5,2) * t246 + t201 * mrSges(5,3);
t142 = m(5) * t151 - t209 * mrSges(5,1) - t191 * mrSges(5,3) - t202 * t181 - t188 * t246 + t146;
t132 = t222 * t136 - t219 * t142;
t207 = (-mrSges(4,1) * t229 + mrSges(4,2) * t226) * qJD(2);
t213 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t247;
t130 = m(4) * t160 - qJDD(3) * mrSges(4,2) + t209 * mrSges(4,3) - qJD(3) * t213 + t207 * t246 + t132;
t138 = t228 * t144 + t225 * t145;
t248 = t229 * t192;
t155 = t238 - t248;
t137 = -m(5) * t155 + t190 * mrSges(5,1) - t191 * mrSges(5,2) + t201 * t188 - t202 * t189 - t138;
t159 = -t169 + t248;
t214 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t246;
t134 = m(4) * t159 + qJDD(3) * mrSges(4,1) - t208 * mrSges(4,3) + qJD(3) * t214 - t207 * t247 + t137;
t123 = t226 * t130 + t229 * t134;
t161 = Ifges(6,5) * t184 + Ifges(6,6) * t183 + Ifges(6,3) * t199;
t163 = Ifges(6,1) * t184 + Ifges(6,4) * t183 + Ifges(6,5) * t199;
t140 = -mrSges(6,1) * t149 + mrSges(6,3) * t148 + Ifges(6,4) * t167 + Ifges(6,2) * t166 + Ifges(6,6) * t187 - t184 * t161 + t199 * t163;
t162 = Ifges(6,4) * t184 + Ifges(6,2) * t183 + Ifges(6,6) * t199;
t141 = mrSges(6,2) * t149 - mrSges(6,3) * t147 + Ifges(6,1) * t167 + Ifges(6,4) * t166 + Ifges(6,5) * t187 + t183 * t161 - t199 * t162;
t177 = Ifges(5,5) * t202 + Ifges(5,6) * t201 - Ifges(5,3) * t246;
t178 = Ifges(5,4) * t202 + Ifges(5,2) * t201 - Ifges(5,6) * t246;
t124 = mrSges(5,2) * t155 - mrSges(5,3) * t151 + Ifges(5,1) * t191 + Ifges(5,4) * t190 - Ifges(5,5) * t209 - pkin(8) * t138 - t225 * t140 + t228 * t141 + t201 * t177 + t178 * t246;
t179 = Ifges(5,1) * t202 + Ifges(5,4) * t201 - Ifges(5,5) * t246;
t235 = mrSges(6,1) * t147 - mrSges(6,2) * t148 + Ifges(6,5) * t167 + Ifges(6,6) * t166 + Ifges(6,3) * t187 + t184 * t162 - t183 * t163;
t125 = -mrSges(5,1) * t155 + mrSges(5,3) * t152 + Ifges(5,4) * t191 + Ifges(5,2) * t190 - Ifges(5,6) * t209 - pkin(4) * t138 - t202 * t177 - t179 * t246 - t235;
t197 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t226 + Ifges(4,2) * t229) * qJD(2);
t198 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t226 + Ifges(4,4) * t229) * qJD(2);
t256 = mrSges(4,1) * t159 - mrSges(4,2) * t160 + Ifges(4,5) * t208 + Ifges(4,6) * t209 + Ifges(4,3) * qJDD(3) + pkin(3) * t137 + qJ(4) * t132 + t219 * t124 + t222 * t125 + (t226 * t197 - t229 * t198) * qJD(2);
t109 = -mrSges(3,1) * t192 + mrSges(3,3) * t176 + t232 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t123 - t256;
t242 = t229 * t130 - t226 * t134;
t121 = m(3) * t176 - t232 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t242;
t131 = t219 * t136 + t222 * t142;
t236 = -m(4) * t171 + t209 * mrSges(4,1) - t208 * mrSges(4,2) - t213 * t247 + t214 * t246 - t131;
t127 = m(3) * t175 + qJDD(2) * mrSges(3,1) - t232 * mrSges(3,2) + t236;
t118 = t230 * t121 - t227 * t127;
t257 = pkin(6) * t118 + t109 * t230;
t252 = t127 * t230;
t122 = m(3) * t192 + t123;
t113 = t121 * t249 - t221 * t122 + t224 * t252;
t196 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t226 + Ifges(4,6) * t229) * qJD(2);
t114 = mrSges(4,2) * t171 - mrSges(4,3) * t159 + Ifges(4,1) * t208 + Ifges(4,4) * t209 + Ifges(4,5) * qJDD(3) - qJ(4) * t131 - qJD(3) * t197 + t222 * t124 - t219 * t125 + t196 * t246;
t233 = -mrSges(5,1) * t151 + mrSges(5,2) * t152 - Ifges(5,5) * t191 - Ifges(5,6) * t190 - pkin(4) * t146 - pkin(8) * t139 - t228 * t140 - t225 * t141 - t202 * t178 + t201 * t179;
t115 = t233 + Ifges(4,6) * qJDD(3) - t196 * t247 + (Ifges(4,2) + Ifges(5,3)) * t209 + Ifges(4,4) * t208 + qJD(3) * t198 - mrSges(4,1) * t171 + mrSges(4,3) * t160 - pkin(3) * t131;
t105 = mrSges(3,1) * t175 - mrSges(3,2) * t176 + Ifges(3,3) * qJDD(2) + pkin(2) * t236 + pkin(7) * t242 + t226 * t114 + t229 * t115;
t107 = mrSges(3,2) * t192 - mrSges(3,3) * t175 + Ifges(3,5) * qJDD(2) - t232 * Ifges(3,6) - pkin(7) * t123 + t229 * t114 - t226 * t115;
t237 = mrSges(2,1) * t211 - mrSges(2,2) * t212 + pkin(1) * t113 + t224 * t105 + t107 * t250 + t257 * t221;
t116 = m(2) * t212 + t118;
t112 = t224 * t122 + (t121 * t227 + t252) * t221;
t110 = m(2) * t211 + t113;
t103 = mrSges(2,2) * t218 - mrSges(2,3) * t211 + t230 * t107 - t227 * t109 + (-t112 * t221 - t113 * t224) * pkin(6);
t102 = -mrSges(2,1) * t218 + mrSges(2,3) * t212 - pkin(1) * t112 - t221 * t105 + (t107 * t227 + t257) * t224;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t223 * t103 - t220 * t102 - qJ(1) * (t223 * t110 + t220 * t116), t103, t107, t114, t124, t141; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t220 * t103 + t223 * t102 + qJ(1) * (-t220 * t110 + t223 * t116), t102, t109, t115, t125, t140; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t237, t237, t105, t256, -Ifges(5,3) * t209 - t233, t235;];
m_new = t1;
