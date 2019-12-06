% Calculate vector of cutting torques with Newton-Euler for
% S5PRRPP3
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRRPP3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP3_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP3_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP3_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP3_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP3_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:12:05
% EndTime: 2019-12-05 16:12:12
% DurationCPUTime: 2.68s
% Computational Cost: add. (27299->257), mult. (55664->310), div. (0->0), fcn. (32659->8), ass. (0->101)
t208 = sin(pkin(8));
t210 = sin(qJ(3));
t232 = qJD(2) * t210;
t239 = cos(pkin(8));
t183 = t208 * qJD(3) + t239 * t232;
t209 = sin(pkin(7));
t240 = cos(pkin(7));
t198 = -t240 * g(1) - t209 * g(2);
t207 = -g(3) + qJDD(1);
t211 = sin(qJ(2));
t213 = cos(qJ(2));
t174 = -t211 * t198 + t213 * t207;
t215 = qJD(2) ^ 2;
t160 = -qJDD(2) * pkin(2) - t215 * pkin(6) - t174;
t212 = cos(qJ(3));
t230 = qJD(2) * qJD(3);
t228 = t212 * t230;
t194 = t210 * qJDD(2) + t228;
t229 = t210 * t230;
t195 = t212 * qJDD(2) - t229;
t142 = (-t194 - t228) * qJ(4) + (-t195 + t229) * pkin(3) + t160;
t175 = t213 * t198 + t211 * t207;
t161 = -t215 * pkin(2) + qJDD(2) * pkin(6) + t175;
t197 = t209 * g(1) - t240 * g(2);
t146 = t212 * t161 - t210 * t197;
t192 = (-pkin(3) * t212 - qJ(4) * t210) * qJD(2);
t214 = qJD(3) ^ 2;
t231 = qJD(2) * t212;
t143 = -t214 * pkin(3) + qJDD(3) * qJ(4) + t192 * t231 + t146;
t222 = t239 * t142 - t208 * t143;
t243 = -2 * qJD(4);
t137 = t183 * t243 + t222;
t182 = -t239 * qJD(3) + t208 * t232;
t138 = t208 * t142 + t239 * t143 + t182 * t243;
t152 = Ifges(5,1) * t183 - Ifges(5,4) * t182 - Ifges(5,5) * t231;
t156 = t182 * mrSges(6,1) - t183 * mrSges(6,3);
t168 = -t182 * mrSges(6,2) - mrSges(6,3) * t231;
t171 = mrSges(6,1) * t231 + t183 * mrSges(6,2);
t172 = -t239 * qJDD(3) + t208 * t194;
t173 = t208 * qJDD(3) + t239 * t194;
t155 = t182 * pkin(4) - t183 * qJ(5);
t238 = t212 ^ 2 * t215;
t242 = -2 * qJD(5);
t133 = -pkin(4) * t238 - t195 * qJ(5) - t182 * t155 + t231 * t242 + t138;
t135 = -qJ(5) * t238 + t195 * pkin(4) + qJDD(5) + ((2 * qJD(4)) + t155) * t183 - t222;
t151 = Ifges(6,1) * t183 - Ifges(6,4) * t231 + Ifges(6,5) * t182;
t221 = mrSges(6,1) * t135 - mrSges(6,3) * t133 - Ifges(6,4) * t173 + Ifges(6,2) * t195 - Ifges(6,6) * t172 - t182 * t151;
t227 = -m(6) * t135 - t195 * mrSges(6,1);
t147 = Ifges(6,5) * t183 - Ifges(6,6) * t231 + Ifges(6,3) * t182;
t234 = -Ifges(5,4) * t183 + Ifges(5,2) * t182 + Ifges(5,6) * t231 + t147;
t236 = m(6) * t133 - t195 * mrSges(6,3);
t245 = t234 * t183 - mrSges(5,1) * t137 + mrSges(5,2) * t138 - Ifges(5,5) * t173 + Ifges(5,6) * t172 - pkin(4) * (-t173 * mrSges(6,2) - t183 * t156 - t168 * t231 + t227) - qJ(5) * (-t172 * mrSges(6,2) - t182 * t156 - t171 * t231 + t236) - t182 * t152 + t221;
t158 = t210 * t161;
t219 = -qJDD(3) * pkin(3) - t214 * qJ(4) + t192 * t232 + qJDD(4) + t158;
t136 = t172 * pkin(4) - t173 * qJ(5) + t183 * t242 + (t197 + (-pkin(4) * t183 - qJ(5) * t182) * qJD(2)) * t212 + t219;
t130 = m(6) * t136 + t172 * mrSges(6,1) - t173 * mrSges(6,3) + t182 * t168 - t183 * t171;
t237 = t212 * t197;
t141 = t219 + t237;
t225 = -mrSges(6,1) * t136 + mrSges(6,2) * t133;
t149 = Ifges(6,4) * t183 - Ifges(6,2) * t231 + Ifges(6,6) * t182;
t235 = -Ifges(5,5) * t183 + Ifges(5,6) * t182 + Ifges(5,3) * t231 - t149;
t119 = -mrSges(5,1) * t141 + mrSges(5,3) * t138 - pkin(4) * t130 + (-Ifges(5,6) + Ifges(6,6)) * t195 + t235 * t183 + (Ifges(5,4) - Ifges(6,5)) * t173 + (-Ifges(5,2) - Ifges(6,3)) * t172 + (-t151 - t152) * t231 + t225;
t223 = mrSges(6,2) * t135 - mrSges(6,3) * t136 + Ifges(6,1) * t173 - Ifges(6,4) * t195 + Ifges(6,5) * t172;
t120 = mrSges(5,2) * t141 - mrSges(5,3) * t137 + Ifges(5,1) * t173 - Ifges(5,4) * t172 - Ifges(5,5) * t195 - qJ(5) * t130 + t235 * t182 - t234 * t231 + t223;
t170 = -mrSges(5,1) * t231 - t183 * mrSges(5,3);
t233 = -t182 * mrSges(5,1) - t183 * mrSges(5,2) - t156;
t241 = -mrSges(5,3) - mrSges(6,2);
t124 = m(5) * t138 + t195 * mrSges(5,2) + t233 * t182 + t241 * t172 + (t170 - t171) * t231 + t236;
t169 = mrSges(5,2) * t231 - t182 * mrSges(5,3);
t125 = m(5) * t137 - t195 * mrSges(5,1) + t233 * t183 + t241 * t173 + (-t168 - t169) * t231 + t227;
t122 = t239 * t124 - t208 * t125;
t127 = -m(5) * t141 - t172 * mrSges(5,1) - t173 * mrSges(5,2) - t182 * t169 - t183 * t170 - t130;
t145 = -t158 - t237;
t179 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t210 + Ifges(4,2) * t212) * qJD(2);
t180 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t210 + Ifges(4,4) * t212) * qJD(2);
t244 = mrSges(4,1) * t145 - mrSges(4,2) * t146 + Ifges(4,5) * t194 + Ifges(4,6) * t195 + Ifges(4,3) * qJDD(3) + pkin(3) * t127 + qJ(4) * t122 + (t179 * t210 - t180 * t212) * qJD(2) + t239 * t119 + t208 * t120;
t193 = (-mrSges(4,1) * t212 + mrSges(4,2) * t210) * qJD(2);
t199 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t232;
t117 = m(4) * t146 - qJDD(3) * mrSges(4,2) + t195 * mrSges(4,3) - qJD(3) * t199 + t193 * t231 + t122;
t200 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t231;
t126 = m(4) * t145 + qJDD(3) * mrSges(4,1) - t194 * mrSges(4,3) + qJD(3) * t200 - t193 * t232 + t127;
t114 = t212 * t117 - t210 * t126;
t110 = m(3) * t175 - t215 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t114;
t121 = t208 * t124 + t239 * t125;
t118 = -m(4) * t160 + t195 * mrSges(4,1) - t194 * mrSges(4,2) - t199 * t232 + t200 * t231 - t121;
t115 = m(3) * t174 + qJDD(2) * mrSges(3,1) - t215 * mrSges(3,2) + t118;
t226 = t213 * t110 - t211 * t115;
t113 = t210 * t117 + t212 * t126;
t178 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t210 + Ifges(4,6) * t212) * qJD(2);
t104 = mrSges(4,2) * t160 - mrSges(4,3) * t145 + Ifges(4,1) * t194 + Ifges(4,4) * t195 + Ifges(4,5) * qJDD(3) - qJ(4) * t121 - qJD(3) * t179 - t208 * t119 + t239 * t120 + t178 * t231;
t108 = Ifges(4,6) * qJDD(3) - t178 * t232 + (Ifges(4,2) + Ifges(5,3)) * t195 + Ifges(4,4) * t194 + qJD(3) * t180 + mrSges(4,3) * t146 - mrSges(4,1) * t160 - pkin(3) * t121 + t245;
t101 = -mrSges(3,2) * t197 - mrSges(3,3) * t174 + Ifges(3,5) * qJDD(2) - t215 * Ifges(3,6) - pkin(6) * t113 + t212 * t104 - t210 * t108;
t103 = mrSges(3,1) * t197 + mrSges(3,3) * t175 + t215 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t113 - t244;
t220 = -mrSges(2,2) * t198 + pkin(5) * t226 + t211 * t101 + t213 * t103 + pkin(1) * (m(3) * t197 - t113) + mrSges(2,1) * t197;
t218 = mrSges(3,1) * t174 - mrSges(3,2) * t175 + Ifges(3,3) * qJDD(2) + pkin(2) * t118 + pkin(6) * t114 + t210 * t104 + t212 * t108;
t111 = (m(2) + m(3)) * t197 - t113;
t107 = t211 * t110 + t213 * t115;
t105 = m(2) * t198 + t226;
t99 = -mrSges(2,1) * t207 + mrSges(2,3) * t198 - pkin(1) * t107 - t218;
t98 = mrSges(2,2) * t207 - mrSges(2,3) * t197 - pkin(5) * t107 + t213 * t101 - t211 * t103;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t240 * t98 - t209 * t99 - qJ(1) * (t209 * t105 + t240 * t111), t98, t101, t104, t120, -t147 * t231 - t182 * t149 + t223; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t209 * t98 + t240 * t99 + qJ(1) * (t240 * t105 - t209 * t111), t99, t103, t108, t119, -t183 * t147 - t221; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t220, t220, t218, t244, -Ifges(5,3) * t195 - t245, Ifges(6,5) * t173 - Ifges(6,6) * t195 + Ifges(6,3) * t172 + t183 * t149 + t151 * t231 - t225;];
m_new = t1;
