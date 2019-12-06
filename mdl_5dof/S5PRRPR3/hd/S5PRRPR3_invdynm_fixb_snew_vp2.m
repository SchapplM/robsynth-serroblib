% Calculate vector of cutting torques with Newton-Euler for
% S5PRRPR3
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
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRRPR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR3_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR3_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR3_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR3_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR3_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:19:12
% EndTime: 2019-12-05 16:19:22
% DurationCPUTime: 6.93s
% Computational Cost: add. (85118->256), mult. (187770->332), div. (0->0), fcn. (125769->10), ass. (0->106)
t229 = sin(pkin(8));
t231 = cos(pkin(8));
t214 = t229 * g(1) - t231 * g(2);
t215 = -t231 * g(1) - t229 * g(2);
t234 = sin(qJ(2));
t237 = cos(qJ(2));
t192 = t234 * t214 + t237 * t215;
t238 = qJD(2) ^ 2;
t187 = -t238 * pkin(2) + qJDD(2) * pkin(6) + t192;
t227 = -g(3) + qJDD(1);
t233 = sin(qJ(3));
t236 = cos(qJ(3));
t173 = -t233 * t187 + t236 * t227;
t254 = qJD(2) * qJD(3);
t253 = t236 * t254;
t211 = t233 * qJDD(2) + t253;
t167 = (-t211 + t253) * qJ(4) + (t233 * t236 * t238 + qJDD(3)) * pkin(3) + t173;
t174 = t236 * t187 + t233 * t227;
t212 = t236 * qJDD(2) - t233 * t254;
t256 = qJD(2) * t233;
t216 = qJD(3) * pkin(3) - qJ(4) * t256;
t226 = t236 ^ 2;
t168 = -t226 * t238 * pkin(3) + t212 * qJ(4) - qJD(3) * t216 + t174;
t228 = sin(pkin(9));
t230 = cos(pkin(9));
t200 = (t228 * t236 + t230 * t233) * qJD(2);
t147 = -0.2e1 * qJD(4) * t200 + t230 * t167 - t228 * t168;
t189 = t230 * t211 + t228 * t212;
t199 = (-t228 * t233 + t230 * t236) * qJD(2);
t144 = (qJD(3) * t199 - t189) * pkin(7) + (t199 * t200 + qJDD(3)) * pkin(4) + t147;
t148 = 0.2e1 * qJD(4) * t199 + t228 * t167 + t230 * t168;
t188 = -t228 * t211 + t230 * t212;
t195 = qJD(3) * pkin(4) - t200 * pkin(7);
t198 = t199 ^ 2;
t145 = -t198 * pkin(4) + t188 * pkin(7) - qJD(3) * t195 + t148;
t232 = sin(qJ(5));
t235 = cos(qJ(5));
t142 = t235 * t144 - t232 * t145;
t178 = t235 * t199 - t232 * t200;
t157 = t178 * qJD(5) + t232 * t188 + t235 * t189;
t179 = t232 * t199 + t235 * t200;
t163 = -t178 * mrSges(6,1) + t179 * mrSges(6,2);
t223 = qJD(3) + qJD(5);
t171 = -t223 * mrSges(6,2) + t178 * mrSges(6,3);
t222 = qJDD(3) + qJDD(5);
t139 = m(6) * t142 + t222 * mrSges(6,1) - t157 * mrSges(6,3) - t179 * t163 + t223 * t171;
t143 = t232 * t144 + t235 * t145;
t156 = -t179 * qJD(5) + t235 * t188 - t232 * t189;
t172 = t223 * mrSges(6,1) - t179 * mrSges(6,3);
t140 = m(6) * t143 - t222 * mrSges(6,2) + t156 * mrSges(6,3) + t178 * t163 - t223 * t172;
t130 = t235 * t139 + t232 * t140;
t181 = -t199 * mrSges(5,1) + t200 * mrSges(5,2);
t193 = -qJD(3) * mrSges(5,2) + t199 * mrSges(5,3);
t127 = m(5) * t147 + qJDD(3) * mrSges(5,1) - t189 * mrSges(5,3) + qJD(3) * t193 - t200 * t181 + t130;
t194 = qJD(3) * mrSges(5,1) - t200 * mrSges(5,3);
t249 = -t232 * t139 + t235 * t140;
t128 = m(5) * t148 - qJDD(3) * mrSges(5,2) + t188 * mrSges(5,3) - qJD(3) * t194 + t199 * t181 + t249;
t123 = t230 * t127 + t228 * t128;
t202 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t233 + Ifges(4,2) * t236) * qJD(2);
t203 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t233 + Ifges(4,4) * t236) * qJD(2);
t176 = Ifges(5,4) * t200 + Ifges(5,2) * t199 + Ifges(5,6) * qJD(3);
t177 = Ifges(5,1) * t200 + Ifges(5,4) * t199 + Ifges(5,5) * qJD(3);
t159 = Ifges(6,4) * t179 + Ifges(6,2) * t178 + Ifges(6,6) * t223;
t160 = Ifges(6,1) * t179 + Ifges(6,4) * t178 + Ifges(6,5) * t223;
t244 = -mrSges(6,1) * t142 + mrSges(6,2) * t143 - Ifges(6,5) * t157 - Ifges(6,6) * t156 - Ifges(6,3) * t222 - t179 * t159 + t178 * t160;
t241 = -mrSges(5,1) * t147 + mrSges(5,2) * t148 - Ifges(5,5) * t189 - Ifges(5,6) * t188 - Ifges(5,3) * qJDD(3) - pkin(4) * t130 - t200 * t176 + t199 * t177 + t244;
t257 = mrSges(4,1) * t173 - mrSges(4,2) * t174 + Ifges(4,5) * t211 + Ifges(4,6) * t212 + Ifges(4,3) * qJDD(3) + pkin(3) * t123 + (t233 * t202 - t236 * t203) * qJD(2) - t241;
t210 = (-mrSges(4,1) * t236 + mrSges(4,2) * t233) * qJD(2);
t255 = qJD(2) * t236;
t218 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t255;
t121 = m(4) * t173 + qJDD(3) * mrSges(4,1) - t211 * mrSges(4,3) + qJD(3) * t218 - t210 * t256 + t123;
t217 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t256;
t250 = -t228 * t127 + t230 * t128;
t122 = m(4) * t174 - qJDD(3) * mrSges(4,2) + t212 * mrSges(4,3) - qJD(3) * t217 + t210 * t255 + t250;
t251 = -t233 * t121 + t236 * t122;
t113 = m(3) * t192 - t238 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t251;
t191 = t237 * t214 - t234 * t215;
t246 = -qJDD(2) * pkin(2) - t191;
t186 = -t238 * pkin(6) + t246;
t169 = -t212 * pkin(3) + qJDD(4) + t216 * t256 + (-qJ(4) * t226 - pkin(6)) * t238 + t246;
t150 = -t188 * pkin(4) - t198 * pkin(7) + t200 * t195 + t169;
t248 = m(6) * t150 - t156 * mrSges(6,1) + t157 * mrSges(6,2) - t178 * t171 + t179 * t172;
t242 = m(5) * t169 - t188 * mrSges(5,1) + t189 * mrSges(5,2) - t199 * t193 + t200 * t194 + t248;
t240 = -m(4) * t186 + t212 * mrSges(4,1) - t211 * mrSges(4,2) - t217 * t256 + t218 * t255 - t242;
t134 = m(3) * t191 + qJDD(2) * mrSges(3,1) - t238 * mrSges(3,2) + t240;
t110 = t234 * t113 + t237 * t134;
t115 = t236 * t121 + t233 * t122;
t252 = t237 * t113 - t234 * t134;
t158 = Ifges(6,5) * t179 + Ifges(6,6) * t178 + Ifges(6,3) * t223;
t131 = -mrSges(6,1) * t150 + mrSges(6,3) * t143 + Ifges(6,4) * t157 + Ifges(6,2) * t156 + Ifges(6,6) * t222 - t179 * t158 + t223 * t160;
t132 = mrSges(6,2) * t150 - mrSges(6,3) * t142 + Ifges(6,1) * t157 + Ifges(6,4) * t156 + Ifges(6,5) * t222 + t178 * t158 - t223 * t159;
t175 = Ifges(5,5) * t200 + Ifges(5,6) * t199 + Ifges(5,3) * qJD(3);
t116 = -mrSges(5,1) * t169 + mrSges(5,3) * t148 + Ifges(5,4) * t189 + Ifges(5,2) * t188 + Ifges(5,6) * qJDD(3) - pkin(4) * t248 + pkin(7) * t249 + qJD(3) * t177 + t235 * t131 + t232 * t132 - t200 * t175;
t117 = mrSges(5,2) * t169 - mrSges(5,3) * t147 + Ifges(5,1) * t189 + Ifges(5,4) * t188 + Ifges(5,5) * qJDD(3) - pkin(7) * t130 - qJD(3) * t176 - t232 * t131 + t235 * t132 + t199 * t175;
t201 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t233 + Ifges(4,6) * t236) * qJD(2);
t104 = -mrSges(4,1) * t186 + mrSges(4,3) * t174 + Ifges(4,4) * t211 + Ifges(4,2) * t212 + Ifges(4,6) * qJDD(3) - pkin(3) * t242 + qJ(4) * t250 + qJD(3) * t203 + t230 * t116 + t228 * t117 - t201 * t256;
t106 = mrSges(4,2) * t186 - mrSges(4,3) * t173 + Ifges(4,1) * t211 + Ifges(4,4) * t212 + Ifges(4,5) * qJDD(3) - qJ(4) * t123 - qJD(3) * t202 - t228 * t116 + t230 * t117 + t201 * t255;
t245 = mrSges(3,1) * t191 - mrSges(3,2) * t192 + Ifges(3,3) * qJDD(2) + pkin(2) * t240 + pkin(6) * t251 + t236 * t104 + t233 * t106;
t243 = mrSges(2,1) * t214 - mrSges(2,2) * t215 + pkin(1) * t110 + t245;
t108 = m(2) * t215 + t252;
t107 = m(2) * t214 + t110;
t102 = -mrSges(3,1) * t227 + mrSges(3,3) * t192 + t238 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t115 - t257;
t101 = mrSges(3,2) * t227 - mrSges(3,3) * t191 + Ifges(3,5) * qJDD(2) - t238 * Ifges(3,6) - pkin(6) * t115 - t233 * t104 + t236 * t106;
t100 = mrSges(2,2) * t227 - mrSges(2,3) * t214 - pkin(5) * t110 + t237 * t101 - t234 * t102;
t99 = -mrSges(2,1) * t227 + mrSges(2,3) * t215 + t234 * t101 + t237 * t102 - pkin(1) * (m(3) * t227 + t115) + pkin(5) * t252;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t231 * t100 - t229 * t99 - qJ(1) * (t231 * t107 + t229 * t108), t100, t101, t106, t117, t132; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t229 * t100 + t231 * t99 + qJ(1) * (-t229 * t107 + t231 * t108), t99, t102, t104, t116, t131; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t243, t243, t245, t257, -t241, -t244;];
m_new = t1;
