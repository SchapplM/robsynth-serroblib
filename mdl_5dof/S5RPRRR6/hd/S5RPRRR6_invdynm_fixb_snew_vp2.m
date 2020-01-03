% Calculate vector of cutting torques with Newton-Euler for
% S5RPRRR6
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRRR6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR6_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR6_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR6_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR6_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR6_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:01:04
% EndTime: 2019-12-31 19:01:12
% DurationCPUTime: 5.79s
% Computational Cost: add. (79124->267), mult. (154787->341), div. (0->0), fcn. (98225->10), ass. (0->110)
t231 = sin(qJ(4));
t232 = sin(qJ(3));
t235 = cos(qJ(4));
t236 = cos(qJ(3));
t201 = (t231 * t232 - t235 * t236) * qJD(1);
t233 = sin(qJ(1));
t237 = cos(qJ(1));
t215 = t233 * g(1) - t237 * g(2);
t206 = qJDD(1) * pkin(1) + t215;
t216 = -t237 * g(1) - t233 * g(2);
t238 = qJD(1) ^ 2;
t208 = -t238 * pkin(1) + t216;
t228 = sin(pkin(9));
t229 = cos(pkin(9));
t189 = t228 * t206 + t229 * t208;
t184 = -t238 * pkin(2) + qJDD(1) * pkin(6) + t189;
t227 = -g(3) + qJDD(2);
t172 = -t232 * t184 + t236 * t227;
t255 = qJD(1) * qJD(3);
t254 = t236 * t255;
t209 = t232 * qJDD(1) + t254;
t159 = (-t209 + t254) * pkin(7) + (t232 * t236 * t238 + qJDD(3)) * pkin(3) + t172;
t173 = t236 * t184 + t232 * t227;
t210 = t236 * qJDD(1) - t232 * t255;
t257 = qJD(1) * t232;
t214 = qJD(3) * pkin(3) - pkin(7) * t257;
t226 = t236 ^ 2;
t160 = -t226 * t238 * pkin(3) + t210 * pkin(7) - qJD(3) * t214 + t173;
t153 = t231 * t159 + t235 * t160;
t202 = (t231 * t236 + t232 * t235) * qJD(1);
t174 = -t202 * qJD(4) - t231 * t209 + t235 * t210;
t185 = t201 * mrSges(5,1) + t202 * mrSges(5,2);
t223 = qJD(3) + qJD(4);
t193 = t223 * mrSges(5,1) - t202 * mrSges(5,3);
t222 = qJDD(3) + qJDD(4);
t186 = t201 * pkin(4) - t202 * pkin(8);
t221 = t223 ^ 2;
t149 = -t221 * pkin(4) + t222 * pkin(8) - t201 * t186 + t153;
t188 = t229 * t206 - t228 * t208;
t247 = -qJDD(1) * pkin(2) - t188;
t165 = -t210 * pkin(3) + t214 * t257 + (-pkin(7) * t226 - pkin(6)) * t238 + t247;
t175 = -t201 * qJD(4) + t235 * t209 + t231 * t210;
t150 = (t201 * t223 - t175) * pkin(8) + (t202 * t223 - t174) * pkin(4) + t165;
t230 = sin(qJ(5));
t234 = cos(qJ(5));
t146 = -t230 * t149 + t234 * t150;
t190 = -t230 * t202 + t234 * t223;
t156 = t190 * qJD(5) + t234 * t175 + t230 * t222;
t191 = t234 * t202 + t230 * t223;
t166 = -t190 * mrSges(6,1) + t191 * mrSges(6,2);
t171 = qJDD(5) - t174;
t194 = qJD(5) + t201;
t176 = -t194 * mrSges(6,2) + t190 * mrSges(6,3);
t142 = m(6) * t146 + t171 * mrSges(6,1) - t156 * mrSges(6,3) - t191 * t166 + t194 * t176;
t147 = t234 * t149 + t230 * t150;
t155 = -t191 * qJD(5) - t230 * t175 + t234 * t222;
t177 = t194 * mrSges(6,1) - t191 * mrSges(6,3);
t143 = m(6) * t147 - t171 * mrSges(6,2) + t155 * mrSges(6,3) + t190 * t166 - t194 * t177;
t250 = -t230 * t142 + t234 * t143;
t129 = m(5) * t153 - t222 * mrSges(5,2) + t174 * mrSges(5,3) - t201 * t185 - t223 * t193 + t250;
t152 = t235 * t159 - t231 * t160;
t192 = -t223 * mrSges(5,2) - t201 * mrSges(5,3);
t148 = -t222 * pkin(4) - t221 * pkin(8) + t202 * t186 - t152;
t245 = -m(6) * t148 + t155 * mrSges(6,1) - t156 * mrSges(6,2) + t190 * t176 - t191 * t177;
t138 = m(5) * t152 + t222 * mrSges(5,1) - t175 * mrSges(5,3) - t202 * t185 + t223 * t192 + t245;
t123 = t231 * t129 + t235 * t138;
t199 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t232 + Ifges(4,2) * t236) * qJD(1);
t200 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t232 + Ifges(4,4) * t236) * qJD(1);
t161 = Ifges(6,5) * t191 + Ifges(6,6) * t190 + Ifges(6,3) * t194;
t163 = Ifges(6,1) * t191 + Ifges(6,4) * t190 + Ifges(6,5) * t194;
t135 = -mrSges(6,1) * t148 + mrSges(6,3) * t147 + Ifges(6,4) * t156 + Ifges(6,2) * t155 + Ifges(6,6) * t171 - t191 * t161 + t194 * t163;
t162 = Ifges(6,4) * t191 + Ifges(6,2) * t190 + Ifges(6,6) * t194;
t136 = mrSges(6,2) * t148 - mrSges(6,3) * t146 + Ifges(6,1) * t156 + Ifges(6,4) * t155 + Ifges(6,5) * t171 + t190 * t161 - t194 * t162;
t180 = Ifges(5,4) * t202 - Ifges(5,2) * t201 + Ifges(5,6) * t223;
t181 = Ifges(5,1) * t202 - Ifges(5,4) * t201 + Ifges(5,5) * t223;
t242 = -mrSges(5,1) * t152 + mrSges(5,2) * t153 - Ifges(5,5) * t175 - Ifges(5,6) * t174 - Ifges(5,3) * t222 - pkin(4) * t245 - pkin(8) * t250 - t234 * t135 - t230 * t136 - t202 * t180 - t201 * t181;
t258 = mrSges(4,1) * t172 - mrSges(4,2) * t173 + Ifges(4,5) * t209 + Ifges(4,6) * t210 + Ifges(4,3) * qJDD(3) + pkin(3) * t123 + (t232 * t199 - t236 * t200) * qJD(1) - t242;
t207 = (-mrSges(4,1) * t236 + mrSges(4,2) * t232) * qJD(1);
t256 = qJD(1) * t236;
t213 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t256;
t121 = m(4) * t172 + qJDD(3) * mrSges(4,1) - t209 * mrSges(4,3) + qJD(3) * t213 - t207 * t257 + t123;
t212 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t257;
t251 = t235 * t129 - t231 * t138;
t122 = m(4) * t173 - qJDD(3) * mrSges(4,2) + t210 * mrSges(4,3) - qJD(3) * t212 + t207 * t256 + t251;
t252 = -t232 * t121 + t236 * t122;
t113 = m(3) * t189 - t238 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t252;
t183 = -t238 * pkin(6) + t247;
t131 = t234 * t142 + t230 * t143;
t244 = m(5) * t165 - t174 * mrSges(5,1) + t175 * mrSges(5,2) + t201 * t192 + t202 * t193 + t131;
t240 = -m(4) * t183 + t210 * mrSges(4,1) - t209 * mrSges(4,2) - t212 * t257 + t213 * t256 - t244;
t125 = m(3) * t188 + qJDD(1) * mrSges(3,1) - t238 * mrSges(3,2) + t240;
t110 = t228 * t113 + t229 * t125;
t115 = t236 * t121 + t232 * t122;
t253 = t229 * t113 - t228 * t125;
t179 = Ifges(5,5) * t202 - Ifges(5,6) * t201 + Ifges(5,3) * t223;
t116 = mrSges(5,2) * t165 - mrSges(5,3) * t152 + Ifges(5,1) * t175 + Ifges(5,4) * t174 + Ifges(5,5) * t222 - pkin(8) * t131 - t230 * t135 + t234 * t136 - t201 * t179 - t223 * t180;
t241 = mrSges(6,1) * t146 - mrSges(6,2) * t147 + Ifges(6,5) * t156 + Ifges(6,6) * t155 + Ifges(6,3) * t171 + t191 * t162 - t190 * t163;
t117 = -mrSges(5,1) * t165 + mrSges(5,3) * t153 + Ifges(5,4) * t175 + Ifges(5,2) * t174 + Ifges(5,6) * t222 - pkin(4) * t131 - t202 * t179 + t223 * t181 - t241;
t198 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t232 + Ifges(4,6) * t236) * qJD(1);
t104 = -mrSges(4,1) * t183 + mrSges(4,3) * t173 + Ifges(4,4) * t209 + Ifges(4,2) * t210 + Ifges(4,6) * qJDD(3) - pkin(3) * t244 + pkin(7) * t251 + qJD(3) * t200 + t231 * t116 + t235 * t117 - t198 * t257;
t106 = mrSges(4,2) * t183 - mrSges(4,3) * t172 + Ifges(4,1) * t209 + Ifges(4,4) * t210 + Ifges(4,5) * qJDD(3) - pkin(7) * t123 - qJD(3) * t199 + t235 * t116 - t231 * t117 + t198 * t256;
t246 = mrSges(3,1) * t188 - mrSges(3,2) * t189 + Ifges(3,3) * qJDD(1) + pkin(2) * t240 + pkin(6) * t252 + t236 * t104 + t232 * t106;
t243 = mrSges(2,1) * t215 - mrSges(2,2) * t216 + Ifges(2,3) * qJDD(1) + pkin(1) * t110 + t246;
t108 = m(2) * t216 - t238 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t253;
t107 = m(2) * t215 + qJDD(1) * mrSges(2,1) - t238 * mrSges(2,2) + t110;
t102 = -mrSges(3,1) * t227 + mrSges(3,3) * t189 + t238 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t115 - t258;
t101 = mrSges(3,2) * t227 - mrSges(3,3) * t188 + Ifges(3,5) * qJDD(1) - t238 * Ifges(3,6) - pkin(6) * t115 - t232 * t104 + t236 * t106;
t100 = -mrSges(2,2) * g(3) - mrSges(2,3) * t215 + Ifges(2,5) * qJDD(1) - t238 * Ifges(2,6) - qJ(2) * t110 + t229 * t101 - t228 * t102;
t99 = Ifges(2,6) * qJDD(1) + t238 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t216 + t228 * t101 + t229 * t102 - pkin(1) * (m(3) * t227 + t115) + qJ(2) * t253;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t237 * t100 - t233 * t99 - pkin(5) * (t237 * t107 + t233 * t108), t100, t101, t106, t116, t136; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t233 * t100 + t237 * t99 + pkin(5) * (-t233 * t107 + t237 * t108), t99, t102, t104, t117, t135; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t243, t243, t246, t258, -t242, t241;];
m_new = t1;
