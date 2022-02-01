% Calculate vector of cutting torques with Newton-Euler for
% S5RPRRP3
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% m [6x1]
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
% Datum: 2022-01-23 09:30
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRRP3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP3_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP3_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP3_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP3_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:29:39
% EndTime: 2022-01-23 09:29:45
% DurationCPUTime: 3.51s
% Computational Cost: add. (40477->262), mult. (80564->323), div. (0->0), fcn. (47650->8), ass. (0->99)
t238 = sin(qJ(4));
t239 = sin(qJ(3));
t241 = cos(qJ(4));
t242 = cos(qJ(3));
t206 = (-t238 * t239 + t241 * t242) * qJD(1);
t264 = qJD(1) * qJD(3);
t261 = t242 * t264;
t214 = t239 * qJDD(1) + t261;
t215 = t242 * qJDD(1) - t239 * t264;
t174 = t206 * qJD(4) + t241 * t214 + t238 * t215;
t207 = (t238 * t242 + t239 * t241) * qJD(1);
t188 = -t206 * mrSges(6,1) + t207 * mrSges(6,2);
t240 = sin(qJ(1));
t243 = cos(qJ(1));
t220 = t240 * g(1) - t243 * g(2);
t211 = qJDD(1) * pkin(1) + t220;
t221 = -t243 * g(1) - t240 * g(2);
t244 = qJD(1) ^ 2;
t213 = -t244 * pkin(1) + t221;
t236 = sin(pkin(8));
t237 = cos(pkin(8));
t192 = t236 * t211 + t237 * t213;
t187 = -t244 * pkin(2) + qJDD(1) * pkin(6) + t192;
t235 = -g(3) + qJDD(2);
t171 = -t239 * t187 + t242 * t235;
t152 = (-t214 + t261) * pkin(7) + (t239 * t242 * t244 + qJDD(3)) * pkin(3) + t171;
t172 = t242 * t187 + t239 * t235;
t266 = qJD(1) * t239;
t219 = qJD(3) * pkin(3) - pkin(7) * t266;
t234 = t242 ^ 2;
t153 = -t234 * t244 * pkin(3) + t215 * pkin(7) - qJD(3) * t219 + t172;
t147 = t241 * t152 - t238 * t153;
t230 = qJDD(3) + qJDD(4);
t231 = qJD(3) + qJD(4);
t139 = -0.2e1 * qJD(5) * t207 + (t206 * t231 - t174) * qJ(5) + (t206 * t207 + t230) * pkin(4) + t147;
t194 = -t231 * mrSges(6,2) + t206 * mrSges(6,3);
t263 = m(6) * t139 + t230 * mrSges(6,1) + t231 * t194;
t136 = -t174 * mrSges(6,3) - t207 * t188 + t263;
t148 = t238 * t152 + t241 * t153;
t173 = -t207 * qJD(4) - t238 * t214 + t241 * t215;
t180 = Ifges(5,4) * t207 + Ifges(5,2) * t206 + Ifges(5,6) * t231;
t181 = Ifges(6,1) * t207 + Ifges(6,4) * t206 + Ifges(6,5) * t231;
t182 = Ifges(5,1) * t207 + Ifges(5,4) * t206 + Ifges(5,5) * t231;
t196 = t231 * pkin(4) - t207 * qJ(5);
t199 = t206 ^ 2;
t142 = -t199 * pkin(4) + t173 * qJ(5) + 0.2e1 * qJD(5) * t206 - t231 * t196 + t148;
t179 = Ifges(6,4) * t207 + Ifges(6,2) * t206 + Ifges(6,6) * t231;
t252 = -mrSges(6,1) * t139 + mrSges(6,2) * t142 - Ifges(6,5) * t174 - Ifges(6,6) * t173 - Ifges(6,3) * t230 - t207 * t179;
t270 = mrSges(5,1) * t147 - mrSges(5,2) * t148 + Ifges(5,5) * t174 + Ifges(5,6) * t173 + Ifges(5,3) * t230 + pkin(4) * t136 + t207 * t180 - t252 + (-t182 - t181) * t206;
t189 = -t206 * mrSges(5,1) + t207 * mrSges(5,2);
t195 = -t231 * mrSges(5,2) + t206 * mrSges(5,3);
t130 = m(5) * t147 + t230 * mrSges(5,1) + t231 * t195 + (-t188 - t189) * t207 + (-mrSges(5,3) - mrSges(6,3)) * t174 + t263;
t197 = t231 * mrSges(6,1) - t207 * mrSges(6,3);
t198 = t231 * mrSges(5,1) - t207 * mrSges(5,3);
t262 = m(6) * t142 + t173 * mrSges(6,3) + t206 * t188;
t133 = m(5) * t148 + t173 * mrSges(5,3) + t206 * t189 + (-t197 - t198) * t231 + (-mrSges(5,2) - mrSges(6,2)) * t230 + t262;
t126 = t241 * t130 + t238 * t133;
t204 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t239 + Ifges(4,2) * t242) * qJD(1);
t205 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t239 + Ifges(4,4) * t242) * qJD(1);
t269 = mrSges(4,1) * t171 - mrSges(4,2) * t172 + Ifges(4,5) * t214 + Ifges(4,6) * t215 + Ifges(4,3) * qJDD(3) + pkin(3) * t126 + (t239 * t204 - t242 * t205) * qJD(1) + t270;
t212 = (-mrSges(4,1) * t242 + mrSges(4,2) * t239) * qJD(1);
t265 = qJD(1) * t242;
t218 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t265;
t124 = m(4) * t171 + qJDD(3) * mrSges(4,1) - t214 * mrSges(4,3) + qJD(3) * t218 - t212 * t266 + t126;
t217 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t266;
t257 = -t238 * t130 + t241 * t133;
t125 = m(4) * t172 - qJDD(3) * mrSges(4,2) + t215 * mrSges(4,3) - qJD(3) * t217 + t212 * t265 + t257;
t258 = -t239 * t124 + t242 * t125;
t116 = m(3) * t192 - t244 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t258;
t191 = t237 * t211 - t236 * t213;
t254 = -qJDD(1) * pkin(2) - t191;
t186 = -t244 * pkin(6) + t254;
t154 = -t215 * pkin(3) + t219 * t266 + (-pkin(7) * t234 - pkin(6)) * t244 + t254;
t145 = -t173 * pkin(4) - t199 * qJ(5) + t207 * t196 + qJDD(5) + t154;
t256 = m(6) * t145 - t173 * mrSges(6,1) + t174 * mrSges(6,2) - t206 * t194 + t207 * t197;
t248 = m(5) * t154 - t173 * mrSges(5,1) + t174 * mrSges(5,2) - t206 * t195 + t207 * t198 + t256;
t246 = -m(4) * t186 + t215 * mrSges(4,1) - t214 * mrSges(4,2) - t217 * t266 + t218 * t265 - t248;
t128 = m(3) * t191 + qJDD(1) * mrSges(3,1) - t244 * mrSges(3,2) + t246;
t113 = t236 * t116 + t237 * t128;
t118 = t242 * t124 + t239 * t125;
t259 = t237 * t116 - t236 * t128;
t253 = -mrSges(6,1) * t145 + mrSges(6,3) * t142 + Ifges(6,4) * t174 + Ifges(6,2) * t173 + Ifges(6,6) * t230 + t231 * t181;
t177 = Ifges(6,5) * t207 + Ifges(6,6) * t206 + Ifges(6,3) * t231;
t251 = mrSges(6,2) * t145 - mrSges(6,3) * t139 + Ifges(6,1) * t174 + Ifges(6,4) * t173 + Ifges(6,5) * t230 + t206 * t177;
t178 = Ifges(5,5) * t207 + Ifges(5,6) * t206 + Ifges(5,3) * t231;
t119 = Ifges(5,4) * t174 + Ifges(5,2) * t173 + Ifges(5,6) * t230 + t231 * t182 - mrSges(5,1) * t154 + mrSges(5,3) * t148 - pkin(4) * t256 + qJ(5) * (-t230 * mrSges(6,2) - t231 * t197 + t262) + (-t178 - t177) * t207 + t253;
t120 = mrSges(5,2) * t154 - mrSges(5,3) * t147 + Ifges(5,1) * t174 + Ifges(5,4) * t173 + Ifges(5,5) * t230 - qJ(5) * t136 + t206 * t178 + (-t179 - t180) * t231 + t251;
t203 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t239 + Ifges(4,6) * t242) * qJD(1);
t107 = -mrSges(4,1) * t186 + mrSges(4,3) * t172 + Ifges(4,4) * t214 + Ifges(4,2) * t215 + Ifges(4,6) * qJDD(3) - pkin(3) * t248 + pkin(7) * t257 + qJD(3) * t205 + t241 * t119 + t238 * t120 - t203 * t266;
t109 = mrSges(4,2) * t186 - mrSges(4,3) * t171 + Ifges(4,1) * t214 + Ifges(4,4) * t215 + Ifges(4,5) * qJDD(3) - pkin(7) * t126 - qJD(3) * t204 - t238 * t119 + t241 * t120 + t203 * t265;
t250 = mrSges(3,1) * t191 - mrSges(3,2) * t192 + Ifges(3,3) * qJDD(1) + pkin(2) * t246 + pkin(6) * t258 + t242 * t107 + t239 * t109;
t249 = mrSges(2,1) * t220 - mrSges(2,2) * t221 + Ifges(2,3) * qJDD(1) + pkin(1) * t113 + t250;
t111 = m(2) * t221 - t244 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t259;
t110 = m(2) * t220 + qJDD(1) * mrSges(2,1) - t244 * mrSges(2,2) + t113;
t105 = -mrSges(3,1) * t235 + mrSges(3,3) * t192 + t244 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t118 - t269;
t104 = mrSges(3,2) * t235 - mrSges(3,3) * t191 + Ifges(3,5) * qJDD(1) - t244 * Ifges(3,6) - pkin(6) * t118 - t239 * t107 + t242 * t109;
t103 = -mrSges(2,2) * g(3) - mrSges(2,3) * t220 + Ifges(2,5) * qJDD(1) - t244 * Ifges(2,6) - qJ(2) * t113 + t237 * t104 - t236 * t105;
t102 = Ifges(2,6) * qJDD(1) + t244 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t221 + t236 * t104 + t237 * t105 - pkin(1) * (m(3) * t235 + t118) + qJ(2) * t259;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t243 * t103 - t240 * t102 - pkin(5) * (t243 * t110 + t240 * t111), t103, t104, t109, t120, -t231 * t179 + t251; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t240 * t103 + t243 * t102 + pkin(5) * (-t240 * t110 + t243 * t111), t102, t105, t107, t119, -t207 * t177 + t253; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t249, t249, t250, t269, t270, -t206 * t181 - t252;];
m_new = t1;
