% Calculate vector of cutting torques with Newton-Euler for
% S5RPRPR1
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRPR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR1_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR1_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR1_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR1_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR1_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:47:08
% EndTime: 2019-12-05 17:47:19
% DurationCPUTime: 4.77s
% Computational Cost: add. (56860->268), mult. (123781->332), div. (0->0), fcn. (77977->8), ass. (0->107)
t236 = sin(qJ(1));
t239 = cos(qJ(1));
t216 = -t239 * g(1) - t236 * g(2);
t253 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t216;
t265 = -pkin(1) - pkin(6);
t264 = mrSges(2,1) - mrSges(3,2);
t263 = Ifges(2,5) - Ifges(3,4);
t262 = (-Ifges(2,6) + Ifges(3,5));
t215 = t236 * g(1) - t239 * g(2);
t240 = qJD(1) ^ 2;
t252 = -t240 * qJ(2) + qJDD(2) - t215;
t188 = t265 * qJDD(1) + t252;
t235 = sin(qJ(3));
t238 = cos(qJ(3));
t179 = t235 * g(3) + t238 * t188;
t259 = qJD(1) * qJD(3);
t257 = t235 * t259;
t210 = t238 * qJDD(1) - t257;
t159 = (-t210 - t257) * qJ(4) + (-t235 * t238 * t240 + qJDD(3)) * pkin(3) + t179;
t180 = -t238 * g(3) + t235 * t188;
t209 = -t235 * qJDD(1) - t238 * t259;
t260 = qJD(1) * t238;
t213 = qJD(3) * pkin(3) - qJ(4) * t260;
t229 = t235 ^ 2;
t160 = -t229 * t240 * pkin(3) + t209 * qJ(4) - qJD(3) * t213 + t180;
t232 = sin(pkin(8));
t233 = cos(pkin(8));
t198 = (-t232 * t235 + t233 * t238) * qJD(1);
t142 = -0.2e1 * qJD(4) * t198 + t233 * t159 - t232 * t160;
t177 = t232 * t209 + t233 * t210;
t197 = (-t232 * t238 - t233 * t235) * qJD(1);
t138 = (qJD(3) * t197 - t177) * pkin(7) + (t197 * t198 + qJDD(3)) * pkin(4) + t142;
t143 = 0.2e1 * qJD(4) * t197 + t232 * t159 + t233 * t160;
t176 = t233 * t209 - t232 * t210;
t187 = qJD(3) * pkin(4) - t198 * pkin(7);
t196 = t197 ^ 2;
t139 = -t196 * pkin(4) + t176 * pkin(7) - qJD(3) * t187 + t143;
t234 = sin(qJ(5));
t237 = cos(qJ(5));
t136 = t237 * t138 - t234 * t139;
t169 = t237 * t197 - t234 * t198;
t150 = t169 * qJD(5) + t234 * t176 + t237 * t177;
t170 = t234 * t197 + t237 * t198;
t155 = -mrSges(6,1) * t169 + mrSges(6,2) * t170;
t222 = qJD(3) + qJD(5);
t164 = -t222 * mrSges(6,2) + t169 * mrSges(6,3);
t221 = qJDD(3) + qJDD(5);
t133 = m(6) * t136 + t221 * mrSges(6,1) - t150 * mrSges(6,3) - t170 * t155 + t222 * t164;
t137 = t234 * t138 + t237 * t139;
t149 = -t170 * qJD(5) + t237 * t176 - t234 * t177;
t165 = t222 * mrSges(6,1) - t170 * mrSges(6,3);
t134 = m(6) * t137 - t221 * mrSges(6,2) + t149 * mrSges(6,3) + t169 * t155 - t222 * t165;
t124 = t237 * t133 + t234 * t134;
t172 = -t197 * mrSges(5,1) + t198 * mrSges(5,2);
t185 = -qJD(3) * mrSges(5,2) + t197 * mrSges(5,3);
t121 = m(5) * t142 + qJDD(3) * mrSges(5,1) - t177 * mrSges(5,3) + qJD(3) * t185 - t198 * t172 + t124;
t186 = qJD(3) * mrSges(5,1) - t198 * mrSges(5,3);
t255 = -t234 * t133 + t237 * t134;
t122 = m(5) * t143 - qJDD(3) * mrSges(5,2) + t176 * mrSges(5,3) - qJD(3) * t186 + t197 * t172 + t255;
t117 = t233 * t121 + t232 * t122;
t261 = qJD(1) * t235;
t208 = (mrSges(4,1) * t235 + mrSges(4,2) * t238) * qJD(1);
t212 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t261;
t114 = m(4) * t179 + qJDD(3) * mrSges(4,1) - t210 * mrSges(4,3) + qJD(3) * t212 - t208 * t260 + t117;
t214 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t260;
t256 = -t232 * t121 + t233 * t122;
t115 = m(4) * t180 - qJDD(3) * mrSges(4,2) + t209 * mrSges(4,3) - qJD(3) * t214 - t208 * t261 + t256;
t110 = -t235 * t114 + t238 * t115;
t109 = t238 * t114 + t235 * t115;
t162 = -t209 * pkin(3) + qJDD(4) + t213 * t260 + (-qJ(4) * t229 + t265) * t240 + t253;
t144 = -t176 * pkin(4) - t196 * pkin(7) + t198 * t187 + t162;
t251 = m(6) * t144 - t149 * mrSges(6,1) + t150 * mrSges(6,2) - t169 * t164 + t170 * t165;
t195 = -qJDD(1) * pkin(1) + t252;
t250 = -m(3) * t195 + (t240 * mrSges(3,3)) - t109;
t152 = Ifges(6,4) * t170 + Ifges(6,2) * t169 + Ifges(6,6) * t222;
t153 = Ifges(6,1) * t170 + Ifges(6,4) * t169 + Ifges(6,5) * t222;
t249 = -mrSges(6,1) * t136 + mrSges(6,2) * t137 - Ifges(6,5) * t150 - Ifges(6,6) * t149 - Ifges(6,3) * t221 - t170 * t152 + t169 * t153;
t151 = Ifges(6,5) * t170 + Ifges(6,6) * t169 + Ifges(6,3) * t222;
t125 = -mrSges(6,1) * t144 + mrSges(6,3) * t137 + Ifges(6,4) * t150 + Ifges(6,2) * t149 + Ifges(6,6) * t221 - t170 * t151 + t222 * t153;
t126 = mrSges(6,2) * t144 - mrSges(6,3) * t136 + Ifges(6,1) * t150 + Ifges(6,4) * t149 + Ifges(6,5) * t221 + t169 * t151 - t222 * t152;
t166 = Ifges(5,5) * t198 + Ifges(5,6) * t197 + (Ifges(5,3) * qJD(3));
t168 = Ifges(5,1) * t198 + Ifges(5,4) * t197 + Ifges(5,5) * qJD(3);
t111 = -mrSges(5,1) * t162 + mrSges(5,3) * t143 + Ifges(5,4) * t177 + Ifges(5,2) * t176 + Ifges(5,6) * qJDD(3) - pkin(4) * t251 + pkin(7) * t255 + qJD(3) * t168 + t237 * t125 + t234 * t126 - t198 * t166;
t167 = Ifges(5,4) * t198 + Ifges(5,2) * t197 + Ifges(5,6) * qJD(3);
t112 = mrSges(5,2) * t162 - mrSges(5,3) * t142 + Ifges(5,1) * t177 + Ifges(5,4) * t176 + Ifges(5,5) * qJDD(3) - pkin(7) * t124 - qJD(3) * t167 - t234 * t125 + t237 * t126 + t197 * t166;
t184 = t265 * t240 + t253;
t199 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t238 - Ifges(4,6) * t235) * qJD(1);
t201 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t238 - Ifges(4,4) * t235) * qJD(1);
t246 = m(5) * t162 - t176 * mrSges(5,1) + t177 * mrSges(5,2) - t197 * t185 + t198 * t186 + t251;
t103 = -mrSges(4,1) * t184 + mrSges(4,3) * t180 + Ifges(4,4) * t210 + Ifges(4,2) * t209 + Ifges(4,6) * qJDD(3) - pkin(3) * t246 + qJ(4) * t256 + qJD(3) * t201 + t233 * t111 + t232 * t112 - t199 * t260;
t200 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t238 - Ifges(4,2) * t235) * qJD(1);
t105 = mrSges(4,2) * t184 - mrSges(4,3) * t179 + Ifges(4,1) * t210 + Ifges(4,4) * t209 + Ifges(4,5) * qJDD(3) - qJ(4) * t117 - qJD(3) * t200 - t232 * t111 + t233 * t112 - t199 * t261;
t191 = t240 * pkin(1) - t253;
t248 = mrSges(3,2) * t195 - mrSges(3,3) * t191 + Ifges(3,1) * qJDD(1) - pkin(6) * t109 - t235 * t103 + t238 * t105;
t129 = -m(4) * t184 + t209 * mrSges(4,1) - t210 * mrSges(4,2) - t212 * t261 - t214 * t260 - t246;
t247 = -mrSges(3,1) * t191 - pkin(2) * t129 - pkin(6) * t110 - t238 * t103 - t235 * t105;
t243 = -m(3) * t191 + t240 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t129;
t245 = -mrSges(2,2) * t216 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t250) + qJ(2) * t243 + mrSges(2,1) * t215 + Ifges(2,3) * qJDD(1) + t248;
t244 = -mrSges(5,1) * t142 + mrSges(5,2) * t143 - Ifges(5,5) * t177 - Ifges(5,6) * t176 - Ifges(5,3) * qJDD(3) - pkin(4) * t124 - t198 * t167 + t197 * t168 + t249;
t242 = -mrSges(4,1) * t179 + mrSges(4,2) * t180 - Ifges(4,5) * t210 - Ifges(4,6) * t209 - Ifges(4,3) * qJDD(3) - pkin(3) * t117 - t200 * t260 - t201 * t261 + t244;
t241 = -mrSges(3,1) * t195 - pkin(2) * t109 + t242;
t127 = m(2) * t216 - t240 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t243;
t108 = -m(3) * g(3) + t110;
t106 = m(2) * t215 - t240 * mrSges(2,2) + t264 * qJDD(1) + t250;
t102 = -qJ(2) * t108 - mrSges(2,3) * t215 - t241 + (-mrSges(2,2) + mrSges(3,3)) * g(3) + t263 * qJDD(1) + (t262 * t240);
t101 = mrSges(2,3) * t216 - pkin(1) * t108 + t264 * g(3) - t262 * qJDD(1) + t263 * t240 + t247;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t239 * t102 - t236 * t101 - pkin(5) * (t239 * t106 + t236 * t127), t102, t248, t105, t112, t126; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t236 * t102 + t239 * t101 + pkin(5) * (-t236 * t106 + t239 * t127), t101, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - (t240 * Ifges(3,5)) + t241, t103, t111, t125; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t245, t245, mrSges(3,2) * g(3) + t240 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t247, -t242, -t244, -t249;];
m_new = t1;
