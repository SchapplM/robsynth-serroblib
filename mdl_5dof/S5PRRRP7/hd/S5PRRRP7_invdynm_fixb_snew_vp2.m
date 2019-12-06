% Calculate vector of cutting torques with Newton-Euler for
% S5PRRRP7
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRRRP7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP7_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP7_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP7_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP7_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP7_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:54:14
% EndTime: 2019-12-05 16:54:26
% DurationCPUTime: 4.62s
% Computational Cost: add. (57059->260), mult. (108089->326), div. (0->0), fcn. (70297->10), ass. (0->107)
t226 = sin(pkin(9));
t228 = cos(pkin(9));
t219 = t226 * g(1) - t228 * g(2);
t220 = -t228 * g(1) - t226 * g(2);
t225 = -g(3) + qJDD(1);
t235 = cos(qJ(2));
t229 = cos(pkin(5));
t232 = sin(qJ(2));
t258 = t229 * t232;
t227 = sin(pkin(5));
t259 = t227 * t232;
t167 = t219 * t258 + t235 * t220 + t225 * t259;
t237 = qJD(2) ^ 2;
t164 = -t237 * pkin(2) + qJDD(2) * pkin(7) + t167;
t196 = -t227 * t219 + t229 * t225;
t231 = sin(qJ(3));
t234 = cos(qJ(3));
t157 = t234 * t164 + t231 * t196;
t215 = (-pkin(3) * t234 - pkin(8) * t231) * qJD(2);
t236 = qJD(3) ^ 2;
t254 = t234 * qJD(2);
t152 = -t236 * pkin(3) + qJDD(3) * pkin(8) + t215 * t254 + t157;
t166 = -t232 * t220 + (t219 * t229 + t225 * t227) * t235;
t163 = -qJDD(2) * pkin(2) - t237 * pkin(7) - t166;
t253 = qJD(2) * qJD(3);
t249 = t234 * t253;
t216 = t231 * qJDD(2) + t249;
t250 = t231 * t253;
t217 = t234 * qJDD(2) - t250;
t155 = (-t216 - t249) * pkin(8) + (-t217 + t250) * pkin(3) + t163;
t230 = sin(qJ(4));
t233 = cos(qJ(4));
t146 = -t230 * t152 + t233 * t155;
t255 = qJD(2) * t231;
t212 = t233 * qJD(3) - t230 * t255;
t186 = t212 * qJD(4) + t230 * qJDD(3) + t233 * t216;
t213 = t230 * qJD(3) + t233 * t255;
t188 = -t212 * mrSges(6,1) + t213 * mrSges(6,2);
t189 = -t212 * mrSges(5,1) + t213 * mrSges(5,2);
t224 = qJD(4) - t254;
t192 = -t224 * mrSges(5,2) + t212 * mrSges(5,3);
t209 = qJDD(4) - t217;
t142 = -0.2e1 * qJD(5) * t213 + (t212 * t224 - t186) * qJ(5) + (t212 * t213 + t209) * pkin(4) + t146;
t191 = -t224 * mrSges(6,2) + t212 * mrSges(6,3);
t252 = m(6) * t142 + t209 * mrSges(6,1) + t224 * t191;
t133 = m(5) * t146 + t209 * mrSges(5,1) + t224 * t192 + (-t188 - t189) * t213 + (-mrSges(5,3) - mrSges(6,3)) * t186 + t252;
t147 = t233 * t152 + t230 * t155;
t185 = -t213 * qJD(4) + t233 * qJDD(3) - t230 * t216;
t193 = t224 * pkin(4) - t213 * qJ(5);
t208 = t212 ^ 2;
t145 = -t208 * pkin(4) + t185 * qJ(5) + 0.2e1 * qJD(5) * t212 - t224 * t193 + t147;
t251 = m(6) * t145 + t185 * mrSges(6,3) + t212 * t188;
t194 = t224 * mrSges(6,1) - t213 * mrSges(6,3);
t256 = -t224 * mrSges(5,1) + t213 * mrSges(5,3) - t194;
t262 = -mrSges(5,2) - mrSges(6,2);
t135 = m(5) * t147 + t185 * mrSges(5,3) + t212 * t189 + t262 * t209 + t256 * t224 + t251;
t132 = -t230 * t133 + t233 * t135;
t214 = (-mrSges(4,1) * t234 + mrSges(4,2) * t231) * qJD(2);
t221 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t255;
t129 = m(4) * t157 - qJDD(3) * mrSges(4,2) + t217 * mrSges(4,3) - qJD(3) * t221 + t214 * t254 + t132;
t156 = -t231 * t164 + t234 * t196;
t151 = -qJDD(3) * pkin(3) - t236 * pkin(8) + t215 * t255 - t156;
t149 = -t185 * pkin(4) - t208 * qJ(5) + t213 * t193 + qJDD(5) + t151;
t247 = -m(6) * t149 + t185 * mrSges(6,1) + t212 * t191;
t138 = -m(5) * t151 + t185 * mrSges(5,1) + t262 * t186 + t212 * t192 + t256 * t213 + t247;
t222 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t254;
t137 = m(4) * t156 + qJDD(3) * mrSges(4,1) - t216 * mrSges(4,3) + qJD(3) * t222 - t214 * t255 + t138;
t123 = t231 * t129 + t234 * t137;
t168 = Ifges(6,5) * t213 + Ifges(6,6) * t212 + Ifges(6,3) * t224;
t169 = Ifges(5,5) * t213 + Ifges(5,6) * t212 + Ifges(5,3) * t224;
t173 = Ifges(5,1) * t213 + Ifges(5,4) * t212 + Ifges(5,5) * t224;
t172 = Ifges(6,1) * t213 + Ifges(6,4) * t212 + Ifges(6,5) * t224;
t244 = -mrSges(6,1) * t149 + mrSges(6,3) * t145 + Ifges(6,4) * t186 + Ifges(6,2) * t185 + Ifges(6,6) * t209 + t224 * t172;
t124 = Ifges(5,4) * t186 + Ifges(5,2) * t185 + Ifges(5,6) * t209 + t224 * t173 - mrSges(5,1) * t151 + mrSges(5,3) * t147 - pkin(4) * (t186 * mrSges(6,2) - t247) + qJ(5) * (-t209 * mrSges(6,2) - t224 * t194 + t251) + (-pkin(4) * t194 - t168 - t169) * t213 + t244;
t139 = -t186 * mrSges(6,3) - t213 * t188 + t252;
t170 = Ifges(6,4) * t213 + Ifges(6,2) * t212 + Ifges(6,6) * t224;
t171 = Ifges(5,4) * t213 + Ifges(5,2) * t212 + Ifges(5,6) * t224;
t242 = mrSges(6,2) * t149 - mrSges(6,3) * t142 + Ifges(6,1) * t186 + Ifges(6,4) * t185 + Ifges(6,5) * t209 + t212 * t168;
t130 = mrSges(5,2) * t151 - mrSges(5,3) * t146 + Ifges(5,1) * t186 + Ifges(5,4) * t185 + Ifges(5,5) * t209 - qJ(5) * t139 + t212 * t169 + (-t170 - t171) * t224 + t242;
t200 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t231 + Ifges(4,2) * t234) * qJD(2);
t201 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t231 + Ifges(4,4) * t234) * qJD(2);
t264 = mrSges(4,1) * t156 - mrSges(4,2) * t157 + Ifges(4,5) * t216 + Ifges(4,6) * t217 + Ifges(4,3) * qJDD(3) + pkin(3) * t138 + pkin(8) * t132 + t233 * t124 + t230 * t130 + (t231 * t200 - t234 * t201) * qJD(2);
t109 = -mrSges(3,1) * t196 + mrSges(3,3) * t167 + t237 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t123 - t264;
t248 = t234 * t129 - t231 * t137;
t121 = m(3) * t167 - t237 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t248;
t131 = t233 * t133 + t230 * t135;
t240 = -m(4) * t163 + t217 * mrSges(4,1) - t216 * mrSges(4,2) - t221 * t255 + t222 * t254 - t131;
t126 = m(3) * t166 + qJDD(2) * mrSges(3,1) - t237 * mrSges(3,2) + t240;
t117 = t235 * t121 - t232 * t126;
t266 = pkin(6) * t117 + t109 * t235;
t243 = -mrSges(6,1) * t142 + mrSges(6,2) * t145 - Ifges(6,5) * t186 - Ifges(6,6) * t185 - Ifges(6,3) * t209 - t213 * t170;
t265 = mrSges(5,1) * t146 - mrSges(5,2) * t147 + Ifges(5,5) * t186 + Ifges(5,6) * t185 + Ifges(5,3) * t209 + pkin(4) * t139 + t213 * t171 - (t173 + t172) * t212 - t243;
t260 = t126 * t235;
t122 = m(3) * t196 + t123;
t113 = t121 * t258 - t227 * t122 + t229 * t260;
t199 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t231 + Ifges(4,6) * t234) * qJD(2);
t114 = mrSges(4,2) * t163 - mrSges(4,3) * t156 + Ifges(4,1) * t216 + Ifges(4,4) * t217 + Ifges(4,5) * qJDD(3) - pkin(8) * t131 - qJD(3) * t200 - t230 * t124 + t233 * t130 + t199 * t254;
t118 = -mrSges(4,1) * t163 + mrSges(4,3) * t157 + Ifges(4,4) * t216 + Ifges(4,2) * t217 + Ifges(4,6) * qJDD(3) - pkin(3) * t131 + qJD(3) * t201 - t199 * t255 - t265;
t105 = mrSges(3,1) * t166 - mrSges(3,2) * t167 + Ifges(3,3) * qJDD(2) + pkin(2) * t240 + pkin(7) * t248 + t231 * t114 + t234 * t118;
t107 = mrSges(3,2) * t196 - mrSges(3,3) * t166 + Ifges(3,5) * qJDD(2) - t237 * Ifges(3,6) - pkin(7) * t123 + t234 * t114 - t231 * t118;
t241 = mrSges(2,1) * t219 - mrSges(2,2) * t220 + pkin(1) * t113 + t229 * t105 + t107 * t259 + t266 * t227;
t115 = m(2) * t220 + t117;
t112 = t229 * t122 + (t121 * t232 + t260) * t227;
t110 = m(2) * t219 + t113;
t103 = mrSges(2,2) * t225 - mrSges(2,3) * t219 + t235 * t107 - t232 * t109 + (-t112 * t227 - t113 * t229) * pkin(6);
t102 = -mrSges(2,1) * t225 + mrSges(2,3) * t220 - pkin(1) * t112 - t227 * t105 + (t107 * t232 + t266) * t229;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t228 * t103 - t226 * t102 - qJ(1) * (t228 * t110 + t226 * t115), t103, t107, t114, t130, -t224 * t170 + t242; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t226 * t103 + t228 * t102 + qJ(1) * (-t226 * t110 + t228 * t115), t102, t109, t118, t124, -t213 * t168 + t244; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t241, t241, t105, t264, t265, -t212 * t172 - t243;];
m_new = t1;
