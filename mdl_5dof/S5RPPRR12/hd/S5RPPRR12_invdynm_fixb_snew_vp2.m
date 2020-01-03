% Calculate vector of cutting torques with Newton-Euler for
% S5RPPRR12
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPPRR12_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR12_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR12_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR12_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR12_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR12_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR12_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:07:01
% EndTime: 2019-12-31 18:07:06
% DurationCPUTime: 3.44s
% Computational Cost: add. (38692->245), mult. (85724->299), div. (0->0), fcn. (55365->8), ass. (0->110)
t233 = qJD(1) ^ 2;
t224 = sin(pkin(8));
t215 = t224 ^ 2;
t225 = cos(pkin(8));
t265 = t225 ^ 2 + t215;
t256 = t265 * mrSges(4,3);
t273 = t233 * t256;
t228 = sin(qJ(1));
t231 = cos(qJ(1));
t204 = g(1) * t228 - t231 * g(2);
t247 = -qJ(2) * t233 + qJDD(2) - t204;
t267 = -pkin(1) - qJ(3);
t272 = -(2 * qJD(1) * qJD(3)) + t267 * qJDD(1) + t247;
t205 = -g(1) * t231 - g(2) * t228;
t271 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t205;
t227 = sin(qJ(4));
t230 = cos(qJ(4));
t251 = t224 * t230 + t225 * t227;
t197 = t251 * qJD(1);
t250 = -t224 * t227 + t225 * t230;
t198 = t250 * qJD(1);
t262 = qJD(4) * t198;
t180 = -t251 * qJDD(1) - t262;
t270 = pkin(3) * t233;
t269 = mrSges(2,1) - mrSges(3,2);
t268 = -Ifges(2,6) + Ifges(3,5);
t266 = Ifges(4,6) * t224;
t177 = t224 * g(3) + t272 * t225;
t161 = (-pkin(6) * qJDD(1) - t224 * t270) * t225 + t177;
t178 = -g(3) * t225 + t272 * t224;
t260 = qJDD(1) * t224;
t162 = -pkin(6) * t260 - t215 * t270 + t178;
t149 = t227 * t161 + t230 * t162;
t172 = mrSges(5,1) * t197 + mrSges(5,2) * t198;
t189 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t198;
t179 = pkin(4) * t197 - pkin(7) * t198;
t232 = qJD(4) ^ 2;
t145 = -pkin(4) * t232 + qJDD(4) * pkin(7) - t179 * t197 + t149;
t246 = qJDD(3) + t271;
t166 = pkin(3) * t260 + (-t265 * pkin(6) + t267) * t233 + t246;
t263 = qJD(4) * t197;
t181 = t250 * qJDD(1) - t263;
t146 = (-t181 + t263) * pkin(7) + (-t180 + t262) * pkin(4) + t166;
t226 = sin(qJ(5));
t229 = cos(qJ(5));
t142 = -t145 * t226 + t146 * t229;
t183 = qJD(4) * t229 - t198 * t226;
t156 = qJD(5) * t183 + qJDD(4) * t226 + t181 * t229;
t184 = qJD(4) * t226 + t198 * t229;
t159 = -mrSges(6,1) * t183 + mrSges(6,2) * t184;
t195 = qJD(5) + t197;
t163 = -mrSges(6,2) * t195 + mrSges(6,3) * t183;
t176 = qJDD(5) - t180;
t138 = m(6) * t142 + mrSges(6,1) * t176 - t156 * mrSges(6,3) - t159 * t184 + t163 * t195;
t143 = t145 * t229 + t146 * t226;
t155 = -qJD(5) * t184 + qJDD(4) * t229 - t181 * t226;
t164 = mrSges(6,1) * t195 - mrSges(6,3) * t184;
t139 = m(6) * t143 - mrSges(6,2) * t176 + t155 * mrSges(6,3) + t159 * t183 - t164 * t195;
t254 = -t138 * t226 + t229 * t139;
t125 = m(5) * t149 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t180 - qJD(4) * t189 - t172 * t197 + t254;
t148 = t161 * t230 - t162 * t227;
t188 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t197;
t144 = -qJDD(4) * pkin(4) - pkin(7) * t232 + t179 * t198 - t148;
t244 = -m(6) * t144 + t155 * mrSges(6,1) - t156 * mrSges(6,2) + t183 * t163 - t164 * t184;
t134 = m(5) * t148 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t181 + qJD(4) * t188 - t172 * t198 + t244;
t119 = t227 * t125 + t230 * t134;
t127 = t229 * t138 + t226 * t139;
t264 = t233 * (Ifges(4,5) * t225 - t266);
t259 = qJDD(1) * t225;
t257 = Ifges(3,4) + t266;
t249 = -qJDD(1) * mrSges(4,3) - t233 * (mrSges(4,1) * t224 + mrSges(4,2) * t225);
t116 = m(4) * t177 + t249 * t225 + t119;
t255 = t230 * t125 - t227 * t134;
t117 = m(4) * t178 + t249 * t224 + t255;
t112 = -t116 * t224 + t225 * t117;
t253 = Ifges(4,1) * t225 - Ifges(4,4) * t224;
t252 = Ifges(4,4) * t225 - Ifges(4,2) * t224;
t111 = t116 * t225 + t117 * t224;
t196 = -qJDD(1) * pkin(1) + t247;
t245 = -m(3) * t196 + t233 * mrSges(3,3) - t111;
t243 = -m(5) * t166 + t180 * mrSges(5,1) - t181 * mrSges(5,2) - t188 * t197 - t198 * t189 - t127;
t150 = Ifges(6,5) * t184 + Ifges(6,6) * t183 + Ifges(6,3) * t195;
t152 = Ifges(6,1) * t184 + Ifges(6,4) * t183 + Ifges(6,5) * t195;
t131 = -mrSges(6,1) * t144 + mrSges(6,3) * t143 + Ifges(6,4) * t156 + Ifges(6,2) * t155 + Ifges(6,6) * t176 - t150 * t184 + t152 * t195;
t151 = Ifges(6,4) * t184 + Ifges(6,2) * t183 + Ifges(6,6) * t195;
t132 = mrSges(6,2) * t144 - mrSges(6,3) * t142 + Ifges(6,1) * t156 + Ifges(6,4) * t155 + Ifges(6,5) * t176 + t150 * t183 - t151 * t195;
t167 = Ifges(5,5) * t198 - Ifges(5,6) * t197 + Ifges(5,3) * qJD(4);
t168 = Ifges(5,4) * t198 - Ifges(5,2) * t197 + Ifges(5,6) * qJD(4);
t113 = mrSges(5,2) * t166 - mrSges(5,3) * t148 + Ifges(5,1) * t181 + Ifges(5,4) * t180 + Ifges(5,5) * qJDD(4) - pkin(7) * t127 - qJD(4) * t168 - t131 * t226 + t132 * t229 - t167 * t197;
t169 = Ifges(5,1) * t198 - Ifges(5,4) * t197 + Ifges(5,5) * qJD(4);
t237 = mrSges(6,1) * t142 - mrSges(6,2) * t143 + Ifges(6,5) * t156 + Ifges(6,6) * t155 + Ifges(6,3) * t176 + t151 * t184 - t152 * t183;
t114 = -mrSges(5,1) * t166 + mrSges(5,3) * t149 + Ifges(5,4) * t181 + Ifges(5,2) * t180 + Ifges(5,6) * qJDD(4) - pkin(4) * t127 + qJD(4) * t169 - t167 * t198 - t237;
t187 = t267 * t233 + t246;
t105 = -mrSges(4,1) * t187 + mrSges(4,3) * t178 + pkin(3) * t243 + pkin(6) * t255 + t252 * qJDD(1) + t227 * t113 + t230 * t114 - t225 * t264;
t107 = mrSges(4,2) * t187 - mrSges(4,3) * t177 - pkin(6) * t119 + t253 * qJDD(1) + t113 * t230 - t114 * t227 - t224 * t264;
t191 = pkin(1) * t233 - t271;
t242 = mrSges(3,2) * t196 - mrSges(3,3) * t191 + Ifges(3,1) * qJDD(1) - qJ(3) * t111 - t105 * t224 + t225 * t107;
t240 = -m(4) * t187 - mrSges(4,1) * t260 - mrSges(4,2) * t259 + t243;
t241 = -mrSges(3,1) * t191 - pkin(2) * (t240 + t273) - qJ(3) * t112 - t105 * t225 - t107 * t224;
t239 = -mrSges(5,1) * t148 + mrSges(5,2) * t149 - Ifges(5,5) * t181 - Ifges(5,6) * t180 - Ifges(5,3) * qJDD(4) - pkin(4) * t244 - pkin(7) * t254 - t229 * t131 - t226 * t132 - t198 * t168 - t197 * t169;
t236 = -m(3) * t191 + t233 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t240;
t238 = -mrSges(2,2) * t205 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t245) + qJ(2) * (t236 - t273) + mrSges(2,1) * t204 + Ifges(2,3) * qJDD(1) + t242;
t235 = -mrSges(4,1) * t177 + mrSges(4,2) * t178 - Ifges(4,5) * t259 - pkin(3) * t119 + t239 + (-t224 * t253 - t225 * t252) * t233;
t234 = -mrSges(3,1) * t196 - pkin(2) * t111 + t235;
t120 = t236 - qJDD(1) * mrSges(2,2) + (-mrSges(2,1) - t256) * t233 + m(2) * t205;
t110 = -m(3) * g(3) + t112;
t108 = m(2) * t204 - mrSges(2,2) * t233 + t269 * qJDD(1) + t245;
t104 = -t234 + t268 * t233 + (Ifges(2,5) - t257) * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) - mrSges(2,3) * t204 - qJ(2) * t110;
t103 = mrSges(2,3) * t205 - pkin(1) * t110 + (-Ifges(3,4) + Ifges(2,5)) * t233 - t268 * qJDD(1) + t269 * g(3) + t241;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t231 * t104 - t228 * t103 - pkin(5) * (t108 * t231 + t120 * t228), t104, t242, t107, t113, t132; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t228 * t104 + t231 * t103 + pkin(5) * (-t108 * t228 + t231 * t120), t103, -mrSges(3,3) * g(3) - t233 * Ifges(3,5) + t257 * qJDD(1) + t234, t105, t114, t131; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t238, t238, mrSges(3,2) * g(3) + Ifges(3,4) * t233 + Ifges(3,5) * qJDD(1) - t241, -Ifges(4,6) * t260 - t235, -t239, t237;];
m_new = t1;
