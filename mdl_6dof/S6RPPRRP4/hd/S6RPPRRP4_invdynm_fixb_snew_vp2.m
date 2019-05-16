% Calculate vector of cutting torques with Newton-Euler for
% S6RPPRRP4
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% m [3x7]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 14:56
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPPRRP4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP4_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP4_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP4_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP4_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP4_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:54:02
% EndTime: 2019-05-05 14:54:08
% DurationCPUTime: 4.05s
% Computational Cost: add. (59186->294), mult. (104839->345), div. (0->0), fcn. (51330->8), ass. (0->113)
t249 = qJD(1) ^ 2;
t243 = sin(qJ(1));
t246 = cos(qJ(1));
t221 = -t246 * g(1) - t243 * g(2);
t261 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t221;
t278 = -pkin(1) - pkin(2);
t194 = t278 * t249 + t261;
t220 = t243 * g(1) - t246 * g(2);
t260 = -t249 * qJ(2) + qJDD(2) - t220;
t197 = t278 * qJDD(1) + t260;
t240 = sin(pkin(9));
t241 = cos(pkin(9));
t163 = -t240 * t194 + t241 * t197;
t159 = qJDD(1) * pkin(3) - t249 * pkin(7) - t163;
t242 = sin(qJ(4));
t245 = cos(qJ(4));
t269 = qJD(1) * qJD(4);
t266 = t245 * t269;
t215 = -t242 * qJDD(1) - t266;
t267 = t242 * t269;
t216 = -t245 * qJDD(1) + t267;
t151 = (-t215 + t266) * pkin(8) + (-t216 - t267) * pkin(4) + t159;
t164 = t241 * t194 + t240 * t197;
t160 = -t249 * pkin(3) - qJDD(1) * pkin(7) + t164;
t236 = g(3) + qJDD(3);
t156 = t245 * t160 + t242 * t236;
t214 = (pkin(4) * t245 + pkin(8) * t242) * qJD(1);
t248 = qJD(4) ^ 2;
t270 = t245 * qJD(1);
t154 = -t248 * pkin(4) + qJDD(4) * pkin(8) - t214 * t270 + t156;
t244 = cos(qJ(5));
t277 = sin(qJ(5));
t148 = t244 * t151 - t277 * t154;
t149 = t277 * t151 + t244 * t154;
t271 = qJD(1) * t242;
t211 = t244 * qJD(4) + t277 * t271;
t212 = -t277 * qJD(4) + t244 * t271;
t223 = qJD(5) + t270;
t166 = -Ifges(7,5) * t212 + Ifges(7,6) * t223 - Ifges(7,3) * t211;
t169 = -Ifges(6,4) * t212 + Ifges(6,2) * t211 + Ifges(6,6) * t223;
t171 = -Ifges(6,1) * t212 + Ifges(6,4) * t211 + Ifges(6,5) * t223;
t180 = -t212 * qJD(5) - t244 * qJDD(4) + t277 * t215;
t181 = t211 * qJD(5) + t277 * qJDD(4) + t244 * t215;
t185 = -t211 * mrSges(7,1) + t212 * mrSges(7,3);
t210 = qJDD(5) - t216;
t184 = -t211 * pkin(5) + t212 * qJ(6);
t222 = t223 ^ 2;
t279 = 2 * qJD(6);
t144 = -t222 * pkin(5) + t210 * qJ(6) + t211 * t184 + t223 * t279 + t149;
t146 = -t210 * pkin(5) - t222 * qJ(6) - t212 * t184 + qJDD(6) - t148;
t170 = -Ifges(7,1) * t212 + Ifges(7,4) * t223 - Ifges(7,5) * t211;
t259 = mrSges(7,1) * t146 - mrSges(7,3) * t144 - Ifges(7,4) * t181 - Ifges(7,2) * t210 - Ifges(7,6) * t180 + t211 * t170;
t193 = t211 * mrSges(7,2) + t223 * mrSges(7,3);
t265 = -m(7) * t146 + t210 * mrSges(7,1) + t223 * t193;
t192 = -t223 * mrSges(7,1) - t212 * mrSges(7,2);
t268 = m(7) * t144 + t210 * mrSges(7,3) + t223 * t192;
t281 = -(t169 - t166) * t212 + mrSges(6,1) * t148 - mrSges(6,2) * t149 + Ifges(6,5) * t181 - Ifges(6,6) * t180 + Ifges(6,3) * t210 + pkin(5) * (-t181 * mrSges(7,2) + t212 * t185 + t265) + qJ(6) * (-t180 * mrSges(7,2) + t211 * t185 + t268) - t211 * t171 - t259;
t155 = -t242 * t160 + t245 * t236;
t153 = -qJDD(4) * pkin(4) - t248 * pkin(8) - t214 * t271 - t155;
t147 = t212 * t279 + (-t211 * t223 - t181) * qJ(6) + (-t212 * t223 + t180) * pkin(5) + t153;
t141 = m(7) * t147 + t180 * mrSges(7,1) - t181 * mrSges(7,3) + t212 * t192 - t211 * t193;
t264 = -mrSges(7,1) * t147 + mrSges(7,2) * t144;
t168 = -Ifges(7,4) * t212 + Ifges(7,2) * t223 - Ifges(7,6) * t211;
t274 = -Ifges(6,5) * t212 + Ifges(6,6) * t211 + Ifges(6,3) * t223 + t168;
t130 = -mrSges(6,1) * t153 + mrSges(6,3) * t149 - pkin(5) * t141 + (t170 + t171) * t223 + t274 * t212 + (Ifges(6,6) - Ifges(7,6)) * t210 + (Ifges(6,4) - Ifges(7,5)) * t181 + (-Ifges(6,2) - Ifges(7,3)) * t180 + t264;
t258 = mrSges(7,2) * t146 - mrSges(7,3) * t147 + Ifges(7,1) * t181 + Ifges(7,4) * t210 + Ifges(7,5) * t180 + t223 * t166;
t131 = mrSges(6,2) * t153 - mrSges(6,3) * t148 + Ifges(6,1) * t181 - Ifges(6,4) * t180 + Ifges(6,5) * t210 - qJ(6) * t141 - t223 * t169 + t274 * t211 + t258;
t191 = t223 * mrSges(6,1) + t212 * mrSges(6,3);
t272 = -t211 * mrSges(6,1) - t212 * mrSges(6,2) + t185;
t275 = -mrSges(6,3) - mrSges(7,2);
t136 = m(6) * t149 - t210 * mrSges(6,2) + t275 * t180 - t223 * t191 + t272 * t211 + t268;
t190 = -t223 * mrSges(6,2) + t211 * mrSges(6,3);
t137 = m(6) * t148 + t210 * mrSges(6,1) + t275 * t181 + t223 * t190 + t272 * t212 + t265;
t133 = t244 * t136 - t277 * t137;
t138 = -m(6) * t153 - t180 * mrSges(6,1) - t181 * mrSges(6,2) + t211 * t190 + t212 * t191 - t141;
t202 = (Ifges(5,6) * qJD(4)) + (-Ifges(5,4) * t242 - Ifges(5,2) * t245) * qJD(1);
t203 = (Ifges(5,5) * qJD(4)) + (-Ifges(5,1) * t242 - Ifges(5,4) * t245) * qJD(1);
t280 = mrSges(5,1) * t155 - mrSges(5,2) * t156 + Ifges(5,5) * t215 + Ifges(5,6) * t216 + Ifges(5,3) * qJDD(4) + pkin(4) * t138 + pkin(8) * t133 - (t202 * t242 - t203 * t245) * qJD(1) + t244 * t130 + t277 * t131;
t276 = mrSges(2,1) + mrSges(3,1);
t213 = (mrSges(5,1) * t245 - mrSges(5,2) * t242) * qJD(1);
t218 = qJD(4) * mrSges(5,1) + mrSges(5,3) * t271;
t128 = m(5) * t156 - qJDD(4) * mrSges(5,2) + t216 * mrSges(5,3) - qJD(4) * t218 - t213 * t270 + t133;
t219 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t270;
t134 = m(5) * t155 + qJDD(4) * mrSges(5,1) - t215 * mrSges(5,3) + qJD(4) * t219 + t213 * t271 + t138;
t125 = t245 * t128 - t242 * t134;
t121 = m(4) * t164 - t249 * mrSges(4,1) + qJDD(1) * mrSges(4,2) + t125;
t132 = t277 * t136 + t244 * t137;
t129 = -m(5) * t159 + t216 * mrSges(5,1) - t215 * mrSges(5,2) + t218 * t271 - t219 * t270 - t132;
t126 = m(4) * t163 - qJDD(1) * mrSges(4,1) - t249 * mrSges(4,2) + t129;
t118 = t241 * t121 - t240 * t126;
t117 = t240 * t121 + t241 * t126;
t124 = t242 * t128 + t245 * t134;
t198 = -t249 * pkin(1) + t261;
t262 = m(3) * t198 + qJDD(1) * mrSges(3,3) + t118;
t123 = -m(4) * t236 - t124;
t200 = -qJDD(1) * pkin(1) + t260;
t257 = -m(3) * t200 + qJDD(1) * mrSges(3,1) + t249 * mrSges(3,3) - t117;
t201 = Ifges(5,3) * qJD(4) + (-Ifges(5,5) * t242 - Ifges(5,6) * t245) * qJD(1);
t112 = mrSges(5,2) * t159 - mrSges(5,3) * t155 + Ifges(5,1) * t215 + Ifges(5,4) * t216 + Ifges(5,5) * qJDD(4) - pkin(8) * t132 - qJD(4) * t202 - t277 * t130 + t244 * t131 - t201 * t270;
t119 = -mrSges(5,1) * t159 + mrSges(5,3) * t156 + Ifges(5,4) * t215 + Ifges(5,2) * t216 + Ifges(5,6) * qJDD(4) - pkin(4) * t132 + qJD(4) * t203 + t201 * t271 - t281;
t110 = mrSges(4,2) * t236 - mrSges(4,3) * t163 - Ifges(4,5) * qJDD(1) - t249 * Ifges(4,6) - pkin(7) * t124 + t245 * t112 - t242 * t119;
t111 = -mrSges(4,1) * t236 + mrSges(4,3) * t164 + t249 * Ifges(4,5) - Ifges(4,6) * qJDD(1) - pkin(3) * t124 - t280;
t256 = mrSges(3,2) * t200 + mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) + t249 * Ifges(3,6) - qJ(3) * t117 + t241 * t110 - t240 * t111;
t255 = mrSges(3,2) * t198 - pkin(2) * t123 - qJ(3) * t118 - t240 * t110 - t241 * t111;
t254 = mrSges(4,1) * t163 - mrSges(4,2) * t164 - Ifges(4,3) * qJDD(1) + pkin(3) * t129 + pkin(7) * t125 + t242 * t112 + t245 * t119;
t253 = -mrSges(3,1) * t200 + mrSges(3,3) * t198 + Ifges(3,2) * qJDD(1) - pkin(2) * t117 - t254;
t250 = -mrSges(2,2) * t221 + qJ(2) * (-t249 * mrSges(3,1) + t262) + pkin(1) * t257 + mrSges(2,1) * t220 + Ifges(2,3) * qJDD(1) + t253;
t122 = -m(3) * g(3) + t123;
t114 = m(2) * t220 + qJDD(1) * mrSges(2,1) - t249 * mrSges(2,2) + t257;
t113 = m(2) * t221 - qJDD(1) * mrSges(2,2) - t276 * t249 + t262;
t108 = -mrSges(2,2) * g(3) - mrSges(2,3) * t220 + Ifges(2,5) * qJDD(1) - t249 * Ifges(2,6) - qJ(2) * t122 + t256;
t107 = mrSges(2,3) * t221 - pkin(1) * t122 + (Ifges(3,4) + Ifges(2,5)) * t249 + (Ifges(2,6) - Ifges(3,6)) * qJDD(1) + t276 * g(3) + t255;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t246 * t108 - t243 * t107 - pkin(6) * (t243 * t113 + t246 * t114), t108, t256, t110, t112, t131, t211 * t168 + t258; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t243 * t108 + t246 * t107 + pkin(6) * (t246 * t113 - t243 * t114), t107, t253, t111, t119, t130, t212 * t166 - t259; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t250, t250, -mrSges(3,1) * g(3) - t249 * Ifges(3,4) + Ifges(3,6) * qJDD(1) - t255, t254, t280, t281, Ifges(7,5) * t181 + Ifges(7,6) * t210 + Ifges(7,3) * t180 - t212 * t168 - t223 * t170 - t264;];
m_new  = t1;
