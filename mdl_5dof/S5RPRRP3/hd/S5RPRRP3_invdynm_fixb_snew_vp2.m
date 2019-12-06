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
% Datum: 2019-12-05 18:04
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:03:16
% EndTime: 2019-12-05 18:03:23
% DurationCPUTime: 3.49s
% Computational Cost: add. (40477->262), mult. (80564->323), div. (0->0), fcn. (47650->8), ass. (0->99)
t242 = sin(qJ(4));
t243 = sin(qJ(3));
t245 = cos(qJ(4));
t246 = cos(qJ(3));
t208 = (-t242 * t243 + t245 * t246) * qJD(1);
t268 = qJD(1) * qJD(3);
t265 = t246 * t268;
t216 = t243 * qJDD(1) + t265;
t217 = t246 * qJDD(1) - t243 * t268;
t176 = t208 * qJD(4) + t245 * t216 + t242 * t217;
t209 = (t242 * t246 + t243 * t245) * qJD(1);
t190 = -t208 * mrSges(6,1) + t209 * mrSges(6,2);
t244 = sin(qJ(1));
t247 = cos(qJ(1));
t223 = t247 * g(2) + t244 * g(3);
t213 = qJDD(1) * pkin(1) + t223;
t222 = t244 * g(2) - t247 * g(3);
t248 = qJD(1) ^ 2;
t215 = -t248 * pkin(1) + t222;
t240 = sin(pkin(8));
t241 = cos(pkin(8));
t194 = t240 * t213 + t241 * t215;
t189 = -t248 * pkin(2) + qJDD(1) * pkin(6) + t194;
t239 = -g(1) + qJDD(2);
t173 = -t243 * t189 + t246 * t239;
t154 = (-t216 + t265) * pkin(7) + (t243 * t246 * t248 + qJDD(3)) * pkin(3) + t173;
t174 = t246 * t189 + t243 * t239;
t270 = qJD(1) * t243;
t221 = qJD(3) * pkin(3) - pkin(7) * t270;
t238 = t246 ^ 2;
t155 = -t238 * t248 * pkin(3) + t217 * pkin(7) - qJD(3) * t221 + t174;
t149 = t245 * t154 - t242 * t155;
t234 = qJDD(3) + qJDD(4);
t235 = qJD(3) + qJD(4);
t141 = -0.2e1 * qJD(5) * t209 + (t208 * t235 - t176) * qJ(5) + (t208 * t209 + t234) * pkin(4) + t149;
t196 = -t235 * mrSges(6,2) + t208 * mrSges(6,3);
t267 = m(6) * t141 + t234 * mrSges(6,1) + t235 * t196;
t138 = -t176 * mrSges(6,3) - t209 * t190 + t267;
t150 = t242 * t154 + t245 * t155;
t175 = -t209 * qJD(4) - t242 * t216 + t245 * t217;
t182 = Ifges(5,4) * t209 + Ifges(5,2) * t208 + Ifges(5,6) * t235;
t183 = Ifges(6,1) * t209 + Ifges(6,4) * t208 + Ifges(6,5) * t235;
t184 = Ifges(5,1) * t209 + Ifges(5,4) * t208 + Ifges(5,5) * t235;
t198 = t235 * pkin(4) - t209 * qJ(5);
t201 = t208 ^ 2;
t144 = -t201 * pkin(4) + t175 * qJ(5) + 0.2e1 * qJD(5) * t208 - t235 * t198 + t150;
t181 = Ifges(6,4) * t209 + Ifges(6,2) * t208 + Ifges(6,6) * t235;
t256 = -mrSges(6,1) * t141 + mrSges(6,2) * t144 - Ifges(6,5) * t176 - Ifges(6,6) * t175 - Ifges(6,3) * t234 - t209 * t181;
t274 = mrSges(5,1) * t149 - mrSges(5,2) * t150 + Ifges(5,5) * t176 + Ifges(5,6) * t175 + Ifges(5,3) * t234 + pkin(4) * t138 + t209 * t182 - t256 + (-t184 - t183) * t208;
t191 = -t208 * mrSges(5,1) + t209 * mrSges(5,2);
t197 = -t235 * mrSges(5,2) + t208 * mrSges(5,3);
t132 = m(5) * t149 + t234 * mrSges(5,1) + t235 * t197 + (-t190 - t191) * t209 + (-mrSges(5,3) - mrSges(6,3)) * t176 + t267;
t199 = t235 * mrSges(6,1) - t209 * mrSges(6,3);
t200 = t235 * mrSges(5,1) - t209 * mrSges(5,3);
t266 = m(6) * t144 + t175 * mrSges(6,3) + t208 * t190;
t135 = m(5) * t150 + t175 * mrSges(5,3) + t208 * t191 + (-t199 - t200) * t235 + (-mrSges(5,2) - mrSges(6,2)) * t234 + t266;
t128 = t245 * t132 + t242 * t135;
t206 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t243 + Ifges(4,2) * t246) * qJD(1);
t207 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t243 + Ifges(4,4) * t246) * qJD(1);
t273 = mrSges(4,1) * t173 - mrSges(4,2) * t174 + Ifges(4,5) * t216 + Ifges(4,6) * t217 + Ifges(4,3) * qJDD(3) + pkin(3) * t128 + (t243 * t206 - t246 * t207) * qJD(1) + t274;
t214 = (-mrSges(4,1) * t246 + mrSges(4,2) * t243) * qJD(1);
t269 = qJD(1) * t246;
t220 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t269;
t126 = m(4) * t173 + qJDD(3) * mrSges(4,1) - t216 * mrSges(4,3) + qJD(3) * t220 - t214 * t270 + t128;
t219 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t270;
t261 = -t242 * t132 + t245 * t135;
t127 = m(4) * t174 - qJDD(3) * mrSges(4,2) + t217 * mrSges(4,3) - qJD(3) * t219 + t214 * t269 + t261;
t262 = -t243 * t126 + t246 * t127;
t118 = m(3) * t194 - t248 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t262;
t193 = t241 * t213 - t240 * t215;
t258 = -qJDD(1) * pkin(2) - t193;
t188 = -t248 * pkin(6) + t258;
t156 = -t217 * pkin(3) + t221 * t270 + (-pkin(7) * t238 - pkin(6)) * t248 + t258;
t147 = -t175 * pkin(4) - t201 * qJ(5) + t209 * t198 + qJDD(5) + t156;
t260 = m(6) * t147 - t175 * mrSges(6,1) + t176 * mrSges(6,2) - t208 * t196 + t209 * t199;
t252 = m(5) * t156 - t175 * mrSges(5,1) + t176 * mrSges(5,2) - t208 * t197 + t209 * t200 + t260;
t250 = -m(4) * t188 + t217 * mrSges(4,1) - t216 * mrSges(4,2) - t219 * t270 + t220 * t269 - t252;
t130 = m(3) * t193 + qJDD(1) * mrSges(3,1) - t248 * mrSges(3,2) + t250;
t115 = t240 * t118 + t241 * t130;
t120 = t246 * t126 + t243 * t127;
t263 = t241 * t118 - t240 * t130;
t257 = -mrSges(6,1) * t147 + mrSges(6,3) * t144 + Ifges(6,4) * t176 + Ifges(6,2) * t175 + Ifges(6,6) * t234 + t235 * t183;
t179 = Ifges(6,5) * t209 + Ifges(6,6) * t208 + Ifges(6,3) * t235;
t255 = mrSges(6,2) * t147 - mrSges(6,3) * t141 + Ifges(6,1) * t176 + Ifges(6,4) * t175 + Ifges(6,5) * t234 + t208 * t179;
t180 = Ifges(5,5) * t209 + Ifges(5,6) * t208 + Ifges(5,3) * t235;
t121 = Ifges(5,4) * t176 + Ifges(5,2) * t175 + Ifges(5,6) * t234 + t235 * t184 - mrSges(5,1) * t156 + mrSges(5,3) * t150 - pkin(4) * t260 + qJ(5) * (-t234 * mrSges(6,2) - t235 * t199 + t266) + (-t180 - t179) * t209 + t257;
t122 = mrSges(5,2) * t156 - mrSges(5,3) * t149 + Ifges(5,1) * t176 + Ifges(5,4) * t175 + Ifges(5,5) * t234 - qJ(5) * t138 + t208 * t180 + (-t181 - t182) * t235 + t255;
t205 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t243 + Ifges(4,6) * t246) * qJD(1);
t109 = -mrSges(4,1) * t188 + mrSges(4,3) * t174 + Ifges(4,4) * t216 + Ifges(4,2) * t217 + Ifges(4,6) * qJDD(3) - pkin(3) * t252 + pkin(7) * t261 + qJD(3) * t207 + t245 * t121 + t242 * t122 - t205 * t270;
t111 = mrSges(4,2) * t188 - mrSges(4,3) * t173 + Ifges(4,1) * t216 + Ifges(4,4) * t217 + Ifges(4,5) * qJDD(3) - pkin(7) * t128 - qJD(3) * t206 - t242 * t121 + t245 * t122 + t205 * t269;
t254 = mrSges(3,1) * t193 - mrSges(3,2) * t194 + Ifges(3,3) * qJDD(1) + pkin(2) * t250 + pkin(6) * t262 + t246 * t109 + t243 * t111;
t253 = mrSges(2,1) * t223 - mrSges(2,2) * t222 + Ifges(2,3) * qJDD(1) + pkin(1) * t115 + t254;
t113 = m(2) * t223 + qJDD(1) * mrSges(2,1) - t248 * mrSges(2,2) + t115;
t112 = m(2) * t222 - t248 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t263;
t107 = -mrSges(3,1) * t239 + mrSges(3,3) * t194 + t248 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t120 - t273;
t106 = mrSges(3,2) * t239 - mrSges(3,3) * t193 + Ifges(3,5) * qJDD(1) - t248 * Ifges(3,6) - pkin(6) * t120 - t243 * t109 + t246 * t111;
t105 = -mrSges(2,2) * g(1) - mrSges(2,3) * t223 + Ifges(2,5) * qJDD(1) - t248 * Ifges(2,6) - qJ(2) * t115 + t241 * t106 - t240 * t107;
t104 = Ifges(2,6) * qJDD(1) + t248 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t222 + t240 * t106 + t241 * t107 - pkin(1) * (m(3) * t239 + t120) + qJ(2) * t263;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t253, t105, t106, t111, t122, -t235 * t181 + t255; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - t244 * t105 - t247 * t104 - pkin(5) * (t247 * t112 - t244 * t113), t104, t107, t109, t121, -t209 * t179 + t257; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t247 * t105 - t244 * t104 + pkin(5) * (-t244 * t112 - t247 * t113), t253, t254, t273, t274, -t208 * t183 - t256;];
m_new = t1;
