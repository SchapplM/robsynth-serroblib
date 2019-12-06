% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRRR6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR6_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR6_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR6_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR6_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:09:39
% EndTime: 2019-12-05 17:09:43
% DurationCPUTime: 2.20s
% Computational Cost: add. (5454->193), mult. (12120->255), div. (0->0), fcn. (12124->8), ass. (0->135)
t177 = cos(qJ(4));
t336 = t177 ^ 2;
t175 = sin(qJ(3));
t285 = sin(qJ(2));
t286 = cos(qJ(3));
t287 = cos(qJ(2));
t146 = t175 * t285 - t286 * t287;
t174 = sin(qJ(4));
t238 = t174 ^ 2 + t336;
t308 = t238 * t146;
t173 = sin(qJ(5));
t176 = cos(qJ(5));
t203 = t173 * t177 + t176 * t174;
t261 = t203 * mrSges(6,3);
t225 = -t261 / 0.2e1;
t289 = -pkin(8) - pkin(7);
t159 = t289 * t174;
t160 = t289 * t177;
t110 = t173 * t159 - t176 * t160;
t282 = t175 * pkin(2);
t166 = pkin(7) + t282;
t274 = pkin(8) + t166;
t141 = t274 * t174;
t142 = t274 * t177;
t94 = -t173 * t141 + t176 * t142;
t311 = t110 + t94;
t333 = t311 * t225;
t307 = qJD(2) + qJD(3);
t216 = -t176 * t141 - t173 * t142;
t332 = -t94 * mrSges(6,1) - t216 * mrSges(6,2);
t215 = t176 * t159 + t173 * t160;
t331 = -t110 * mrSges(6,1) - t215 * mrSges(6,2);
t316 = -mrSges(6,2) / 0.2e1;
t317 = mrSges(6,1) / 0.2e1;
t87 = t203 * t146;
t299 = -t173 * t174 + t176 * t177;
t89 = t299 * t146;
t273 = -t89 * t316 + t87 * t317;
t291 = m(6) / 0.2e1;
t323 = pkin(4) * t291;
t190 = (-t173 * t89 + t176 * t87) * t323 + t273;
t96 = mrSges(6,1) * t203 + mrSges(6,2) * t299;
t259 = t146 * t96;
t226 = t261 / 0.2e1;
t147 = -t175 * t287 - t286 * t285;
t86 = t299 * t147;
t300 = (t225 + t226) * t86;
t283 = pkin(4) * t174;
t314 = m(6) * t283;
t315 = t146 * t314;
t319 = -t315 / 0.2e1 - t259 / 0.2e1 + t190 + t300;
t330 = qJD(1) * t319;
t237 = t259 / 0.2e1 + t300;
t211 = t237 + t273;
t106 = t146 * t147;
t88 = t203 * t147;
t21 = m(6) * (t86 * t89 + t88 * t87 - t106) + m(5) * (t238 - 0.1e1) * t106;
t242 = t21 * qJD(1);
t253 = t177 * mrSges(5,2);
t256 = t174 * mrSges(5,1);
t207 = t253 + t256;
t318 = t315 / 0.2e1 + t190 + t237 + (t207 / 0.2e1 + t253 / 0.2e1 + t256 / 0.2e1) * t146;
t327 = t318 * qJD(4) + t211 * qJD(5) + t242;
t210 = t237 - t273;
t326 = -qJD(4) * t319 + qJD(5) * t210 - t242;
t239 = Ifges(6,5) * t299 - Ifges(6,6) * t203;
t25 = t239 + t332;
t322 = t25 * qJD(5);
t32 = t239 + t331;
t321 = t32 * qJD(5);
t219 = t86 * mrSges(6,1) - t88 * mrSges(6,2);
t320 = t219 * qJD(5);
t208 = t177 * mrSges(5,1) - t174 * mrSges(5,2);
t97 = -mrSges(6,1) * t299 + mrSges(6,2) * t203;
t306 = mrSges(4,1) + t208 - t97;
t281 = t177 * pkin(4);
t168 = -pkin(3) - t281;
t257 = t168 * t96;
t232 = t286 * pkin(2);
t167 = -t232 - pkin(3);
t151 = t167 - t281;
t258 = t151 * t96;
t305 = t258 / 0.2e1 + t257 / 0.2e1 + t333;
t304 = Ifges(5,4) * t336 + (-Ifges(5,4) * t174 + (Ifges(5,1) - Ifges(5,2)) * t177) * t174;
t218 = Ifges(6,4) * t299 ^ 2 + (-Ifges(6,4) * t203 + (-Ifges(6,2) + Ifges(6,1)) * t299) * t203;
t303 = (t151 + t168) * t283;
t298 = qJD(1) * t210;
t222 = t177 * t286;
t223 = t174 * t286;
t118 = (-t173 * t222 - t176 * t223) * pkin(2);
t119 = (-t173 * t223 + t176 * t222) * pkin(2);
t201 = t118 * t317 + t119 * t316;
t295 = (-mrSges(5,1) * t223 / 0.2e1 - mrSges(5,2) * t222 / 0.2e1) * pkin(2) + t201 + (t118 * t176 + t119 * t173) * t323;
t294 = (-t168 / 0.2e1 - t151 / 0.2e1) * t96 - t218;
t293 = t218 + t304;
t292 = m(5) / 0.2e1;
t290 = m(5) * pkin(2);
t284 = pkin(4) * t173;
t268 = pkin(4) * qJD(4);
t263 = t299 * mrSges(6,3);
t131 = t146 * t282;
t233 = -t151 * t147 + t216 * t87 - t94 * t89;
t231 = t176 * t263;
t224 = -t167 * t147 - t308 * t166;
t95 = t97 * t283;
t191 = t95 + t293;
t197 = t167 * t207;
t8 = t151 * t314 + t191 + t197 + t258;
t206 = t8 * qJD(2) - t330;
t179 = -t118 * t261 + t119 * t263 - t306 * t282 + (t238 * mrSges(5,3) - mrSges(4,2)) * t232;
t195 = t238 * t286;
t23 = m(6) * (t118 * t216 + t94 * t119 + t151 * t282) + (t195 * t166 + t167 * t175) * t290 + t179;
t178 = (-t195 * t147 * pkin(2) + t131 + t224) * t292 + (t118 * t88 - t119 * t86 + t131 + t233) * t291;
t181 = (-t110 * t89 - t168 * t147 + t215 * t87) * t291 + (pkin(3) * t147 - pkin(7) * t308) * t292;
t9 = t178 - t181;
t205 = t9 * qJD(1) + t23 * qJD(2);
t17 = t218 + t258;
t204 = t17 * qJD(2) + t298;
t200 = pkin(3) * t207;
t196 = Ifges(5,5) * t177 - Ifges(5,6) * t174 - t261 * t284 + t239;
t12 = t168 * t314 + t191 - t200 + t257;
t183 = -t95 + t200 / 0.2e1 - t197 / 0.2e1;
t2 = t183 - t303 * m(6) / 0.2e1 + t294 + t295 - t304;
t194 = -t2 * qJD(2) + t12 * qJD(3) - t330;
t18 = t218 + t257;
t4 = t201 + t294;
t193 = -t4 * qJD(2) + t18 * qJD(3) + t298;
t188 = t146 * mrSges(4,2) - t308 * mrSges(5,3) + t306 * t147 - t87 * t261 - t89 * t263;
t150 = (mrSges(6,1) * t173 + mrSges(6,2) * t176) * pkin(4);
t187 = t150 * qJD(4);
t143 = t150 * qJD(5);
t5 = t311 * t226 + t201 + t218 + t305;
t3 = t178 + t181 + t188;
t1 = t303 * t291 - t183 + t293 + t295 + t305 - t333;
t6 = [t307 * t21, t3 * qJD(3) + (-t285 * mrSges(3,1) - t287 * mrSges(3,2) + t188 + 0.2e1 * t233 * t291 + 0.2e1 * t224 * t292 + m(4) * (t147 * t232 - t131)) * qJD(2) + t327, t3 * qJD(2) + (0.2e1 * t181 + t188) * qJD(3) + t327, (m(6) * (t176 * pkin(4) * t86 + t88 * t284) + t208 * t147 + t219) * qJD(4) + t320 + t307 * t318, qJD(4) * t219 + t307 * t211 + t320; t9 * qJD(3) + t326, t23 * qJD(3) + t8 * qJD(4) + t17 * qJD(5), ((-pkin(3) * t175 + t195 * pkin(7)) * t290 + m(6) * (t110 * t119 + t118 * t215 + t168 * t282) + t179) * qJD(3) + t1 * qJD(4) + t5 * qJD(5) + t205, t1 * qJD(3) + (-t208 * t166 + t196 + t332) * qJD(4) + t322 + (-t231 + m(6) * (t173 * t216 - t176 * t94)) * t268 + t206, t5 * qJD(3) + t25 * qJD(4) + t204 + t322; -t9 * qJD(2) + t326, -t2 * qJD(4) - t4 * qJD(5) - t205, t12 * qJD(4) + t18 * qJD(5), (-t208 * pkin(7) + t196 + t331) * qJD(4) + t321 + (-t231 + m(6) * (-t110 * t176 + t173 * t215)) * t268 + t194, t32 * qJD(4) + t193 + t321; t307 * t319, t2 * qJD(3) - t206, -t194, -t143, -t143 - t187; -t307 * t210, t4 * qJD(3) - t204, -t193, t187, 0;];
Cq = t6;
