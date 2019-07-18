% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
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
% Datum: 2019-07-18 17:22
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(4,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR1_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_coriolismatJ_fixb_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR1_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:20:47
% EndTime: 2019-07-18 17:20:56
% DurationCPUTime: 2.74s
% Computational Cost: add. (5242->283), mult. (11043->380), div. (0->0), fcn. (10577->6), ass. (0->183)
t178 = sin(qJ(2));
t171 = t178 * pkin(1);
t159 = t178 * pkin(2) + t171;
t339 = m(5) * t159;
t180 = cos(qJ(2));
t305 = sin(qJ(4));
t306 = cos(qJ(4));
t150 = t178 * t305 - t180 * t306;
t177 = sin(qJ(5));
t173 = t177 ^ 2;
t179 = cos(qJ(5));
t175 = t179 ^ 2;
t322 = t173 + t175;
t225 = t322 * t150;
t181 = pkin(2) + pkin(1);
t234 = t305 * t181;
t161 = t234 + pkin(4);
t152 = -t178 * t306 - t180 * t305;
t236 = t306 * t152;
t155 = -mrSges(6,1) * t179 + mrSges(6,2) * t177;
t267 = t152 * t155;
t241 = -t267 / 0.2e1 - mrSges(6,3) * t225 / 0.2e1;
t304 = m(5) * t181;
t317 = m(6) / 0.2e1;
t338 = (-t161 * t225 + t181 * t236) * t317 + (-t150 * t305 + t236) * t304 / 0.2e1 + t241;
t218 = Ifges(6,5) * t177 + Ifges(6,6) * t179;
t205 = t152 * t218;
t307 = t179 / 0.2e1;
t310 = t177 / 0.2e1;
t167 = Ifges(6,4) * t179;
t158 = Ifges(6,1) * t177 + t167;
t258 = t179 * t158;
t297 = Ifges(6,4) * t177;
t157 = Ifges(6,2) * t179 + t297;
t262 = t177 * t157;
t327 = -t258 / 0.2e1 + t262 / 0.2e1;
t219 = Ifges(6,2) * t177 - t167;
t65 = -Ifges(6,6) * t152 + t150 * t219;
t220 = Ifges(6,1) * t179 - t297;
t67 = -Ifges(6,5) * t152 - t150 * t220;
t203 = t65 * t307 + t67 * t310 - t205 / 0.2e1 + Ifges(5,6) * t152 + (-Ifges(5,5) + t327) * t150;
t229 = (-pkin(3) - qJ(3)) * t178;
t165 = t180 * qJ(3);
t257 = t180 * pkin(3) + t165;
t319 = t305 * t229 + t257 * t306;
t328 = t319 * t155;
t330 = t319 * mrSges(5,1);
t117 = -t306 * t229 + t257 * t305;
t336 = t117 * mrSges(5,2);
t337 = t203 + t328 - t330 + t336;
t335 = t117 * t177;
t334 = t117 * t179;
t333 = t117 * t305;
t274 = t319 * t117;
t331 = t328 / 0.2e1 - t330 / 0.2e1;
t311 = -t177 / 0.2e1;
t166 = Ifges(6,5) * t179;
t295 = Ifges(6,6) * t177;
t209 = -t166 / 0.2e1 + t295 / 0.2e1;
t329 = Ifges(5,4) + t209;
t237 = t306 * t319;
t125 = pkin(4) * t150 + t159;
t51 = t125 * t177 - t334;
t284 = t179 * t51;
t49 = t125 * t179 + t335;
t287 = t177 * t49;
t217 = t284 - t287;
t269 = t150 * t177;
t101 = mrSges(6,2) * t152 + mrSges(6,3) * t269;
t268 = t150 * t179;
t103 = -mrSges(6,1) * t152 + mrSges(6,3) * t268;
t302 = t101 * t310 + t103 * t307;
t69 = pkin(4) * t268 + t335;
t70 = pkin(4) * t269 - t334;
t326 = t302 + (t177 * t70 + t179 * t69) * t317;
t285 = t179 * mrSges(6,2);
t288 = t177 * mrSges(6,1);
t156 = t285 + t288;
t100 = t156 * t152;
t160 = t181 * t180;
t228 = t152 * mrSges(5,1) + t150 * mrSges(5,2);
t126 = pkin(4) * t152 - t160;
t50 = t126 * t179 - t177 * t319;
t52 = t126 * t177 + t179 * t319;
t99 = t156 * t150;
t325 = -t319 * t100 + t52 * t101 + t50 * t103 - t117 * t99 + t160 * t228;
t324 = t178 ^ 2;
t323 = Ifges(3,4) + Ifges(4,4);
t321 = 0.3e1 / 0.4e1 * t297 + t157 / 0.4e1;
t320 = t219 * t307 + t220 * t311;
t261 = t179 * t101;
t264 = t177 * t103;
t318 = t336 / 0.2e1 + (t284 / 0.2e1 - t287 / 0.2e1) * mrSges(6,3) + (t217 * t317 + t261 / 0.2e1 - t264 / 0.2e1) * pkin(4);
t316 = m(4) * pkin(1);
t315 = -mrSges(6,1) / 0.2e1;
t314 = mrSges(6,2) / 0.2e1;
t313 = Ifges(6,3) / 0.2e1;
t308 = -t179 / 0.2e1;
t303 = m(6) * t181;
t301 = mrSges(5,3) * t150;
t300 = mrSges(5,3) * t152;
t299 = mrSges(6,3) * t152;
t253 = t177 * t299;
t102 = -t150 * mrSges(6,2) + t253;
t104 = t150 * mrSges(6,1) + t179 * t299;
t188 = t67 * t308 + t65 * t310 - t329 * t152 + (-Ifges(6,3) + Ifges(5,1) - Ifges(5,2)) * t150;
t68 = Ifges(6,5) * t150 - t152 * t220;
t282 = t179 * t68;
t66 = Ifges(6,6) * t150 + t152 * t219;
t286 = t177 * t66;
t191 = -t282 / 0.2e1 + t286 / 0.2e1 + t329 * t150;
t281 = t180 * mrSges(4,2);
t1 = t51 * t102 + t49 * t104 - t160 * t339 + m(6) * (t49 * t50 + t51 * t52 + t274) + (t159 * mrSges(5,1) + t191) * t150 + (-t159 * mrSges(5,2) + t188) * t152 + (pkin(1) * mrSges(4,2) - t323) * t324 + (-pkin(1) * t281 + (Ifges(3,1) + Ifges(4,1) - Ifges(3,2) - Ifges(4,2) + (-(2 * mrSges(4,1)) - t316) * pkin(1)) * t178 + t323 * t180) * t180 + t325;
t294 = t1 * qJD(1);
t283 = t179 * t52;
t2 = m(6) * (t50 * t69 + t52 * t70 + t274) + t70 * t102 + t69 * t104 + t191 * t150 + t188 * t152 + t325;
t280 = t2 * qJD(1);
t279 = t69 * t177;
t143 = -t262 / 0.2e1;
t7 = -t117 * t267 + t52 * t104 + (t68 * t311 + t66 * t308 - mrSges(6,3) * t283 - t150 * t218 / 0.2e1 + (t158 * t307 + t143) * t152) * t152 + (-t102 + t253) * t50;
t278 = t7 * qJD(1);
t277 = t70 * t179;
t260 = t179 * t102;
t263 = t177 * t104;
t272 = t117 * t152;
t10 = (t100 + t300) * t152 + (-t260 + t263 + t301) * t150 + m(6) * (-t272 + (t177 * t50 - t283) * t150) + m(5) * (-t150 * t319 - t272) + (m(4) * qJ(3) + mrSges(4,3)) * (t180 ^ 2 + t324);
t276 = qJD(1) * t10;
t192 = t339 / 0.2e1 + (t177 * t51 + t179 * t49) * t317 + t302;
t11 = m(4) * t171 + t178 * mrSges(4,1) + t192 - t228 + t281 - t338;
t275 = t11 * qJD(1);
t195 = -pkin(4) * t225 * t317 + t241;
t14 = t195 + t228 - t326;
t270 = t14 * qJD(1);
t199 = (t285 / 0.2e1 + t288 / 0.2e1) * t150;
t206 = t263 / 0.2e1 - t260 / 0.2e1;
t17 = t199 + t206;
t266 = t17 * qJD(1);
t255 = mrSges(6,3) * t279;
t254 = mrSges(6,3) * t277;
t252 = -Ifges(6,2) / 0.4e1 + Ifges(6,1) / 0.4e1;
t251 = t161 * t264;
t250 = t161 * t261;
t240 = t177 * t306;
t239 = t179 * t306;
t238 = t181 * t306;
t235 = t306 * t156;
t226 = t166 - t295;
t53 = t258 / 0.2e1 + t143 - t320;
t223 = -t238 / 0.2e1;
t221 = mrSges(6,3) * (t173 / 0.2e1 + t175 / 0.2e1);
t216 = t277 - t279;
t185 = (t155 - mrSges(5,1)) * t234 + (mrSges(6,3) * t322 - mrSges(5,2)) * t238;
t201 = t322 * t306;
t25 = (t161 * t201 - t234 * t306) * t303 + t185;
t182 = -m(6) * (t216 * t161 + (t239 * t52 - t240 * t50 - t237 + t333) * t181) / 0.2e1 - t336 / 0.2e1 + t251 / 0.2e1 - t250 / 0.2e1 + t255 / 0.2e1 - t254 / 0.2e1 + t100 * t234 / 0.2e1 + t223 * t260 + (-t99 + t263) * t238 / 0.2e1 - t331;
t4 = t182 + t318 + t331;
t214 = -t4 * qJD(1) + t25 * qJD(2);
t187 = t320 + t327;
t32 = t181 * t235 + t187;
t183 = t155 * t223 + Ifges(6,2) * t175 / 0.4e1 + (-t175 / 0.4e1 + t173 / 0.4e1) * Ifges(6,1) + t161 * t221 + (t158 - t219) * t177 / 0.4e1 + t321 * t179;
t196 = t117 * t156 / 0.2e1 - t286 / 0.4e1 + t282 / 0.4e1;
t189 = (-0.3e1 / 0.4e1 * t295 + 0.3e1 / 0.4e1 * t166) * t150 + t196;
t207 = t102 * t311 + t104 * t308;
t197 = t207 * t161;
t211 = t314 * t51 + t315 * t49;
t6 = t197 + (t313 + t183) * t152 + t189 + t211;
t213 = t6 * qJD(1) - t32 * qJD(2);
t210 = t314 * t70 + t315 * t69;
t200 = t207 * pkin(4);
t193 = -mrSges(6,2) * t239 / 0.2e1 + t240 * t315;
t26 = (t235 / 0.2e1 + t193) * t181 + t187;
t184 = pkin(4) * t221 + (t158 / 0.4e1 + t167 / 0.4e1 + t252 * t177) * t177 + (-t179 * t252 + t321) * t179;
t9 = t200 + (t313 + t184) * t152 + t189 + t210;
t198 = t9 * qJD(1) - t26 * qJD(2) + t53 * qJD(4);
t190 = -Ifges(6,3) * t152 / 0.2e1 + t196 + (t226 / 0.4e1 + t209) * t150;
t27 = t156 * t223 + t181 * t193 + t53;
t18 = t199 - t206;
t15 = t195 + t326;
t12 = t192 + t338;
t8 = t152 * t184 + t190 + t200 - t210;
t5 = t152 * t183 + t190 + t197 - t211;
t3 = (t155 / 0.2e1 - mrSges(5,1) / 0.2e1) * t319 - t182 + t203 + t318;
t13 = [qJD(2) * t1 + qJD(3) * t10 + qJD(4) * t2 - qJD(5) * t7, t12 * qJD(3) + t3 * qJD(4) + t5 * qJD(5) + t294 + ((-t237 - t333) * t304 + m(6) * (t161 * t217 - t181 * t237) + t234 * t300 + t250 - t251 + (-mrSges(4,3) * pkin(1) + Ifges(3,5) + Ifges(4,5)) * t180 + (mrSges(4,2) * qJ(3) - Ifges(3,6) - Ifges(4,6)) * t178 + (t301 + t99) * t238 + (-t316 - mrSges(4,1)) * t165 + t217 * mrSges(6,3) + t337) * qJD(2), qJD(2) * t12 + qJD(4) * t15 + qJD(5) * t18 + t276, t3 * qJD(2) + t15 * qJD(3) + t8 * qJD(5) + t280 + (t254 - t255 + (m(6) * t216 + t261 - t264) * pkin(4) + t337) * qJD(4), -t278 + t5 * qJD(2) + t18 * qJD(3) + t8 * qJD(4) + (-mrSges(6,1) * t52 - mrSges(6,2) * t50 + t205) * qJD(5); -qJD(3) * t11 - qJD(4) * t4 + qJD(5) * t6 - t294, qJD(4) * t25 - qJD(5) * t32, -t275, (pkin(4) * t201 * t303 + t185) * qJD(4) + t27 * qJD(5) + t214, t27 * qJD(4) + (t155 * t161 + t226) * qJD(5) + t213; qJD(2) * t11 - qJD(4) * t14 - qJD(5) * t17 - t276, t275, 0, -t270, -qJD(5) * t156 - t266; qJD(2) * t4 + qJD(3) * t14 + qJD(5) * t9 - t280, -qJD(5) * t26 - t214, t270, t53 * qJD(5), (pkin(4) * t155 + t226) * qJD(5) + t198; -qJD(2) * t6 + qJD(3) * t17 - qJD(4) * t9 + t278, qJD(4) * t26 - t213, t266, -t198, 0;];
Cq  = t13;
