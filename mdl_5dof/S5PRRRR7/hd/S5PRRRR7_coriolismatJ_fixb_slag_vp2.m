% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRRR7
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
% Datum: 2019-12-05 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRRR7_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR7_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR7_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR7_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR7_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:11:53
% EndTime: 2019-12-05 17:12:01
% DurationCPUTime: 3.95s
% Computational Cost: add. (8117->235), mult. (18354->335), div. (0->0), fcn. (19714->8), ass. (0->143)
t326 = qJD(3) + qJD(4);
t189 = sin(qJ(4));
t190 = sin(qJ(3));
t193 = cos(qJ(4));
t194 = cos(qJ(3));
t239 = t193 * t194;
t168 = -t189 * t190 + t239;
t169 = -t189 * t194 - t193 * t190;
t188 = sin(qJ(5));
t192 = cos(qJ(5));
t128 = t168 * t188 - t169 * t192;
t227 = t192 * t168 + t169 * t188;
t334 = Ifges(6,5) * t227 - Ifges(6,6) * t128;
t298 = -pkin(7) - pkin(6);
t174 = t298 * t190;
t175 = t298 * t194;
t138 = t174 * t189 - t193 * t175;
t110 = pkin(8) * t168 + t138;
t312 = t193 * t174 + t189 * t175;
t327 = t169 * pkin(8) + t312;
t346 = -t110 * t188 + t192 * t327;
t352 = t346 * mrSges(6,2);
t62 = t110 * t192 + t188 * t327;
t356 = t62 * mrSges(6,1);
t360 = -t356 / 0.2e1 - t352 / 0.2e1;
t363 = 0.2e1 * t360 + t334;
t364 = t363 * qJD(5);
t129 = -mrSges(5,1) * t169 + mrSges(5,2) * t168;
t184 = -pkin(3) * t194 - pkin(2);
t143 = -pkin(4) * t168 + t184;
t275 = Ifges(6,4) * t128;
t319 = t227 / 0.2e1;
t323 = mrSges(6,1) * t128 + mrSges(6,2) * t227;
t331 = -t128 / 0.2e1;
t9 = t143 * t323 + (Ifges(6,1) * t227 - t275) * t128 / 0.2e1 + (Ifges(6,2) * t227 + t275) * t331 + (0.2e1 * Ifges(6,4) * t227 + (Ifges(6,1) - Ifges(6,2)) * t128) * t319;
t362 = -Ifges(5,4) * t169 ^ 2 + (Ifges(5,4) * t168 + (-Ifges(5,1) + Ifges(5,2)) * t169) * t168 + t184 * t129 + t9;
t343 = -t138 * mrSges(5,1) - t312 * mrSges(5,2) + Ifges(5,5) * t168 + Ifges(5,6) * t169 + t334;
t347 = -t352 - t356;
t361 = t343 + t347;
t288 = pkin(4) * t169;
t290 = pkin(3) * t190;
t149 = -t288 + t290;
t357 = m(6) * t149;
t355 = t188 * t346 - t192 * t62;
t302 = m(6) * pkin(4);
t349 = -t302 / 0.2e1;
t195 = cos(qJ(2));
t152 = t169 * t195;
t153 = t168 * t195;
t102 = t152 * t192 - t153 * t188;
t105 = t152 * t188 + t153 * t192;
t276 = t102 * mrSges(6,1) / 0.2e1 - t105 * mrSges(6,2) / 0.2e1;
t226 = t276 - t153 * mrSges(5,2) / 0.2e1 + t152 * mrSges(5,1) / 0.2e1;
t342 = (t102 * t192 + t105 * t188) * t349 - t226;
t191 = sin(qJ(2));
t242 = t190 * t191;
t150 = t189 * t242 - t191 * t239;
t151 = t169 * t191;
t100 = t150 * t192 - t151 * t188;
t329 = t100 * mrSges(6,1);
t299 = t329 / 0.2e1;
t341 = m(5) * t290;
t228 = t150 * t188 + t192 * t151;
t314 = t228 * mrSges(6,2);
t333 = t329 - t314;
t325 = t349 * t355 - t360;
t318 = pkin(4) * t192;
t256 = t194 * mrSges(4,1);
t311 = t190 * mrSges(4,2) - t256;
t205 = t299 - t329 / 0.2e1;
t310 = qJD(1) * t205;
t243 = t189 * t192;
t160 = (-t188 * t193 - t243) * pkin(3);
t158 = t160 * mrSges(6,1);
t244 = t188 * t189;
t161 = (t192 * t193 - t244) * pkin(3);
t260 = t161 * mrSges(6,2);
t307 = (mrSges(5,1) * t189 + mrSges(5,2) * t193) * pkin(3) - t158 + t260;
t306 = qJD(5) * t205;
t214 = 0.2e1 * t299 - t314;
t305 = t214 * qJD(5);
t304 = m(5) / 0.2e1;
t303 = m(6) / 0.2e1;
t292 = t161 / 0.2e1;
t291 = -t195 / 0.2e1;
t289 = pkin(3) * t193;
t287 = pkin(4) * t188;
t183 = pkin(4) + t289;
t156 = -pkin(3) * t244 + t183 * t192;
t157 = pkin(3) * t243 + t183 * t188;
t277 = t156 * t100 + t157 * t228;
t258 = t190 * mrSges(4,1);
t255 = t194 * mrSges(4,2);
t43 = -t157 * mrSges(6,1) - t156 * mrSges(6,2);
t253 = qJD(5) * t43;
t248 = t156 * t227;
t247 = t157 * t128;
t245 = t188 * t128;
t241 = t191 * t195;
t240 = t192 * t227;
t234 = t190 ^ 2 + t194 ^ 2;
t25 = m(6) * (-t100 * t105 + t102 * t228 - t241) + m(5) * (-t150 * t153 + t151 * t152 - t241) + m(4) * (-0.1e1 + t234) * t241;
t238 = t25 * qJD(1);
t232 = mrSges(6,3) * t331;
t231 = -t227 * mrSges(6,3) / 0.2e1;
t230 = t323 * t291;
t224 = -t318 / 0.2e1 - t156 / 0.2e1;
t130 = -mrSges(5,1) * t168 - mrSges(5,2) * t169;
t73 = -mrSges(6,1) * t227 + mrSges(6,2) * t128;
t1 = t149 * t73 + (-pkin(2) * mrSges(4,1) - Ifges(4,4) * t190 + pkin(3) * t130) * t190 + t184 * t341 + t143 * t357 + (-pkin(2) * mrSges(4,2) + Ifges(4,4) * t194 + (Ifges(4,1) - Ifges(4,2)) * t190) * t194 + t362;
t196 = (t341 + t357) * t195 / 0.2e1;
t199 = (t156 * t102 + t157 * t105) * t303 + (t152 * t193 + t153 * t189) * pkin(3) * t304 + t226;
t3 = (t129 / 0.2e1 + t323 / 0.2e1) * t195 + t196 + t199;
t221 = -t3 * qJD(1) + t1 * qJD(2);
t2 = (-m(6) * t143 - t73) * t288 + t362;
t209 = -t100 * t232 + t228 * t231 + (t129 + t323) * t291 + (t100 * t331 + t228 * t319) * mrSges(6,3);
t198 = t195 * t288 * t303 + t209;
t8 = t198 + t342;
t220 = t8 * qJD(1) + t2 * qJD(2);
t13 = t230 - t276;
t218 = t13 * qJD(1) + t9 * qJD(2);
t217 = (-t287 / 0.2e1 - t157 / 0.2e1) * mrSges(6,1);
t213 = t150 * mrSges(5,1) - t151 * mrSges(5,2) + t333;
t211 = (t100 * t192 + t188 * t228) * t302;
t200 = (-t100 * t161 + t160 * t228 + t277) * t303;
t24 = -t211 / 0.2e1 + t200;
t42 = -m(6) * (t156 * t160 + t157 * t161) + t307;
t201 = ((t157 + t160) * t346 + (-t156 + t161) * t62) * t303 + t360;
t203 = t227 * t292 + t160 * t331 - t247 / 0.2e1 - t248 / 0.2e1;
t6 = ((t245 / 0.2e1 + t240 / 0.2e1) * pkin(4) + t203) * mrSges(6,3) + t201 + t325;
t208 = t24 * qJD(1) + t6 * qJD(2) - t42 * qJD(3);
t207 = t43 * qJD(3) - t310;
t173 = (t188 * mrSges(6,1) + t192 * mrSges(6,2)) * pkin(4);
t38 = -t158 / 0.2e1 + t217 + (t292 + t224) * mrSges(6,2);
t204 = -qJD(3) * t38 + qJD(4) * t173 + t310;
t167 = t173 * qJD(5);
t39 = -t260 / 0.2e1 + t158 / 0.2e1 + t224 * mrSges(6,2) + t217;
t19 = t211 / 0.2e1 + t200 + t213;
t12 = t230 + t276;
t7 = t198 - t342;
t5 = t203 * mrSges(6,3) + t231 * t318 + t232 * t287 + t201 - t325 + t343;
t4 = (t255 + t258) * t291 + (-t255 / 0.2e1 - t258 / 0.2e1) * t195 - t196 + t199 + t209;
t10 = [t25 * qJD(2), t4 * qJD(3) + t7 * qJD(4) + t12 * qJD(5) + t238 + ((t234 * mrSges(4,3) - mrSges(3,2)) * t195 + (-mrSges(3,1) + t130 + t73 + t311) * t191 + 0.2e1 * (t102 * t346 + t105 * t62 + t143 * t191) * t303 + 0.2e1 * (t138 * t153 + t152 * t312 + t184 * t191) * t304 + m(4) * (t234 * t195 * pkin(6) - t191 * pkin(2)) + (-t102 * t128 + t105 * t227) * mrSges(6,3) + (t152 * t169 + t153 * t168) * mrSges(5,3)) * qJD(2), t4 * qJD(2) + (mrSges(4,2) * t242 - t191 * t256 + t213) * qJD(3) + t19 * qJD(4) + 0.2e1 * (t277 * t303 + (pkin(3) * t151 * t189 + t150 * t289) * t304) * qJD(3) + t305, t7 * qJD(2) + t19 * qJD(3) + (t211 + t213) * qJD(4) + t305, t12 * qJD(2) + t333 * qJD(5) + t326 * t214; -qJD(3) * t3 + qJD(4) * t8 + qJD(5) * t13 - t238, qJD(3) * t1 + qJD(4) * t2 + qJD(5) * t9, t5 * qJD(4) + t364 + t221 + (m(6) * (-t156 * t62 + t157 * t346) + Ifges(4,5) * t194 - Ifges(4,6) * t190 + (m(5) * (-t138 * t193 + t189 * t312) + (-t168 * t193 + t169 * t189) * mrSges(5,3)) * pkin(3) + t311 * pkin(6) + (-t247 - t248) * mrSges(6,3) + t361) * qJD(3), t5 * qJD(3) + t364 + t220 + ((m(6) * t355 + (-t240 - t245) * mrSges(6,3)) * pkin(4) + t361) * qJD(4), (t334 + t347) * qJD(5) + t218 + t326 * t363; qJD(2) * t3 + qJD(4) * t24 - t306, qJD(4) * t6 - t221, -qJD(4) * t42 + t253, ((t160 * t192 + t161 * t188) * t302 - t307) * qJD(4) + t39 * qJD(5) + t208, t39 * qJD(4) + t207 + t253; -qJD(2) * t8 - qJD(3) * t24 - t306, -qJD(3) * t6 - t220, qJD(5) * t38 - t208, -t167, -t167 - t204; -t13 * qJD(2) + t205 * t326, -t218, -qJD(4) * t38 - t207, t204, 0;];
Cq = t10;
