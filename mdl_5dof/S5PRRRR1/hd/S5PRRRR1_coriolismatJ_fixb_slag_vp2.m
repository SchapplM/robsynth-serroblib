% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2019-07-18 13:29
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRRR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(2,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR1_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_coriolismatJ_fixb_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR1_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:28:22
% EndTime: 2019-07-18 13:28:28
% DurationCPUTime: 2.29s
% Computational Cost: add. (3280->266), mult. (9057->412), div. (0->0), fcn. (9410->8), ass. (0->170)
t168 = sin(qJ(5));
t280 = t168 / 0.2e1;
t172 = cos(qJ(5));
t278 = -t172 / 0.2e1;
t266 = Ifges(6,4) * t168;
t150 = Ifges(6,1) * t172 - t266;
t277 = t172 / 0.2e1;
t162 = Ifges(6,4) * t172;
t295 = -Ifges(6,2) * t168 + t162;
t147 = Ifges(6,2) * t172 + t266;
t149 = Ifges(6,1) * t168 + t162;
t303 = t147 * t280 + t149 * t278;
t46 = t150 * t280 + t277 * t295 - t303;
t164 = t168 ^ 2;
t166 = t172 ^ 2;
t221 = t164 + t166;
t298 = t221 * mrSges(6,3);
t308 = mrSges(5,2) - t298;
t169 = sin(qJ(4));
t171 = sin(qJ(2));
t173 = cos(qJ(4));
t174 = cos(qJ(3));
t224 = t173 * t174;
t170 = sin(qJ(3));
t227 = t170 * t171;
t131 = -t169 * t227 + t171 * t224;
t175 = cos(qJ(2));
t103 = t131 * t172 - t168 * t175;
t307 = -t103 / 0.2e1;
t306 = t175 / 0.2e1;
t251 = t172 * mrSges(6,1);
t257 = t168 * mrSges(6,2);
t205 = t251 - t257;
t305 = mrSges(5,1) + t205;
t143 = t169 * t170 - t224;
t161 = Ifges(6,5) * t172;
t264 = Ifges(6,6) * t168;
t296 = t161 - t264;
t235 = t143 * t296;
t101 = -t131 * t168 - t172 * t175;
t241 = t103 * t168;
t274 = pkin(2) * t170;
t304 = (m(5) * t306 - m(6) * (t101 * t172 + t241) / 0.2e1) * t274;
t302 = qJD(3) + qJD(4);
t301 = t147 / 0.4e1 - t150 / 0.4e1;
t299 = Ifges(5,1) - Ifges(5,2);
t297 = t175 * t171;
t144 = t169 * t174 + t170 * t173;
t232 = t144 * t168;
t91 = -mrSges(6,2) * t143 - mrSges(6,3) * t232;
t231 = t144 * t172;
t93 = mrSges(6,1) * t143 - mrSges(6,3) * t231;
t294 = t168 * t91 + t172 * t93;
t129 = t144 * t171;
t250 = t172 * mrSges(6,2);
t258 = t168 * mrSges(6,1);
t193 = t250 / 0.2e1 + t258 / 0.2e1;
t204 = t250 + t258;
t281 = -t168 / 0.2e1;
t216 = mrSges(6,3) * t281;
t285 = t129 / 0.2e1;
t286 = mrSges(6,3) / 0.2e1;
t179 = t103 * t216 + t129 * t193 + t204 * t285 + t241 * t286;
t240 = t103 * t172;
t243 = t101 * t168;
t191 = t131 - t240 + t243;
t20 = m(6) * t191 * t129;
t223 = t20 * qJD(1);
t291 = qJD(5) * t179 + t223;
t80 = t204 * t144;
t234 = t143 * t168;
t90 = -mrSges(6,2) * t144 + mrSges(6,3) * t234;
t233 = t143 * t172;
t92 = mrSges(6,1) * t144 + mrSges(6,3) * t233;
t96 = mrSges(5,1) * t144 - mrSges(5,2) * t143;
t290 = t96 * t306 - t131 * t80 / 0.2e1 + t90 * t307 - t101 * t92 / 0.2e1;
t289 = m(6) / 0.2e1;
t288 = mrSges(6,1) / 0.2e1;
t287 = -mrSges(6,2) / 0.2e1;
t284 = -t143 / 0.2e1;
t283 = -t144 / 0.2e1;
t282 = t144 / 0.2e1;
t279 = t169 / 0.2e1;
t275 = pkin(2) * t169;
t273 = pkin(2) * t173;
t272 = pkin(2) * t174;
t271 = mrSges(5,3) * t143;
t270 = mrSges(5,3) * t144;
t267 = Ifges(5,4) * t144;
t263 = Ifges(6,3) * t144;
t132 = t143 * t175;
t259 = t132 * mrSges(5,2);
t59 = Ifges(6,6) * t143 + t144 * t295;
t256 = t168 * t59;
t254 = t168 * t92;
t253 = t168 * t93;
t61 = Ifges(6,5) * t143 + t144 * t150;
t249 = t172 * t61;
t248 = t172 * t90;
t247 = t172 * t91;
t102 = t132 * t168 + t171 * t172;
t242 = t102 * t168;
t104 = -t132 * t172 + t168 * t171;
t239 = t104 * t172;
t130 = t144 * t175;
t238 = t129 * t130;
t237 = t129 * t169;
t236 = t130 * t173;
t220 = t170 ^ 2 + t174 ^ 2;
t17 = m(6) * (t101 * t102 + t103 * t104 + t238) + m(4) * (-0.1e1 + t220) * t297 + (-t131 * t132 + t238 - t297) * m(5);
t228 = t17 * qJD(1);
t219 = pkin(2) * t237;
t222 = t221 * t219;
t119 = t131 * t273;
t213 = -t233 / 0.2e1;
t214 = t234 / 0.2e1;
t215 = Ifges(6,5) * t213 + Ifges(6,6) * t214 + t263 / 0.2e1;
t209 = t305 * t169;
t208 = t102 * t216 + t239 * t286 + t259 / 0.2e1 - t305 * t130 / 0.2e1;
t201 = Ifges(6,5) * t168 + Ifges(6,6) * t172;
t137 = Ifges(5,4) * t143;
t176 = pkin(2) ^ 2;
t196 = t168 * t90 + t172 * t92 + t96;
t58 = Ifges(6,6) * t144 - t143 * t295;
t60 = Ifges(6,5) * t144 - t143 * t150;
t97 = mrSges(5,1) * t143 + mrSges(5,2) * t144;
t1 = t59 * t214 - t58 * t232 / 0.2e1 + t61 * t213 + t60 * t231 / 0.2e1 + t143 * (t263 - t235) / 0.2e1 + (-Ifges(5,2) * t143 + t267) * t283 + (Ifges(4,4) * t174 - t196 * pkin(2)) * t174 + (-Ifges(4,4) * t170 + (t97 + t294) * pkin(2) + (Ifges(4,1) - Ifges(4,2) + (-m(6) * t221 - m(5)) * t176) * t174) * t170 + (t144 * t299 - 0.2e1 * t137) * t284 + (t296 * t144 - t267 + (-Ifges(5,1) + Ifges(6,3)) * t143) * t282;
t177 = 0.2e1 * (m(5) * (-t132 * t169 - t236) / 0.4e1 + m(6) * (-t236 + (t239 - t242) * t169) / 0.4e1) * pkin(2) + t208;
t79 = t204 * t143;
t185 = t247 / 0.2e1 - t253 / 0.2e1 + t79 / 0.2e1;
t2 = t129 * t185 + t177 + t290 + t304;
t200 = -t2 * qJD(1) + t1 * qJD(2);
t195 = t161 / 0.2e1 - t264 / 0.2e1;
t184 = t195 * t143;
t4 = t196 * t272 + (-t256 / 0.2e1 + t249 / 0.2e1 - t137 + t184) * t143 + (t58 * t280 + t60 * t278 + (Ifges(5,4) - t195) * t144 + (-Ifges(6,3) + t299) * t143) * t144;
t183 = -t129 * t247 / 0.2e1 - t290 + (t253 - t79) * t285;
t5 = -t259 / 0.2e1 + (mrSges(5,1) / 0.2e1 + t205 / 0.2e1) * t130 + (-t239 / 0.2e1 + t242 / 0.2e1) * mrSges(6,3) + t183;
t199 = qJD(1) * t5 - qJD(2) * t4;
t192 = pkin(2) * t204;
t186 = t173 * t192;
t35 = t186 - t46;
t189 = t173 * t205;
t81 = t144 * t147;
t82 = t144 * t149;
t7 = -t235 / 0.4e1 + (-t61 / 0.4e1 + t81 / 0.4e1) * t172 + (t82 / 0.4e1 + t59 / 0.4e1) * t168 + ((t170 * t288 + t279 * t93) * t172 + (t170 * t287 + t279 * t91) * t168) * pkin(2) + (t301 * t172 + (t295 / 0.4e1 + t149 / 0.4e1) * t168 + (t189 / 0.2e1 + (t166 / 0.2e1 + t164 / 0.2e1) * t169 * mrSges(6,3)) * pkin(2)) * t144 + t215;
t198 = -t7 * qJD(2) - t35 * qJD(3);
t178 = (t205 * t285 + (-t240 / 0.2e1 + t243 / 0.2e1) * mrSges(6,3)) * t144 + t101 * t91 / 0.2e1 + t93 * t307;
t194 = t102 * t288 + t104 * t287;
t10 = t178 - t194;
t14 = (-t247 + t253) * t272 + (t201 * t284 - t82 * t277 + t59 * t278 + (t61 - t81) * t281) * t144;
t197 = t10 * qJD(1) + t14 * qJD(2);
t188 = t129 * t308 - t305 * t131;
t187 = -Ifges(5,6) * t144 + t201 * t282 + t58 * t277 + t60 * t280 + (-Ifges(5,5) + t303) * t143;
t16 = (t185 * t173 + (t248 / 0.2e1 - t254 / 0.2e1 + t80 / 0.2e1) * t169) * pkin(2);
t21 = ((-t173 * t191 + t237) * pkin(2) - t222) * t289;
t34 = -pkin(2) * t209 + (m(6) * (-0.1e1 + t221) * t176 * t169 - t308 * pkin(2)) * t173;
t182 = t21 * qJD(1) + t16 * qJD(2) + t34 * qJD(3);
t180 = t235 / 0.4e1 - t256 / 0.4e1 + t249 / 0.4e1 - t168 * t82 / 0.4e1 - t172 * t81 / 0.4e1 - t301 * t231 - (t295 + t149) * t232 / 0.4e1;
t13 = -t263 / 0.2e1 + t184 + t180;
t29 = (-t149 / 0.2e1 - t295 / 0.2e1) * t172 + (-t150 / 0.2e1 + t147 / 0.2e1) * t168;
t181 = -qJD(2) * t13 + qJD(3) * t29 - qJD(4) * t46;
t12 = t180 + t215;
t30 = -t186 / 0.2e1 - t193 * t273 + t46;
t15 = t188 + t21;
t11 = t178 + t194;
t9 = t16 + t187;
t8 = t12 + pkin(2) * t189 * t283 + (t251 / 0.2e1 - t257 / 0.2e1) * t274 - (t144 * t298 + t294) * t275 / 0.2e1;
t6 = t183 + t208;
t3 = t177 + t183 - t304 + (-t170 * mrSges(4,1) - t174 * mrSges(4,2)) * t175;
t18 = [t17 * qJD(2) + t20 * t302, t3 * qJD(3) + t6 * qJD(4) + t11 * qJD(5) + t228 + (t132 * t271 + m(6) * (-t102 * t172 - t104 * t168) * t272 + t104 * t91 + t102 * t93 + (t270 + t80) * t130 + (mrSges(4,3) * t220 - mrSges(3,2)) * t175 + (t170 * mrSges(4,2) - mrSges(3,1) + t97 + (-m(5) * pkin(2) - mrSges(4,1)) * t174) * t171) * qJD(2), t3 * qJD(2) + t15 * qJD(4) + (-t171 * t174 * mrSges(4,1) + mrSges(4,2) * t227 + t188 + m(5) * (-t119 - t219) + 0.2e1 * (-t119 - t222) * t289) * qJD(3) + t291, qJD(2) * t6 + qJD(3) * t15 + qJD(4) * t188 + t291, t11 * qJD(2) + (-mrSges(6,1) * t103 - mrSges(6,2) * t101) * qJD(5) + t302 * t179; -qJD(3) * t2 + qJD(4) * t5 + qJD(5) * t10 - t228, qJD(3) * t1 - qJD(4) * t4 + qJD(5) * t14, t9 * qJD(4) + t8 * qJD(5) + t200 + (Ifges(4,5) * t174 - Ifges(4,6) * t170 + t187 + ((t79 + t271) * t173 + (t248 - t254 - t270) * t169) * pkin(2)) * qJD(3), qJD(3) * t9 + qJD(4) * t187 + qJD(5) * t12 + t199, t8 * qJD(3) + t12 * qJD(4) + (-t144 * t201 + t174 * t192) * qJD(5) + t197; qJD(2) * t2 + qJD(4) * t21 - t223, qJD(4) * t16 - qJD(5) * t7 - t200, qJD(4) * t34 - qJD(5) * t35, t30 * qJD(5) + (-t173 * t308 - t209) * qJD(4) * pkin(2) + t182, t30 * qJD(4) + (-t205 * t275 + t296) * qJD(5) + t198; -qJD(2) * t5 - qJD(3) * t21 - t223, -qJD(3) * t16 + qJD(5) * t13 - t199, -qJD(5) * t29 - t182, t46 * qJD(5), qJD(5) * t296 - t181; -t10 * qJD(2), qJD(3) * t7 - qJD(4) * t13 - t197, qJD(4) * t29 - t198, t181, 0;];
Cq  = t18;
