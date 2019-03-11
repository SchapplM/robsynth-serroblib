% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPPPRR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR3_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR3_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:33:30
% EndTime: 2019-03-09 01:33:34
% DurationCPUTime: 2.54s
% Computational Cost: add. (8340->294), mult. (15285->421), div. (0->0), fcn. (16068->8), ass. (0->184)
t207 = cos(pkin(9));
t346 = t207 * qJ(2);
t205 = sin(pkin(9));
t345 = -t205 / 0.2e1;
t204 = sin(pkin(10));
t206 = cos(pkin(10));
t210 = sin(qJ(5));
t314 = cos(qJ(5));
t182 = t314 * t204 + t210 * t206;
t209 = sin(qJ(6));
t202 = t209 ^ 2;
t211 = cos(qJ(6));
t203 = t211 ^ 2;
t332 = t202 + t203;
t256 = (-0.1e1 + t332) * t182;
t344 = m(7) * t256;
t325 = pkin(1) + pkin(2);
t232 = -t205 * t325 - qJ(4) + t346;
t226 = -pkin(7) + t232;
t159 = t226 * t206;
t224 = t204 * t226;
t107 = t314 * t159 - t210 * t224;
t235 = -t210 * t204 + t314 * t206;
t300 = t209 * mrSges(7,3);
t266 = t182 * t300;
t134 = -mrSges(7,2) * t235 + t266;
t272 = t211 * t134;
t296 = t211 * mrSges(7,3);
t136 = mrSges(7,1) * t235 + t182 * t296;
t276 = t209 * t136;
t236 = t276 / 0.2e1 - t272 / 0.2e1;
t194 = t205 * qJ(2);
t243 = t207 * t325 + pkin(3) + t194;
t234 = t206 * pkin(4) + t243;
t110 = pkin(5) * t235 + t182 * pkin(8) + t234;
t51 = t107 * t211 + t209 * t110;
t295 = t211 * t51;
t50 = -t209 * t107 + t110 * t211;
t251 = -t209 * t50 + t295;
t326 = m(7) / 0.2e1;
t343 = t236 + (t107 - t251) * t326;
t305 = Ifges(7,2) * t211;
t306 = Ifges(7,4) * t209;
t185 = t305 + t306;
t199 = Ifges(7,4) * t211;
t307 = Ifges(7,1) * t209;
t187 = t199 + t307;
t315 = t211 / 0.2e1;
t318 = -t209 / 0.2e1;
t342 = t185 * t318 + t187 * t315;
t269 = t204 ^ 2 + t206 ^ 2;
t341 = (-m(5) * t232 + mrSges(5,3)) * t269;
t183 = -mrSges(7,1) * t211 + t209 * mrSges(7,2);
t121 = t182 * t183;
t340 = -t121 / 0.2e1;
t339 = m(5) * t345;
t198 = Ifges(7,5) * t211;
t265 = -t198 / 0.2e1;
t304 = Ifges(7,6) * t209;
t337 = Ifges(6,4) + t265 + t304 / 0.2e1;
t297 = t211 * mrSges(7,2);
t301 = t209 * mrSges(7,1);
t184 = t297 + t301;
t123 = t184 * t182;
t334 = mrSges(6,3) * t182 + t123;
t178 = t235 * mrSges(6,2);
t141 = -t182 * mrSges(6,1) - t178;
t333 = -t209 * Ifges(7,2) + t199;
t239 = t297 / 0.2e1 + t301 / 0.2e1;
t320 = t184 / 0.2e1;
t331 = t239 + t320;
t106 = t159 * t210 + t314 * t224;
t311 = pkin(8) * t235;
t312 = pkin(5) * t182;
t144 = t311 - t312;
t58 = t209 * t106 + t144 * t211;
t59 = -t106 * t211 + t209 * t144;
t250 = -t209 * t58 + t211 * t59;
t283 = t235 * t209;
t133 = mrSges(7,2) * t182 + mrSges(7,3) * t283;
t282 = t235 * t211;
t135 = -t182 * mrSges(7,1) + mrSges(7,3) * t282;
t164 = t235 * t205;
t145 = -t209 * t164 - t207 * t211;
t146 = t164 * t211 - t209 * t207;
t323 = t145 / 0.2e1;
t324 = -t123 / 0.2e1;
t329 = (t164 * t106 + t145 * t58 + t146 * t59) * t326 + t135 * t323 + t146 * t133 / 0.2e1 + t164 * t324 - t207 * t141 / 0.2e1;
t328 = -m(7) * pkin(5) - mrSges(6,1) + t183;
t327 = m(6) / 0.2e1;
t162 = t182 * t205;
t322 = t162 / 0.2e1;
t321 = -t235 / 0.2e1;
t317 = t209 / 0.2e1;
t316 = -t211 / 0.2e1;
t288 = t146 * t211;
t289 = t145 * t209;
t246 = t288 - t289;
t238 = t164 - t246;
t21 = (-t162 * t256 - t235 * t238) * t326;
t310 = t21 * qJD(5);
t309 = mrSges(6,3) * t235;
t280 = t182 * t162;
t214 = (-t164 * t235 - t280) * t327 + (-t235 * t246 - t280) * t326 + t269 * t339;
t222 = t339 + (m(7) * t332 + m(6)) * t345;
t26 = t214 + t222;
t293 = qJD(1) * t26;
t165 = t235 * t207;
t279 = t182 * t165;
t163 = t182 * t207;
t284 = t235 * t163;
t36 = (t332 * t279 - t284) * t326 + (t279 - t284) * t327;
t292 = qJD(1) * t36;
t216 = (t288 / 0.2e1 - t289 / 0.2e1) * t182 * mrSges(7,3) + t134 * t323 - t146 * t136 / 0.2e1 + t121 * t322;
t244 = -t209 * t165 + t205 * t211;
t245 = t165 * t211 + t209 * t205;
t219 = -t245 * mrSges(7,2) / 0.2e1 + t244 * mrSges(7,1) / 0.2e1;
t10 = t216 - t219;
t291 = t10 * qJD(1);
t290 = t106 * t163;
t237 = t134 * t318 + t136 * t316;
t255 = mrSges(7,3) * (t203 / 0.2e1 + t202 / 0.2e1);
t16 = t121 * t321 + (t182 * t255 + t237) * t182;
t287 = t16 * qJD(1);
t286 = t162 * t163;
t218 = -t235 * t255 + (-t311 * t332 + t312) * t326;
t220 = (t209 * t59 + t211 * t58) * t326 + t133 * t317 + t135 * t315;
t18 = t178 + (-t183 / 0.2e1 + mrSges(6,1)) * t182 + t218 - t220;
t285 = t18 * qJD(1);
t281 = t182 * t106;
t103 = Ifges(7,6) * t235 - t182 * t333;
t278 = t209 * t103;
t277 = t209 * t135;
t188 = t211 * Ifges(7,1) - t306;
t105 = Ifges(7,5) * t235 - t182 * t188;
t274 = t211 * t105;
t273 = t211 * t133;
t230 = t239 * t235;
t23 = t230 + t236;
t271 = t23 * qJD(1);
t39 = t321 * t344;
t270 = t39 * qJD(1);
t259 = t198 - t304;
t258 = t235 * t332;
t253 = Ifges(7,5) * t209 + Ifges(7,6) * t211;
t102 = -Ifges(7,6) * t182 - t235 * t333;
t104 = -Ifges(7,5) * t182 - t188 * t235;
t122 = t235 * t184;
t1 = -t107 * t123 - t106 * t122 + t59 * t134 + t51 * t133 + t58 * t136 + t50 * t135 + m(7) * (t106 * t107 + t50 * t58 + t51 * t59) + t234 * t141 - (t274 / 0.2e1 - t278 / 0.2e1 - t337 * t235) * t235 + (t104 * t316 + t102 * t317 - t337 * t182 - (Ifges(7,3) - Ifges(6,1) + Ifges(6,2)) * t235) * t182;
t217 = -t122 / 0.2e1 + t343;
t8 = (t324 + (t106 + t250) * t326 + t273 / 0.2e1 - t277 / 0.2e1) * t182 - t217 * t235;
t252 = t1 * qJD(1) + t8 * qJD(3);
t5 = -t106 * t121 + t51 * t136 + (-mrSges(7,3) * t295 + t103 * t316 + t105 * t318 + t342 * t182 + t253 * t321) * t182 + (-t134 + t266) * t50;
t249 = -t5 * qJD(1) + t16 * qJD(3);
t9 = t165 * t309 - t245 * t134 - t244 * t136 - m(7) * (t50 * t244 + t51 * t245 + t290) - m(6) * (t107 * t165 + t290) - mrSges(3,3) + t334 * t163 - m(3) * qJ(2) + (-m(4) * t194 - m(5) * t243 - m(6) * t234 - mrSges(5,1) * t206 - mrSges(6,1) * t235 + mrSges(5,2) * t204 + mrSges(6,2) * t182 - mrSges(4,1)) * t205 + (-m(4) * t346 - mrSges(4,2) + t341) * t207;
t248 = -t9 * qJD(1) + t36 * qJD(3);
t12 = t334 * t182 - (t272 - t276 - t309) * t235 + m(7) * (-t235 * t251 - t281) + m(6) * (-t107 * t235 - t281) + t341;
t247 = -qJD(1) * t12 - qJD(3) * t39;
t242 = -t58 * mrSges(7,1) / 0.2e1 + t59 * mrSges(7,2) / 0.2e1;
t233 = t182 * t253;
t229 = -t184 / 0.2e1 + t239;
t31 = m(7) * t238 * t162;
t212 = t245 * t296 / 0.2e1 - t244 * t300 / 0.2e1 + (pkin(8) * t332 * t326 - mrSges(6,2) / 0.2e1) * t165 + (-pkin(5) * t326 - mrSges(6,1) / 0.2e1 + t183 / 0.2e1) * t163;
t4 = -t122 * t322 + t343 * t162 - t212 + t329;
t228 = t4 * qJD(1) + t31 * qJD(2) + t21 * qJD(3);
t40 = t235 * t344;
t227 = t8 * qJD(1) + t21 * qJD(2) + t40 * qJD(3);
t34 = t229 * t162;
t69 = -pkin(5) * t184 + (t187 / 0.2e1 + t333 / 0.2e1) * t211 + (t188 / 0.2e1 - t185 / 0.2e1) * t209;
t213 = pkin(8) * t255 + (-t188 / 0.4e1 + t185 / 0.4e1 + t305 / 0.4e1) * t211 + (t187 / 0.4e1 + t333 / 0.4e1 + t199 / 0.2e1 + t307 / 0.4e1) * t209;
t215 = t237 * pkin(8) + pkin(5) * t340 + t106 * t320 - t278 / 0.4e1 + t274 / 0.4e1;
t7 = -(0.3e1 / 0.4e1 * t304 - t198 / 0.4e1 + t265) * t235 + (Ifges(7,3) / 0.2e1 + t213) * t182 + t215 + t242;
t77 = t229 * t235;
t223 = t7 * qJD(1) - t34 * qJD(2) + t77 * qJD(3) + t69 * qJD(5);
t78 = t331 * t235;
t35 = t331 * t162;
t25 = t214 - t222;
t24 = t230 - t236;
t19 = t340 + t218 + t220;
t11 = t216 + t219;
t6 = t235 * t259 / 0.4e1 - Ifges(7,5) * t282 / 0.2e1 + Ifges(7,6) * t283 / 0.2e1 + t215 - t242 + (-Ifges(7,3) / 0.2e1 + t213) * t182;
t3 = t217 * t162 + t212 + t329;
t2 = qJD(2) * t36 + qJD(4) * t39 + qJD(5) * t8 + qJD(6) * t16;
t13 = [-qJD(2) * t9 + qJD(4) * t12 + qJD(5) * t1 - qJD(6) * t5, 0.2e1 * ((t145 * t244 + t146 * t245 + t286) * t326 + (t164 * t165 + t286) * t327 + (-m(6) / 0.2e1 + m(5) * (-0.1e1 + t269) / 0.2e1) * t205 * t207) * qJD(2) + t25 * qJD(4) + t3 * qJD(5) + t11 * qJD(6) + t248, t2, qJD(2) * t25 + qJD(5) * t19 + qJD(6) * t24 - t247, t3 * qJD(2) + t19 * qJD(4) + t6 * qJD(6) + t252 + (t104 * t317 + t102 * t315 + pkin(5) * t122 - t233 / 0.2e1 + Ifges(6,6) * t182 + t106 * mrSges(6,2) - (Ifges(6,5) + t342) * t235 + t328 * t107 + (m(7) * t250 + t273 - t277) * pkin(8) + t250 * mrSges(7,3)) * qJD(5), t11 * qJD(2) + t24 * qJD(4) + t6 * qJD(5) + (-t51 * mrSges(7,1) - t50 * mrSges(7,2) + t233) * qJD(6) + t249; qJD(4) * t26 + qJD(5) * t4 + qJD(6) * t10 - t248, t31 * qJD(5), -t292 + t310, t293 (t328 * t164 + (mrSges(6,2) - (m(7) * pkin(8) + mrSges(7,3)) * t332) * t162) * qJD(5) + t35 * qJD(6) + t228, t291 + t35 * qJD(5) + (-mrSges(7,1) * t146 - mrSges(7,2) * t145) * qJD(6); t2, t292 + t310, t40 * qJD(5), t270 (t121 + m(7) * (pkin(8) * t258 - t312) + mrSges(7,3) * t258 + t141) * qJD(5) - t78 * qJD(6) + t227, -t78 * qJD(5) + qJD(6) * t121 + t287; -qJD(2) * t26 - qJD(5) * t18 - qJD(6) * t23 + t247, -t293, -t270, 0, -t285, -t184 * qJD(6) - t271; -qJD(2) * t4 + qJD(4) * t18 + qJD(6) * t7 - t252, -qJD(6) * t34 - t228, qJD(6) * t77 - t227, t285, t69 * qJD(6) (pkin(8) * t183 + t259) * qJD(6) + t223; -qJD(2) * t10 + qJD(4) * t23 - qJD(5) * t7 - t249, t34 * qJD(5) - t291, -t77 * qJD(5) - t287, t271, -t223, 0;];
Cq  = t13;
