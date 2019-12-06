% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-05 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRPR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR1_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR1_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:37:52
% EndTime: 2019-12-05 18:38:03
% DurationCPUTime: 5.00s
% Computational Cost: add. (15214->236), mult. (29122->317), div. (0->0), fcn. (33520->8), ass. (0->142)
t212 = sin(qJ(3));
t213 = sin(qJ(2));
t215 = cos(qJ(3));
t216 = cos(qJ(2));
t194 = -t212 * t213 + t215 * t216;
t195 = -t212 * t216 - t215 * t213;
t209 = sin(pkin(9));
t210 = cos(pkin(9));
t162 = t194 * t209 - t195 * t210;
t211 = sin(qJ(5));
t214 = cos(qJ(5));
t235 = t210 * t194 + t195 * t209;
t110 = t162 * t214 + t211 * t235;
t318 = -t162 * t211 + t214 * t235;
t338 = Ifges(6,5) * t318 - Ifges(6,6) * t110;
t305 = -pkin(7) - pkin(6);
t200 = t305 * t213;
t201 = t305 * t216;
t171 = t200 * t212 - t215 * t201;
t151 = qJ(4) * t194 + t171;
t314 = t215 * t200 + t212 * t201;
t321 = t195 * qJ(4) + t314;
t337 = -t209 * t151 + t210 * t321;
t343 = -pkin(8) * t162 + t337;
t93 = t210 * t151 + t209 * t321;
t67 = pkin(8) * t235 + t93;
t358 = -t211 * t67 + t214 * t343;
t363 = t358 * mrSges(6,2);
t41 = t211 * t343 + t214 * t67;
t366 = t41 * mrSges(6,1);
t369 = -t366 / 0.2e1 - t363 / 0.2e1;
t373 = t338 + 0.2e1 * t369;
t374 = t373 * qJD(5);
t112 = mrSges(5,1) * t162 + mrSges(5,2) * t235;
t206 = -pkin(2) * t216 - pkin(1);
t175 = -t194 * pkin(3) + t206;
t255 = t235 ^ 2;
t257 = t162 ^ 2;
t260 = t318 ^ 2;
t261 = t110 ^ 2;
t117 = -pkin(4) * t235 + t175;
t329 = t318 * mrSges(6,2);
t58 = t110 * mrSges(6,1) + t329;
t340 = t117 * t58;
t344 = (Ifges(6,1) - Ifges(6,2)) * t318;
t372 = t340 + t175 * t112 + t206 * (-mrSges(4,1) * t195 + mrSges(4,2) * t194) + (t260 - t261) * Ifges(6,4) + (t255 - t257) * Ifges(5,4) + t110 * t344 + (Ifges(5,1) - Ifges(5,2)) * t235 * t162;
t355 = -t171 * mrSges(4,1) - t93 * mrSges(5,1) - t314 * mrSges(4,2) - t337 * mrSges(5,2) + Ifges(4,5) * t194 + Ifges(5,5) * t235 + Ifges(4,6) * t195 - Ifges(5,6) * t162 + t338;
t359 = -t363 - t366;
t371 = t355 + t359;
t370 = m(6) * t117 - mrSges(6,1) * t318 + mrSges(6,2) * t110;
t202 = pkin(3) * t210 + pkin(4);
t298 = pkin(3) * t209;
t182 = t202 * t214 - t211 * t298;
t184 = t202 * t211 + t214 * t298;
t364 = m(6) * (-t182 * t41 + t184 * t358);
t353 = m(5) * t175;
t347 = t209 * t337 - t210 * t93;
t302 = -t110 / 0.2e1;
t303 = -t329 / 0.2e1;
t245 = t210 * t235;
t248 = t209 * t162;
t324 = (t248 / 0.2e1 + t245 / 0.2e1) * pkin(3);
t308 = m(5) * pkin(3);
t241 = -t308 / 0.2e1;
t322 = -t364 / 0.2e1 + t347 * t241 - t369;
t244 = t210 * t212;
t188 = (-t209 * t215 - t244) * pkin(2);
t247 = t209 * t212;
t189 = (t210 * t215 - t247) * pkin(2);
t153 = t188 * t214 - t189 * t211;
t152 = t153 * mrSges(6,1);
t154 = t188 * t211 + t189 * t214;
t272 = t154 * mrSges(6,2);
t320 = t188 * mrSges(5,1) - t189 * mrSges(5,2) - (mrSges(4,1) * t212 + mrSges(4,2) * t215) * pkin(2) + t152 - t272;
t310 = m(5) / 0.2e1;
t309 = m(6) / 0.2e1;
t205 = pkin(2) * t215 + pkin(3);
t181 = -pkin(2) * t247 + t210 * t205;
t180 = pkin(4) + t181;
t183 = pkin(2) * t244 + t205 * t209;
t141 = t180 * t214 - t183 * t211;
t301 = -t141 / 0.2e1;
t300 = t154 / 0.2e1;
t299 = pkin(3) * t195;
t208 = t213 * pkin(2);
t281 = Ifges(4,4) * t195;
t113 = -mrSges(5,1) * t235 + mrSges(5,2) * t162;
t124 = pkin(4) * t162 - t299;
t119 = t124 + t208;
t178 = t208 - t299;
t234 = Ifges(4,4) * t194 + (-Ifges(4,1) + Ifges(4,2)) * t195;
t1 = m(4) * t206 * t208 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t213) * t213 + (-mrSges(4,1) * t208 + t234) * t194 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t216 + (Ifges(3,1) - Ifges(3,2)) * t213) * t216 + (-mrSges(4,2) * t208 - t281) * t195 + (t353 + t113) * t178 + t370 * t119 + t372;
t278 = t1 * qJD(1);
t2 = (-pkin(3) * t113 - t281) * t195 + t234 * t194 - t299 * t353 + t370 * t124 + t372;
t266 = t2 * qJD(1);
t5 = -Ifges(6,4) * t260 - t340 + (Ifges(6,4) * t110 - t344) * t110;
t265 = t5 * qJD(1);
t10 = -m(6) * (-t110 * t358 + t318 * t41) - m(5) * (-t162 * t337 + t235 * t93) + (-t261 - t260) * mrSges(6,3) + (-t257 - t255) * mrSges(5,3);
t263 = qJD(1) * t10;
t142 = t180 * t211 + t183 * t214;
t127 = t142 * mrSges(6,1);
t29 = t141 * mrSges(6,2) + t127;
t262 = qJD(5) * t29;
t259 = t141 * t318;
t258 = t142 * t110;
t220 = (-t110 * t141 + t142 * t318) * t309 + (-t162 * t181 + t183 * t235) * t310;
t227 = -t112 - t58;
t231 = t119 * t309 + t178 * t310;
t16 = t220 + t227 - t231;
t256 = t16 * qJD(1);
t240 = t308 / 0.2e1;
t224 = (-t162 * t210 + t209 * t235) * t240;
t250 = t184 * t318;
t252 = t182 * t110;
t17 = t195 * t240 + 0.2e1 * (-t252 / 0.4e1 + t250 / 0.4e1 - t124 / 0.4e1) * m(6) + t224 + t227;
t254 = t17 * qJD(1);
t253 = t181 * t235;
t251 = t183 * t162;
t27 = 0.2e1 * t302 * mrSges(6,1) + 0.2e1 * t303;
t243 = t27 * qJD(1);
t179 = t184 * mrSges(6,1);
t72 = t182 * mrSges(6,2) + t179;
t242 = t72 * qJD(5);
t239 = t182 * t318 * mrSges(6,3);
t238 = t184 * t110 * mrSges(6,3);
t237 = -t127 / 0.2e1 - t179 / 0.2e1;
t236 = t301 - t182 / 0.2e1;
t23 = -m(5) * (t181 * t188 + t183 * t189) - m(6) * (t141 * t153 + t142 * t154) - t320;
t218 = ((t183 + t188) * t337 + (-t181 + t189) * t93) * t310 + ((t142 + t153) * t358 + (-t141 + t154) * t41) * t309 + t369;
t222 = -t188 * t162 / 0.2e1 + t189 * t235 / 0.2e1 - t253 / 0.2e1 - t251 / 0.2e1;
t228 = t153 * t302 + t300 * t318;
t4 = ((t182 / 0.2e1 + t301) * t318 - (-t184 / 0.2e1 + t142 / 0.2e1) * t110 + t228) * mrSges(6,3) + (t222 + t324) * mrSges(5,3) + t218 + t322;
t233 = t4 * qJD(1) - t23 * qJD(2);
t232 = t29 * qJD(2);
t19 = -t152 / 0.2e1 + (t300 + t236) * mrSges(6,2) + t237;
t225 = t19 * qJD(2) - t72 * qJD(3);
t28 = t329 / 0.2e1 + t303;
t24 = t195 * t241 + t224 + (t124 + t250 - t252) * t309;
t21 = t220 + t231;
t20 = -t272 / 0.2e1 + t152 / 0.2e1 + t236 * mrSges(6,2) + t237;
t3 = t218 - t239 / 0.2e1 - t238 / 0.2e1 + (-t259 / 0.2e1 - t258 / 0.2e1 + t228) * mrSges(6,3) + (t222 - t324) * mrSges(5,3) - t322 + t355;
t6 = [qJD(2) * t1 + qJD(3) * t2 - qJD(4) * t10 - qJD(5) * t5, t3 * qJD(3) + t21 * qJD(4) + t374 + t278 + (Ifges(3,5) * t216 - Ifges(3,6) * t213 + 0.2e1 * (-t181 * t93 + t183 * t337) * t310 + 0.2e1 * (-t141 * t41 + t142 * t358) * t309 + (m(4) * (-t171 * t215 + t212 * t314) + (-t194 * t215 + t195 * t212) * mrSges(4,3)) * pkin(2) + (-mrSges(3,1) * t216 + mrSges(3,2) * t213) * pkin(6) + (-t258 - t259) * mrSges(6,3) + (-t251 - t253) * mrSges(5,3) + t371) * qJD(2), t3 * qJD(2) + t24 * qJD(4) + t374 + t266 + (-t239 - t238 + t364 + (m(5) * t347 + (-t245 - t248) * mrSges(5,3)) * pkin(3) + t371) * qJD(3), qJD(2) * t21 + qJD(3) * t24 + qJD(5) * t28 - t263, -t265 + t28 * qJD(4) + (t338 + t359) * qJD(5) + (qJD(2) + qJD(3)) * t373; qJD(3) * t4 + qJD(4) * t16 - t278, -qJD(3) * t23 - t262, t320 * qJD(3) + t20 * qJD(5) + 0.2e1 * ((t153 * t182 + t154 * t184) * t309 + (t188 * t210 + t189 * t209) * t240) * qJD(3) + t233, t256, t20 * qJD(3) - t232 - t262; -qJD(2) * t4 + qJD(4) * t17 - t266, qJD(5) * t19 - t233, -t242, t254, t225 - t242; -qJD(2) * t16 - qJD(3) * t17 - qJD(5) * t27 + t263, -t256, -t254, 0, -t243; qJD(4) * t27 + t265, -qJD(3) * t19 + t232, -t225, t243, 0;];
Cq = t6;
