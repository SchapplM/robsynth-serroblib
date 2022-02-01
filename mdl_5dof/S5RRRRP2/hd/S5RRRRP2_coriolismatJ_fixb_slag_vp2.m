% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% m [6x1]
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
% Datum: 2022-01-20 11:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRRP2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP2_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP2_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP2_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:48:47
% EndTime: 2022-01-20 11:48:54
% DurationCPUTime: 3.33s
% Computational Cost: add. (7940->269), mult. (15550->323), div. (0->0), fcn. (14836->6), ass. (0->179)
t187 = sin(qJ(4));
t188 = sin(qJ(3));
t190 = cos(qJ(3));
t304 = cos(qJ(4));
t162 = -t187 * t188 + t304 * t190;
t163 = -t187 * t190 - t304 * t188;
t110 = -mrSges(5,1) * t162 - mrSges(5,2) * t163;
t336 = Ifges(5,2) + Ifges(6,2);
t338 = Ifges(5,1) + Ifges(6,1);
t221 = t336 - t338;
t369 = Ifges(6,4) + Ifges(5,4);
t224 = t369 * t163;
t350 = -Ifges(4,4) * t188 + (Ifges(4,1) - Ifges(4,2)) * t190;
t372 = t190 ^ 2;
t356 = t372 * Ifges(4,4);
t381 = (pkin(3) * t110 + t350) * t188 + t356 + (t221 * t162 - t224) * t163;
t309 = -pkin(8) - pkin(7);
t174 = t309 * t188;
t175 = t309 * t190;
t324 = t304 * t174 + t187 * t175;
t344 = t324 * mrSges(5,2);
t145 = t163 * qJ(5);
t342 = t324 + t145;
t357 = t342 * mrSges(6,2);
t118 = t187 * t174 - t304 * t175;
t359 = t118 * mrSges(5,1);
t252 = t162 * qJ(5);
t96 = t118 + t252;
t375 = t96 * mrSges(6,1);
t380 = -t344 - t357 - t359 - t375;
t189 = sin(qJ(2));
t296 = t189 * pkin(1);
t179 = pkin(7) + t296;
t287 = pkin(8) + t179;
t159 = t287 * t188;
t160 = t287 * t190;
t326 = -t304 * t159 - t187 * t160;
t343 = t326 * mrSges(5,2);
t341 = t326 + t145;
t358 = t341 * mrSges(6,2);
t106 = -t187 * t159 + t304 * t160;
t360 = t106 * mrSges(5,1);
t65 = t106 + t252;
t376 = t65 * mrSges(6,1);
t379 = -t343 - t358 - t360 - t376;
t378 = t344 / 0.2e1 + t357 / 0.2e1 + t359 / 0.2e1 + t375 / 0.2e1;
t377 = t343 / 0.2e1 + t358 / 0.2e1 + t360 / 0.2e1 + t376 / 0.2e1;
t109 = -mrSges(6,1) * t162 - mrSges(6,2) * t163;
t297 = pkin(4) * t163;
t298 = pkin(3) * t188;
t136 = -t297 + t298;
t253 = t136 * t109;
t373 = t253 + t381;
t371 = mrSges(5,1) + mrSges(6,1);
t370 = mrSges(5,2) + mrSges(6,2);
t362 = mrSges(6,3) * t341;
t361 = mrSges(6,3) * t342;
t355 = t188 ^ 2 + t372;
t191 = cos(qJ(2));
t301 = pkin(1) * t191;
t129 = t163 * t301;
t171 = mrSges(4,1) * t188 + mrSges(4,2) * t190;
t236 = t304 * pkin(3);
t180 = t236 + pkin(4);
t130 = t162 * t301;
t254 = t130 * t187;
t311 = m(5) * pkin(3);
t312 = m(6) / 0.2e1;
t315 = t371 * t129 / 0.2e1 - t370 * t130 / 0.2e1;
t354 = (pkin(3) * t254 + t129 * t180) * t312 + (t304 * t129 + t254) * t311 / 0.2e1 - t171 * t301 / 0.2e1 + t315;
t353 = t236 - t180;
t352 = t342 + t341;
t351 = t353 * t312;
t349 = (t342 / 0.2e1 + t341 / 0.2e1) * mrSges(6,3);
t346 = m(6) * t297;
t242 = t369 * t162;
t340 = (mrSges(6,3) + mrSges(5,3)) * (t129 * t163 + t130 * t162);
t241 = qJD(1) + qJD(2);
t107 = -t163 * mrSges(6,1) + t162 * mrSges(6,2);
t299 = pkin(3) * t187;
t238 = t162 * t299;
t250 = t180 * t163;
t27 = 0.2e1 * (t136 / 0.4e1 - t238 / 0.4e1 - t250 / 0.4e1) * m(6) + t107;
t333 = t241 * t27;
t68 = t107 - t346;
t332 = t241 * t68;
t331 = t221 * t163;
t216 = (Ifges(5,6) + Ifges(6,6)) * t163 + (Ifges(5,5) + Ifges(6,5)) * t162;
t279 = mrSges(6,3) * t162;
t230 = -t279 / 0.2e1;
t325 = pkin(4) * t230 + t216;
t108 = -mrSges(5,1) * t163 + mrSges(5,2) * t162;
t182 = -t190 * pkin(3) - pkin(2);
t128 = -t162 * pkin(4) + t182;
t119 = t128 - t301;
t168 = t182 - t301;
t322 = (t182 / 0.2e1 + t168 / 0.2e1) * t108 + (t119 / 0.2e1 + t128 / 0.2e1) * t107;
t321 = t370 * t304;
t320 = t355 * t191;
t317 = -mrSges(4,1) * t190 + mrSges(4,2) * t188;
t314 = t163 ^ 2;
t313 = m(5) / 0.2e1;
t310 = m(6) * pkin(4);
t303 = m(5) * t189;
t302 = m(6) * t189;
t300 = pkin(2) * t171;
t286 = t27 * qJD(3) + t68 * qJD(4);
t158 = t168 * t298;
t251 = t168 * t108;
t259 = t119 * t107;
t39 = t341 * t279;
t208 = t251 - t39 + t259;
t181 = -pkin(2) - t301;
t249 = t181 * t171;
t78 = t119 * t136;
t5 = t249 + m(5) * t158 + m(6) * t78 + (t242 + t362) * t162 + t208 + t373;
t262 = t5 * qJD(1);
t206 = t242 + t331;
t79 = t109 * t297;
t209 = -t314 * t369 - t79;
t7 = -t119 * t346 + (t206 + t362) * t162 + t208 + t209;
t261 = t7 * qJD(1);
t220 = (t162 ^ 2 + t314) * mrSges(6,3);
t23 = m(6) * (t162 * t65 + t163 * t341) + t220;
t260 = qJD(1) * t23;
t258 = t128 * t107;
t202 = (t355 * mrSges(4,3) - mrSges(3,2)) * t191 + (-mrSges(3,1) + t109 + t110 + t317) * t189;
t13 = m(6) * (t129 * t341 + t130 * t65) + m(5) * (t106 * t130 + t129 * t326) + (t119 * t302 + t168 * t303 + m(4) * (t320 * t179 + t181 * t189) + t202) * pkin(1) + t340;
t256 = t13 * qJD(1);
t248 = t182 * t108;
t240 = -t310 / 0.2e1;
t239 = t310 / 0.2e1;
t237 = t163 * t299;
t234 = t296 / 0.2e1;
t48 = t342 * t279;
t233 = -t39 / 0.2e1 - t48 / 0.2e1;
t229 = t279 / 0.2e1;
t223 = t352 * t163;
t222 = (t65 + t96) * t162;
t219 = pkin(4) * t229;
t218 = t162 * t236;
t217 = -t237 / 0.2e1;
t169 = t182 * t298;
t97 = t128 * t136;
t195 = t253 + (t158 + t169) * t313 + (t78 + t97) * t312 + t233;
t1 = t195 + t258 / 0.2e1 + t259 / 0.2e1 + t251 / 0.2e1 + t248 / 0.2e1 + t249 / 0.2e1 + t356 - t300 / 0.2e1 + t110 * t298 + t352 * t229 + t350 * t188 - t224 * t163 - t354 + (t331 / 0.2e1 + t242 + (-t338 / 0.2e1 + t336 / 0.2e1) * t163) * t162;
t207 = t248 - t48 + t258;
t6 = -t300 + m(5) * t169 + m(6) * t97 + (t242 + t361) * t162 + t207 + t373;
t213 = t1 * qJD(1) + t6 * qJD(2);
t192 = (-t224 + (-t119 - t128) * t239) * t163 + (t349 + t206) * t162 - t79 + t233 + t322;
t3 = (mrSges(5,2) / 0.2e1 + mrSges(6,2) / 0.2e1) * t130 + (t240 - mrSges(5,1) / 0.2e1 - mrSges(6,1) / 0.2e1) * t129 + t192;
t8 = -t128 * t346 + (t206 + t361) * t162 + t207 + t209;
t212 = t3 * qJD(1) + t8 * qJD(2);
t16 = (t234 - t223 / 0.2e1 - t222 / 0.2e1) * m(6) - t220;
t24 = m(6) * (t162 * t96 + t163 * t342) + t220;
t211 = -qJD(1) * t16 + qJD(2) * t24;
t137 = mrSges(6,3) * t237;
t138 = mrSges(5,3) * t237;
t197 = t137 / 0.2e1 + t138 / 0.2e1 + t180 * t230 + mrSges(5,3) * t217 + (t217 + t218 / 0.2e1) * mrSges(6,3);
t194 = t65 * t351 + t197 - t377;
t201 = -t65 * t240 + t377;
t10 = t194 + t219 + t201;
t193 = t96 * t351 + t197 - t378;
t200 = -t96 * t240 + t378;
t12 = t193 + t219 + t200;
t80 = t321 * pkin(3) + (-m(6) * t353 + t371) * t299;
t205 = t10 * qJD(1) + t12 * qJD(2) - t80 * qJD(3);
t198 = -mrSges(5,3) * t218 + Ifges(4,5) * t190 - Ifges(4,6) * t188 - t180 * t279 + t137 + t138 + t216;
t67 = t68 * qJD(5);
t38 = (t238 + t250 + t136) * t312;
t30 = t38 * qJD(3);
t29 = t38 * qJD(5);
t26 = t27 * qJD(5);
t17 = (t223 + t222) * t312 + m(6) * t234 + t220;
t11 = t193 - t200 + t325;
t9 = t194 - t201 + t325;
t4 = t129 * t239 + t192 + t315;
t2 = t195 + (t349 + t242) * t162 + (-pkin(2) / 0.2e1 + t181 / 0.2e1) * t171 + t322 + t354 + t381;
t14 = [qJD(2) * t13 + qJD(3) * t5 + qJD(4) * t7 + qJD(5) * t23, t2 * qJD(3) + t4 * qJD(4) + t17 * qJD(5) + t256 + (0.2e1 * (t129 * t342 + t130 * t96) * t312 + 0.2e1 * (t118 * t130 + t129 * t324) * t313 + (t128 * t302 + t182 * t303 + m(4) * (-pkin(2) * t189 + t320 * pkin(7)) + t202) * pkin(1) + t340) * qJD(2), t262 + t2 * qJD(2) + (t198 + (-t106 * t304 + t187 * t326) * t311 + m(6) * (-t180 * t65 + t299 * t341) + t317 * t179 + t379) * qJD(3) + t9 * qJD(4) + t29, t261 + t4 * qJD(2) + t9 * qJD(3) + ((-m(6) * t65 - t279) * pkin(4) + t216 + t379) * qJD(4), qJD(2) * t17 + t260 + t30; qJD(3) * t1 + qJD(4) * t3 - qJD(5) * t16 - t256, qJD(3) * t6 + qJD(4) * t8 + qJD(5) * t24, (t198 + m(6) * (-t180 * t96 + t299 * t342) + (-t118 * t304 + t187 * t324) * t311 + t317 * pkin(7) + t380) * qJD(3) + t11 * qJD(4) + t29 + t213, t11 * qJD(3) + ((-m(6) * t96 - t279) * pkin(4) + t216 + t380) * qJD(4) + t212, t211 + t30; -qJD(2) * t1 + qJD(4) * t10 - t26 - t262, qJD(4) * t12 - t213 - t26, -t80 * qJD(4), ((-t310 - t371) * t187 - t321) * qJD(4) * pkin(3) + t205, -t333; -qJD(2) * t3 - qJD(3) * t10 - t261 - t67, -qJD(3) * t12 - t212 - t67, -t205, 0, -t332; qJD(2) * t16 - t260 + t286, -t211 + t286, t333, t332, 0;];
Cq = t14;
