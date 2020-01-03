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
% Datum: 2020-01-03 12:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:11:08
% EndTime: 2020-01-03 12:11:14
% DurationCPUTime: 3.35s
% Computational Cost: add. (7940->274), mult. (15550->327), div. (0->0), fcn. (14836->6), ass. (0->181)
t186 = sin(qJ(4));
t187 = sin(qJ(3));
t189 = cos(qJ(3));
t305 = cos(qJ(4));
t161 = -t186 * t187 + t189 * t305;
t162 = -t186 * t189 - t187 * t305;
t110 = -mrSges(5,1) * t161 - mrSges(5,2) * t162;
t337 = Ifges(5,2) + Ifges(6,2);
t339 = Ifges(5,1) + Ifges(6,1);
t222 = t337 - t339;
t371 = Ifges(5,4) + Ifges(6,4);
t225 = t371 * t162;
t351 = -Ifges(4,4) * t187 + (Ifges(4,1) - Ifges(4,2)) * t189;
t374 = t189 ^ 2;
t358 = t374 * Ifges(4,4);
t382 = (pkin(3) * t110 + t351) * t187 + t358 + (t161 * t222 - t225) * t162;
t310 = -pkin(8) - pkin(7);
t173 = t310 * t187;
t174 = t310 * t189;
t326 = t173 * t305 + t174 * t186;
t344 = t326 * mrSges(5,2);
t144 = t162 * qJ(5);
t342 = t326 + t144;
t359 = t342 * mrSges(6,2);
t118 = t173 * t186 - t174 * t305;
t361 = t118 * mrSges(5,1);
t253 = t161 * qJ(5);
t96 = t118 + t253;
t377 = t96 * mrSges(6,1);
t381 = -t344 - t359 - t361 - t377;
t188 = sin(qJ(2));
t297 = t188 * pkin(1);
t178 = pkin(7) + t297;
t288 = pkin(8) + t178;
t158 = t288 * t187;
t159 = t288 * t189;
t327 = -t158 * t305 - t159 * t186;
t343 = t327 * mrSges(5,2);
t341 = t327 + t144;
t360 = t341 * mrSges(6,2);
t106 = -t158 * t186 + t159 * t305;
t362 = t106 * mrSges(5,1);
t65 = t106 + t253;
t378 = t65 * mrSges(6,1);
t380 = -t343 - t360 - t362 - t378;
t379 = t344 / 0.2e1 + t359 / 0.2e1 + t361 / 0.2e1 + t377 / 0.2e1;
t109 = -mrSges(6,1) * t161 - mrSges(6,2) * t162;
t298 = pkin(4) * t162;
t299 = pkin(3) * t187;
t136 = -t298 + t299;
t254 = t136 * t109;
t375 = t254 + t382;
t311 = m(6) * pkin(4);
t241 = t311 / 0.2e1;
t373 = mrSges(5,1) + mrSges(6,1);
t372 = mrSges(5,2) + mrSges(6,2);
t368 = -t343 / 0.2e1 - t360 / 0.2e1 - t362 / 0.2e1;
t364 = mrSges(6,3) * t341;
t363 = mrSges(6,3) * t342;
t190 = cos(qJ(2));
t302 = pkin(1) * t190;
t129 = t162 * t302;
t130 = t161 * t302;
t340 = mrSges(6,3) + mrSges(5,3);
t357 = t340 * (t129 * t162 + t130 * t161);
t356 = t187 ^ 2 + t374;
t170 = mrSges(4,1) * t187 + mrSges(4,2) * t189;
t238 = t305 * pkin(3);
t179 = t238 + pkin(4);
t255 = t130 * t186;
t312 = m(5) * pkin(3);
t313 = m(6) / 0.2e1;
t317 = t373 * t129 / 0.2e1 - t372 * t130 / 0.2e1;
t355 = (pkin(3) * t255 + t129 * t179) * t313 + (t129 * t305 + t255) * t312 / 0.2e1 - t170 * t302 / 0.2e1 + t317;
t354 = t238 - t179;
t353 = t342 + t341;
t352 = m(6) * t354;
t300 = pkin(3) * t186;
t239 = t162 * t300;
t137 = mrSges(6,3) * t239;
t138 = mrSges(5,3) * t239;
t220 = t161 * t238;
t350 = t137 / 0.2e1 + t138 / 0.2e1 + mrSges(6,3) * t220 / 0.2e1 - t340 * t239 / 0.2e1;
t349 = (t342 / 0.2e1 + t341 / 0.2e1) * mrSges(6,3);
t346 = m(6) * t298;
t243 = t371 * t161;
t242 = qJD(1) + qJD(2);
t146 = t161 * mrSges(6,2);
t107 = -t162 * mrSges(6,1) + t146;
t240 = t161 * t300;
t251 = t179 * t162;
t27 = 0.2e1 * (t240 / 0.4e1 + t251 / 0.4e1 - t136 / 0.4e1) * m(6) - t107;
t334 = t242 * t27;
t234 = mrSges(6,1) + t311;
t68 = t162 * t234 - t146;
t333 = t242 * t68;
t332 = t162 * t222;
t108 = -mrSges(5,1) * t162 + mrSges(5,2) * t161;
t181 = -pkin(3) * t189 - pkin(2);
t128 = -pkin(4) * t161 + t181;
t119 = t128 - t302;
t167 = t181 - t302;
t324 = (t181 / 0.2e1 + t167 / 0.2e1) * t108 + (t119 / 0.2e1 + t128 / 0.2e1) * t107;
t323 = t372 * t305;
t322 = t356 * t190;
t319 = -mrSges(4,1) * t189 + mrSges(4,2) * t187;
t315 = t162 ^ 2;
t314 = m(5) / 0.2e1;
t304 = m(5) * t188;
t303 = m(6) * t188;
t301 = pkin(2) * t170;
t287 = -qJD(3) * t27 - qJD(4) * t68;
t280 = mrSges(6,3) * t161;
t157 = t167 * t299;
t252 = t167 * t108;
t260 = t119 * t107;
t39 = t341 * t280;
t206 = t252 - t39 + t260;
t180 = -pkin(2) - t302;
t250 = t180 * t170;
t78 = t119 * t136;
t5 = t250 + m(5) * t157 + m(6) * t78 + (t243 + t364) * t161 + t206 + t375;
t263 = t5 * qJD(1);
t204 = t243 + t332;
t79 = t109 * t298;
t208 = -t315 * t371 - t79;
t7 = -t119 * t346 + (t204 + t364) * t161 + t206 + t208;
t262 = t7 * qJD(1);
t221 = (t161 ^ 2 + t315) * mrSges(6,3);
t23 = m(6) * (t161 * t65 + t162 * t341) + t221;
t261 = qJD(1) * t23;
t259 = t128 * t107;
t199 = (mrSges(4,3) * t356 - mrSges(3,2)) * t190 + (-mrSges(3,1) + t109 + t110 + t319) * t188;
t13 = m(6) * (t129 * t341 + t130 * t65) + m(5) * (t106 * t130 + t129 * t327) + (t119 * t303 + t167 * t304 + m(4) * (t178 * t322 + t180 * t188) + t199) * pkin(1) + t357;
t257 = t13 * qJD(1);
t249 = t181 * t108;
t236 = t297 / 0.2e1;
t48 = t342 * t280;
t235 = -t39 / 0.2e1 - t48 / 0.2e1;
t231 = -t280 / 0.2e1;
t230 = t280 / 0.2e1;
t224 = t353 * t162;
t223 = (t65 + t96) * t161;
t217 = mrSges(6,1) / 0.2e1 + t241;
t216 = (Ifges(5,6) + Ifges(6,6)) * t162 + (Ifges(5,5) + Ifges(6,5)) * t161;
t168 = t181 * t299;
t97 = t128 * t136;
t192 = t254 + (t157 + t168) * t314 + (t78 + t97) * t313 + t235;
t2 = -t301 / 0.2e1 + t192 + t259 / 0.2e1 + t260 / 0.2e1 + t252 / 0.2e1 + t249 / 0.2e1 + t250 / 0.2e1 + t110 * t299 + t358 + t353 * t230 + t351 * t187 - t225 * t162 - t355 + (t332 / 0.2e1 + t243 + (-t339 / 0.2e1 + t337 / 0.2e1) * t162) * t161;
t205 = t249 - t48 + t259;
t6 = -t301 + m(5) * t168 + m(6) * t97 + (t243 + t363) * t161 + t205 + t375;
t212 = t2 * qJD(1) + t6 * qJD(2);
t191 = (-t225 + (-t119 - t128) * t241) * t162 + (t349 + t204) * t161 - t79 + t235 + t324;
t4 = (mrSges(5,2) / 0.2e1 + mrSges(6,2) / 0.2e1) * t130 + (-mrSges(5,1) / 0.2e1 - t217) * t129 + t191;
t8 = -t128 * t346 + (t204 + t363) * t161 + t205 + t208;
t211 = t4 * qJD(1) + t8 * qJD(2);
t16 = (t236 - t224 / 0.2e1 - t223 / 0.2e1) * m(6) - t221;
t24 = m(6) * (t161 * t96 + t162 * t342) + t221;
t210 = -qJD(1) * t16 + qJD(2) * t24;
t207 = t179 * t231 + t350;
t195 = t378 / 0.2e1 - t65 * t352 / 0.2e1 - t368;
t143 = pkin(4) * t231;
t201 = t143 + t368;
t10 = -t378 / 0.2e1 - t65 * t241 + t179 * t230 + t195 + t201 - t350;
t193 = t354 * t313 * t96 + t207 - t379;
t198 = t96 * t241 + t379;
t12 = pkin(4) * t230 + t193 + t198;
t80 = t323 * pkin(3) + (-t352 + t373) * t300;
t203 = qJD(1) * t10 - qJD(2) * t12 + qJD(3) * t80;
t196 = -mrSges(5,3) * t220 + Ifges(4,5) * t189 - Ifges(4,6) * t187 - t179 * t280 + t137 + t138 + t216;
t66 = t68 * qJD(5);
t38 = (t240 + t251 + t136) * t313;
t30 = t38 * qJD(3);
t29 = t38 * qJD(5);
t25 = t27 * qJD(5);
t17 = (t224 + t223) * t313 + m(6) * t236 + t221;
t11 = t193 + t143 - t198 + t216;
t9 = -t217 * t65 - t195 + t201 + t207 + t216;
t3 = t129 * t241 + t191 + t317;
t1 = t192 + (-pkin(2) / 0.2e1 + t180 / 0.2e1) * t170 + (t349 + t243) * t161 + t324 + t355 + t382;
t14 = [qJD(2) * t13 + qJD(3) * t5 + qJD(4) * t7 + qJD(5) * t23, t1 * qJD(3) + t3 * qJD(4) + t17 * qJD(5) + t257 + (0.2e1 * (t129 * t342 + t130 * t96) * t313 + 0.2e1 * (t118 * t130 + t129 * t326) * t314 + (t128 * t303 + t181 * t304 + m(4) * (-pkin(2) * t188 + pkin(7) * t322) + t199) * pkin(1) + t357) * qJD(2), t263 + t1 * qJD(2) + (t196 + (-t106 * t305 + t186 * t327) * t312 + m(6) * (-t179 * t65 + t300 * t341) + t319 * t178 + t380) * qJD(3) + t9 * qJD(4) + t29, t262 + t3 * qJD(2) + t9 * qJD(3) + ((-m(6) * t65 - t280) * pkin(4) + t216 + t380) * qJD(4), qJD(2) * t17 + t261 + t30; qJD(3) * t2 + qJD(4) * t4 - qJD(5) * t16 - t257, qJD(3) * t6 + qJD(4) * t8 + qJD(5) * t24, (t196 + (-t118 * t305 + t186 * t326) * t312 + m(6) * (-t179 * t96 + t300 * t342) + t319 * pkin(7) + t381) * qJD(3) + t11 * qJD(4) + t29 + t212, t11 * qJD(3) + ((-m(6) * t96 - t280) * pkin(4) + t216 + t381) * qJD(4) + t211, t210 + t30; -qJD(2) * t2 - qJD(4) * t10 + t25 - t263, qJD(4) * t12 - t212 + t25, -t80 * qJD(4), ((-mrSges(5,1) - t234) * t186 - t323) * qJD(4) * pkin(3) - t203, t334; -qJD(2) * t4 + qJD(3) * t10 - t262 + t66, -qJD(3) * t12 - t211 + t66, t203, 0, t333; qJD(2) * t16 - t261 + t287, -t210 + t287, -t334, -t333, 0;];
Cq = t14;
