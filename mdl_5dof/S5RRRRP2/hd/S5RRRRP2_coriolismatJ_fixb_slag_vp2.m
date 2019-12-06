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
% Datum: 2019-12-05 18:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:47:36
% EndTime: 2019-12-05 18:47:44
% DurationCPUTime: 3.38s
% Computational Cost: add. (7940->274), mult. (15550->326), div. (0->0), fcn. (14836->6), ass. (0->179)
t187 = sin(qJ(4));
t188 = sin(qJ(3));
t190 = cos(qJ(3));
t305 = cos(qJ(4));
t162 = -t187 * t188 + t305 * t190;
t163 = -t187 * t190 - t305 * t188;
t110 = -mrSges(5,1) * t162 - mrSges(5,2) * t163;
t337 = Ifges(5,2) + Ifges(6,2);
t339 = Ifges(5,1) + Ifges(6,1);
t223 = t337 - t339;
t372 = Ifges(6,4) + Ifges(5,4);
t226 = t372 * t163;
t352 = -Ifges(4,4) * t188 + (Ifges(4,1) - Ifges(4,2)) * t190;
t375 = t190 ^ 2;
t359 = t375 * Ifges(4,4);
t383 = (pkin(3) * t110 + t352) * t188 + t359 + (t223 * t162 - t226) * t163;
t310 = -pkin(8) - pkin(7);
t174 = t310 * t188;
t175 = t310 * t190;
t326 = t305 * t174 + t187 * t175;
t345 = t326 * mrSges(5,2);
t145 = t163 * qJ(5);
t343 = t326 + t145;
t360 = t343 * mrSges(6,2);
t118 = t187 * t174 - t305 * t175;
t362 = t118 * mrSges(5,1);
t253 = t162 * qJ(5);
t96 = t118 + t253;
t378 = t96 * mrSges(6,1);
t382 = -t345 - t360 - t362 - t378;
t189 = sin(qJ(2));
t297 = t189 * pkin(1);
t179 = pkin(7) + t297;
t288 = pkin(8) + t179;
t159 = t288 * t188;
t160 = t288 * t190;
t327 = -t305 * t159 - t187 * t160;
t344 = t327 * mrSges(5,2);
t342 = t327 + t145;
t361 = t342 * mrSges(6,2);
t106 = -t187 * t159 + t305 * t160;
t363 = t106 * mrSges(5,1);
t65 = t106 + t253;
t379 = t65 * mrSges(6,1);
t381 = -t344 - t361 - t363 - t379;
t380 = t344 / 0.2e1 + t361 / 0.2e1 + t363 / 0.2e1 + t379 / 0.2e1;
t109 = -mrSges(6,1) * t162 - mrSges(6,2) * t163;
t298 = pkin(4) * t163;
t299 = pkin(3) * t188;
t136 = -t298 + t299;
t254 = t136 * t109;
t376 = t254 + t383;
t311 = m(6) * pkin(4);
t241 = t311 / 0.2e1;
t374 = mrSges(5,1) + mrSges(6,1);
t373 = mrSges(5,2) + mrSges(6,2);
t368 = -t345 / 0.2e1 - t360 / 0.2e1 - t362 / 0.2e1;
t365 = mrSges(6,3) * t342;
t364 = mrSges(6,3) * t343;
t191 = cos(qJ(2));
t302 = pkin(1) * t191;
t129 = t163 * t302;
t130 = t162 * t302;
t340 = mrSges(6,3) + mrSges(5,3);
t358 = t340 * (t129 * t163 + t130 * t162);
t357 = t188 ^ 2 + t375;
t171 = mrSges(4,1) * t188 + mrSges(4,2) * t190;
t238 = t305 * pkin(3);
t180 = t238 + pkin(4);
t255 = t130 * t187;
t312 = m(5) * pkin(3);
t313 = m(6) / 0.2e1;
t317 = t374 * t129 / 0.2e1 - t373 * t130 / 0.2e1;
t356 = (pkin(3) * t255 + t129 * t180) * t313 + (t305 * t129 + t255) * t312 / 0.2e1 - t171 * t302 / 0.2e1 + t317;
t355 = t238 - t180;
t354 = t343 + t342;
t353 = m(6) * t355;
t300 = pkin(3) * t187;
t239 = t163 * t300;
t137 = mrSges(6,3) * t239;
t138 = mrSges(5,3) * t239;
t221 = t162 * t238;
t351 = t137 / 0.2e1 + t138 / 0.2e1 + mrSges(6,3) * t221 / 0.2e1 - t340 * t239 / 0.2e1;
t350 = (t343 / 0.2e1 + t342 / 0.2e1) * mrSges(6,3);
t347 = m(6) * t298;
t243 = t372 * t162;
t242 = qJD(1) + qJD(2);
t107 = -t163 * mrSges(6,1) + t162 * mrSges(6,2);
t240 = t162 * t300;
t251 = t180 * t163;
t27 = 0.2e1 * (t136 / 0.4e1 - t240 / 0.4e1 - t251 / 0.4e1) * m(6) + t107;
t334 = t242 * t27;
t68 = t107 - t347;
t333 = t242 * t68;
t332 = t223 * t163;
t108 = -mrSges(5,1) * t163 + mrSges(5,2) * t162;
t182 = -t190 * pkin(3) - pkin(2);
t128 = -t162 * pkin(4) + t182;
t119 = t128 - t302;
t168 = t182 - t302;
t324 = (t182 / 0.2e1 + t168 / 0.2e1) * t108 + (t119 / 0.2e1 + t128 / 0.2e1) * t107;
t323 = t373 * t305;
t322 = t357 * t191;
t319 = -mrSges(4,1) * t190 + mrSges(4,2) * t188;
t315 = t163 ^ 2;
t314 = m(5) / 0.2e1;
t304 = m(5) * t189;
t303 = m(6) * t189;
t301 = pkin(2) * t171;
t287 = t27 * qJD(3) + t68 * qJD(4);
t280 = mrSges(6,3) * t162;
t158 = t168 * t299;
t252 = t168 * t108;
t260 = t119 * t107;
t39 = t342 * t280;
t207 = t252 - t39 + t260;
t181 = -pkin(2) - t302;
t250 = t181 * t171;
t78 = t119 * t136;
t5 = t250 + m(5) * t158 + m(6) * t78 + (t243 + t365) * t162 + t207 + t376;
t263 = t5 * qJD(1);
t205 = t243 + t332;
t79 = t109 * t298;
t209 = -t315 * t372 - t79;
t7 = -t119 * t347 + (t205 + t365) * t162 + t207 + t209;
t262 = t7 * qJD(1);
t222 = (t162 ^ 2 + t315) * mrSges(6,3);
t23 = m(6) * (t162 * t65 + t163 * t342) + t222;
t261 = qJD(1) * t23;
t259 = t128 * t107;
t200 = (t357 * mrSges(4,3) - mrSges(3,2)) * t191 + (-mrSges(3,1) + t109 + t110 + t319) * t189;
t13 = m(6) * (t129 * t342 + t130 * t65) + m(5) * (t106 * t130 + t129 * t327) + (t119 * t303 + t168 * t304 + m(4) * (t322 * t179 + t181 * t189) + t200) * pkin(1) + t358;
t257 = t13 * qJD(1);
t249 = t182 * t108;
t236 = t297 / 0.2e1;
t48 = t343 * t280;
t235 = -t39 / 0.2e1 - t48 / 0.2e1;
t232 = -t280 / 0.2e1;
t231 = t280 / 0.2e1;
t225 = t354 * t163;
t224 = (t65 + t96) * t162;
t218 = mrSges(6,1) / 0.2e1 + t241;
t217 = (Ifges(5,6) + Ifges(6,6)) * t163 + (Ifges(5,5) + Ifges(6,5)) * t162;
t169 = t182 * t299;
t97 = t128 * t136;
t193 = t254 + (t158 + t169) * t314 + (t78 + t97) * t313 + t235;
t1 = t193 - t301 / 0.2e1 + t259 / 0.2e1 + t260 / 0.2e1 + t252 / 0.2e1 + t250 / 0.2e1 + t249 / 0.2e1 + t110 * t299 + t359 + t354 * t231 + t352 * t188 - t226 * t163 - t356 + (t332 / 0.2e1 + t243 + (-t339 / 0.2e1 + t337 / 0.2e1) * t163) * t162;
t206 = t249 - t48 + t259;
t6 = -t301 + m(5) * t169 + m(6) * t97 + (t243 + t364) * t162 + t206 + t376;
t213 = t1 * qJD(1) + t6 * qJD(2);
t192 = (-t226 + (-t119 - t128) * t241) * t163 + (t350 + t205) * t162 - t79 + t235 + t324;
t3 = (mrSges(5,2) / 0.2e1 + mrSges(6,2) / 0.2e1) * t130 + (-mrSges(5,1) / 0.2e1 - t218) * t129 + t192;
t8 = -t128 * t347 + (t205 + t364) * t162 + t206 + t209;
t212 = t3 * qJD(1) + t8 * qJD(2);
t16 = (t236 - t225 / 0.2e1 - t224 / 0.2e1) * m(6) - t222;
t24 = m(6) * (t162 * t96 + t163 * t343) + t222;
t211 = -qJD(1) * t16 + qJD(2) * t24;
t208 = t180 * t232 + t351;
t194 = t355 * t313 * t65 + t208 - t380;
t199 = t65 * t241 + t380;
t10 = pkin(4) * t231 + t194 + t199;
t196 = t378 / 0.2e1 - t96 * t353 / 0.2e1 - t368;
t143 = pkin(4) * t232;
t202 = t143 + t368;
t12 = -t378 / 0.2e1 - t96 * t241 + t180 * t231 + t196 + t202 - t351;
t80 = t323 * pkin(3) + (-t353 + t374) * t300;
t204 = t10 * qJD(1) - t12 * qJD(2) - t80 * qJD(3);
t197 = -mrSges(5,3) * t221 + Ifges(4,5) * t190 - Ifges(4,6) * t188 - t180 * t280 + t137 + t138 + t217;
t67 = t68 * qJD(5);
t38 = (t240 + t251 + t136) * t313;
t30 = t38 * qJD(3);
t29 = t38 * qJD(5);
t26 = t27 * qJD(5);
t17 = (t225 + t224) * t313 + m(6) * t236 + t222;
t11 = -t218 * t96 - t196 + t202 + t208 + t217;
t9 = t194 + t143 - t199 + t217;
t4 = t129 * t241 + t192 + t317;
t2 = t193 + (-pkin(2) / 0.2e1 + t181 / 0.2e1) * t171 + (t350 + t243) * t162 + t324 + t356 + t383;
t14 = [qJD(2) * t13 + qJD(3) * t5 + qJD(4) * t7 + qJD(5) * t23, t2 * qJD(3) + t4 * qJD(4) + t17 * qJD(5) + t257 + (0.2e1 * (t129 * t343 + t130 * t96) * t313 + 0.2e1 * (t118 * t130 + t129 * t326) * t314 + (t128 * t303 + t182 * t304 + m(4) * (-pkin(2) * t189 + t322 * pkin(7)) + t200) * pkin(1) + t358) * qJD(2), t263 + t2 * qJD(2) + (t197 + (-t106 * t305 + t187 * t327) * t312 + m(6) * (-t180 * t65 + t300 * t342) + t319 * t179 + t381) * qJD(3) + t9 * qJD(4) + t29, t262 + t4 * qJD(2) + t9 * qJD(3) + ((-m(6) * t65 - t280) * pkin(4) + t217 + t381) * qJD(4), qJD(2) * t17 + t261 + t30; qJD(3) * t1 + qJD(4) * t3 - qJD(5) * t16 - t257, qJD(3) * t6 + qJD(4) * t8 + qJD(5) * t24, (t197 + (-t118 * t305 + t187 * t326) * t312 + m(6) * (-t180 * t96 + t300 * t343) + t319 * pkin(7) + t382) * qJD(3) + t11 * qJD(4) + t29 + t213, t11 * qJD(3) + ((-m(6) * t96 - t280) * pkin(4) + t217 + t382) * qJD(4) + t212, t211 + t30; -qJD(2) * t1 + qJD(4) * t10 - t26 - t263, -qJD(4) * t12 - t213 - t26, -t80 * qJD(4), ((-t311 - t374) * t187 - t323) * qJD(4) * pkin(3) + t204, -t334; -qJD(2) * t3 - qJD(3) * t10 - t262 - t67, qJD(3) * t12 - t212 - t67, -t204, 0, -t333; qJD(2) * t16 - t261 + t287, -t211 + t287, t334, t333, 0;];
Cq = t14;
