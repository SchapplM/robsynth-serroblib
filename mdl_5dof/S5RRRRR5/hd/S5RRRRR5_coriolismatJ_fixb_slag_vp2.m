% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2022-01-20 12:02
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRRR5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR5_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR5_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR5_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:01:45
% EndTime: 2022-01-20 12:01:54
% DurationCPUTime: 4.01s
% Computational Cost: add. (9984->318), mult. (20155->404), div. (0->0), fcn. (18011->8), ass. (0->195)
t225 = cos(qJ(4));
t219 = t225 ^ 2;
t430 = -m(6) / 0.2e1;
t427 = pkin(4) * t430;
t226 = cos(qJ(3));
t221 = sin(qJ(4));
t218 = t221 ^ 2;
t287 = t218 + t219;
t383 = t287 * t226;
t227 = cos(qJ(2));
t344 = pkin(1) * t227;
t214 = pkin(2) + t344;
t222 = sin(qJ(3));
t223 = sin(qJ(2));
t291 = t223 * t226;
t168 = pkin(1) * t291 + t214 * t222;
t162 = pkin(8) + t168;
t330 = pkin(9) + t162;
t137 = t330 * t221;
t138 = t330 * t225;
t220 = sin(qJ(5));
t224 = cos(qJ(5));
t265 = -t224 * t137 - t138 * t220;
t66 = -t137 * t220 + t138 * t224;
t426 = -t66 * mrSges(6,1) - t265 * mrSges(6,2);
t343 = pkin(2) * t222;
t212 = pkin(8) + t343;
t329 = pkin(9) + t212;
t185 = t329 * t221;
t186 = t329 * t225;
t117 = -t185 * t220 + t186 * t224;
t264 = -t224 * t185 - t186 * t220;
t425 = -t117 * mrSges(6,1) - t264 * mrSges(6,2);
t359 = -pkin(9) - pkin(8);
t205 = t359 * t221;
t206 = t359 * t225;
t142 = t205 * t220 - t206 * t224;
t263 = t224 * t205 + t206 * t220;
t424 = -t142 * mrSges(6,1) - t263 * mrSges(6,2);
t207 = t222 * t223 * pkin(1);
t167 = t214 * t226 - t207;
t306 = t225 * mrSges(5,2);
t307 = t221 * mrSges(5,1);
t254 = t307 / 0.2e1 + t306 / 0.2e1;
t189 = -t220 * t225 - t224 * t221;
t84 = t189 * t167;
t417 = t84 / 0.2e1;
t418 = -mrSges(6,2) / 0.2e1;
t188 = -t220 * t221 + t224 * t225;
t85 = t188 * t167;
t394 = mrSges(6,1) * t417 + t85 * t418;
t421 = t254 * t167 + (t220 * t85 + t224 * t84) * t427 - t394;
t172 = t226 * t344 - t207;
t99 = t188 * t172;
t416 = -t99 / 0.2e1;
t98 = t189 * t172;
t393 = t98 * mrSges(6,1) / 0.2e1 + mrSges(6,2) * t416;
t420 = t254 * t172 + (t220 * t99 + t224 * t98) * t427 - t393;
t342 = pkin(2) * t226;
t156 = t189 * t342;
t157 = t188 * t342;
t415 = t156 / 0.2e1;
t388 = mrSges(6,1) * t415 + t157 * t418;
t419 = t254 * t342 + (t156 * t224 + t157 * t220) * t427 - t388;
t288 = Ifges(6,5) * t188 + Ifges(6,6) * t189;
t22 = t288 + t426;
t412 = t22 * qJD(5);
t31 = t288 + t425;
t411 = t31 * qJD(5);
t36 = t288 + t424;
t410 = t36 * qJD(5);
t404 = -t142 / 0.2e1;
t405 = -t117 / 0.2e1;
t409 = t405 + t404;
t406 = -t66 / 0.2e1;
t408 = t406 + t404;
t407 = t406 + t405;
t171 = (t222 * t227 + t291) * pkin(1);
t339 = pkin(4) * t225;
t215 = -pkin(3) - t339;
t122 = -mrSges(6,1) * t188 - mrSges(6,2) * t189;
t195 = -mrSges(5,1) * t225 + mrSges(5,2) * t221;
t260 = t195 / 0.2e1 + t122 / 0.2e1 - mrSges(4,1) / 0.2e1;
t267 = t287 * t172;
t400 = -m(5) * (-pkin(3) * t171 + pkin(8) * t267) / 0.2e1 + (t142 * t99 + t171 * t215 + t263 * t98) * t430 - t260 * t171;
t375 = Ifges(6,4) * t188 ^ 2 + (-Ifges(6,4) * t189 + (Ifges(6,2) - Ifges(6,1)) * t188) * t189;
t340 = pkin(4) * t221;
t374 = t122 * t340 + Ifges(5,4) * t219 + (-Ifges(5,4) * t221 + (-Ifges(5,2) + Ifges(5,1)) * t225) * t221;
t371 = m(6) / 0.2e1;
t396 = 0.2e1 * t371;
t372 = m(5) / 0.2e1;
t395 = 0.2e1 * t372;
t377 = t287 * mrSges(5,3);
t196 = t306 + t307;
t213 = -pkin(3) - t342;
t169 = t213 * t196;
t311 = t188 * mrSges(6,3);
t78 = t264 * t311;
t382 = t169 / 0.2e1 + t78 / 0.2e1;
t161 = -pkin(3) - t167;
t134 = t161 * t196;
t54 = t265 * t311;
t381 = t134 / 0.2e1 + t54 / 0.2e1;
t193 = t215 - t342;
t258 = -mrSges(6,1) * t189 + mrSges(6,2) * t188;
t105 = t193 * t258;
t309 = t189 * mrSges(6,3);
t79 = t117 * t309;
t380 = t105 / 0.2e1 + t79 / 0.2e1;
t55 = t66 * t309;
t149 = t161 - t339;
t72 = t149 * t258;
t379 = t72 / 0.2e1 + t55 / 0.2e1;
t259 = -mrSges(4,2) + t377;
t279 = -mrSges(4,1) + t122 + t195;
t378 = t259 * t167 + (t85 * t188 + t84 * t189) * mrSges(6,3) + t279 * t168;
t376 = (t156 * t189 + t157 * t188) * mrSges(6,3);
t373 = t259 * t172 + t279 * t171 + (-mrSges(3,1) * t223 - mrSges(3,2) * t227) * pkin(1);
t368 = -mrSges(4,2) / 0.2e1;
t364 = -t265 / 0.2e1;
t356 = -t264 / 0.2e1;
t351 = -t263 / 0.2e1;
t349 = -t167 / 0.2e1;
t347 = t172 / 0.2e1;
t345 = m(6) * t222;
t341 = pkin(3) * t196;
t322 = pkin(4) * qJD(4);
t310 = t188 * t99;
t308 = t189 * t98;
t136 = t149 * t340;
t240 = -t375 - t374;
t13 = (t188 * t265 + t66 * t189) * mrSges(6,3) + t240 - t134 - t72 - t55 - t54 - m(6) * t136;
t305 = t13 * qJD(1);
t268 = t287 * t167;
t16 = m(6) * (t149 * t168 + t265 * t84 + t66 * t85) + m(5) * (t161 * t168 + t162 * t268) + t378;
t302 = t16 * qJD(1);
t17 = (t308 + t310) * mrSges(6,3) + m(6) * (t149 * t171 + t265 * t98 + t66 * t99) + m(5) * (t161 * t171 + t162 * t267) + m(4) * (-t167 * t171 + t168 * t172) + t373;
t301 = t17 * qJD(1);
t19 = -t375 - t72;
t296 = t19 * qJD(1);
t286 = mrSges(6,3) * t308;
t285 = mrSges(6,3) * t310;
t107 = t142 * t309;
t284 = t224 * t311;
t273 = t415 + t417;
t272 = t157 / 0.2e1 + t85 / 0.2e1;
t271 = t349 + t347;
t266 = t287 * t212;
t184 = t193 * t340;
t18 = -t105 + (t117 * t189 + t188 * t264) * mrSges(6,3) + t240 - t169 - t79 - t78 - m(6) * t184;
t239 = t375 + t379 + t380;
t231 = t239 + (t136 + t184) * t371 + (t407 * t189 + (t356 + t364) * t188) * mrSges(6,3) + t374 + t381 + t382;
t2 = t231 + t420;
t257 = t2 * qJD(1) - t18 * qJD(2);
t21 = -t105 - t375;
t236 = t407 * t309 + t239;
t9 = t236 - t393;
t256 = t9 * qJD(1) - t21 * qJD(2);
t237 = t279 * t222 + t259 * t226;
t30 = m(6) * (t117 * t157 + t156 * t264) + t376 + (t193 * t345 + m(5) * (t383 * t212 + t213 * t222) + t237) * pkin(2);
t228 = (t260 * t222 + t226 * t368) * pkin(2) + t260 * t168 + (t168 * t213 + t167 * t266 + (t161 * t222 + t162 * t383) * pkin(2)) * t372 + (t117 * t85 + t149 * t343 + t156 * t265 + t157 * t66 + t168 * t193 + t264 * t84) * t371;
t4 = t271 * mrSges(4,2) + ((-t98 / 0.2e1 + t273) * t189 + (t416 + t272) * t188) * mrSges(6,3) + t228 + (t342 / 0.2e1 - t271) * t377 + t400;
t255 = t4 * qJD(1) + t30 * qJD(2);
t248 = t220 * pkin(4) * t309 + Ifges(5,5) * t225 - Ifges(5,6) * t221 + t288;
t113 = t215 * t258;
t247 = t107 / 0.2e1 + t113 / 0.2e1 + t375;
t106 = t263 * t311;
t194 = t215 * t340;
t20 = -t113 - t106 - t107 + (t142 * t189 + t188 * t263) * mrSges(6,3) + t240 - m(6) * t194 + t341;
t233 = t106 / 0.2e1 - t341 / 0.2e1 + t247 + t374;
t230 = t233 + (t408 * t189 + (t351 + t364) * t188) * mrSges(6,3) + (t136 + t194) * t371 + t379 + t381;
t6 = t230 + t421;
t229 = t233 + (t409 * t189 + (t351 + t356) * t188) * mrSges(6,3) + (t184 + t194) * t371 + t380 + t382;
t8 = t229 + t419;
t243 = t6 * qJD(1) + t8 * qJD(2) - t20 * qJD(3);
t235 = t408 * t309 + t247 + t379;
t11 = t235 - t394;
t234 = t409 * t309 + t247 + t380;
t14 = t234 - t388;
t25 = -t113 - t375;
t242 = t11 * qJD(1) + t14 * qJD(2) - t25 * qJD(3);
t192 = (mrSges(6,1) * t220 + mrSges(6,2) * t224) * pkin(4);
t24 = (t364 + t265 / 0.2e1) * mrSges(6,2) + (t406 + t66 / 0.2e1) * mrSges(6,1);
t33 = (t356 + t264 / 0.2e1) * mrSges(6,2) + (t405 + t117 / 0.2e1) * mrSges(6,1);
t38 = (t351 + t263 / 0.2e1) * mrSges(6,2) + (t404 + t142 / 0.2e1) * mrSges(6,1);
t238 = -qJD(1) * t24 - qJD(2) * t33 - qJD(3) * t38 + qJD(4) * t192;
t187 = t192 * qJD(5);
t15 = t234 + t388;
t12 = t235 + t394;
t10 = t236 + t393;
t7 = t229 - t419;
t5 = t230 - t421;
t3 = t172 * t368 + t285 / 0.2e1 + t286 / 0.2e1 + mrSges(4,2) * t349 + (t272 * t188 + t273 * t189) * mrSges(6,3) + t228 + t347 * t377 + (t167 + t342) * mrSges(5,3) * (t218 / 0.2e1 + t219 / 0.2e1) - t400;
t1 = t231 - t420;
t23 = [qJD(2) * t17 + qJD(3) * t16 - qJD(4) * t13 - qJD(5) * t19, t3 * qJD(3) + t1 * qJD(4) + t10 * qJD(5) + t301 + (t285 + t286 + (t117 * t99 + t171 * t193 + t264 * t98) * t396 + (t171 * t213 + t172 * t266) * t395 + m(4) * (-t171 * t226 + t172 * t222) * pkin(2) + t373) * qJD(2), t3 * qJD(2) + t5 * qJD(4) + t12 * qJD(5) + t302 + ((t142 * t85 + t168 * t215 + t263 * t84) * t396 + (-pkin(3) * t168 + pkin(8) * t268) * t395 + t378) * qJD(3), -t305 + t1 * qJD(2) + t5 * qJD(3) + (t195 * t162 + t248 + t426) * qJD(4) + t412 + (-t284 + m(6) * (t220 * t265 - t224 * t66)) * t322, t10 * qJD(2) + t12 * qJD(3) + t22 * qJD(4) - t296 + t412; qJD(3) * t4 + qJD(4) * t2 + qJD(5) * t9 - t301, qJD(3) * t30 - qJD(4) * t18 - qJD(5) * t21, t7 * qJD(4) + t15 * qJD(5) + t255 + (m(6) * (t142 * t157 + t156 * t263) + t376 + (t215 * t345 + m(5) * (-pkin(3) * t222 + t383 * pkin(8)) + t237) * pkin(2)) * qJD(3), t7 * qJD(3) + (t195 * t212 + t248 + t425) * qJD(4) + t411 + (-t284 + m(6) * (-t117 * t224 + t220 * t264)) * t322 + t257, t15 * qJD(3) + t31 * qJD(4) + t256 + t411; -qJD(2) * t4 + qJD(4) * t6 + qJD(5) * t11 - t302, qJD(4) * t8 + qJD(5) * t14 - t255, -qJD(4) * t20 - qJD(5) * t25, (t195 * pkin(8) + t248 + t424) * qJD(4) + t410 + (-t284 + m(6) * (-t142 * t224 + t220 * t263)) * t322 + t243, t36 * qJD(4) + t242 + t410; -qJD(2) * t2 - qJD(3) * t6 + qJD(5) * t24 + t305, -qJD(3) * t8 + qJD(5) * t33 - t257, qJD(5) * t38 - t243, -t187, -t187 - t238; -qJD(2) * t9 - qJD(3) * t11 - qJD(4) * t24 + t296, -qJD(3) * t14 - qJD(4) * t33 - t256, -qJD(4) * t38 - t242, t238, 0;];
Cq = t23;
