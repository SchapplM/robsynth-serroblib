% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR2_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR2_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:27:11
% EndTime: 2019-12-05 18:27:24
% DurationCPUTime: 4.73s
% Computational Cost: add. (14278->222), mult. (27105->311), div. (0->0), fcn. (31732->8), ass. (0->144)
t185 = sin(pkin(9));
t186 = cos(pkin(9));
t189 = sin(qJ(2));
t192 = cos(qJ(2));
t167 = -t185 * t189 + t186 * t192;
t168 = -t185 * t192 - t186 * t189;
t188 = sin(qJ(4));
t191 = cos(qJ(4));
t152 = t167 * t188 - t168 * t191;
t187 = sin(qJ(5));
t190 = cos(qJ(5));
t212 = t191 * t167 + t168 * t188;
t316 = -t152 * t187 + t190 * t212;
t335 = Ifges(6,5) * t316;
t106 = t152 * t190 + t187 * t212;
t348 = Ifges(6,6) * t106;
t272 = -qJ(3) - pkin(6);
t176 = t272 * t189;
t177 = t272 * t192;
t230 = t185 * t176 - t186 * t177;
t130 = -pkin(7) * t167 - t230;
t305 = t186 * t176 + t185 * t177;
t317 = pkin(7) * t168 + t305;
t329 = t188 * t130 + t191 * t317;
t344 = -pkin(8) * t152 + t329;
t244 = t190 * t344;
t78 = t191 * t130 - t188 * t317;
t54 = -pkin(8) * t212 + t78;
t246 = t187 * t54;
t365 = t244 + t246;
t368 = t365 * mrSges(6,2);
t37 = -t187 * t344 + t190 * t54;
t369 = t37 * mrSges(6,1);
t371 = t369 / 0.2e1 - t368 / 0.2e1 + t335 / 0.2e1 - t348 / 0.2e1;
t375 = 0.2e1 * t371;
t376 = t375 * qJD(5);
t310 = Ifges(5,5) * t212;
t323 = Ifges(5,6) * t152;
t343 = t369 - t368 + t335 - t348;
t346 = t329 * mrSges(5,2);
t362 = t78 * mrSges(5,1);
t374 = t343 + t310 - t323 - t346 + t362;
t373 = t310 / 0.2e1 - t323 / 0.2e1 - t346 / 0.2e1 + t362 / 0.2e1 + t371;
t225 = -pkin(2) * t192 - pkin(1);
t205 = -t167 * pkin(3) + t225;
t108 = -pkin(4) * t212 + t205;
t372 = m(6) * t108 - mrSges(6,1) * t316 + mrSges(6,2) * t106;
t256 = t106 * mrSges(6,3);
t367 = t187 * t365 + t190 * t37;
t293 = t316 / 0.2e1;
t327 = -t106 / 0.2e1;
t338 = -t316 / 0.2e1;
t351 = t106 / 0.2e1;
t366 = (t293 + t338) * Ifges(6,5) + (t327 + t351) * Ifges(6,6);
t361 = 0.2e1 * Ifges(6,4) * t316 + (Ifges(6,1) - Ifges(6,2)) * t106;
t334 = t316 * mrSges(6,2);
t347 = t106 * mrSges(6,1);
t219 = t334 + t347;
t360 = t108 * t219;
t270 = Ifges(6,4) * t106;
t358 = Ifges(6,1) * t316 - t270;
t271 = -t347 / 0.2e1 - t334 / 0.2e1;
t226 = t334 / 0.2e1;
t184 = pkin(2) * t186 + pkin(3);
t286 = pkin(2) * t185;
t159 = t191 * t184 - t188 * t286;
t158 = pkin(4) + t159;
t160 = t184 * t188 + t191 * t286;
t233 = t160 * t190;
t125 = t158 * t187 + t233;
t126 = -t159 * t187 - t233;
t345 = t125 + t126;
t339 = -Ifges(5,4) * t212 - (Ifges(5,1) - Ifges(5,2)) * t152;
t50 = Ifges(6,2) * t316 + t270;
t5 = t361 * t293 + t50 * t327 + t358 * t351 + t360;
t328 = t152 ^ 2;
t287 = -t152 / 0.2e1;
t282 = t152 * pkin(4);
t320 = t152 * mrSges(5,1);
t227 = -pkin(4) * t190 / 0.2e1;
t298 = m(6) * pkin(4);
t228 = t298 / 0.2e1;
t253 = t316 * mrSges(6,3);
t284 = pkin(4) * t187;
t318 = -t367 * t228 + t256 * t284 / 0.2e1 - t227 * t253 - t373;
t146 = t212 * mrSges(5,2);
t234 = t160 * t187;
t124 = t158 * t190 - t234;
t127 = t159 * t190 - t234;
t288 = t127 / 0.2e1;
t289 = -t124 / 0.2e1;
t299 = m(6) / 0.2e1;
t303 = (t345 * t365 - (-t124 + t127) * t37) * t299 + ((t288 + t289) * t316 + t345 * t327) * mrSges(6,3) + t373;
t302 = t167 ^ 2;
t301 = t168 ^ 2;
t300 = m(5) / 0.2e1;
t297 = -t365 / 0.2e1;
t285 = pkin(2) * t189;
t210 = -pkin(3) * t168 + t285;
t201 = t210 + t282;
t214 = t146 + t320;
t215 = -t168 * mrSges(4,1) + t167 * mrSges(4,2);
t249 = t152 * mrSges(5,3);
t1 = (t189 ^ 2 - t192 ^ 2) * Ifges(3,4) + pkin(1) * mrSges(3,2) * t192 + (t301 - t302) * Ifges(4,4) - t78 * t249 - t360 + (pkin(1) * mrSges(3,1) + (-Ifges(3,1) + Ifges(3,2)) * t192) * t189 + t358 * t327 + t50 * t351 + mrSges(4,1) * t167 * t285 + (-t210 * mrSges(5,2) + mrSges(5,3) * t78 + Ifges(5,4) * t152) * t152 + (t210 * mrSges(5,1) + t339) * t212 + (mrSges(4,2) * t285 + (Ifges(4,1) - Ifges(4,2)) * t167) * t168 + (-m(4) * t285 - t215) * t225 + t361 * t338 + (-m(5) * t210 - t214) * t205 - t372 * t201;
t259 = t1 * qJD(1);
t251 = t127 * mrSges(6,2);
t2 = t205 * t320 - t328 * Ifges(5,4) + (t205 * mrSges(5,2) - t339) * t212 + t5 + t372 * t282;
t242 = t2 * qJD(1);
t241 = t5 * qJD(1);
t10 = (t106 ^ 2 + t316 ^ 2) * mrSges(6,3) + (t212 ^ 2 + t328) * mrSges(5,3) + (t301 + t302) * mrSges(4,3) + m(6) * (-t106 * t365 - t316 * t37) + m(5) * (-t152 * t329 - t212 * t78) + m(4) * (t167 * t230 + t168 * t305);
t240 = qJD(1) * t10;
t29 = -t125 * mrSges(6,1) - t124 * mrSges(6,2);
t239 = qJD(5) * t29;
t238 = t106 * t190;
t237 = t316 * t187;
t229 = m(4) * pkin(2) / 0.2e1;
t195 = (-t106 * t124 + t125 * t316) * t299 + (-t152 * t159 + t160 * t212) * t300 + (t167 * t185 + t168 * t186) * t229;
t196 = t189 * t229 + t201 * t299 + t210 * t300;
t14 = t195 - t196 - t214 - t215 - t219;
t236 = t14 * qJD(1);
t15 = 0.2e1 * t287 * mrSges(5,1) + (-t238 / 0.2e1 + t237 / 0.2e1 + t287) * t298 - t146 + 0.2e1 * t271;
t235 = t15 * qJD(1);
t25 = 0.2e1 * t351 * mrSges(6,1) + 0.2e1 * t226;
t231 = t25 * qJD(1);
t211 = t227 + t289;
t123 = t126 * mrSges(6,1);
t200 = -t160 * mrSges(5,1) - t159 * mrSges(5,2) + t123 - t251;
t22 = -m(6) * (t124 * t126 + t125 * t127) - t200;
t4 = t303 + t318;
t207 = t4 * qJD(1) - t22 * qJD(2);
t7 = (t297 + t246 / 0.2e1 + t244 / 0.2e1) * mrSges(6,2) + t366;
t206 = t7 * qJD(1) + t29 * qJD(2);
t204 = (-t284 / 0.2e1 - t125 / 0.2e1) * mrSges(6,1);
t175 = (t187 * mrSges(6,1) + t190 * mrSges(6,2)) * pkin(4);
t27 = -t123 / 0.2e1 + t204 + (t288 + t211) * mrSges(6,2);
t9 = (t297 + t365 / 0.2e1) * mrSges(6,2) + t366;
t199 = -qJD(1) * t9 - qJD(2) * t27 + qJD(4) * t175;
t169 = t175 * qJD(5);
t28 = -t251 / 0.2e1 + t123 / 0.2e1 + t211 * mrSges(6,2) + t204;
t26 = t226 + t347 / 0.2e1 + t271;
t19 = t195 + t196;
t16 = (t152 + t237 - t238) * t228;
t3 = t303 - t318;
t6 = [-qJD(2) * t1 + qJD(3) * t10 + qJD(4) * t2 + qJD(5) * t5, t19 * qJD(3) + t3 * qJD(4) + t376 - t259 + (-t230 * mrSges(4,1) + (m(4) * (t185 * t305 - t186 * t230) + (-t167 * t186 + t168 * t185) * mrSges(4,3)) * pkin(2) - t305 * mrSges(4,2) - t124 * t253 - t159 * t212 * mrSges(5,3) + (-mrSges(3,1) * t192 + mrSges(3,2) * t189) * pkin(6) - t160 * t249 - t125 * t256 + 0.2e1 * (t124 * t37 + t125 * t365) * t299 + 0.2e1 * (t159 * t78 + t160 * t329) * t300 + Ifges(4,5) * t167 + Ifges(4,6) * t168 - Ifges(3,6) * t189 + Ifges(3,5) * t192 + t374) * qJD(2), qJD(2) * t19 + qJD(4) * t16 + qJD(5) * t26 + t240, t3 * qJD(2) + t16 * qJD(3) + t376 + t242 + ((m(6) * t367 + (-t106 * t187 - t190 * t316) * mrSges(6,3)) * pkin(4) + t374) * qJD(4), t26 * qJD(3) + t343 * qJD(5) + t241 + (qJD(2) + qJD(4)) * t375; qJD(3) * t14 + qJD(4) * t4 + qJD(5) * t7 + t259, -qJD(4) * t22 + t239, t236, ((t126 * t190 + t127 * t187) * t298 + t200) * qJD(4) + t28 * qJD(5) + t207, t28 * qJD(4) + t206 + t239; -qJD(2) * t14 - qJD(4) * t15 + qJD(5) * t25 - t240, -t236, 0, -t235, t231; -qJD(2) * t4 + qJD(3) * t15 + qJD(5) * t9 - t242, qJD(5) * t27 - t207, t235, -t169, -t169 - t199; -qJD(2) * t7 - qJD(3) * t25 - qJD(4) * t9 - t241, -qJD(4) * t27 - t206, -t231, t199, 0;];
Cq = t6;
