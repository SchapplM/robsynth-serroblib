% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-05 18:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR2_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR2_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR2_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR2_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:11:23
% EndTime: 2019-12-05 18:11:33
% DurationCPUTime: 3.58s
% Computational Cost: add. (13136->204), mult. (25230->284), div. (0->0), fcn. (29919->8), ass. (0->128)
t333 = qJD(3) + qJD(4);
t149 = sin(pkin(9));
t150 = cos(pkin(9));
t153 = sin(qJ(3));
t156 = cos(qJ(3));
t138 = -t153 * t149 + t150 * t156;
t139 = t149 * t156 + t153 * t150;
t152 = sin(qJ(4));
t155 = cos(qJ(4));
t119 = t138 * t152 + t139 * t155;
t151 = sin(qJ(5));
t154 = cos(qJ(5));
t170 = t155 * t138 - t139 * t152;
t271 = -t119 * t151 + t154 * t170;
t294 = Ifges(6,5) * t271;
t86 = t119 * t154 + t151 * t170;
t307 = Ifges(6,6) * t86;
t212 = pkin(6) + qJ(2);
t142 = t212 * t149;
t143 = t212 * t150;
t121 = -t153 * t142 + t143 * t156;
t101 = t138 * pkin(7) + t121;
t120 = -t156 * t142 - t143 * t153;
t163 = -pkin(7) * t139 + t120;
t257 = -t152 * t101 + t155 * t163;
t273 = -pkin(8) * t119 + t257;
t61 = t155 * t101 + t152 * t163;
t45 = pkin(8) * t170 + t61;
t311 = -t151 * t45 + t154 * t273;
t317 = t311 * mrSges(6,2);
t32 = t151 * t273 + t154 * t45;
t324 = t32 * mrSges(6,1);
t326 = -t324 / 0.2e1 + t294 / 0.2e1 - t307 / 0.2e1 - t317 / 0.2e1;
t331 = 0.2e1 * t326;
t332 = t331 * qJD(5);
t302 = -t324 + t294 - t307 - t317;
t303 = -t61 * mrSges(5,1) - t257 * mrSges(5,2) + Ifges(5,5) * t170 - Ifges(5,6) * t119;
t330 = t302 + t303;
t178 = -pkin(2) * t150 - pkin(1);
t122 = -pkin(3) * t138 + t178;
t91 = -pkin(4) * t170 + t122;
t329 = m(6) * t91 - mrSges(6,1) * t271 + mrSges(6,2) * t86;
t244 = t86 / 0.2e1;
t284 = -t86 / 0.2e1;
t299 = t271 / 0.2e1;
t319 = (t299 - t271 / 0.2e1) * Ifges(6,5) + (t284 + t244) * Ifges(6,6);
t328 = qJD(1) * t319;
t327 = qJD(5) * t319;
t279 = t119 * mrSges(5,1);
t174 = t170 * mrSges(5,2) + t279;
t210 = Ifges(5,4) * t119;
t240 = t119 / 0.2e1;
t241 = -t119 / 0.2e1;
t269 = t170 / 0.2e1;
t306 = t86 * mrSges(6,1);
t177 = t271 * mrSges(6,2) + t306;
t237 = Ifges(6,4) * t86;
t5 = t91 * t177 + (0.2e1 * Ifges(6,4) * t271 + (Ifges(6,1) - Ifges(6,2)) * t86) * t299 + (Ifges(6,2) * t271 + t237) * t284 + (Ifges(6,1) * t271 - t237) * t244;
t323 = (Ifges(5,2) * t170 + t210) * t241 + (Ifges(5,1) * t170 - t210) * t240 + (0.2e1 * Ifges(5,4) * t170 + (Ifges(5,1) - Ifges(5,2)) * t119) * t269 + t122 * t174 + t5;
t320 = t151 * t311 - t154 * t32;
t298 = t306 / 0.2e1;
t196 = t151 * t86;
t146 = pkin(3) * t155 + pkin(4);
t185 = t152 * t154;
t130 = pkin(3) * t185 + t146 * t151;
t135 = (-t151 * t155 - t185) * pkin(3);
t305 = t130 + t135;
t19 = 0.2e1 * t299 * mrSges(6,2) + 0.2e1 * t298;
t20 = t298 - t306 / 0.2e1;
t254 = -pkin(4) / 0.2e1;
t229 = t139 * pkin(3);
t297 = m(5) * t229;
t230 = t119 * pkin(4);
t253 = m(6) * pkin(4);
t276 = -t320 * t253 / 0.2e1 - t326;
t132 = t135 * mrSges(6,1);
t186 = t151 * t152;
t136 = (t154 * t155 - t186) * pkin(3);
t199 = t136 * mrSges(6,2);
t258 = (mrSges(5,1) * t152 + mrSges(5,2) * t155) * pkin(3) - t132 + t199;
t255 = m(6) / 0.2e1;
t129 = -pkin(3) * t186 + t146 * t154;
t239 = -t129 / 0.2e1;
t238 = t136 / 0.2e1;
t214 = t271 * mrSges(6,3);
t133 = t139 * mrSges(4,1);
t93 = t229 + t230;
t1 = t178 * t133 + (pkin(3) * (-mrSges(5,1) * t170 + mrSges(5,2) * t119) - Ifges(4,4) * t139) * t139 + t122 * t297 + (Ifges(4,4) * t138 + t178 * mrSges(4,2) + (Ifges(4,1) - Ifges(4,2)) * t139) * t138 + t323 + t329 * t93;
t205 = t1 * qJD(1);
t201 = t129 * t86;
t200 = t130 * t271;
t198 = t151 * mrSges(6,1);
t195 = t151 * t271;
t193 = t154 * t86;
t2 = t329 * t230 + t323;
t192 = t2 * qJD(1);
t191 = t5 * qJD(1);
t10 = (t271 ^ 2 + t86 ^ 2) * mrSges(6,3) + (t119 ^ 2 + t170 ^ 2) * mrSges(5,3) + (t138 ^ 2 + t139 ^ 2) * mrSges(4,3) + m(6) * (t271 * t32 - t311 * t86) + m(5) * (-t119 * t257 + t170 * t61) + m(4) * (-t120 * t139 + t121 * t138) + (m(3) * qJ(2) + mrSges(3,3)) * (t149 ^ 2 + t150 ^ 2);
t190 = qJD(1) * t10;
t125 = t130 * mrSges(6,1);
t47 = t129 * mrSges(6,2) + t125;
t189 = qJD(5) * t47;
t13 = t279 + 0.2e1 * t269 * mrSges(5,2) + (t240 - t195 / 0.2e1 + t193 / 0.2e1) * t253 + t19;
t188 = t13 * qJD(1);
t161 = m(5) * (-t119 * t155 + t152 * t170) * pkin(3);
t172 = t297 / 0.2e1;
t15 = t172 + t138 * mrSges(4,2) + t133 + 0.2e1 * (t93 / 0.4e1 + t201 / 0.4e1 - t200 / 0.4e1) * m(6) - t161 / 0.2e1 + t174 + t177;
t187 = t15 * qJD(1);
t184 = t19 * qJD(1);
t182 = t253 / 0.2e1;
t181 = t154 * t254;
t169 = t239 + t181;
t157 = (t305 * t311 + (-t129 + t136) * t32) * t255 + t326;
t158 = (t238 + t239) * t271 + t305 * t284;
t3 = (t240 + t241) * Ifges(5,6) + (t269 - t170 / 0.2e1) * Ifges(5,5) + ((t196 / 0.2e1 + t154 * t299) * pkin(4) + t158) * mrSges(6,3) + t157 + t276;
t46 = -m(6) * (t129 * t135 + t130 * t136) + t258;
t167 = t3 * qJD(1) - t46 * qJD(3);
t166 = -t47 * qJD(3) + t328;
t165 = -t125 / 0.2e1 + t198 * t254;
t141 = (t154 * mrSges(6,2) + t198) * pkin(4);
t39 = -t132 / 0.2e1 + (t238 + t169) * mrSges(6,2) + t165;
t160 = -qJD(3) * t39 + qJD(4) * t141 - t328;
t140 = t141 * qJD(5);
t40 = -t199 / 0.2e1 + t132 / 0.2e1 + t169 * mrSges(6,2) + t165;
t18 = t172 + t161 / 0.2e1 + (t93 + t200 - t201) * t255;
t14 = -t279 / 0.2e1 + (-t193 + t195) * t182 + (t182 + mrSges(5,1) / 0.2e1) * t119 + t20;
t4 = t181 * t214 + t157 + (t196 * t254 + t158) * mrSges(6,3) - t276 + t303;
t6 = [qJD(2) * t10 + qJD(3) * t1 + qJD(4) * t2 + qJD(5) * t5, qJD(3) * t18 + qJD(4) * t14 + qJD(5) * t20 + t190, t18 * qJD(2) + t4 * qJD(4) + t332 + t205 + (-t130 * t86 * mrSges(6,3) - t129 * t214 + m(6) * (-t129 * t32 + t130 * t311) + Ifges(4,5) * t138 - Ifges(4,6) * t139 - t120 * mrSges(4,2) - t121 * mrSges(4,1) + (m(5) * (t152 * t257 - t155 * t61) + (-t119 * t152 - t155 * t170) * mrSges(5,3)) * pkin(3) + t330) * qJD(3), t14 * qJD(2) + t4 * qJD(3) + t332 + t192 + ((m(6) * t320 + (-t154 * t271 - t196) * mrSges(6,3)) * pkin(4) + t330) * qJD(4), t20 * qJD(2) + t302 * qJD(5) + t333 * t331 + t191; qJD(3) * t15 + qJD(4) * t13 + qJD(5) * t19 - t190, 0, t187, t188, t184; -qJD(2) * t15 + qJD(4) * t3 - t205 + t327, -t187, -qJD(4) * t46 - t189, ((t135 * t154 + t136 * t151) * t253 - t258) * qJD(4) + t40 * qJD(5) + t167, t40 * qJD(4) + t166 - t189; -qJD(2) * t13 - qJD(3) * t3 - t192 + t327, -t188, qJD(5) * t39 - t167, -t140, -t140 - t160; -qJD(2) * t19 - t333 * t319 - t191, -t184, -qJD(4) * t39 - t166, t160, 0;];
Cq = t6;
