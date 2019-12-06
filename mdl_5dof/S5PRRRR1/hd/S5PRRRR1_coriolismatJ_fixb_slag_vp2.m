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
% Datum: 2019-12-05 17:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:02:54
% EndTime: 2019-12-05 17:03:00
% DurationCPUTime: 2.43s
% Computational Cost: add. (3280->266), mult. (9057->409), div. (0->0), fcn. (9410->8), ass. (0->168)
t174 = sin(qJ(5));
t286 = t174 / 0.2e1;
t178 = cos(qJ(5));
t284 = -t178 / 0.2e1;
t272 = Ifges(6,4) * t174;
t156 = Ifges(6,1) * t178 - t272;
t283 = t178 / 0.2e1;
t168 = Ifges(6,4) * t178;
t301 = -Ifges(6,2) * t174 + t168;
t153 = Ifges(6,2) * t178 + t272;
t155 = Ifges(6,1) * t174 + t168;
t309 = t153 * t286 + t155 * t284;
t48 = t156 * t286 + t301 * t283 - t309;
t170 = t174 ^ 2;
t172 = t178 ^ 2;
t227 = t170 + t172;
t304 = t227 * mrSges(6,3);
t313 = mrSges(5,2) - t304;
t295 = m(6) / 0.2e1;
t175 = sin(qJ(4));
t177 = sin(qJ(2));
t179 = cos(qJ(4));
t180 = cos(qJ(3));
t231 = t179 * t180;
t176 = sin(qJ(3));
t234 = t176 * t177;
t137 = -t175 * t234 + t177 * t231;
t181 = cos(qJ(2));
t109 = t137 * t178 - t181 * t174;
t312 = -t109 / 0.2e1;
t150 = t175 * t180 + t176 * t179;
t288 = -t150 / 0.2e1;
t257 = t178 * mrSges(6,1);
t263 = t174 * mrSges(6,2);
t211 = t257 - t263;
t311 = t211 + mrSges(5,1);
t107 = -t137 * t174 - t178 * t181;
t248 = t109 * t174;
t280 = pkin(2) * t176;
t310 = (-m(5) * t181 / 0.2e1 + (t107 * t178 + t248) * t295) * t280;
t308 = qJD(3) + qJD(4);
t307 = -t156 / 0.4e1 + t153 / 0.4e1;
t305 = Ifges(5,1) - Ifges(5,2);
t303 = t181 * t177;
t167 = Ifges(6,5) * t178;
t270 = Ifges(6,6) * t174;
t302 = t167 - t270;
t149 = t175 * t176 - t231;
t239 = t150 * t174;
t97 = -mrSges(6,2) * t149 - mrSges(6,3) * t239;
t238 = t150 * t178;
t99 = mrSges(6,1) * t149 - mrSges(6,3) * t238;
t300 = t174 * t97 + t178 * t99;
t135 = t150 * t177;
t256 = t178 * mrSges(6,2);
t264 = t174 * mrSges(6,1);
t199 = t256 / 0.2e1 + t264 / 0.2e1;
t210 = t256 + t264;
t222 = -mrSges(6,3) * t174 / 0.2e1;
t291 = t135 / 0.2e1;
t292 = mrSges(6,3) / 0.2e1;
t185 = t109 * t222 + t199 * t135 + t210 * t291 + t248 * t292;
t247 = t109 * t178;
t250 = t107 * t174;
t197 = t137 - t247 + t250;
t20 = m(6) * t197 * t135;
t229 = t20 * qJD(1);
t297 = t185 * qJD(5) + t229;
t102 = mrSges(5,1) * t150 - mrSges(5,2) * t149;
t86 = t210 * t150;
t241 = t149 * t174;
t96 = -mrSges(6,2) * t150 + mrSges(6,3) * t241;
t240 = t149 * t178;
t98 = mrSges(6,1) * t150 + mrSges(6,3) * t240;
t296 = t181 * t102 / 0.2e1 - t137 * t86 / 0.2e1 + t96 * t312 - t107 * t98 / 0.2e1;
t294 = mrSges(6,1) / 0.2e1;
t293 = -mrSges(6,2) / 0.2e1;
t289 = t149 / 0.2e1;
t287 = t150 / 0.2e1;
t285 = t175 / 0.2e1;
t281 = pkin(2) * t175;
t279 = pkin(2) * t179;
t278 = pkin(2) * t180;
t277 = mrSges(5,3) * t149;
t276 = mrSges(5,3) * t150;
t273 = Ifges(5,4) * t150;
t269 = Ifges(6,3) * t150;
t138 = t149 * t181;
t265 = t138 * mrSges(5,2);
t64 = Ifges(6,6) * t149 + t150 * t301;
t262 = t174 * t64;
t260 = t174 * t98;
t259 = t174 * t99;
t66 = Ifges(6,5) * t149 + t150 * t156;
t255 = t178 * t66;
t254 = t178 * t96;
t253 = t178 * t97;
t108 = t138 * t174 + t177 * t178;
t249 = t108 * t174;
t110 = -t138 * t178 + t174 * t177;
t246 = t110 * t178;
t136 = t150 * t181;
t245 = t135 * t136;
t244 = t135 * t175;
t243 = t136 * t179;
t242 = t149 * t302;
t226 = t176 ^ 2 + t180 ^ 2;
t17 = m(6) * (t107 * t108 + t109 * t110 + t245) + m(4) * (-0.1e1 + t226) * t303 + (-t137 * t138 + t245 - t303) * m(5);
t237 = t17 * qJD(1);
t225 = pkin(2) * t244;
t228 = t227 * t225;
t125 = t137 * t279;
t219 = -t240 / 0.2e1;
t220 = t241 / 0.2e1;
t221 = Ifges(6,5) * t219 + Ifges(6,6) * t220 + t269 / 0.2e1;
t215 = t311 * t175;
t213 = t108 * t222 + t246 * t292 + t265 / 0.2e1 - t311 * t136 / 0.2e1;
t207 = Ifges(6,5) * t174 + Ifges(6,6) * t178;
t103 = mrSges(5,1) * t149 + mrSges(5,2) * t150;
t143 = Ifges(5,4) * t149;
t182 = pkin(2) ^ 2;
t202 = t174 * t96 + t178 * t98 + t102;
t63 = Ifges(6,6) * t150 - t149 * t301;
t65 = Ifges(6,5) * t150 - t156 * t149;
t1 = t64 * t220 - t63 * t239 / 0.2e1 + t66 * t219 + t65 * t238 / 0.2e1 + t269 * t289 + t273 * t288 + (Ifges(4,4) * t180 - t202 * pkin(2)) * t180 + (-Ifges(4,4) * t176 + (t103 + t300) * pkin(2) + (Ifges(4,1) - Ifges(4,2) + (-m(6) * t227 - m(5)) * t182) * t180) * t176 + (t150 * t302 - t273) * t287 + (-t302 * t289 + t143 + (-Ifges(5,1) + Ifges(6,3)) * t287 + (-Ifges(5,2) + t305) * t288) * t149;
t183 = 0.2e1 * (m(5) * (-t138 * t175 - t243) / 0.4e1 + m(6) * (-t243 + (t246 - t249) * t175) / 0.4e1) * pkin(2) + t213;
t85 = t210 * t149;
t191 = t253 / 0.2e1 - t259 / 0.2e1 + t85 / 0.2e1;
t2 = t191 * t135 + t183 + t296 - t310;
t206 = -t2 * qJD(1) + t1 * qJD(2);
t201 = t167 / 0.2e1 - t270 / 0.2e1;
t190 = t201 * t149;
t4 = t202 * t278 + (t255 / 0.2e1 - t262 / 0.2e1 - t143 + t190) * t149 + (t65 * t284 + t63 * t286 + (Ifges(5,4) - t201) * t150 + (-Ifges(6,3) + t305) * t149) * t150;
t189 = -t135 * t253 / 0.2e1 - t296 + (t259 - t85) * t291;
t5 = -t265 / 0.2e1 + (mrSges(5,1) / 0.2e1 + t211 / 0.2e1) * t136 + (-t246 / 0.2e1 + t249 / 0.2e1) * mrSges(6,3) + t189;
t205 = t5 * qJD(1) - t4 * qJD(2);
t198 = pkin(2) * t210;
t192 = t179 * t198;
t35 = t192 - t48;
t195 = t179 * t211;
t87 = t150 * t153;
t88 = t150 * t155;
t7 = -t242 / 0.4e1 + (-t66 / 0.4e1 + t87 / 0.4e1) * t178 + (t88 / 0.4e1 + t64 / 0.4e1) * t174 + ((t176 * t294 + t99 * t285) * t178 + (t176 * t293 + t97 * t285) * t174) * pkin(2) + (t307 * t178 + (t155 / 0.4e1 + t301 / 0.4e1) * t174 + (t195 / 0.2e1 + (t172 / 0.2e1 + t170 / 0.2e1) * t175 * mrSges(6,3)) * pkin(2)) * t150 + t221;
t204 = -t7 * qJD(2) - t35 * qJD(3);
t184 = (t211 * t291 + (-t247 / 0.2e1 + t250 / 0.2e1) * mrSges(6,3)) * t150 + t107 * t97 / 0.2e1 + t99 * t312;
t200 = t108 * t294 + t110 * t293;
t10 = t184 - t200;
t14 = (t207 * t289 + t64 * t283 - t88 * t284 + (t66 - t87) * t286) * t150 + (t253 - t259) * t278;
t203 = t10 * qJD(1) - t14 * qJD(2);
t194 = t313 * t135 - t311 * t137;
t193 = -Ifges(5,6) * t150 + t207 * t287 + t63 * t283 + t65 * t286 + (-Ifges(5,5) + t309) * t149;
t16 = (t191 * t179 + (t254 / 0.2e1 - t260 / 0.2e1 + t86 / 0.2e1) * t175) * pkin(2);
t21 = ((-t197 * t179 + t244) * pkin(2) - t228) * t295;
t34 = -pkin(2) * t215 + (m(6) * (-0.1e1 + t227) * t182 * t175 - t313 * pkin(2)) * t179;
t188 = t21 * qJD(1) + t16 * qJD(2) + t34 * qJD(3);
t186 = t242 / 0.4e1 - t262 / 0.4e1 + t255 / 0.4e1 - t174 * t88 / 0.4e1 - t178 * t87 / 0.4e1 - t307 * t238 - (t301 + t155) * t239 / 0.4e1;
t13 = -t269 / 0.2e1 + t190 + t186;
t29 = (-t155 / 0.2e1 - t301 / 0.2e1) * t178 + (-t156 / 0.2e1 + t153 / 0.2e1) * t174;
t187 = -qJD(2) * t13 + qJD(3) * t29 - qJD(4) * t48;
t12 = t186 + t221;
t30 = -t192 / 0.2e1 - t199 * t279 + t48;
t15 = t194 + t21;
t11 = t184 + t200;
t9 = t16 + t193;
t8 = t12 + pkin(2) * t195 * t288 + (-t263 / 0.2e1 + t257 / 0.2e1) * t280 - (t150 * t304 + t300) * t281 / 0.2e1;
t6 = t189 + t213;
t3 = t183 + t189 + t310 + (-t176 * mrSges(4,1) - t180 * mrSges(4,2)) * t181;
t18 = [t17 * qJD(2) + t308 * t20, t3 * qJD(3) + t6 * qJD(4) + t11 * qJD(5) + t237 + (t138 * t277 + m(6) * (-t108 * t178 - t110 * t174) * t278 + t110 * t97 + t108 * t99 + (t276 + t86) * t136 + (mrSges(4,3) * t226 - mrSges(3,2)) * t181 + (t176 * mrSges(4,2) - mrSges(3,1) + t103 + (-m(5) * pkin(2) - mrSges(4,1)) * t180) * t177) * qJD(2), t3 * qJD(2) + t15 * qJD(4) + (-t177 * t180 * mrSges(4,1) + mrSges(4,2) * t234 + t194 + 0.2e1 * (-t125 - t228) * t295 + m(5) * (-t125 - t225)) * qJD(3) + t297, t6 * qJD(2) + t15 * qJD(3) + qJD(4) * t194 + t297, t11 * qJD(2) + (-mrSges(6,1) * t109 - mrSges(6,2) * t107) * qJD(5) + t308 * t185; -qJD(3) * t2 + qJD(4) * t5 + qJD(5) * t10 - t237, qJD(3) * t1 - qJD(4) * t4 - qJD(5) * t14, t9 * qJD(4) + t8 * qJD(5) + t206 + (Ifges(4,5) * t180 - Ifges(4,6) * t176 + t193 + ((t85 + t277) * t179 + (t254 - t260 - t276) * t175) * pkin(2)) * qJD(3), t9 * qJD(3) + qJD(4) * t193 + t12 * qJD(5) + t205, t8 * qJD(3) + t12 * qJD(4) + (-t150 * t207 + t180 * t198) * qJD(5) + t203; qJD(2) * t2 + qJD(4) * t21 - t229, qJD(4) * t16 - qJD(5) * t7 - t206, qJD(4) * t34 - qJD(5) * t35, t30 * qJD(5) + (-t179 * t313 - t215) * qJD(4) * pkin(2) + t188, t30 * qJD(4) + (-t211 * t281 + t302) * qJD(5) + t204; -qJD(2) * t5 - qJD(3) * t21 - t229, -qJD(3) * t16 + qJD(5) * t13 - t205, -qJD(5) * t29 - t188, t48 * qJD(5), qJD(5) * t302 - t187; -t10 * qJD(2), qJD(3) * t7 - qJD(4) * t13 - t203, qJD(4) * t29 - t204, t187, 0;];
Cq = t18;
