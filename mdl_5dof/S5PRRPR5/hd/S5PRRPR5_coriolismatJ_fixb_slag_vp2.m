% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRPR5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR5_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR5_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR5_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR5_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:25:28
% EndTime: 2019-12-05 16:25:37
% DurationCPUTime: 2.74s
% Computational Cost: add. (5312->279), mult. (12547->415), div. (0->0), fcn. (13683->10), ass. (0->173)
t187 = cos(qJ(3));
t170 = -pkin(3) * t187 - pkin(2);
t305 = m(5) * t170;
t184 = sin(qJ(3));
t250 = cos(pkin(10));
t220 = t250 * t184;
t303 = qJ(4) + pkin(7);
t231 = t303 * t187;
t249 = sin(pkin(10));
t124 = t220 * t303 + t249 * t231;
t182 = cos(pkin(5));
t181 = sin(pkin(5));
t185 = sin(qJ(2));
t237 = t181 * t185;
t145 = t182 * t187 - t184 * t237;
t236 = t181 * t187;
t146 = t182 * t184 + t185 * t236;
t198 = t249 * t145 + t250 * t146;
t304 = t124 * t198;
t188 = cos(qJ(2));
t302 = m(4) * t188;
t186 = cos(qJ(5));
t172 = Ifges(6,5) * t186;
t183 = sin(qJ(5));
t269 = Ifges(6,6) * t183;
t301 = Ifges(5,4) - t172 / 0.2e1 + t269 / 0.2e1;
t224 = t249 * pkin(3);
t168 = t224 + pkin(8);
t177 = t183 ^ 2;
t179 = t186 ^ 2;
t230 = t177 + t179;
t300 = t168 * t230;
t151 = -t249 * t187 - t220;
t260 = t186 * mrSges(6,2);
t265 = t183 * mrSges(6,1);
t157 = t260 + t265;
t108 = t157 * t151;
t267 = t151 * mrSges(5,3);
t299 = t108 + t267;
t173 = Ifges(6,4) * t186;
t273 = Ifges(6,1) * t183;
t160 = t173 + t273;
t298 = t184 ^ 2 + t187 ^ 2;
t219 = t249 * t184;
t149 = -t250 * t187 + t219;
t278 = t184 * pkin(3);
t111 = -pkin(4) * t151 + pkin(8) * t149 + t278;
t59 = t111 * t186 + t124 * t183;
t60 = t111 * t183 - t124 * t186;
t297 = -t59 * t183 + t60 * t186;
t296 = t184 * mrSges(4,1) + t187 * mrSges(4,2);
t295 = -t187 * mrSges(4,1) + t184 * mrSges(4,2);
t216 = Ifges(6,2) * t183 - t173;
t274 = mrSges(6,3) * t151;
t228 = t183 * t274;
t113 = -t149 * mrSges(6,2) + t228;
t233 = t186 * t113;
t115 = t149 * mrSges(6,1) + t186 * t274;
t234 = t183 * t115;
t203 = -t233 / 0.2e1 + t234 / 0.2e1;
t294 = mrSges(6,3) * t230;
t288 = m(5) * pkin(3);
t293 = t249 * t288 - mrSges(5,2);
t292 = -t219 * t303 + t250 * t231;
t156 = -mrSges(6,1) * t186 + mrSges(6,2) * t183;
t225 = t250 * pkin(3);
t169 = -t225 - pkin(4);
t291 = m(6) * t169 - t250 * t288 - mrSges(5,1) + t156;
t290 = m(5) / 0.2e1;
t289 = m(6) / 0.2e1;
t235 = t181 * t188;
t73 = -t183 * t198 - t186 * t235;
t287 = t73 / 0.2e1;
t106 = t151 * t156;
t286 = t106 / 0.2e1;
t285 = -t149 / 0.2e1;
t284 = -t151 / 0.2e1;
t283 = -t183 / 0.2e1;
t282 = t183 / 0.2e1;
t281 = -t186 / 0.2e1;
t280 = t186 / 0.2e1;
t272 = Ifges(6,4) * t183;
t270 = Ifges(6,2) * t186;
t268 = t149 * mrSges(5,3);
t95 = -t250 * t145 + t249 * t146;
t266 = t151 * t95;
t264 = t183 * t73;
t78 = Ifges(6,6) * t149 + t216 * t151;
t263 = t183 * t78;
t110 = pkin(4) * t149 + pkin(8) * t151 + t170;
t56 = t110 * t183 + t186 * t292;
t259 = t186 * t56;
t74 = -t183 * t235 + t186 * t198;
t258 = t186 * t74;
t161 = Ifges(6,1) * t186 - t272;
t80 = Ifges(6,5) * t149 - t151 * t161;
t257 = t186 * t80;
t131 = t151 * t235;
t252 = t95 * t131;
t251 = t95 * t157;
t210 = -t258 + t264;
t10 = m(6) * (t198 + t210) * t95;
t248 = t10 * qJD(1);
t132 = t149 * t235;
t104 = t132 * t183 + t186 * t237;
t247 = t104 * t183;
t105 = -t132 * t186 + t183 * t237;
t246 = t105 * t186;
t244 = t124 * t131;
t243 = t124 * t151;
t238 = t181 ^ 2 * t185;
t13 = m(6) * (t104 * t73 + t105 * t74 - t252) + (-t145 * t181 * t184 + t146 * t236 - t238) * t302 + (-t132 * t198 - t238 * t188 - t252) * m(5);
t242 = t13 * qJD(1);
t229 = t288 / 0.2e1;
t191 = (-t149 * t300 - t169 * t151) * t289 + t156 * t284 + (-t249 * t149 + t250 * t151) * t229 + t285 * t294;
t240 = t149 * t183;
t112 = mrSges(6,2) * t151 + mrSges(6,3) * t240;
t239 = t149 * t186;
t114 = -mrSges(6,1) * t151 + mrSges(6,3) * t239;
t196 = (t183 * t60 + t186 * t59) * t289 + t112 * t282 + t114 * t280 + t184 * t229;
t147 = t149 * mrSges(5,2);
t218 = -t151 * mrSges(5,1) - t147;
t14 = t191 - t196 - t218;
t241 = t14 * qJD(2);
t204 = t260 / 0.2e1 + t265 / 0.2e1;
t200 = t204 * t149;
t22 = t200 + t203;
t232 = t22 * qJD(2);
t223 = t240 / 0.2e1;
t222 = -t239 / 0.2e1;
t221 = -t235 / 0.2e1;
t217 = t172 - t269;
t158 = t270 + t272;
t215 = Ifges(6,5) * t183 + Ifges(6,6) * t186;
t107 = t157 * t149;
t55 = t110 * t186 - t183 * t292;
t211 = t183 * t55 - t259;
t212 = t292 * t95 + t304;
t189 = (-t235 * t278 + t212 - t304) * t290 + (t59 * t73 + t60 * t74 + t212) * t289 + t114 * t287 + t74 * t112 / 0.2e1 - t198 * t108 / 0.2e1 + (t218 + t296) * t221 + (-t107 / 0.2e1 + t211 * t289 - t290 * t292 + t203) * t95;
t190 = (t246 - t247) * t168 * t289 + t296 * t221 - (-mrSges(5,2) / 0.2e1 + t249 * t229) * t132 - (t169 * t289 - mrSges(5,1) / 0.2e1 + t156 / 0.2e1 - t250 * t229) * t131 + (-t247 / 0.2e1 + t246 / 0.2e1) * mrSges(6,3);
t1 = -t189 + t190;
t116 = mrSges(5,1) * t149 - mrSges(5,2) * t151;
t77 = -Ifges(6,6) * t151 + t216 * t149;
t79 = -Ifges(6,5) * t151 - t161 * t149;
t3 = -t124 * t107 - t292 * t108 + t56 * t112 + t60 * t113 + t55 * t114 + t59 * t115 - t170 * t147 + (-pkin(2) * mrSges(4,1) - Ifges(4,4) * t184 + pkin(3) * t116) * t184 + m(6) * (t292 * t124 + t55 * t59 + t56 * t60) + t278 * t305 + (t263 / 0.2e1 - t257 / 0.2e1 + t301 * t149) * t149 + (-pkin(2) * mrSges(4,2) + Ifges(4,4) * t187 + (Ifges(4,1) - Ifges(4,2)) * t184) * t187 + (-t170 * mrSges(5,1) + t77 * t282 + t79 * t281 - t301 * t151 + (-Ifges(6,3) + Ifges(5,1) - Ifges(5,2)) * t149) * t151;
t214 = -t1 * qJD(1) + t3 * qJD(2);
t6 = -t124 * t106 + t56 * t115 + (t80 * t283 + t78 * t281 - mrSges(6,3) * t259 + t215 * t285 + (t158 * t283 + t160 * t280) * t151) * t151 + (-t113 + t228) * t55;
t195 = (t258 / 0.2e1 - t264 / 0.2e1) * t274 + t113 * t287 - t74 * t115 / 0.2e1 + t95 * t286;
t205 = t104 * mrSges(6,1) / 0.2e1 - t105 * mrSges(6,2) / 0.2e1;
t7 = t195 - t205;
t213 = t7 * qJD(1) - t6 * qJD(2);
t11 = t299 * t151 + (-t233 + t234 + t268) * t149 + m(6) * (t211 * t149 - t243) + m(5) * (-t149 * t292 - t243);
t194 = (t210 * t149 - t266) * t289 + (-t149 * t198 - t266) * t290;
t197 = (t186 * t104 + t183 * t105) * t289 + t237 * t290;
t16 = t194 - t197;
t209 = -qJD(1) * t16 - qJD(2) * t11;
t207 = -t59 * mrSges(6,1) / 0.2e1 + t60 * mrSges(6,2) / 0.2e1;
t202 = t151 * t215;
t201 = t204 * t95;
t18 = -t251 / 0.2e1 + t201;
t192 = (t179 / 0.2e1 + t177 / 0.2e1) * t168 * mrSges(6,3) + (-t161 / 0.4e1 + t158 / 0.4e1 + t270 / 0.4e1) * t186 + (t160 / 0.4e1 - t216 / 0.4e1 + t173 / 0.2e1 + t273 / 0.4e1) * t183;
t193 = (t113 * t283 + t115 * t281) * t168 + t124 * t157 / 0.2e1 + t169 * t286 - t263 / 0.4e1 + t257 / 0.4e1;
t5 = (-0.3e1 / 0.4e1 * t269 + 0.3e1 / 0.4e1 * t172) * t149 + (Ifges(6,3) / 0.2e1 + t192) * t151 + t193 + t207;
t54 = t169 * t157 + (t160 / 0.2e1 - t216 / 0.2e1) * t186 + (t161 / 0.2e1 - t158 / 0.2e1) * t183;
t199 = t18 * qJD(1) - t5 * qJD(2) - t54 * qJD(3);
t23 = t200 - t203;
t19 = t251 / 0.2e1 + t201;
t17 = t191 + t196;
t15 = t194 + t197;
t8 = t195 + t205;
t4 = t149 * t217 / 0.4e1 + Ifges(6,5) * t222 + Ifges(6,6) * t223 + Ifges(6,3) * t284 + t192 * t151 + t193 - t207;
t2 = t189 + t190;
t9 = [t13 * qJD(2) + t10 * qJD(3), t2 * qJD(3) + t15 * qJD(4) + t8 * qJD(5) + t242 + (t104 * t115 + t105 * t113 + t299 * t131 + t132 * t268 + 0.2e1 * (t104 * t55 + t105 * t56 - t244) * t289 + 0.2e1 * (-t132 * t292 - t244) * t290 + ((t298 * mrSges(4,3) - mrSges(3,2)) * t188 + t298 * pkin(7) * t302 + (-m(4) * pkin(2) - mrSges(3,1) + t116 + t295 + t305) * t185) * t181) * qJD(2), t248 + t2 * qJD(2) + (-t146 * mrSges(4,1) - t145 * mrSges(4,2) + t291 * t198 + (-m(6) * t300 - t293 - t294) * t95) * qJD(3) + t19 * qJD(5), qJD(2) * t15, t8 * qJD(2) + t19 * qJD(3) + (-mrSges(6,1) * t74 - mrSges(6,2) * t73) * qJD(5); -qJD(3) * t1 + qJD(4) * t16 + qJD(5) * t7 - t242, qJD(3) * t3 + qJD(4) * t11 - qJD(5) * t6, (t224 * t267 + t225 * t268 + Ifges(4,5) * t187 - Ifges(4,6) * t184 - Ifges(5,5) * t149 + Ifges(5,6) * t151 + t160 * t222 + t158 * t223 - t202 / 0.2e1 + t79 * t282 + t77 * t280 - t169 * t107 - t293 * t124 + (m(6) * t297 + t186 * t112 - t183 * t114) * t168 + t291 * t292 + t295 * pkin(7) + t297 * mrSges(6,3)) * qJD(3) + t17 * qJD(4) + t4 * qJD(5) + t214, qJD(3) * t17 + qJD(5) * t23 - t209, t4 * qJD(3) + t23 * qJD(4) + (-mrSges(6,1) * t56 - mrSges(6,2) * t55 + t202) * qJD(5) + t213; qJD(2) * t1 - qJD(5) * t18 - t248, qJD(4) * t14 + qJD(5) * t5 - t214, t54 * qJD(5), t241, (t156 * t168 + t217) * qJD(5) - t199; -qJD(2) * t16, -qJD(3) * t14 - qJD(5) * t22 + t209, -t241, 0, -t157 * qJD(5) - t232; -t7 * qJD(2) + t18 * qJD(3), -qJD(3) * t5 + qJD(4) * t22 - t213, t199, t232, 0;];
Cq = t9;
