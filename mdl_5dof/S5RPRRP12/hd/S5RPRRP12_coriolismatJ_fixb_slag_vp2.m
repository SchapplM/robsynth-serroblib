% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 18:57
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRP12_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP12_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP12_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP12_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP12_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:56:18
% EndTime: 2019-12-31 18:56:25
% DurationCPUTime: 3.02s
% Computational Cost: add. (3484->382), mult. (7372->537), div. (0->0), fcn. (5672->4), ass. (0->189)
t197 = cos(qJ(4));
t321 = Ifges(5,6) + Ifges(6,6);
t325 = t197 * t321;
t195 = sin(qJ(4));
t193 = t195 ^ 2;
t194 = t197 ^ 2;
t256 = t193 + t194;
t324 = mrSges(5,3) + mrSges(6,3);
t198 = cos(qJ(3));
t289 = pkin(7) * t198;
t196 = sin(qJ(3));
t292 = pkin(3) * t196;
t152 = qJ(2) - t289 + t292;
t135 = t197 * t152;
t199 = -pkin(1) - pkin(6);
t232 = -t195 * t199 + pkin(4);
t258 = t197 * t198;
t244 = qJ(5) * t258;
t33 = t196 * t232 + t135 - t244;
t262 = t195 * t198;
t259 = t196 * t199;
t53 = t195 * t152 + t197 * t259;
t38 = -qJ(5) * t262 + t53;
t323 = m(6) * (t195 * t38 + t33 * t197);
t322 = Ifges(5,5) + Ifges(6,5);
t320 = Ifges(5,3) + Ifges(6,3);
t310 = m(6) * pkin(4);
t246 = mrSges(6,1) + t310;
t273 = t197 * mrSges(5,2);
t280 = t195 * mrSges(5,1);
t157 = t273 + t280;
t319 = t196 * t157;
t275 = t196 * mrSges(6,2);
t142 = -mrSges(6,3) * t262 - t275;
t143 = -t196 * mrSges(5,2) - mrSges(5,3) * t262;
t318 = t142 + t143;
t276 = t196 * mrSges(6,1);
t146 = -mrSges(6,3) * t258 + t276;
t147 = t196 * mrSges(5,1) - mrSges(5,3) * t258;
t317 = t146 + t147;
t191 = Ifges(6,4) * t197;
t160 = Ifges(6,1) * t195 + t191;
t192 = Ifges(5,4) * t197;
t161 = Ifges(5,1) * t195 + t192;
t164 = pkin(3) * t198 + pkin(7) * t196;
t139 = t197 * t164;
t257 = t198 * t199;
t58 = -t195 * t257 + t139;
t59 = t195 * t164 + t197 * t257;
t218 = -t195 * t58 + t197 * t59;
t261 = t196 * t197;
t34 = qJ(5) * t261 + t198 * t232 + t139;
t263 = t195 * t196;
t39 = qJ(5) * t263 + t59;
t315 = -t195 * t34 + t197 * t39;
t248 = Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1;
t221 = -Ifges(6,2) * t195 + t191;
t84 = Ifges(6,6) * t196 + t198 * t221;
t222 = -Ifges(5,2) * t195 + t192;
t86 = Ifges(5,6) * t196 + t198 * t222;
t314 = t248 * t196 + t84 / 0.2e1 + t86 / 0.2e1;
t313 = m(5) / 0.2e1;
t312 = -m(6) / 0.2e1;
t311 = m(6) / 0.2e1;
t309 = -t34 / 0.2e1;
t112 = mrSges(6,1) * t258 - mrSges(6,2) * t262;
t308 = t112 / 0.2e1;
t272 = t197 * mrSges(6,2);
t279 = t195 * mrSges(6,1);
t225 = t272 + t279;
t116 = t198 * t225;
t307 = -t116 / 0.2e1;
t291 = pkin(4) * t195;
t231 = -t199 + t291;
t137 = t231 * t198;
t306 = t137 / 0.2e1;
t305 = -t143 / 0.2e1;
t144 = mrSges(6,1) * t198 + mrSges(6,3) * t261;
t304 = t144 / 0.2e1;
t303 = -t146 / 0.2e1;
t302 = t146 / 0.2e1;
t286 = -qJ(5) - pkin(7);
t153 = t286 * t195;
t301 = t153 / 0.2e1;
t300 = t195 / 0.2e1;
t299 = t196 / 0.2e1;
t298 = t197 / 0.2e1;
t297 = -t198 / 0.2e1;
t296 = t198 / 0.2e1;
t295 = m(6) * t137;
t288 = t197 * pkin(4);
t181 = -pkin(3) - t288;
t294 = m(6) * t181;
t293 = m(6) * t196;
t290 = pkin(4) * t198 ^ 2;
t52 = -t195 * t259 + t135;
t37 = t52 - t244;
t285 = -t33 + t37;
t282 = Ifges(5,4) * t195;
t281 = Ifges(6,4) * t195;
t274 = t197 * mrSges(5,1);
t114 = t225 * t196;
t136 = t231 * t196;
t140 = -t198 * mrSges(6,2) + mrSges(6,3) * t263;
t141 = -t198 * mrSges(5,2) + mrSges(5,3) * t263;
t145 = t198 * mrSges(5,1) + mrSges(5,3) * t261;
t249 = Ifges(6,5) / 0.2e1 + Ifges(5,5) / 0.2e1;
t223 = Ifges(6,1) * t197 - t281;
t88 = Ifges(6,5) * t196 + t198 * t223;
t224 = Ifges(5,1) * t197 - t282;
t90 = Ifges(5,5) * t196 + t198 * t224;
t205 = -t88 / 0.2e1 - t90 / 0.2e1 - t249 * t196;
t83 = Ifges(6,6) * t198 - t196 * t221;
t85 = Ifges(5,6) * t198 - t196 * t222;
t87 = Ifges(6,5) * t198 - t196 * t223;
t89 = Ifges(5,5) * t198 - t196 * t224;
t3 = -t137 * t114 - t136 * t116 + t38 * t140 + t53 * t141 + t39 * t142 + t59 * t143 + t33 * t144 + t52 * t145 + t34 * t146 + t58 * t147 + m(6) * (-t136 * t137 + t33 * t34 + t38 * t39) + m(5) * (t52 * t58 + t53 * t59) + (-qJ(2) * mrSges(4,2) + Ifges(4,4) * t196 + t314 * t195 + t205 * t197) * t196 + (qJ(2) * mrSges(4,1) - Ifges(4,4) * t198 + t199 * t319 + (t87 / 0.2e1 + t89 / 0.2e1 + t249 * t198) * t197 + (-t83 / 0.2e1 - t85 / 0.2e1 - t248 * t198) * t195 + (-Ifges(4,1) + Ifges(4,2) + (-m(5) * t199 + t157) * t199 + t320) * t196) * t198;
t269 = t3 * qJD(1);
t226 = -t195 * mrSges(5,2) + t274;
t113 = t226 * t198;
t158 = Ifges(6,2) * t197 + t281;
t117 = t198 * t158;
t159 = Ifges(5,2) * t197 + t282;
t118 = t198 * t159;
t119 = t198 * t160;
t120 = t198 * t161;
t245 = m(6) * t285;
t4 = t137 * t112 + t37 * t142 + t52 * t143 - t53 * t147 + (t245 - t146) * t38 + (-t199 * t113 + (t33 * mrSges(6,3) + t52 * mrSges(5,3) + t117 / 0.2e1 + t118 / 0.2e1 + t205) * t195 + (-t38 * mrSges(6,3) - t53 * mrSges(5,3) - t119 / 0.2e1 - t120 / 0.2e1 + (t116 + t295) * pkin(4) - t314) * t197) * t198;
t268 = t4 * qJD(1);
t250 = -mrSges(5,2) / 0.2e1 - mrSges(6,2) / 0.2e1;
t227 = t250 * t195;
t233 = t194 / 0.2e1 + t193 / 0.2e1;
t251 = mrSges(5,1) / 0.2e1 + mrSges(6,1) / 0.2e1;
t7 = (t113 / 0.2e1 + t308) * t198 + t227 + ((t142 / 0.2e1 + t143 / 0.2e1) * t195 + t324 * t198 * t233) * t196 + ((t302 + t147 / 0.2e1) * t196 + (pkin(4) / 0.2e1 + t290 / 0.2e1 + t33 * t299 - t196 * t37 / 0.2e1) * m(6) + t251) * t197;
t267 = t7 * qJD(1);
t266 = -t226 - mrSges(4,1);
t10 = t196 * mrSges(4,1) + t198 * mrSges(4,2) + mrSges(3,3) + t317 * t197 + t318 * t195 + (m(4) + m(3)) * qJ(2) + t323 + m(5) * (t195 * t53 + t197 * t52);
t265 = qJD(1) * t10;
t17 = (-t195 * t142 - t197 * t146 - t323) * t198;
t264 = qJD(1) * t17;
t260 = t196 * t198;
t255 = qJD(3) * t196;
t254 = t310 / 0.2e1;
t253 = t136 * t311;
t243 = -t263 / 0.2e1;
t241 = -t262 / 0.2e1;
t240 = -t261 / 0.2e1;
t236 = t158 / 0.2e1 + t159 / 0.2e1;
t235 = -t160 / 0.2e1 - t161 / 0.2e1;
t189 = Ifges(6,5) * t197;
t190 = Ifges(5,5) * t197;
t234 = t190 / 0.4e1 + t189 / 0.4e1;
t230 = t256 * mrSges(6,3);
t219 = -t195 * t52 + t197 * t53;
t23 = 0.4e1 * (m(5) / 0.4e1 + m(6) / 0.4e1) * (-0.1e1 + t256) * t260;
t156 = t286 * t197;
t207 = (t153 * t197 - t156 * t195) * t311;
t210 = t142 * t298 + t195 * t303;
t211 = m(5) * t218;
t212 = t137 + t315;
t213 = -t195 * t33 + t197 * t38 + t136;
t5 = t207 + (t307 + (-t140 / 0.2e1 - t141 / 0.2e1) * t197 + (t304 + t145 / 0.2e1) * t195 + t212 * t312 - t211 / 0.2e1) * t196 + (-t319 - t114 / 0.2e1 + t197 * t305 + t147 * t300 + t213 * t312 - m(5) * (t219 - 0.2e1 * t259) / 0.2e1 - t210) * t198;
t217 = -t5 * qJD(1) + t23 * qJD(2);
t216 = -t153 * t195 - t156 * t197;
t121 = t195 * t246 + t272;
t55 = -t258 * t310 - t112;
t215 = qJD(1) * t55 - qJD(3) * t121;
t214 = t254 + t251;
t18 = (-t273 / 0.2e1 + t157 / 0.2e1 - t280 / 0.2e1) * t198;
t154 = -t197 * mrSges(6,1) + t195 * mrSges(6,2);
t200 = -t199 * t157 / 0.2e1 - t233 * pkin(7) * mrSges(5,3) + (-t161 / 0.4e1 - t160 / 0.4e1 - t192 / 0.4e1 - t191 / 0.4e1 + mrSges(6,3) * t301 + (Ifges(5,2) / 0.4e1 + Ifges(6,2) / 0.4e1) * t195) * t195 + (-t159 / 0.4e1 - t158 / 0.4e1 + t156 * mrSges(6,3) / 0.2e1 + (Ifges(5,1) / 0.4e1 + Ifges(6,1) / 0.4e1) * t197 + (-Ifges(5,4) / 0.4e1 - Ifges(6,4) / 0.4e1) * t195 + (t294 / 0.2e1 + t154 / 0.2e1) * pkin(4)) * t197;
t201 = -(t245 / 0.2e1 + t303) * t156 + (-pkin(7) * t147 / 0.2e1 - t118 / 0.4e1 - t117 / 0.4e1 + t90 / 0.4e1 + t88 / 0.4e1 + mrSges(6,2) * t306 + (t37 / 0.2e1 - t33 / 0.2e1) * mrSges(6,3)) * t197 - pkin(3) * t113 / 0.2e1 + t142 * t301 + t181 * t308;
t202 = pkin(7) * t305 - t120 / 0.4e1 - t119 / 0.4e1 - t86 / 0.4e1 - t84 / 0.4e1 + mrSges(6,1) * t306 + (t295 / 0.2e1 + t116 / 0.2e1) * pkin(4);
t204 = mrSges(6,1) * t309 + t39 * mrSges(6,2) / 0.2e1 - t58 * mrSges(5,1) / 0.2e1 + t59 * mrSges(5,2) / 0.2e1;
t2 = (m(6) * t309 - t144 / 0.2e1) * pkin(4) + t202 * t195 + (t249 * t197 + (-0.3e1 / 0.4e1 * Ifges(5,6) - 0.3e1 / 0.4e1 * Ifges(6,6)) * t195 + t234) * t196 + (-Ifges(5,3) / 0.2e1 - Ifges(6,3) / 0.2e1 + t200) * t198 + t201 + t204;
t9 = pkin(3) * t157 - t181 * t225 - t154 * t291 + (-t192 / 0.2e1 - t191 / 0.2e1 + t235) * t197 + (-pkin(4) * t294 + (Ifges(6,4) / 0.2e1 + Ifges(5,4) / 0.2e1) * t195 + (Ifges(5,2) / 0.2e1 - Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t197 + t236) * t195;
t208 = t2 * qJD(1) - t18 * qJD(2) - t9 * qJD(3);
t203 = m(6) * ((-t153 * t198 + t38) * t197 + (t156 * t198 - t33) * t195);
t13 = (-t275 / 0.2e1 - t142 / 0.2e1) * t197 + (-t276 / 0.2e1 + t302) * t195 - t253 - t203 / 0.2e1;
t26 = m(6) * t216 + t230;
t57 = (-0.1e1 / 0.2e1 + t233) * t293;
t206 = -qJD(1) * t13 + qJD(2) * t57 + qJD(3) * t26;
t56 = (t256 + 0.1e1) * t293 / 0.2e1;
t19 = t157 * t297 + t307 + t241 * t310 + (-t195 * t214 + t197 * t250) * t198;
t14 = t203 / 0.2e1 - t253 + (-t272 / 0.2e1 - t279 / 0.2e1) * t196 + t210;
t8 = m(6) * (t196 * t285 - t290) * t298 + t227 + t214 * t197 + (t113 + t112) * t297 + t318 * t243 + t317 * t240 - t324 * t256 * t260 / 0.2e1;
t6 = t319 * t296 + t116 * t299 + (t196 * t212 + t198 * t213) * t311 + (t219 * t198 + (t218 - 0.2e1 * t257) * t196) * t313 + t207 + (-t319 - t114) * t297 + (t144 + t145) * t243 + t317 * t241 + (t140 + t141) * t261 / 0.2e1 + t318 * t258 / 0.2e1;
t1 = t200 * t198 + t201 + t234 * t196 + t34 * t254 + pkin(4) * t304 + ((-Ifges(5,6) / 0.4e1 - Ifges(6,6) / 0.4e1) * t196 + t202) * t195 - t204 + t320 * t296 + t321 * t263 / 0.2e1 + t322 * t240;
t11 = [qJD(2) * t10 + qJD(3) * t3 + qJD(4) * t4 + qJD(5) * t17, qJD(3) * t6 + qJD(4) * t8 + t265, t269 + t6 * qJD(2) + t1 * qJD(4) + t14 * qJD(5) + (-Ifges(4,5) + t235 * t197 + t236 * t195 + (-m(5) * pkin(3) + t266) * t199) * t255 + (m(6) * (-t136 * t181 + t153 * t34 - t156 * t39) - mrSges(4,2) * t257 + pkin(3) * t319 + t153 * t144 - t136 * t154 - t156 * t140 - t181 * t114 - Ifges(4,6) * t198 + (t197 * t141 - t195 * t145 + t211) * pkin(7) + (t87 + t89) * t300 + (t83 + t85) * t298 + (t195 * t322 + t325) * t296 + t315 * mrSges(6,3) + t218 * mrSges(5,3)) * qJD(3), t8 * qJD(2) + t1 * qJD(3) + t268 + (-mrSges(5,1) * t53 - mrSges(5,2) * t52 - mrSges(6,2) * t37 + (-t325 + (mrSges(6,3) * pkin(4) - t322) * t195) * t198 - t246 * t38) * qJD(4), qJD(3) * t14 + t264; -qJD(3) * t5 - qJD(4) * t7 - t265, t23 * qJD(3), t19 * qJD(4) + t56 * qJD(5) + (t154 + t266) * t255 + (mrSges(5,3) * t256 - mrSges(4,2) + t230) * qJD(3) * t198 + 0.2e1 * ((t256 * t289 - t292) * t313 + (t181 * t196 + t198 * t216) * t311) * qJD(3) + t217, -t267 + t19 * qJD(3) + ((mrSges(5,2) + mrSges(6,2)) * t195 + (-mrSges(5,1) - t246) * t197) * qJD(4) * t196, t56 * qJD(3); qJD(2) * t5 + qJD(4) * t2 - qJD(5) * t13 - t269, -qJD(4) * t18 + qJD(5) * t57 - t217, -qJD(4) * t9 + qJD(5) * t26, t208 + (-mrSges(6,2) * t153 - mrSges(6,3) * t288 - pkin(7) * t274 + t189 + t190 + (mrSges(5,2) * pkin(7) - t321) * t195 + t246 * t156) * qJD(4), t206; qJD(2) * t7 - qJD(3) * t2 + qJD(5) * t55 - t268, qJD(3) * t18 + t267, -qJD(5) * t121 - t208, 0, t215; qJD(3) * t13 - qJD(4) * t55 - t264, -t57 * qJD(3), qJD(4) * t121 - t206, -t215, 0;];
Cq = t11;
