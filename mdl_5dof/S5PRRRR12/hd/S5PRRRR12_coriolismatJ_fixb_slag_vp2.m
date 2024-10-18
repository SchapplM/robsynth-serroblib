% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
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
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRRR12_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR12_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR12_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR12_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR12_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-28 18:07:25
% EndTime: 2024-09-28 18:07:28
% DurationCPUTime: 1.91s
% Computational Cost: add. (9139->302), mult. (25527->435), div. (0->0), fcn. (27132->12), ass. (0->200)
t193 = sin(qJ(3));
t196 = cos(qJ(3));
t278 = sin(pkin(5));
t234 = t278 * sin(qJ(2));
t235 = cos(qJ(2)) * t278;
t154 = -t193 * t234 + t196 * t235;
t155 = -t193 * t235 - t196 * t234;
t192 = sin(qJ(4));
t195 = cos(qJ(4));
t109 = t192 * t154 - t195 * t155;
t188 = sin(pkin(6));
t187 = t188 ^ 2;
t329 = t109 * t187;
t191 = sin(qJ(5));
t194 = cos(qJ(5));
t165 = (-mrSges(6,1) * t194 + mrSges(6,2) * t191) * t188;
t265 = t188 * t165;
t239 = t195 * t154 + t192 * t155;
t314 = t239 * mrSges(5,2);
t326 = t109 * mrSges(5,1);
t328 = t109 * t265 - t314 - t326;
t327 = -m(6) / 0.2e1;
t299 = mrSges(6,1) / 0.2e1;
t298 = -mrSges(6,2) / 0.2e1;
t292 = t195 * pkin(3);
t184 = pkin(4) + t292;
t266 = t184 * t187;
t325 = t109 * t266;
t291 = t196 * pkin(2);
t185 = pkin(3) + t291;
t257 = t192 * t193;
t162 = -pkin(2) * t257 + t195 * t185;
t157 = pkin(4) + t162;
t324 = t157 * t329;
t189 = cos(pkin(6));
t190 = cos(pkin(5));
t90 = -t188 * t239 + t190 * t189;
t283 = t188 * t90;
t323 = t109 * t283;
t322 = pkin(4) * t329;
t262 = t188 * t194;
t263 = t188 * t191;
t321 = 0.2e1 * Ifges(6,4) * t262 + Ifges(6,5) * t189 + (Ifges(6,1) - Ifges(6,2)) * t263;
t320 = qJD(2) + qJD(3);
t175 = -t189 * mrSges(6,2) + mrSges(6,3) * t262;
t256 = t193 * t195;
t163 = pkin(2) * t256 + t192 * t185;
t186 = t188 * pkin(10);
t144 = t163 + t186;
t259 = t189 * t194;
t96 = -t191 * t144 + t157 * t259;
t80 = t96 * t175;
t174 = t189 * mrSges(6,1) - mrSges(6,3) * t263;
t260 = t189 * t191;
t97 = t194 * t144 + t157 * t260;
t81 = t97 * t174;
t319 = t80 / 0.2e1 - t81 / 0.2e1;
t293 = t192 * pkin(3);
t176 = t186 + t293;
t133 = t194 * t176 + t184 * t260;
t276 = t133 * t174;
t132 = -t191 * t176 + t184 * t259;
t277 = t132 * t175;
t318 = -t276 / 0.2e1 + t277 / 0.2e1;
t172 = (-t192 * t196 - t256) * pkin(2);
t164 = t172 * mrSges(5,1);
t264 = t188 * t172;
t173 = (t195 * t196 - t257) * pkin(2);
t284 = t173 * mrSges(5,2);
t116 = t172 * t259 - t191 * t173;
t93 = t116 * t174;
t117 = t172 * t260 + t194 * t173;
t94 = t117 * t175;
t317 = -(t193 * mrSges(4,1) + t196 * mrSges(4,2)) * pkin(2) - t165 * t264 + t164 - t284 + t93 + t94;
t170 = pkin(4) * t259 - pkin(10) * t263;
t171 = pkin(4) * t260 + pkin(10) * t262;
t312 = t194 * t239;
t221 = -t109 * t260 + t312;
t313 = t191 * t239;
t222 = -t109 * t259 - t313;
t244 = t265 / 0.2e1;
t60 = t222 * t174;
t61 = t221 * t175;
t316 = -t109 * t244 + t326 / 0.2e1 + (t170 * t222 + t171 * t221 - t322) * t327 - t60 / 0.2e1 - t61 / 0.2e1;
t166 = (t191 * mrSges(6,1) + t194 * mrSges(6,2)) * t188;
t228 = t188 * t190 + t189 * t239;
t255 = t194 * t109;
t63 = t228 * t191 + t255;
t296 = -t63 / 0.2e1;
t258 = t191 * t109;
t62 = t228 * t194 - t258;
t210 = (t191 * t296 - t62 * t194 / 0.2e1) * t188 * mrSges(6,3) + t62 * t175 / 0.2e1 + t174 * t296 + t90 * t166 / 0.2e1;
t309 = t221 * t298 + t222 * t299;
t198 = t210 - t309;
t310 = qJD(1) * t198;
t21 = m(6) * (t63 * t221 + t62 * t222 + t323);
t254 = t21 * qJD(1);
t308 = qJD(5) * t198 - t254;
t208 = t210 + t309;
t307 = t208 * qJD(5) + t254;
t241 = -t262 / 0.2e1;
t236 = mrSges(6,3) * t241;
t243 = -t263 / 0.2e1;
t237 = mrSges(6,3) * t243;
t246 = -t166 * t188 / 0.2e1;
t306 = t157 * t246 + t96 * t236 + t97 * t237 + t319;
t305 = t132 * t236 + t133 * t237 + t184 * t246 + t318;
t303 = 2 * qJD(3);
t302 = m(5) / 0.2e1;
t301 = m(6) / 0.2e1;
t300 = -mrSges(6,1) / 0.2e1;
t297 = mrSges(6,2) / 0.2e1;
t295 = pkin(4) * t166;
t290 = -t109 * t162 + t163 * t239;
t289 = Ifges(6,4) * t191;
t285 = t162 * mrSges(5,2);
t70 = -t189 * t255 - t313;
t280 = t70 * t174;
t71 = -t189 * t258 + t312;
t279 = t71 * t175;
t160 = (-t191 * t195 - t192 * t259) * pkin(3);
t274 = t160 * t174;
t161 = (-t192 * t260 + t194 * t195) * pkin(3);
t273 = t161 * t175;
t272 = t163 * t187;
t271 = t170 * t175;
t270 = t170 * t194;
t269 = t171 * t174;
t268 = t171 * t191;
t267 = t172 * t187;
t167 = Ifges(6,5) * t262 - Ifges(6,6) * t263;
t261 = t189 * t167;
t24 = m(6) * (t62 * t70 + t63 * t71 + t323);
t253 = t24 * qJD(1);
t252 = mrSges(5,1) * t293;
t251 = mrSges(5,2) * t292;
t250 = pkin(3) * t302;
t249 = t97 * t221 + t96 * t222 - t324;
t248 = t187 * t293;
t247 = t188 * t293;
t245 = -t265 / 0.2e1;
t242 = t263 / 0.2e1;
t238 = t165 * t247;
t150 = Ifges(6,6) * t189 + (Ifges(6,2) * t194 + t289) * t188;
t169 = (Ifges(6,1) * t194 - t289) * t188;
t233 = t150 * t242 + t169 * t243 - t261 / 0.2e1 + t321 * t241;
t232 = t150 * t243 + t169 * t242 + t261 / 0.2e1 + t321 * t262 / 0.2e1;
t223 = t280 / 0.2e1 + t279 / 0.2e1;
t197 = (t244 - mrSges(5,1) / 0.2e1) * t109 + t223 + t316;
t110 = -t191 * t162 - t163 * t259;
t111 = t194 * t162 - t163 * t260;
t209 = (t110 * t62 + t111 * t63 + t163 * t283 + t96 * t70 + t97 * t71 - t324) * t301;
t2 = t209 + t197;
t156 = t163 * mrSges(5,1);
t88 = t110 * t174;
t89 = t111 * t175;
t216 = t163 * t265 - t156 - t285 + t88 + t89;
t32 = -m(6) * (t96 * t110 + t97 * t111 - t157 * t272) - t216;
t231 = t2 * qJD(1) - t32 * qJD(2);
t29 = -t80 + t81 + (t157 * t166 + (t97 * t191 + t96 * t194) * mrSges(6,3)) * t188 + t233;
t230 = -t29 * qJD(2) + t310;
t33 = -m(6) * (t96 * t116 + t97 * t117 + t157 * t267) - m(5) * (t162 * t172 + t163 * t173) - t317;
t200 = (t132 * t222 + t133 * t221 - t325) * t301 + (-t109 * t195 + t192 * t239) * t250;
t202 = (t116 * t62 + t117 * t63 - t90 * t264 + t249) * t301 + (t173 * t109 + t172 * t239 + t290) * t302;
t6 = -t200 + t202;
t229 = t6 * qJD(1) - t33 * qJD(2);
t227 = t71 * t298 + t70 * t299;
t226 = t110 * t300 + t111 * t297;
t225 = t116 * t299 + t117 * t298;
t224 = t160 * t300 + t161 * t297;
t199 = t156 / 0.2e1 - t88 / 0.2e1 - t89 / 0.2e1 + (t132 * t110 + t133 * t111 + t160 * t96 + t161 * t97 + (-t157 * t293 - t163 * t184) * t187) * t327 - t274 / 0.2e1 - t273 / 0.2e1 + t252 / 0.2e1;
t212 = t164 / 0.2e1 + t93 / 0.2e1 + t94 / 0.2e1 + (pkin(4) * t267 + t170 * t116 + t171 * t117) * t301;
t13 = (-t172 / 0.2e1 - t293 / 0.2e1 - t163 / 0.2e1) * t265 + (-t173 / 0.2e1 + t292 / 0.2e1 + t162 / 0.2e1) * mrSges(5,2) + t199 + t212;
t211 = t238 - t251 - t252 + t273 + t274;
t47 = m(6) * (t132 * t160 + t133 * t161 - t184 * t248) + t211;
t207 = (t132 * t70 + t133 * t71 + t160 * t62 + t161 * t63 + t90 * t247 - t325) * t301;
t5 = t207 + t197;
t218 = t5 * qJD(1) - t13 * qJD(2) + t47 * qJD(3);
t15 = ((t184 / 0.2e1 + t157 / 0.2e1) * t166 + ((t132 / 0.2e1 + t96 / 0.2e1) * t194 + (t133 / 0.2e1 + t97 / 0.2e1) * t191) * mrSges(6,3)) * t188 + t225 + t233 - t318 - t319;
t35 = t277 - t276 + (-t184 * t166 + (-t132 * t194 - t133 * t191) * mrSges(6,3)) * t188 + t232;
t217 = -t15 * qJD(2) + t35 * qJD(3) + t310;
t11 = t210 - t227;
t204 = (-t295 / 0.2e1 + (-t268 / 0.2e1 - t270 / 0.2e1) * mrSges(6,3)) * t188 + t271 / 0.2e1 - t269 / 0.2e1 + t232;
t203 = t204 + t306;
t18 = t203 + t226;
t201 = t204 + t305;
t27 = t201 + t224;
t38 = t271 - t269 + (-t295 + (-t268 - t270) * mrSges(6,3)) * t188 + t232;
t214 = t11 * qJD(1) + t18 * qJD(2) + t27 * qJD(3) + t38 * qJD(4);
t213 = t155 * mrSges(4,1) - t154 * mrSges(4,2) + t328 + t60 + t61;
t205 = -t314 - t109 * t245 - t326 / 0.2e1 + t223 - t316;
t26 = t201 - t224;
t17 = t203 - t226;
t16 = t225 + t232 + t305 + t306;
t14 = t238 / 0.2e1 + t172 * t245 + t163 * t244 - t251 / 0.2e1 - t285 / 0.2e1 - t284 / 0.2e1 - t199 + t212;
t12 = t210 + t227;
t4 = t207 + t205;
t3 = t200 + t202 + t213;
t1 = t209 + t205;
t7 = [t24 * qJD(4) + t320 * t21, t3 * qJD(3) + t1 * qJD(4) + (-mrSges(3,1) * t234 - mrSges(3,2) * t235 + t213 + 0.2e1 * t249 * t301 + 0.2e1 * t290 * t302 + m(4) * (t193 * pkin(2) * t154 + t155 * t291)) * qJD(2) + t307, t3 * qJD(2) + t213 * qJD(3) + t4 * qJD(4) + t200 * t303 + t307, t253 + t1 * qJD(2) + t4 * qJD(3) + (t279 + t280 + m(6) * (t170 * t70 + t171 * t71 - t322) + t328) * qJD(4) + t12 * qJD(5), t12 * qJD(4) + (-t63 * mrSges(6,1) - t62 * mrSges(6,2)) * qJD(5) + t320 * t208; t6 * qJD(3) + t2 * qJD(4) + t308, -t33 * qJD(3) - t32 * qJD(4) - t29 * qJD(5), t317 * qJD(3) + t14 * qJD(4) + t16 * qJD(5) + ((t132 * t116 + t133 * t117 + t172 * t266) * t301 + (t172 * t195 + t173 * t192) * t250) * t303 + t229, t14 * qJD(3) + (m(6) * (-pkin(4) * t272 + t170 * t110 + t171 * t111) + t216) * qJD(4) + t17 * qJD(5) + t231, t16 * qJD(3) + t17 * qJD(4) + (-t97 * mrSges(6,1) - t96 * mrSges(6,2) + t167) * qJD(5) + t230; -t6 * qJD(2) + t5 * qJD(4) + t308, -t13 * qJD(4) - t15 * qJD(5) - t229, t47 * qJD(4) + t35 * qJD(5), (m(6) * (-pkin(4) * t248 + t170 * t160 + t171 * t161) + t211) * qJD(4) + t26 * qJD(5) + t218, t26 * qJD(4) + (-t133 * mrSges(6,1) - t132 * mrSges(6,2) + t167) * qJD(5) + t217; -t2 * qJD(2) - t5 * qJD(3) + t11 * qJD(5) - t253, t13 * qJD(3) + t18 * qJD(5) - t231, t27 * qJD(5) - t218, t38 * qJD(5), (-t171 * mrSges(6,1) - t170 * mrSges(6,2) + t167) * qJD(5) + t214; -t11 * qJD(4) - t320 * t198, t15 * qJD(3) - t18 * qJD(4) - t230, -t27 * qJD(4) - t217, -t214, 0;];
Cq = t7;
