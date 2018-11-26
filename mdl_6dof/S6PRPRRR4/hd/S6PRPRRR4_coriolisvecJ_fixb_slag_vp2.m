% Calculate vector of centrifugal and coriolis load on the joints for
% S6PRPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:05
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6PRPRRR4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR4_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR4_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:05:22
% EndTime: 2018-11-23 15:05:31
% DurationCPUTime: 8.64s
% Computational Cost: add. (8966->511), mult. (22992->725), div. (0->0), fcn. (18037->12), ass. (0->242)
t193 = cos(pkin(12));
t201 = cos(qJ(4));
t191 = sin(pkin(12));
t197 = sin(qJ(4));
t253 = t191 * t197;
t169 = -t201 * t193 + t253;
t192 = sin(pkin(6));
t202 = cos(qJ(2));
t251 = t192 * t202;
t207 = t169 * t251;
t278 = pkin(8) + qJ(3);
t178 = t278 * t191;
t179 = t278 * t193;
t310 = -t201 * t178 - t179 * t197;
t328 = qJD(1) * t207 - t169 * qJD(3) + qJD(4) * t310;
t164 = t169 * qJD(4);
t170 = t191 * t201 + t193 * t197;
t165 = t170 * qJD(4);
t198 = sin(qJ(2));
t246 = qJD(1) * t192;
t237 = t198 * t246;
t341 = pkin(4) * t165 + pkin(9) * t164 - t237;
t186 = -pkin(3) * t193 - pkin(2);
t129 = pkin(4) * t169 - pkin(9) * t170 + t186;
t138 = -t178 * t197 + t179 * t201;
t196 = sin(qJ(5));
t200 = cos(qJ(5));
t241 = qJD(5) * t200;
t242 = qJD(5) * t196;
t330 = t129 * t241 - t138 * t242 + t341 * t196 + t328 * t200;
t340 = -t328 * t196 + t341 * t200;
t131 = t200 * t138;
t281 = pkin(10) * t200;
t339 = t164 * t281 + pkin(5) * t165 + (-t131 + (pkin(10) * t170 - t129) * t196) * qJD(5) + t340;
t211 = -t196 * t164 + t170 * t241;
t338 = pkin(10) * t211 - t330;
t296 = -pkin(10) - pkin(9);
t239 = qJD(5) * t296;
t243 = qJD(2) * t193;
t162 = -qJD(2) * t253 + t201 * t243;
t256 = t162 * t196;
t163 = t170 * qJD(2);
t127 = pkin(4) * t163 - pkin(9) * t162;
t177 = qJD(2) * qJ(3) + t237;
t194 = cos(pkin(6));
t245 = qJD(1) * t194;
t183 = t193 * t245;
t134 = t183 + (-pkin(8) * qJD(2) - t177) * t191;
t147 = t193 * t177 + t191 * t245;
t135 = pkin(8) * t243 + t147;
t80 = t134 * t201 - t197 * t135;
t55 = t196 * t127 + t200 * t80;
t337 = pkin(10) * t256 + t196 * t239 - t55;
t54 = t200 * t127 - t196 * t80;
t336 = -pkin(5) * t163 + t162 * t281 + t200 * t239 - t54;
t144 = qJD(4) * t200 - t163 * t196;
t159 = qJD(5) - t162;
t81 = t134 * t197 + t135 * t201;
t77 = qJD(4) * pkin(9) + t81;
t244 = qJD(1) * t202;
t236 = t192 * t244;
t221 = qJD(3) - t236;
t157 = qJD(2) * t186 + t221;
t91 = -pkin(4) * t162 - pkin(9) * t163 + t157;
t43 = -t196 * t77 + t200 * t91;
t44 = t196 * t91 + t200 * t77;
t219 = t44 * t196 + t43 * t200;
t223 = Ifges(6,5) * t200 - Ifges(6,6) * t196;
t274 = Ifges(6,4) * t200;
t225 = -Ifges(6,2) * t196 + t274;
t275 = Ifges(6,4) * t196;
t227 = Ifges(6,1) * t200 - t275;
t228 = mrSges(6,1) * t196 + mrSges(6,2) * t200;
t283 = t200 / 0.2e1;
t284 = -t196 / 0.2e1;
t145 = qJD(4) * t196 + t163 * t200;
t289 = t145 / 0.2e1;
t276 = Ifges(6,4) * t145;
t69 = t144 * Ifges(6,2) + t159 * Ifges(6,6) + t276;
t141 = Ifges(6,4) * t144;
t70 = t145 * Ifges(6,1) + t159 * Ifges(6,5) + t141;
t76 = -qJD(4) * pkin(4) - t80;
t335 = t144 * t225 / 0.2e1 + t227 * t289 + t159 * t223 / 0.2e1 + t76 * t228 + t69 * t284 + t70 * t283 - t219 * mrSges(6,3);
t199 = cos(qJ(6));
t195 = sin(qJ(6));
t36 = pkin(10) * t144 + t44;
t263 = t195 * t36;
t35 = -pkin(10) * t145 + t43;
t27 = pkin(5) * t159 + t35;
t12 = t199 * t27 - t263;
t262 = t199 * t36;
t13 = t195 * t27 + t262;
t156 = qJD(2) * t165;
t151 = Ifges(7,3) * t156;
t232 = t199 * t144 - t145 * t195;
t90 = t144 * t195 + t145 * t199;
t282 = Ifges(7,4) * t90;
t154 = qJD(6) + t159;
t288 = -t154 / 0.2e1;
t298 = -t90 / 0.2e1;
t300 = -t232 / 0.2e1;
t85 = Ifges(7,4) * t232;
t42 = Ifges(7,1) * t90 + Ifges(7,5) * t154 + t85;
t60 = -pkin(5) * t144 + t76;
t334 = t151 + (Ifges(7,5) * t232 - Ifges(7,6) * t90) * t288 + (t12 * t232 + t13 * t90) * mrSges(7,3) + (-Ifges(7,2) * t90 + t42 + t85) * t300 - t60 * (mrSges(7,1) * t90 + mrSges(7,2) * t232) + (Ifges(7,1) * t232 - t282) * t298;
t333 = -t162 / 0.2e1;
t73 = t200 * t129 - t138 * t196;
t59 = pkin(5) * t169 - t170 * t281 + t73;
t255 = t170 * t196;
t74 = t196 * t129 + t131;
t63 = -pkin(10) * t255 + t74;
t23 = -t195 * t63 + t199 * t59;
t332 = qJD(6) * t23 + t195 * t339 - t338 * t199;
t24 = t195 * t59 + t199 * t63;
t331 = -qJD(6) * t24 + t338 * t195 + t199 * t339;
t329 = -qJD(5) * t74 + t340;
t106 = qJD(3) * t170 + qJD(4) * t138;
t208 = t170 * t251;
t139 = qJD(1) * t208;
t327 = t106 - t139;
t158 = Ifges(5,4) * t162;
t322 = t163 * Ifges(5,1) / 0.2e1;
t326 = t157 * mrSges(5,2) + t158 / 0.2e1 + Ifges(5,5) * qJD(4) + t322 + t335;
t155 = qJD(2) * t164;
t101 = -qJD(5) * t145 + t155 * t196;
t168 = (qJD(3) + t236) * qJD(2);
t56 = qJD(4) * t80 - t169 * t168;
t252 = t192 * t198;
t235 = qJD(2) * t252;
t231 = qJD(1) * t235;
t99 = pkin(4) * t156 + pkin(9) * t155 + t231;
t17 = t196 * t99 + t200 * t56 + t91 * t241 - t242 * t77;
t11 = pkin(10) * t101 + t17;
t100 = qJD(5) * t144 - t155 * t200;
t18 = -qJD(5) * t44 - t196 * t56 + t200 * t99;
t8 = pkin(5) * t156 - pkin(10) * t100 + t18;
t2 = qJD(6) * t12 + t11 * t199 + t195 * t8;
t3 = -qJD(6) * t13 - t11 * t195 + t199 * t8;
t30 = qJD(6) * t232 + t100 * t199 + t101 * t195;
t31 = -qJD(6) * t90 - t100 * t195 + t101 * t199;
t325 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t30 + Ifges(7,6) * t31;
t302 = t30 / 0.2e1;
t301 = t31 / 0.2e1;
t41 = Ifges(7,2) * t232 + Ifges(7,6) * t154 + t282;
t323 = t41 / 0.2e1;
t286 = t156 / 0.2e1;
t321 = Ifges(5,2) * t333;
t180 = t296 * t196;
t181 = t296 * t200;
t143 = t180 * t195 - t181 * t199;
t320 = -qJD(6) * t143 - t195 * t337 + t336 * t199;
t142 = t180 * t199 + t181 * t195;
t319 = qJD(6) * t142 + t336 * t195 + t199 * t337;
t247 = t191 ^ 2 + t193 ^ 2;
t317 = mrSges(4,3) * t247;
t50 = -mrSges(7,1) * t232 + mrSges(7,2) * t90;
t313 = m(7) * t60 + t50;
t214 = t195 * t196 - t199 * t200;
t123 = t214 * t170;
t174 = t195 * t200 + t196 * t199;
t108 = t174 * t162;
t308 = qJD(5) + qJD(6);
t133 = t308 * t174;
t312 = t108 - t133;
t109 = t214 * t162;
t132 = t308 * t214;
t311 = t109 - t132;
t309 = t17 * t200 - t18 * t196;
t306 = t18 * mrSges(6,1) - t17 * mrSges(6,2) + Ifges(6,5) * t100 + Ifges(6,6) * t101 + t325;
t305 = m(4) / 0.2e1;
t304 = Ifges(7,4) * t302 + Ifges(7,2) * t301 + Ifges(7,6) * t286;
t303 = Ifges(7,1) * t302 + Ifges(7,4) * t301 + Ifges(7,5) * t286;
t299 = t232 / 0.2e1;
t297 = t90 / 0.2e1;
t293 = t100 / 0.2e1;
t292 = t101 / 0.2e1;
t291 = -t144 / 0.2e1;
t290 = -t145 / 0.2e1;
t287 = t154 / 0.2e1;
t285 = -t159 / 0.2e1;
t277 = Ifges(5,4) * t163;
t160 = -t191 * t252 + t193 * t194;
t161 = t191 * t194 + t193 * t252;
t215 = t201 * t160 - t161 * t197;
t57 = qJD(4) * t81 + t168 * t170;
t273 = t215 * t57;
t272 = t310 * t57;
t261 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t144 + mrSges(6,2) * t145 + mrSges(5,3) * t163;
t230 = -mrSges(4,1) * t193 + mrSges(4,2) * t191;
t248 = -mrSges(5,1) * t162 + mrSges(5,2) * t163 + qJD(2) * t230;
t240 = t50 + t261;
t238 = t192 ^ 2 * t244;
t110 = t156 * mrSges(5,1) - t155 * mrSges(5,2);
t234 = t247 * t168;
t229 = mrSges(6,1) * t200 - mrSges(6,2) * t196;
t226 = Ifges(6,1) * t196 + t274;
t224 = Ifges(6,2) * t200 + t275;
t222 = Ifges(6,5) * t196 + Ifges(6,6) * t200;
t220 = -t17 * t196 - t18 * t200;
t114 = t160 * t197 + t161 * t201;
t213 = -t114 * t200 + t196 * t251;
t95 = -t114 * t196 - t200 * t251;
t51 = t195 * t213 + t199 * t95;
t52 = t195 * t95 - t199 * t213;
t218 = t43 * t196 - t44 * t200;
t103 = -mrSges(6,2) * t159 + mrSges(6,3) * t144;
t104 = mrSges(6,1) * t159 - mrSges(6,3) * t145;
t217 = t103 * t200 - t104 * t196;
t216 = -(-t177 * t191 + t183) * t191 + t147 * t193;
t210 = t216 * t202;
t205 = t12 * mrSges(7,1) + t157 * mrSges(5,1) + t43 * mrSges(6,1) - Ifges(5,6) * qJD(4) + t321 - t277 / 0.2e1 + t154 * Ifges(7,3) + t90 * Ifges(7,5) + t232 * Ifges(7,6) + t159 * Ifges(6,3) + t145 * Ifges(6,5) + t144 * Ifges(6,6) - t13 * mrSges(7,2) - t44 * mrSges(6,2);
t203 = qJD(2) ^ 2;
t187 = -pkin(5) * t200 - pkin(4);
t172 = -qJD(2) * pkin(2) + t221;
t152 = Ifges(6,3) * t156;
t148 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t162;
t122 = t174 * t170;
t107 = pkin(5) * t255 - t310;
t79 = qJD(2) * t208 + qJD(4) * t114;
t78 = -qJD(2) * t207 + qJD(4) * t215;
t72 = -mrSges(6,2) * t156 + mrSges(6,3) * t101;
t71 = mrSges(6,1) * t156 - mrSges(6,3) * t100;
t67 = mrSges(7,1) * t154 - mrSges(7,3) * t90;
t66 = -mrSges(7,2) * t154 + mrSges(7,3) * t232;
t65 = pkin(5) * t211 + t106;
t64 = pkin(5) * t256 + t81;
t58 = -mrSges(6,1) * t101 + mrSges(6,2) * t100;
t49 = t100 * Ifges(6,1) + t101 * Ifges(6,4) + t156 * Ifges(6,5);
t48 = t100 * Ifges(6,4) + t101 * Ifges(6,2) + t156 * Ifges(6,6);
t47 = t123 * t308 + t174 * t164;
t46 = -t133 * t170 + t164 * t214;
t39 = qJD(5) * t213 - t196 * t78 + t200 * t235;
t38 = qJD(5) * t95 + t196 * t235 + t200 * t78;
t34 = -pkin(5) * t101 + t57;
t26 = -mrSges(7,2) * t156 + mrSges(7,3) * t31;
t25 = mrSges(7,1) * t156 - mrSges(7,3) * t30;
t16 = t199 * t35 - t263;
t15 = -t195 * t35 - t262;
t14 = -mrSges(7,1) * t31 + mrSges(7,2) * t30;
t7 = -qJD(6) * t52 - t195 * t38 + t199 * t39;
t6 = qJD(6) * t51 + t195 * t39 + t199 * t38;
t1 = [-t114 * t156 * mrSges(5,3) + t38 * t103 + t39 * t104 + t78 * t148 + t51 * t25 + t52 * t26 + t6 * t66 + t7 * t67 + t95 * t71 - t213 * t72 + t240 * t79 - (-mrSges(5,3) * t155 + t14 + t58) * t215 + ((-mrSges(3,1) * t203 + qJD(2) * t248) * t198 + (-t110 + (-mrSges(3,2) + t317) * t203) * t202) * t192 + m(6) * (-t17 * t213 + t18 * t95 + t38 * t44 + t39 * t43 + t76 * t79 - t273) + m(7) * (t12 * t7 + t13 * t6 + t2 * t52 - t215 * t34 + t3 * t51 + t60 * t79) + m(5) * (t114 * t56 + t78 * t81 - t79 * t80 - t273) + m(4) * (-t160 * t191 + t161 * t193) * t168 + 0.2e1 * (t192 * t210 * t305 + (m(5) * (t157 * t192 - t238) / 0.2e1 + (t172 * t192 - t238) * t305) * t198) * qJD(2); (-t138 * t156 + t155 * t310 + t164 * t80 - t165 * t81 - t169 * t56 + t170 * t57) * mrSges(5,3) - t310 * t58 + (-Ifges(7,5) * t123 - Ifges(7,6) * t122) * t286 + (-Ifges(7,4) * t123 - Ifges(7,2) * t122) * t301 + (-Ifges(7,1) * t123 - Ifges(7,4) * t122) * t302 + (-t12 * t46 - t122 * t2 + t123 * t3 + t13 * t47) * mrSges(7,3) + t34 * (mrSges(7,1) * t122 - mrSges(7,2) * t123) + t329 * t104 + (t17 * t74 + t18 * t73 + t327 * t76 + t329 * t43 + t330 * t44 - t272) * m(6) + t330 * t103 + (-(t172 * t198 + t210) * t246 - pkin(2) * t231 + qJ(3) * t234 + qJD(3) * t216) * m(4) - (t322 + t326) * t164 + (t138 * t56 - t157 * t237 + t186 * t231 - t327 * t80 + t328 * t81 - t272) * m(5) + t328 * t148 + t331 * t67 + (t107 * t34 + t2 * t24 + t23 * t3 + (-t139 + t65) * t60 + t332 * t13 + t331 * t12) * m(7) + t332 * t66 + (t155 * t169 - t156 * t170 + t164 * t333 - t163 * t165 / 0.2e1) * Ifges(5,4) + (t223 * t286 + t57 * t228 + t227 * t293 + t225 * t292 + mrSges(5,2) * t231 - Ifges(5,1) * t155 + t48 * t284 + t49 * t283 + t220 * mrSges(6,3) + (-t200 * t69 / 0.2e1 + t70 * t284 + t76 * t229 + t224 * t291 + t226 * t290 + t222 * t285 + t218 * mrSges(6,3)) * qJD(5)) * t170 + ((mrSges(5,1) * t169 + t230) * qJD(2) - t248) * t237 - t240 * t139 + t261 * t106 + (t321 + t205) * t165 + (qJD(2) * t221 * t247 + t234) * mrSges(4,3) + (t151 / 0.2e1 + t152 / 0.2e1 + (Ifges(6,3) / 0.2e1 + Ifges(5,2) + Ifges(7,3) / 0.2e1) * t156 + t306) * t169 + t186 * t110 + t47 * t323 + (Ifges(7,5) * t46 + Ifges(7,6) * t47) * t287 + (Ifges(7,1) * t46 + Ifges(7,4) * t47) * t297 + (Ifges(7,4) * t46 + Ifges(7,2) * t47) * t299 - t123 * t303 - t122 * t304 + t23 * t25 + t24 * t26 + t46 * t42 / 0.2e1 + t60 * (-mrSges(7,1) * t47 + mrSges(7,2) * t46) + t65 * t50 + t73 * t71 + t74 * t72 + t107 * t14; -t214 * t25 + t174 * t26 + t196 * t72 + t200 * t71 + t312 * t67 + t311 * t66 + t217 * qJD(5) + (m(4) + m(5)) * t231 - t203 * t317 - t240 * t163 + (-t148 - t217) * t162 - m(5) * (t162 * t81 - t163 * t80) - m(4) * t216 * qJD(2) + t110 + (t12 * t312 + t13 * t311 - t163 * t60 + t174 * t2 - t214 * t3) * m(7) + (-t159 * t218 - t163 * t76 - t220) * m(6); (-mrSges(7,1) * t312 + mrSges(7,2) * t311) * t60 + (-t12 * t311 + t13 * t312 - t174 * t3 - t2 * t214) * mrSges(7,3) + (-Ifges(7,4) * t109 - Ifges(7,2) * t108) * t300 + (-Ifges(7,5) * t109 - Ifges(7,6) * t108) * t288 + (-Ifges(7,1) * t109 - Ifges(7,4) * t108) * t298 + (t80 * mrSges(5,3) - t158 / 0.2e1 + (Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t163 - t326) * t162 + (-Ifges(7,4) * t132 - Ifges(7,2) * t133) * t299 + (-mrSges(5,1) - t229) * t57 + (t277 / 0.2e1 + t81 * mrSges(5,3) - t205) * t163 + t48 * t283 + (Ifges(7,5) * t174 - Ifges(7,6) * t214 + t222) * t286 + t34 * (mrSges(7,1) * t214 + mrSges(7,2) * t174) + (Ifges(7,4) * t174 - Ifges(7,2) * t214) * t301 + (Ifges(7,1) * t174 - Ifges(7,4) * t214) * t302 - t214 * t304 + (t109 / 0.2e1 - t132 / 0.2e1) * t42 + (t108 / 0.2e1 - t133 / 0.2e1) * t41 + (-Ifges(7,5) * t132 - Ifges(7,6) * t133) * t287 + (-Ifges(7,1) * t132 - Ifges(7,4) * t133) * t297 - t261 * t81 + (t313 * t196 * pkin(5) + t335) * qJD(5) + t319 * t66 + t320 * t67 + (t12 * t320 + t13 * t319 + t142 * t3 + t143 * t2 + t187 * t34 - t60 * t64) * m(7) + (m(6) * t309 + (-m(6) * t219 - t196 * t103 - t200 * t104) * qJD(5) - t196 * t71 + t200 * t72) * pkin(9) + t309 * mrSges(6,3) + t196 * t49 / 0.2e1 + t187 * t14 + (-pkin(4) * t57 - t43 * t54 - t44 * t55 - t76 * t81) * m(6) + t224 * t292 + t226 * t293 + t174 * t303 - t56 * mrSges(5,2) - pkin(4) * t58 - t64 * t50 - t55 * t103 - t54 * t104 + t142 * t25 + t143 * t26 - t80 * t148 - Ifges(5,5) * t155 - Ifges(5,6) * t156; t152 - m(7) * (t12 * t15 + t13 * t16) + (Ifges(6,5) * t144 - Ifges(6,6) * t145) * t285 + t90 * t323 + t306 + (t195 * t26 + t199 * t25 + m(7) * (t195 * t2 + t199 * t3) - t313 * t145 + (-t195 * t67 + t199 * t66 + m(7) * (-t12 * t195 + t13 * t199)) * qJD(6)) * pkin(5) + (-Ifges(6,2) * t145 + t141 + t70) * t291 + t69 * t289 + (Ifges(6,1) * t144 - t276) * t290 + (t144 * t43 + t145 * t44) * mrSges(6,3) - t16 * t66 - t15 * t67 - t43 * t103 + t44 * t104 - t76 * (mrSges(6,1) * t145 + mrSges(6,2) * t144) + t334; -t12 * t66 + t13 * t67 + t41 * t297 + t325 + t334;];
tauc  = t1(:);
