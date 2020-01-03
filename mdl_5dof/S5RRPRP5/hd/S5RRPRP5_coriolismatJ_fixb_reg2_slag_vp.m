% Calculate inertial parameters regressor of coriolis matrix for
% S5RRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% cmat_reg [(5*5)x(5*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRPRP5_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP5_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:55:02
% EndTime: 2019-12-31 19:55:09
% DurationCPUTime: 3.64s
% Computational Cost: add. (5649->264), mult. (10871->306), div. (0->0), fcn. (12324->6), ass. (0->213)
t182 = qJD(2) + qJD(4);
t185 = sin(qJ(2));
t186 = cos(qJ(2));
t301 = cos(pkin(8));
t225 = t301 * t186;
t300 = sin(pkin(8));
t170 = -t185 * t300 + t225;
t224 = t300 * t186;
t171 = -t185 * t301 - t224;
t184 = sin(qJ(4));
t319 = cos(qJ(4));
t140 = t170 * t319 + t184 * t171;
t294 = t140 ^ 2;
t166 = t319 * t171;
t288 = t184 * t170;
t333 = -t166 + t288;
t344 = t333 ^ 2;
t43 = t344 - t294;
t346 = t43 * qJD(1);
t335 = t333 * qJD(1);
t345 = t140 * t335;
t341 = t140 * qJ(5);
t343 = t333 * pkin(4);
t76 = -t341 + t343;
t322 = -t140 / 0.2e1;
t323 = t333 / 0.2e1;
t178 = -pkin(2) * t186 - pkin(1);
t148 = -pkin(3) * t170 + t178;
t216 = -pkin(4) * t140 - qJ(5) * t333;
t65 = t148 + t216;
t342 = t65 * t140;
t307 = t65 * t333;
t255 = t140 * qJD(4);
t338 = qJD(2) * t140 + t255;
t337 = t140 * qJD(1);
t336 = t140 * qJD(3);
t254 = t140 * qJD(5);
t334 = t182 * t43;
t221 = pkin(2) * t301 + pkin(3);
t234 = t300 * pkin(2);
t164 = t184 * t221 + t234 * t319;
t313 = -qJ(3) - pkin(6);
t218 = t313 * t300;
t206 = t185 * t218;
t147 = -t225 * t313 + t206;
t316 = t170 * pkin(7);
t194 = t147 + t316;
t189 = t319 * t194;
t219 = t313 * t301;
t207 = t185 * t219;
t145 = t224 * t313 + t207;
t314 = t171 * pkin(7);
t193 = t145 + t314;
t191 = t184 * t193;
t187 = t191 / 0.2e1 + t189 / 0.2e1;
t144 = t186 * t219 - t206;
t200 = t144 - t316;
t106 = t319 * t200;
t146 = t186 * t218 + t207;
t111 = t146 + t314;
t289 = t184 * t111;
t217 = t106 / 0.2e1 - t289 / 0.2e1;
t38 = t187 + t217;
t283 = t38 * qJD(1);
t24 = qJD(2) * t164 + t283;
t163 = t184 * t234 - t221 * t319;
t107 = t319 * t193;
t192 = t184 * t194;
t188 = -t107 / 0.2e1 + t192 / 0.2e1;
t195 = t184 * t200;
t235 = t319 * t111;
t190 = -t195 / 0.2e1 - t235 / 0.2e1;
t35 = t188 - t190;
t302 = -qJD(1) * t35 - qJD(2) * t163;
t68 = t192 - t107;
t70 = t189 + t191;
t215 = t140 * t70 + t333 * t68;
t332 = qJD(3) * t215;
t67 = -t106 + t289;
t69 = t235 + t195;
t199 = (t67 - t70) * t333 + (t68 + t69) * t140;
t331 = t199 * qJD(1);
t210 = t344 + t294;
t330 = t210 * qJD(1);
t329 = t215 * qJD(1);
t222 = -t166 / 0.2e1;
t135 = t222 + t166 / 0.2e1;
t37 = t187 - t217;
t239 = -qJD(2) * t37 + qJD(3) * t135 - qJD(4) * t70;
t328 = qJD(2) * t199 + qJD(3) * t210;
t327 = t171 ^ 2;
t326 = pkin(4) / 0.2e1;
t325 = -qJ(5) / 0.2e1;
t324 = -t333 / 0.2e1;
t158 = qJ(5) + t164;
t321 = -t158 / 0.2e1;
t159 = -pkin(4) + t163;
t320 = -t159 / 0.2e1;
t315 = t171 * pkin(3);
t181 = t185 * pkin(2);
t97 = 0.2e1 * t222 + t288;
t312 = qJD(2) * t38 - qJD(3) * t97;
t311 = -qJD(2) * t67 - qJD(4) * t37;
t310 = qJD(2) * pkin(2);
t220 = t67 * t68 + t69 * t70;
t150 = t181 - t315;
t66 = t150 + t76;
t3 = t65 * t66 + t220;
t309 = t3 * qJD(1);
t4 = t65 * t76;
t308 = t4 * qJD(1);
t304 = -qJD(3) * t333 - qJD(4) * t38;
t11 = t148 * t150 + t220;
t298 = t11 * qJD(1);
t226 = t163 / 0.2e1 + t320;
t227 = t321 + t164 / 0.2e1;
t196 = -t140 * t226 + t227 * t333;
t204 = pkin(4) * t322 + qJ(5) * t324;
t14 = t196 - t204;
t296 = t14 * qJD(1);
t15 = -t140 * t66 + t307;
t292 = t15 * qJD(1);
t16 = -t333 * t66 - t342;
t291 = t16 * qJD(1);
t19 = -t140 * t76 + t307;
t287 = t19 * qJD(1);
t20 = -t333 * t76 - t342;
t286 = t20 * qJD(1);
t179 = t181 / 0.2e1;
t223 = t179 - t315 / 0.2e1;
t25 = (t325 + t321) * t140 + (t326 + t320) * t333 + t223;
t285 = t25 * qJD(1);
t39 = (t144 + t147) * t171 + (-t145 + t146) * t170;
t282 = t39 * qJD(1);
t41 = 0.2e1 * pkin(4) * t323 + 0.2e1 * qJ(5) * t322;
t281 = t41 * qJD(1);
t202 = t163 * t324 + t164 * t322;
t47 = t202 + t223;
t276 = t47 * qJD(1);
t49 = t144 * t145 + t146 * t147 + t178 * t181;
t275 = t49 * qJD(1);
t52 = -t140 * t150 + t148 * t333;
t274 = t52 * qJD(1);
t53 = t140 * t148 + t150 * t333;
t273 = t53 * qJD(1);
t74 = t145 * t171 + t147 * t170;
t270 = t74 * qJD(1);
t268 = t97 * qJD(1);
t266 = qJD(1) * t148;
t265 = qJD(1) * t186;
t197 = (t300 * t170 / 0.2e1 + t301 * t171 / 0.2e1) * pkin(2);
t109 = -t181 / 0.2e1 + t197;
t264 = t109 * qJD(1);
t169 = t170 ^ 2;
t110 = t169 - t327;
t263 = t110 * qJD(1);
t112 = -t170 * t181 - t171 * t178;
t262 = t112 * qJD(1);
t113 = t170 * t178 - t171 * t181;
t261 = t113 * qJD(1);
t260 = t135 * qJD(1);
t259 = t135 * qJD(2);
t118 = t135 * qJD(4);
t258 = t344 * qJD(1);
t256 = t333 * qJD(2);
t251 = t333 * qJD(4);
t143 = t169 + t327;
t250 = t143 * qJD(1);
t248 = t163 * qJD(4);
t246 = t170 * qJD(1);
t245 = t171 * qJD(1);
t244 = t171 * qJD(2);
t175 = -t185 ^ 2 + t186 ^ 2;
t243 = t175 * qJD(1);
t242 = t185 * qJD(2);
t241 = t186 * qJD(2);
t240 = -t248 + qJD(5);
t238 = pkin(1) * t185 * qJD(1);
t237 = pkin(1) * t265;
t236 = t65 * t335;
t232 = t140 * t266;
t231 = t333 * t266;
t230 = t170 * t245;
t229 = t170 * t244;
t228 = t185 * t241;
t198 = -t226 * t70 + t227 * t68;
t205 = t325 * t69 + t326 * t67;
t2 = t198 + t205;
t75 = -t158 * t163 + t159 * t164;
t214 = qJD(1) * t2 + qJD(2) * t75;
t213 = qJD(2) * t35 + t336;
t36 = t188 + t190;
t212 = -qJD(2) * t36 - qJD(4) * t68;
t211 = qJD(2) * t69 - qJD(4) * t36;
t72 = qJD(4) * t97 + t256;
t209 = qJD(2) * t97 + t251;
t208 = -qJD(4) * t35 + t336;
t201 = (t251 + t256) * t140;
t183 = qJ(5) * qJD(5);
t176 = t185 * t265;
t174 = t182 * qJ(5);
t168 = t170 * qJD(2);
t155 = t164 * qJD(4);
t149 = t158 * qJD(5);
t108 = t179 + t197;
t79 = t333 * t337;
t48 = -t202 + t223;
t42 = t338 * t333;
t26 = t140 * t158 / 0.2e1 + t159 * t323 + t343 / 0.2e1 - t341 / 0.2e1 + t223;
t21 = -t155 - t24;
t13 = t196 + t204;
t1 = t198 - t205;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t228, t175 * qJD(2), 0, -t228, 0, 0, -pkin(1) * t242, -pkin(1) * t241, 0, 0, -t229, t110 * qJD(2), 0, t229, 0, 0, t112 * qJD(2), t113 * qJD(2), qJD(2) * t39 + qJD(3) * t143, qJD(2) * t49 + qJD(3) * t74, t42, -t334, 0, -t201, 0, 0, qJD(2) * t52 + t148 * t251, qJD(2) * t53 + t148 * t255, t328, qJD(2) * t11 + t332, t42, 0, t334, 0, 0, -t201, qJD(2) * t15 + qJD(4) * t19 + t254 * t333, t328, qJD(2) * t16 + qJD(4) * t20 + qJD(5) * t344, qJD(2) * t3 + qJD(4) * t4 - qJD(5) * t307 + t332; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t176, t243, t241, -t176, -t242, 0, -pkin(6) * t241 - t238, pkin(6) * t242 - t237, 0, 0, -t230, t263, t168, t230, t244, 0, qJD(2) * t144 + t262, -qJD(2) * t146 + t261, t282 + (-t170 * t301 + t171 * t300) * t310, t275 + (t144 * t301 + t146 * t300) * t310 + t108 * qJD(3), t79, -t346, t338, -t345, -t72, 0, t274 + t311, -t211 + t273, t331 + (t140 * t163 - t164 * t333) * qJD(2), t298 + (t163 * t67 + t164 * t69) * qJD(2) + t48 * qJD(3), t79, t338, t346, 0, t72, -t345, t292 + t311, t331 + (t140 * t159 - t158 * t333) * qJD(2) + t13 * qJD(4) + t254, t211 + t291, t309 + (t158 * t69 + t159 * t67) * qJD(2) + t26 * qJD(3) + t1 * qJD(4) + t37 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t250, qJD(2) * t108 + t270, 0, 0, 0, 0, 0, 0, t118, 0, t330, qJD(2) * t48 + t329, 0, 0, 0, 0, 0, 0, t118, t330, 0, qJD(2) * t26 - qJD(5) * t135 + t329; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t345, -t346, t338, -t345, -t209, 0, t231 + t239, -t212 + t232, 0, 0, t345, t338, t346, 0, t209, -t345, t239 + t287, t13 * qJD(2) + qJD(4) * t216 + t254, t212 + t286, t308 + t1 * qJD(2) + (-pkin(4) * t70 - qJ(5) * t68) * qJD(4) + t70 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t345, t338, t258, -t236 - t239; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t176, -t243, 0, t176, 0, 0, t238, t237, 0, 0, t230, -t263, 0, -t230, 0, 0, qJD(3) * t171 - t262, -qJD(3) * t170 - t261, -t282, qJD(3) * t109 - t275, -t79, t346, 0, t345, -t118, 0, -t274 + t304, -t208 - t273, -t331, -qJD(3) * t47 - t298, -t79, 0, -t346, 0, t118, t345, -t292 + t304, qJD(4) * t14 - t331, t208 - t291, -qJD(3) * t25 + qJD(4) * t2 + qJD(5) * t38 - t309; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t155, t248, 0, 0, 0, 0, 0, 0, 0, 0, -t155, 0, t240, qJD(4) * t75 + t149; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t245, -t246, 0, t264, 0, 0, 0, 0, 0, 0, -t335, -t337, 0, -t276, 0, 0, 0, 0, 0, 0, -t335, 0, t337, -t285; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t260, 0, t21, t248 - t302, 0, 0, 0, 0, 0, 0, t260, 0, t21, t296, t240 + t302, (-pkin(4) * t164 - qJ(5) * t163) * qJD(4) + t149 + t214; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t182, t158 * t182 + t283; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t244, t168, -t250, -qJD(2) * t109 - t270, 0, 0, 0, 0, 0, 0, t72, t338, -t330, qJD(2) * t47 - t329, 0, 0, 0, 0, 0, 0, t72, -t330, -t338, qJD(2) * t25 + qJD(4) * t41 - qJD(5) * t97 - t329; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t245, t246, 0, -t264, 0, 0, 0, 0, 0, 0, t335, t337, 0, t276, 0, 0, 0, 0, 0, 0, t335, 0, -t337, t285; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t268, t337, 0, 0, 0, 0, 0, 0, 0, 0, t268, 0, -t337, t281; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t268; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t345, t346, 0, t345, t259, 0, -t231 + t312, -t213 - t232, 0, 0, -t345, 0, -t346, 0, -t259, t345, -t287 + t312, -qJD(2) * t14, t213 - t286, -qJD(2) * t2 - qJD(3) * t41 - t308; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t260, 0, t24, t302, 0, 0, 0, 0, 0, 0, -t260, 0, t24, -t296, qJD(5) - t302, t183 - t214; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t268, -t337, 0, 0, 0, 0, 0, 0, 0, 0, -t268, 0, t337, -t281; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t182, t174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t345, 0, -t258, t236 - t312; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t182, -qJ(5) * qJD(4) - qJD(2) * t158 - t283; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t268; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t182, -t174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t5;