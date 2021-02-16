% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% tauc_reg [5x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:25
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPR9_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR9_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:23:55
% EndTime: 2021-01-15 23:24:16
% DurationCPUTime: 4.66s
% Computational Cost: add. (4055->371), mult. (10373->539), div. (0->0), fcn. (7310->8), ass. (0->189)
t188 = sin(qJ(5));
t191 = cos(qJ(5));
t189 = sin(qJ(3));
t190 = sin(qJ(2));
t241 = qJD(1) * t190;
t223 = t189 * t241;
t192 = cos(qJ(3));
t233 = t192 * qJD(2);
t150 = t223 - t233;
t234 = t189 * qJD(2);
t152 = t192 * t241 + t234;
t186 = sin(pkin(9));
t187 = cos(pkin(9));
t205 = -t186 * t150 + t187 * t152;
t235 = qJD(5) * t188;
t95 = t187 * t150 + t186 * t152;
t267 = t191 * t95;
t193 = cos(qJ(2));
t230 = qJD(1) * qJD(2);
t220 = t193 * t230;
t229 = qJD(2) * qJD(3);
t110 = -qJD(3) * t223 + (t220 + t229) * t192;
t224 = t193 * t234;
t236 = qJD(3) * t192;
t225 = t190 * t236;
t280 = t224 + t225;
t111 = t280 * qJD(1) + t189 * t229;
t60 = t186 * t110 + t187 * t111;
t61 = t187 * t110 - t186 * t111;
t10 = -qJD(5) * t267 - t188 * t60 + t191 * t61 - t205 * t235;
t231 = t193 * qJD(1);
t171 = -qJD(3) + t231;
t164 = -qJD(5) + t171;
t44 = t188 * t205 + t267;
t265 = t44 * t164;
t291 = t10 - t265;
t206 = t188 * t95 - t191 * t205;
t290 = t206 * t44;
t197 = t206 * qJD(5) - t188 * t61 - t191 * t60;
t266 = t164 * t206;
t289 = t197 + t266;
t288 = t206 ^ 2 - t44 ^ 2;
t180 = pkin(6) * t231;
t162 = qJD(2) * pkin(7) + t180;
t157 = -t193 * pkin(2) - t190 * pkin(7) - pkin(1);
t139 = t157 * qJD(1);
t257 = t189 * t139;
t103 = t192 * t162 + t257;
t71 = -t150 * qJ(4) + t103;
t269 = t187 * t71;
t102 = t192 * t139 - t189 * t162;
t70 = -t152 * qJ(4) + t102;
t62 = -t171 * pkin(3) + t70;
t25 = t186 * t62 + t269;
t282 = t95 * pkin(8);
t16 = t25 - t282;
t15 = t16 * t235;
t176 = t190 * t230;
t208 = pkin(2) * t190 - pkin(7) * t193;
t154 = t208 * qJD(2);
t140 = qJD(1) * t154;
t212 = pkin(6) * t176;
t247 = -t192 * t140 - t189 * t212;
t198 = -t103 * qJD(3) - t247;
t23 = pkin(3) * t176 - t110 * qJ(4) - t152 * qJD(4) + t198;
t238 = qJD(3) * t189;
t203 = t139 * t236 + t189 * t140 - t162 * t238;
t196 = -t192 * t212 + t203;
t28 = -t111 * qJ(4) - t150 * qJD(4) + t196;
t6 = -t186 * t28 + t187 * t23;
t4 = pkin(4) * t176 - t61 * pkin(8) + t6;
t161 = -qJD(2) * pkin(2) + pkin(6) * t241;
t108 = t150 * pkin(3) + qJD(4) + t161;
t52 = t95 * pkin(4) + t108;
t287 = -t188 * t4 + t52 * t44 + t15;
t7 = t186 * t23 + t187 * t28;
t5 = -t60 * pkin(8) + t7;
t227 = -t188 * t5 + t191 * t4;
t285 = t52 * t206 + t227;
t153 = t208 * qJD(1);
t135 = t189 * t153;
t274 = qJ(4) + pkin(7);
t216 = qJD(3) * t274;
t232 = t192 * qJD(4);
t254 = t190 * t192;
t255 = t189 * t193;
t284 = -t189 * t216 + t232 - t135 - (-pkin(6) * t254 - qJ(4) * t255) * qJD(1);
t253 = t192 * t193;
t202 = pkin(3) * t190 - qJ(4) * t253;
t244 = pkin(6) * t223 + t192 * t153;
t283 = t202 * qJD(1) + t189 * qJD(4) + t192 * t216 + t244;
t144 = t186 * t192 + t187 * t189;
t133 = t144 * qJD(3);
t249 = t144 * t231 - t133;
t143 = t186 * t189 - t187 * t192;
t281 = t171 * t143;
t279 = -0.2e1 * t230;
t278 = pkin(8) * t205;
t273 = t284 * t186 + t283 * t187;
t272 = t283 * t186 - t284 * t187;
t209 = -t180 + (-t189 * t231 + t238) * pkin(3);
t275 = pkin(3) * t186;
t183 = t190 * pkin(6);
t91 = t191 * t143 + t188 * t144;
t271 = -t91 * qJD(5) + t249 * t188 + t281 * t191;
t92 = -t188 * t143 + t191 * t144;
t270 = t92 * qJD(5) + t281 * t188 - t249 * t191;
t173 = pkin(6) * t253;
t245 = t192 * t154 + t234 * t183;
t39 = -t190 * t232 + t202 * qJD(2) + (-t173 + (qJ(4) * t190 - t157) * t189) * qJD(3) + t245;
t246 = t189 * t154 + t157 * t236;
t48 = (-pkin(6) * qJD(2) - qJ(4) * qJD(3)) * t254 + (-qJD(4) * t190 + (-pkin(6) * qJD(3) - qJ(4) * qJD(2)) * t193) * t189 + t246;
t13 = t186 * t39 + t187 * t48;
t66 = t186 * t71;
t31 = t187 * t70 - t66;
t24 = t187 * t62 - t66;
t14 = -t171 * pkin(4) + t24 - t278;
t268 = t191 * t14;
t243 = t189 * t157 + t173;
t256 = t189 * t190;
t104 = -qJ(4) * t256 + t243;
t146 = t192 * t157;
t99 = -qJ(4) * t254 + t146 + (-pkin(6) * t189 - pkin(3)) * t193;
t50 = t187 * t104 + t186 * t99;
t264 = -t249 * pkin(4) + t209;
t263 = t110 * t189;
t262 = t150 * t171;
t261 = t152 * t171;
t260 = t161 * t189;
t259 = t161 * t192;
t258 = t171 * t192;
t195 = qJD(1) ^ 2;
t252 = t193 * t195;
t194 = qJD(2) ^ 2;
t251 = t194 * t190;
t250 = t194 * t193;
t159 = t274 * t189;
t160 = t274 * t192;
t107 = -t186 * t159 + t187 * t160;
t155 = pkin(3) * t256 + t183;
t184 = t190 ^ 2;
t242 = -t193 ^ 2 + t184;
t240 = qJD(2) * t190;
t239 = qJD(2) * t193;
t237 = qJD(3) * t190;
t113 = t280 * pkin(3) + pkin(6) * t239;
t178 = -t192 * pkin(3) - pkin(2);
t226 = t189 * t237;
t222 = t193 * t233;
t93 = t111 * pkin(3) + pkin(6) * t220;
t218 = qJD(5) * t14 + t5;
t12 = -t186 * t48 + t187 * t39;
t30 = -t186 * t70 - t269;
t49 = -t186 * t104 + t187 * t99;
t215 = pkin(1) * t279;
t106 = -t187 * t159 - t186 * t160;
t214 = t150 + t233;
t213 = -t152 + t234;
t82 = -t143 * pkin(8) + t107;
t211 = pkin(4) * t241 + t281 * pkin(8) + qJD(5) * t82 + t273;
t81 = -t144 * pkin(8) + t106;
t210 = t249 * pkin(8) + qJD(5) * t81 - t272;
t2 = t188 * t14 + t191 * t16;
t124 = t143 * t190;
t32 = -t193 * pkin(4) + t124 * pkin(8) + t49;
t123 = t144 * t190;
t34 = -t123 * pkin(8) + t50;
t207 = t188 * t32 + t191 * t34;
t73 = t191 * t123 - t188 * t124;
t74 = -t188 * t123 - t191 * t124;
t204 = qJD(1) * t184 - t171 * t193;
t175 = t187 * pkin(3) + pkin(4);
t201 = t188 * t175 + t191 * t275;
t200 = t191 * t175 - t188 * t275;
t117 = t143 * pkin(4) + t178;
t100 = t123 * pkin(4) + t155;
t76 = t190 * t133 + t186 * t224 - t187 * t222;
t75 = t143 * t237 - t144 * t239;
t63 = t152 * pkin(3) + pkin(4) * t205;
t51 = -t75 * pkin(4) + t113;
t33 = t60 * pkin(4) + t93;
t20 = t31 - t278;
t19 = t30 + t282;
t18 = t74 * qJD(5) - t188 * t76 - t191 * t75;
t17 = -t73 * qJD(5) + t188 * t75 - t191 * t76;
t9 = t75 * pkin(8) + t13;
t8 = pkin(4) * t240 + t76 * pkin(8) + t12;
t1 = -t188 * t16 + t268;
t3 = [0, 0, 0, 0.2e1 * t193 * t176, t242 * t279, t250, -t251, 0, -pkin(6) * t250 + t190 * t215, pkin(6) * t251 + t193 * t215, t110 * t254 + (t222 - t226) * t152, (-t150 * t192 - t152 * t189) * t239 + (-t263 - t111 * t192 + (t150 * t189 - t152 * t192) * qJD(3)) * t190, t171 * t226 - t110 * t193 + (t152 * t190 + t192 * t204) * qJD(2), t171 * t225 + t111 * t193 + (-t150 * t190 - t189 * t204) * qJD(2), (-t171 - t231) * t240, -(-t157 * t238 + t245) * t171 + (t161 * t236 + pkin(6) * t111 + (qJD(1) * t146 + t102) * qJD(2)) * t190 + ((pkin(6) * t150 + t260) * qJD(2) + (t257 + (pkin(6) * t171 + t162) * t192) * qJD(3) + t247) * t193, (-t193 * pkin(6) * t238 + t246) * t171 + t203 * t193 + (pkin(6) * t110 - t161 * t238) * t190 + ((pkin(6) * t152 + t259) * t193 + (-pkin(6) * t258 - t243 * qJD(1) - t103) * t190) * qJD(2), -t108 * t75 + t113 * t95 - t12 * t171 + t93 * t123 + t155 * t60 - t6 * t193 + (qJD(1) * t49 + t24) * t240, -t108 * t76 + t113 * t205 - t93 * t124 + t13 * t171 + t155 * t61 + t7 * t193 + (-qJD(1) * t50 - t25) * t240, -t12 * t205 - t7 * t123 + t6 * t124 - t13 * t95 + t24 * t76 + t25 * t75 - t49 * t61 - t50 * t60, t108 * t113 + t24 * t12 + t25 * t13 + t93 * t155 + t6 * t49 + t7 * t50, t10 * t74 - t17 * t206, -t10 * t73 - t17 * t44 + t18 * t206 + t197 * t74, -t10 * t193 - t17 * t164 + (qJD(1) * t74 - t206) * t240, -t197 * t193 + t18 * t164 + (-qJD(1) * t73 - t44) * t240, (-t164 - t231) * t240, -(-t188 * t9 + t191 * t8) * t164 - t227 * t193 + t51 * t44 - t100 * t197 + t33 * t73 + t52 * t18 + (t164 * t207 + t193 * t2) * qJD(5) + ((-t188 * t34 + t191 * t32) * qJD(1) + t1) * t240, t100 * t10 - t15 * t193 + t52 * t17 + t33 * t74 - t51 * t206 + ((-qJD(5) * t34 + t8) * t164 + t4 * t193) * t188 + ((qJD(5) * t32 + t9) * t164 + t218 * t193) * t191 + (-qJD(1) * t207 - t2) * t240; 0, 0, 0, -t190 * t252, t242 * t195, 0, 0, 0, t195 * pkin(1) * t190, pkin(1) * t252, -t152 * t258 + t263, (t110 + t262) * t192 + (-t111 + t261) * t189, -t171 * t236 + (t171 * t253 + t190 * t213) * qJD(1), t171 * t238 + (-t171 * t255 + t190 * t214) * qJD(1), t171 * t241, -pkin(2) * t111 + t244 * t171 + (pkin(7) * t258 + t260) * qJD(3) + ((-pkin(7) * t234 - t102) * t190 + (-pkin(6) * t214 - t260) * t193) * qJD(1), -pkin(2) * t110 - t135 * t171 + (-t189 * pkin(7) * t171 + t259) * qJD(3) + (-t161 * t253 + (-pkin(7) * t233 + t103) * t190 + (t171 * t254 + t193 * t213) * pkin(6)) * qJD(1), t93 * t143 + t178 * t60 + t209 * t95 + t273 * t171 - t249 * t108 + (qJD(2) * t106 - t24) * t241, t93 * t144 + t178 * t61 + t209 * t205 - t272 * t171 + t281 * t108 + (-qJD(2) * t107 + t25) * t241, -t106 * t61 - t107 * t60 - t7 * t143 - t6 * t144 + t205 * t273 - t24 * t281 + t249 * t25 + t272 * t95, t6 * t106 + t7 * t107 + t209 * t108 + t93 * t178 - t273 * t24 - t272 * t25, t10 * t92 - t206 * t271, -t10 * t91 + t197 * t92 + t206 * t270 - t271 * t44, -t271 * t164 + (qJD(2) * t92 + t206) * t241, t270 * t164 + (-qJD(2) * t91 + t44) * t241, t164 * t241, -t117 * t197 + t33 * t91 + t270 * t52 + t264 * t44 + (t188 * t210 + t191 * t211) * t164 + ((-t188 * t82 + t191 * t81) * qJD(2) - t1) * t241, t117 * t10 + t33 * t92 + t271 * t52 - t264 * t206 + (-t188 * t211 + t191 * t210) * t164 + (-(t188 * t81 + t191 * t82) * qJD(2) + t2) * t241; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t152 * t150, -t150 ^ 2 + t152 ^ 2, t110 - t262, -t111 - t261, t176, -t103 * t171 - t161 * t152 + t198, -t102 * t171 + t161 * t150 - t196, -t108 * t205 + t30 * t171 + (-t152 * t95 + t176 * t187) * pkin(3) + t6, t108 * t95 - t31 * t171 + (-t152 * t205 - t176 * t186) * pkin(3) - t7, (-t186 * t60 - t187 * t61) * pkin(3) + (t25 + t30) * t205 + (-t24 + t31) * t95, -t24 * t30 - t25 * t31 + (-t108 * t152 + t186 * t7 + t187 * t6) * pkin(3), -t290, t288, t291, t289, t176, t200 * t176 + (-t188 * t20 + t191 * t19) * t164 - t63 * t44 + (t164 * t201 - t2) * qJD(5) + t285, -t201 * t176 - t191 * t5 - (t188 * t19 + t191 * t20) * t164 + t63 * t206 + (t164 * t200 - t268) * qJD(5) + t287; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t171 * t205 + t60, t171 * t95 + t61, -t205 ^ 2 - t95 ^ 2, t205 * t24 + t25 * t95 + t93, 0, 0, 0, 0, 0, -t197 + t266, t10 + t265; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t290, t288, t291, t289, t176, (-qJD(5) - t164) * t2 + t285, -t1 * t164 - t191 * t218 + t287;];
tauc_reg = t3;
