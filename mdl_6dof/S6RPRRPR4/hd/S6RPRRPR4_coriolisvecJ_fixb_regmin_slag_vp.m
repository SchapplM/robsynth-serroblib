% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% tauc_reg [6x32]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRPR4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:10:16
% EndTime: 2019-03-09 05:10:27
% DurationCPUTime: 4.85s
% Computational Cost: add. (8591->362), mult. (22916->500), div. (0->0), fcn. (18570->10), ass. (0->202)
t188 = cos(pkin(10));
t193 = cos(qJ(3));
t244 = t188 * t193;
t186 = sin(pkin(10));
t191 = sin(qJ(3));
t245 = t186 * t191;
t211 = -t244 + t245;
t149 = t211 * qJD(1);
t161 = t186 * t193 + t188 * t191;
t150 = t161 * qJD(1);
t190 = sin(qJ(4));
t272 = cos(qJ(4));
t131 = t272 * t149 + t150 * t190;
t129 = qJD(6) + t131;
t185 = sin(pkin(11));
t187 = cos(pkin(11));
t189 = sin(qJ(6));
t192 = cos(qJ(6));
t160 = t185 * t192 + t187 * t189;
t243 = t192 * t187;
t246 = t185 * t189;
t158 = -t243 + t246;
t311 = t129 * t158;
t234 = qJD(1) * qJD(3);
t227 = t193 * t234;
t172 = t188 * t227;
t228 = t191 * t234;
t146 = -t186 * t228 + t172;
t207 = -t190 * t149 + t272 * t150;
t242 = t186 * t227 + t188 * t228;
t96 = t207 * qJD(4) + t190 * t146 + t272 * t242;
t312 = -t311 * t129 + t160 * t96;
t310 = t129 * t160;
t297 = t131 * t185;
t309 = pkin(5) * t297;
t308 = pkin(9) * t297;
t221 = -t310 * t129 - t158 * t96;
t184 = qJD(3) + qJD(4);
t121 = t184 * t185 + t187 * t207;
t119 = -t187 * t184 + t185 * t207;
t293 = t192 * t119;
t72 = t121 * t189 + t293;
t266 = t207 * t72;
t307 = t221 + t266;
t296 = t131 * t187;
t270 = pkin(7) + qJ(2);
t167 = t270 * t186;
t162 = qJD(1) * t167;
t169 = t270 * t188;
t163 = qJD(1) * t169;
t213 = t162 * t191 - t163 * t193;
t124 = -pkin(8) * t149 - t213;
t117 = t272 * t124;
t280 = -t193 * t162 - t163 * t191;
t123 = -pkin(8) * t150 + t280;
t118 = qJD(3) * pkin(3) + t123;
t70 = t190 * t118 + t117;
t64 = qJ(5) * t184 + t70;
t179 = -pkin(2) * t188 - pkin(1);
t164 = t179 * qJD(1) + qJD(2);
t137 = t149 * pkin(3) + t164;
t71 = t131 * pkin(4) - qJ(5) * t207 + t137;
t38 = -t185 * t64 + t187 * t71;
t104 = -t242 * pkin(8) - qJD(2) * t149 + t280 * qJD(3);
t203 = t161 * qJD(2);
t201 = qJD(1) * t203;
t105 = -t146 * pkin(8) + t213 * qJD(3) - t201;
t229 = qJD(4) * t272;
t239 = qJD(4) * t190;
t196 = t272 * t104 + t190 * t105 + t118 * t229 - t124 * t239;
t27 = t184 * qJD(5) + t196;
t95 = t272 * t146 - t149 * t229 - t150 * t239 - t190 * t242;
t37 = t242 * pkin(3) + t96 * pkin(4) - t95 * qJ(5) - qJD(5) * t207;
t10 = t185 * t37 + t187 * t27;
t7 = t10 * t187;
t306 = -t38 * t296 + t7;
t305 = pkin(5) * t207 + pkin(9) * t296;
t265 = t185 * t95;
t29 = t190 * t104 - t272 * t105 + t118 * t239 + t124 * t229;
t15 = pkin(5) * t265 + t29;
t16 = pkin(5) * t131 - pkin(9) * t121 + t38;
t39 = t185 * t71 + t187 * t64;
t24 = -pkin(9) * t119 + t39;
t3 = t16 * t192 - t189 * t24;
t116 = t190 * t124;
t69 = t272 * t118 - t116;
t63 = -t184 * pkin(4) + qJD(5) - t69;
t55 = t119 * pkin(5) + t63;
t304 = t15 * t158 - t3 * t207 + t310 * t55;
t4 = t16 * t189 + t192 * t24;
t303 = t15 * t160 + t4 * t207 - t311 * t55;
t214 = t119 * t189 - t121 * t192;
t262 = t214 * t207;
t302 = t262 + t312;
t215 = -qJD(6) * t121 - t265;
t236 = qJD(6) * t192;
t259 = -t119 * t236 + t95 * t243;
t22 = t215 * t189 + t259;
t301 = t22 * t160 + t311 * t214;
t23 = -t214 * qJD(6) + t160 * t95;
t300 = -t22 * t158 - t160 * t23 + t310 * t214 + t311 * t72;
t254 = t131 * t184;
t299 = t95 + t254;
t298 = t129 * t214;
t292 = t207 * t131;
t257 = t207 * t184;
t291 = -t96 + t257;
t289 = -t131 ^ 2 + t207 ^ 2;
t9 = -t185 * t27 + t187 * t37;
t288 = -t131 * t39 - t9;
t97 = pkin(4) * t207 + qJ(5) * t131;
t287 = t137 * t131 - t196;
t173 = pkin(3) * t229 + qJD(5);
t76 = t272 * t123 - t116;
t79 = pkin(3) * t150 + t97;
t45 = t185 * t79 + t187 * t76;
t283 = -t173 * t187 + t45;
t50 = t185 * t97 + t187 * t69;
t282 = -qJD(5) * t187 + t50;
t75 = t190 * t123 + t117;
t222 = pkin(3) * t239 - t75;
t281 = t129 * t207;
t248 = t167 * t193;
t126 = -pkin(8) * t161 - t169 * t191 - t248;
t212 = t167 * t191 - t169 * t193;
t127 = -pkin(8) * t211 - t212;
t279 = t272 * t126 - t190 * t127;
t278 = -qJD(6) + t129;
t275 = t29 * t185 + t207 * t39;
t274 = -t137 * t207 - t29;
t273 = -t29 * t187 - t207 * t38;
t271 = t187 * pkin(5);
t181 = t187 * pkin(9);
t154 = t161 * qJD(3);
t199 = -qJD(3) * t248 + qJD(2) * t244 + (-qJD(2) * t186 - qJD(3) * t169) * t191;
t106 = -t154 * pkin(8) + t199;
t153 = t211 * qJD(3);
t195 = t212 * qJD(3) - t203;
t107 = t153 * pkin(8) + t195;
t42 = t279 * qJD(4) + t272 * t106 + t190 * t107;
t206 = -t190 * t161 - t211 * t272;
t102 = t206 * qJD(4) - t272 * t153 - t190 * t154;
t136 = t272 * t161 - t190 * t211;
t103 = t136 * qJD(4) - t190 * t153 + t272 * t154;
t48 = pkin(3) * t154 + pkin(4) * t103 - qJ(5) * t102 - qJD(5) * t136;
t12 = t185 * t48 + t187 * t42;
t141 = pkin(3) * t211 + t179;
t86 = -pkin(4) * t206 - qJ(5) * t136 + t141;
t92 = t190 * t126 + t272 * t127;
t52 = t185 * t86 + t187 * t92;
t264 = t187 * t95;
t258 = t102 * t185;
t251 = t136 * t185;
t250 = t136 * t187;
t241 = t186 ^ 2 + t188 ^ 2;
t240 = qJD(3) * t150;
t237 = qJD(6) * t136;
t235 = qJD(1) * qJD(2);
t2 = pkin(5) * t96 - pkin(9) * t264 + t9;
t5 = -pkin(9) * t265 + t10;
t232 = -t189 * t5 + t192 * t2;
t230 = qJD(1) * t245;
t11 = -t185 * t42 + t187 * t48;
t49 = -t185 * t69 + t187 * t97;
t44 = -t185 * t76 + t187 * t79;
t51 = -t185 * t92 + t187 * t86;
t225 = t241 * qJD(1) ^ 2;
t180 = -t272 * pkin(3) - pkin(4);
t223 = t222 + t309;
t220 = t189 * t2 + t192 * t5;
t219 = -t9 * t185 + t7;
t218 = t185 * t38 - t187 * t39;
t30 = -pkin(5) * t206 - pkin(9) * t250 + t51;
t41 = -pkin(9) * t251 + t52;
t217 = -t189 * t41 + t192 * t30;
t216 = t189 * t30 + t192 * t41;
t210 = 0.2e1 * t241 * t235;
t177 = pkin(3) * t190 + qJ(5);
t155 = (-pkin(9) - t177) * t185;
t209 = -qJD(6) * t155 + t283 + t308;
t156 = t177 * t187 + t181;
t208 = qJD(6) * t156 + t173 * t185 + t305 + t44;
t166 = (-pkin(9) - qJ(5)) * t185;
t205 = -qJD(6) * t166 + t282 + t308;
t168 = qJ(5) * t187 + t181;
t204 = qJD(5) * t185 + qJD(6) * t168 + t305 + t49;
t202 = t102 * t63 + t136 * t29 - t279 * t95;
t200 = -pkin(4) * t95 - qJ(5) * t96 + (-qJD(5) + t63) * t131;
t198 = -t177 * t96 + t180 * t95 + (-t173 + t63) * t131;
t43 = t92 * qJD(4) + t190 * t106 - t272 * t107;
t178 = -pkin(4) - t271;
t165 = t180 - t271;
t99 = t158 * t136;
t98 = t160 * t136;
t58 = pkin(5) * t251 - t279;
t56 = t70 - t309;
t33 = t160 * t102 + t236 * t250 - t237 * t246;
t32 = -t158 * t102 - t160 * t237;
t20 = pkin(5) * t258 + t43;
t8 = -pkin(9) * t258 + t12;
t6 = pkin(5) * t103 - t102 * t181 + t11;
t1 = [0, 0, 0, 0, 0, t210, qJ(2) * t210, t146 * t161 - t150 * t153, -t146 * t211 + t153 * t149 - t150 * t154 - t161 * t242, -t153 * qJD(3), -t154 * qJD(3), 0, t195 * qJD(3) + t164 * t154 + t179 * t242, -t199 * qJD(3) + t179 * t146 - t164 * t153, t102 * t207 + t136 * t95, -t102 * t131 - t103 * t207 - t136 * t96 + t206 * t95, t102 * t184, -t103 * t184, 0, t137 * t103 + t141 * t96 - t43 * t184 + (t154 * t131 - t206 * t242) * pkin(3), t137 * t102 + t141 * t95 - t42 * t184 + (t136 * t242 + t154 * t207) * pkin(3), t103 * t38 + t11 * t131 + t119 * t43 + t185 * t202 - t206 * t9 + t51 * t96, t10 * t206 - t103 * t39 - t12 * t131 + t121 * t43 + t187 * t202 - t52 * t96, -t11 * t121 - t12 * t119 + (-t102 * t38 - t136 * t9 - t51 * t95) * t187 + (-t10 * t136 - t102 * t39 - t52 * t95) * t185, t10 * t52 + t11 * t38 + t12 * t39 - t279 * t29 + t43 * t63 + t51 * t9, -t214 * t32 - t22 * t99, t214 * t33 - t22 * t98 + t23 * t99 - t32 * t72, -t103 * t214 + t129 * t32 - t206 * t22 - t96 * t99, -t103 * t72 - t129 * t33 + t206 * t23 - t96 * t98, t103 * t129 - t206 * t96 (-t189 * t8 + t192 * t6) * t129 + t217 * t96 - t232 * t206 + t3 * t103 + t20 * t72 + t58 * t23 + t15 * t98 + t55 * t33 + (-t129 * t216 + t206 * t4) * qJD(6) -(t189 * t6 + t192 * t8) * t129 - t216 * t96 + t220 * t206 - t4 * t103 - t20 * t214 + t58 * t22 - t15 * t99 + t55 * t32 + (-t129 * t217 + t206 * t3) * qJD(6); 0, 0, 0, 0, 0, -t225, -qJ(2) * t225, 0, 0, 0, 0, 0, t240 + t242, t172 + (-t149 - t230) * qJD(3), 0, 0, 0, 0, 0, t96 + t257, t95 - t254, -t119 * t207 - t131 * t297 + t187 * t96, -t121 * t207 - t131 * t296 - t185 * t96 (-t185 ^ 2 - t187 ^ 2) * t95 - (t119 * t187 - t121 * t185) * t131, t10 * t185 - t131 * t218 + t9 * t187 - t207 * t63, 0, 0, 0, 0, 0, t221 - t266, t262 - t312; 0, 0, 0, 0, 0, 0, 0, t150 * t149, -t149 ^ 2 + t150 ^ 2, t172 + (t149 - t230) * qJD(3), t240 - t242, 0, -t164 * t150 - t201, t164 * t149 + t211 * t235, t292, t289, t299, t291, 0, t75 * t184 + (-t131 * t150 - t184 * t239) * pkin(3) + t274, t76 * t184 + (-t150 * t207 - t184 * t229) * pkin(3) + t287, t119 * t222 - t131 * t44 + t185 * t198 + t273, t121 * t222 + t131 * t45 + t187 * t198 + t275, t44 * t121 + t283 * t119 + (t121 * t173 + t288) * t185 + t306, -t173 * t218 + t177 * t219 + t180 * t29 + t222 * t63 - t38 * t44 - t39 * t45, t301, t300, t302, t307, -t281 (t155 * t192 - t156 * t189) * t96 + t165 * t23 + t223 * t72 + (t189 * t209 - t192 * t208) * t129 + t304 -(t155 * t189 + t156 * t192) * t96 + t165 * t22 - t223 * t214 + (t189 * t208 + t192 * t209) * t129 + t303; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t292, t289, t299, t291, 0, t184 * t70 + t274, t69 * t184 + t287, -t70 * t119 - t49 * t131 + t185 * t200 + t273, -t70 * t121 + t50 * t131 + t187 * t200 + t275, t49 * t121 + t282 * t119 + (qJD(5) * t121 + t288) * t185 + t306, -t29 * pkin(4) + qJ(5) * t219 - qJD(5) * t218 - t38 * t49 - t39 * t50 - t63 * t70, t301, t300, t302, t307, -t281 (t166 * t192 - t168 * t189) * t96 + t178 * t23 - t56 * t72 + (t189 * t205 - t192 * t204) * t129 + t304 -(t166 * t189 + t168 * t192) * t96 + t178 * t22 + t56 * t214 + (t189 * t204 + t192 * t205) * t129 + t303; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121 * t131 + t265, -t119 * t131 + t264, -t119 ^ 2 - t121 ^ 2, t119 * t39 + t121 * t38 + t29, 0, 0, 0, 0, 0, t23 - t298, -t129 * t293 + (-t121 * t129 + t215) * t189 + t259; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t214 * t72, t214 ^ 2 - t72 ^ 2, t72 * t129 + t22, -t23 - t298, t96, t214 * t55 + t278 * t4 + t232, t278 * t3 + t55 * t72 - t220;];
tauc_reg  = t1;
