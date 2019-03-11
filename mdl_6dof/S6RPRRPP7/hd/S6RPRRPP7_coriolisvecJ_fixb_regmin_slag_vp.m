% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
% 
% Output:
% tauc_reg [6x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRPP7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:52:22
% EndTime: 2019-03-09 04:52:28
% DurationCPUTime: 3.06s
% Computational Cost: add. (2925->366), mult. (5969->473), div. (0->0), fcn. (3331->4), ass. (0->184)
t113 = sin(qJ(3));
t196 = qJD(1) * t113;
t103 = qJD(4) + t196;
t112 = sin(qJ(4));
t115 = cos(qJ(3));
t195 = qJD(1) * t115;
t173 = t112 * t195;
t114 = cos(qJ(4));
t185 = t114 * qJD(3);
t78 = t173 - t185;
t218 = t103 * t78;
t189 = qJD(4) * t115;
t168 = t112 * t189;
t40 = (t113 * t185 + t168) * qJD(1) - qJD(4) * t185;
t251 = -t40 - t218;
t250 = t40 - t218;
t117 = -pkin(1) - pkin(7);
t213 = qJ(5) * t114;
t233 = pkin(4) + pkin(5);
t136 = t233 * t112 - t213;
t249 = t117 - t136;
t194 = qJD(3) * t112;
t80 = t114 * t195 + t194;
t212 = qJD(4) * t80;
t174 = t112 * t196;
t95 = qJD(3) * t174;
t124 = t95 - t212;
t217 = t103 * t80;
t248 = t124 - t217;
t183 = qJD(1) * qJD(3);
t165 = t115 * t183;
t247 = qJ(5) * t165 + t103 * qJD(5);
t246 = t233 * t165;
t159 = t124 + t212;
t160 = qJD(4) * t78 - t40;
t192 = qJD(3) * t115;
t245 = (qJD(1) * t78 + t160 * t113 + t80 * t192) * t112 + (qJD(1) * t80 + t159 * t113 - t78 * t192) * t114;
t179 = 0.2e1 * qJD(1);
t86 = pkin(3) * t113 - pkin(8) * t115 + qJ(2);
t61 = t86 * qJD(1);
t99 = t117 * qJD(1) + qJD(2);
t85 = t113 * t99;
t63 = qJD(3) * pkin(8) + t85;
t27 = t112 * t61 + t114 * t63;
t96 = t103 * qJ(5);
t20 = t96 + t27;
t155 = pkin(4) * t165;
t172 = t112 * t192;
t190 = qJD(4) * t114;
t191 = qJD(4) * t112;
t153 = pkin(3) * t115 + pkin(8) * t113;
t77 = qJD(3) * t153 + qJD(2);
t53 = t77 * qJD(1);
t161 = -t114 * t53 + t99 * t172 + t63 * t190 + t61 * t191;
t7 = -t155 + t161;
t244 = -t103 * t20 + t7;
t243 = -qJ(6) * t124 + t78 * qJD(6);
t216 = t115 * t99;
t64 = -qJD(3) * pkin(3) - t216;
t129 = qJ(5) * t80 - t64;
t13 = -t233 * t78 + qJD(6) + t129;
t242 = (qJD(6) + t13) * t80;
t241 = 0.2e1 * t247;
t76 = t80 ^ 2;
t240 = -t103 ^ 2 - t76;
t239 = t112 * qJD(5) + t85;
t26 = -t112 * t63 + t114 * t61;
t200 = qJD(5) - t26;
t238 = qJ(5) * t192 + t113 * qJD(5);
t209 = t112 * qJ(5);
t235 = -t233 * t114 - t209;
t234 = t78 ^ 2;
t193 = qJD(3) * t113;
t8 = -pkin(4) * t124 + t40 * qJ(5) - t80 * qJD(5) + t99 * t193;
t232 = t112 * t8;
t231 = t114 * t8;
t23 = pkin(4) * t78 - t129;
t230 = t23 * t80;
t229 = t80 * t78;
t228 = pkin(8) - qJ(6);
t207 = t112 * t115;
t82 = t153 * qJD(1);
t164 = t114 * t82 - t99 * t207;
t206 = t113 * t114;
t93 = t228 * t114;
t227 = (qJ(6) * t206 - t233 * t115) * qJD(1) - t164 - qJD(4) * t93 + t112 * qJD(6);
t184 = t114 * qJD(6);
t223 = t112 * t82 + t114 * t216;
t29 = qJ(5) * t195 + t223;
t226 = qJ(6) * t174 - t228 * t191 - t184 - t29;
t225 = t103 * t136 - t239;
t149 = pkin(4) * t112 - t213;
t224 = -t103 * t149 + t239;
t205 = t113 * t117;
t222 = t112 * t86 + t114 * t205;
t221 = qJ(5) * t124;
t16 = qJ(6) * t78 + t27;
t12 = t16 + t96;
t220 = t103 * t12;
t215 = t40 * t112;
t214 = t78 * qJ(5);
t211 = t103 * t112;
t210 = t103 * t114;
t208 = t112 * t113;
t118 = qJD(3) ^ 2;
t204 = t118 * t113;
t203 = t118 * t115;
t119 = qJD(1) ^ 2;
t202 = t119 * qJ(2);
t15 = qJ(6) * t80 + t26;
t201 = qJD(5) - t15;
t111 = t115 ^ 2;
t198 = t113 ^ 2 - t111;
t197 = -t118 - t119;
t188 = qJD(4) * t117;
t187 = qJD(5) * t114;
t170 = t117 * t192;
t182 = t112 * t77 + t114 * t170 + t86 * t190;
t167 = t113 * t188;
t181 = t112 * t170 + t114 * t167 + t86 * t191;
t169 = t103 * t190;
t180 = -qJD(1) * t210 - t103 * t172 - t113 * t169;
t178 = pkin(8) * t211;
t177 = pkin(8) * t210;
t35 = t113 * qJ(5) + t222;
t176 = pkin(8) * t192;
t175 = qJD(2) * t179;
t171 = t115 * t185;
t166 = t114 * t189;
t97 = t112 * t205;
t163 = t114 * t86 - t97;
t162 = -t64 + t216;
t158 = -t80 + t194;
t157 = t78 + t185;
t156 = qJD(4) * t113 + qJD(1);
t154 = t113 * t165;
t9 = -t233 * t103 + t201;
t152 = t112 * t9 + t114 * t12;
t151 = -t112 * t12 + t114 * t9;
t150 = pkin(4) * t114 + t209;
t19 = -pkin(4) * t103 + t200;
t148 = t112 * t20 - t114 * t19;
t147 = t112 * t19 + t114 * t20;
t146 = t114 * t77 - t181;
t145 = qJD(1) * t111 - t103 * t113;
t144 = -t113 * t23 + t176;
t143 = t113 * t64 - t176;
t3 = pkin(5) * t124 - t8;
t142 = -t112 * t3 - t13 * t190;
t141 = t114 * t3 - t13 * t191;
t140 = -t117 + t149;
t139 = qJ(6) * t40 + t161;
t138 = t103 * t27 - t161;
t137 = -t115 * t40 - t80 * t193;
t132 = -t165 + t229;
t131 = -t112 * t53 - t99 * t171 - t61 * t190 + t63 * t191;
t130 = -t114 * t154 - t137;
t128 = -t112 * t167 + t182;
t4 = -t131 + t247;
t126 = t103 * t26 + t131;
t125 = t78 * t193 + (t124 - t95) * t115;
t122 = -qJD(4) * t148 + t112 * t7 + t114 * t4;
t121 = t139 - t246;
t120 = t103 * t171 - t156 * t211 - t130;
t92 = t228 * t112;
t87 = -pkin(3) - t150;
t69 = pkin(3) - t235;
t45 = t140 * t115;
t36 = -pkin(4) * t113 - t163;
t34 = t249 * t115;
t32 = pkin(4) * t80 + t214;
t31 = qJ(6) * t207 + t35;
t30 = -pkin(4) * t195 - t164;
t25 = t97 + (-qJ(6) * t115 - t86) * t114 - t233 * t113;
t24 = -t233 * t80 - t214;
t18 = (qJD(4) * t150 - t187) * t115 - t140 * t193;
t14 = -pkin(4) * t192 - t146;
t11 = t128 + t238;
t10 = (t235 * qJD(4) + t187) * t115 - t249 * t193;
t6 = qJ(6) * t166 + (qJD(6) * t115 + (-qJ(6) * qJD(3) - t188) * t113) * t112 + t182 + t238;
t5 = (qJ(6) * t193 - t77) * t114 + (qJ(6) * t191 - t233 * qJD(3) - t184) * t115 + t181;
t2 = -qJD(6) * t80 + t121;
t1 = t4 + t243;
t17 = [0, 0, 0, 0, t175, qJ(2) * t175, -0.2e1 * t154, 0.2e1 * t198 * t183, -t204, -t203, 0, -t117 * t204 + (qJ(2) * t192 + qJD(2) * t113) * t179, -t117 * t203 + (-qJ(2) * t193 + qJD(2) * t115) * t179, t114 * t137 - t80 * t168 (t112 * t80 + t114 * t78) * t193 + (t215 + t114 * t124 + (t112 * t78 - t114 * t80) * qJD(4)) * t115, -t103 * t168 - t113 * t40 + (t114 * t145 + t115 * t80) * qJD(3), -t103 * t166 + t113 * t124 + (-t112 * t145 - t115 * t78) * qJD(3) (t103 + t196) * t192, t146 * t103 - t161 * t113 + (t117 * t124 + t64 * t190) * t115 + ((t163 * qJD(1) + t26) * t115 + (t162 * t112 + t117 * t78) * t113) * qJD(3), -t128 * t103 + t131 * t113 + (t117 * t40 - t64 * t191) * t115 + ((-t222 * qJD(1) - t27) * t115 + (t162 * t114 + t117 * t80) * t113) * qJD(3), -t103 * t14 + t18 * t78 - t124 * t45 + (-t23 * t194 - t7) * t113 + (t23 * t190 + t232 + (-qJD(1) * t36 - t19) * qJD(3)) * t115, -t11 * t78 + t14 * t80 + t35 * t124 - t36 * t40 + t148 * t193 + (-qJD(4) * t147 - t112 * t4 + t114 * t7) * t115, t103 * t11 - t18 * t80 + t40 * t45 + (t23 * t185 + t4) * t113 + (t23 * t191 - t231 + (qJD(1) * t35 + t20) * qJD(3)) * t115, t11 * t20 + t14 * t19 + t18 * t23 + t35 * t4 + t36 * t7 + t45 * t8, -t10 * t78 - t103 * t5 + t34 * t124 + (t13 * t194 - t2) * t113 + ((-qJD(1) * t25 - t9) * qJD(3) + t142) * t115, t10 * t80 + t103 * t6 - t34 * t40 + (-t13 * t185 + t1) * t113 + ((qJD(1) * t31 + t12) * qJD(3) + t141) * t115, t25 * t40 - t31 * t124 - t5 * t80 + t6 * t78 + t151 * t193 + (qJD(4) * t152 + t1 * t112 - t114 * t2) * t115, t1 * t31 + t10 * t13 + t12 * t6 + t2 * t25 + t3 * t34 + t5 * t9; 0, 0, 0, 0, -t119, -t202, 0, 0, 0, 0, 0, t197 * t113, t197 * t115, 0, 0, 0, 0, 0 (-t156 * t114 - t172) * t103 + t125 (t156 * t112 - t171) * t103 + t130, t125 + t180, t245, t120, -t148 * qJD(1) + (qJD(3) * t147 - t8) * t115 + (qJD(3) * t23 + t122) * t113, t115 * t124 + (t78 - t173) * t193 + t180, t120, -t245, t151 * qJD(1) + (qJD(3) * t152 + t3) * t115 + (-qJD(3) * t13 + qJD(4) * t151 + t1 * t114 + t112 * t2) * t113; 0, 0, 0, 0, 0, 0, t115 * t119 * t113, -t198 * t119, 0, 0, 0, -t115 * t202, t113 * t202, t80 * t210 - t215, t248 * t112 + t114 * t251, t169 + (t103 * t206 + t158 * t115) * qJD(1), -t103 * t191 + (-t103 * t208 + t157 * t115) * qJD(1), -t103 * t195, pkin(3) * t124 - t164 * t103 - t157 * t85 + (t112 * t64 - t177) * qJD(4) + (t112 * t143 - t26 * t115) * qJD(1), pkin(3) * t40 + t223 * t103 + t158 * t85 + (t114 * t64 + t178) * qJD(4) + (t114 * t143 + t27 * t115) * qJD(1), t103 * t30 - t231 - t124 * t87 - t224 * t78 + (t112 * t23 - t177) * qJD(4) + (-t112 * t144 + t115 * t19) * qJD(1), t29 * t78 - t30 * t80 + (t159 * pkin(8) + t103 * t19 + t4) * t114 + (t160 * pkin(8) + t244) * t112, -t103 * t29 - t232 + t40 * t87 + t224 * t80 + (-t114 * t23 - t178) * qJD(4) + (t114 * t144 - t115 * t20) * qJD(1), t122 * pkin(8) - t19 * t30 - t20 * t29 - t224 * t23 + t8 * t87, t124 * t69 + t225 * t78 + t227 * t103 + (-t13 * t208 + (-qJD(3) * t92 + t9) * t115) * qJD(1) + t141, -t40 * t69 - t225 * t80 + t226 * t103 + (t13 * t206 + (qJD(3) * t93 - t12) * t115) * qJD(1) - t142, t40 * t92 - t124 * t93 + t227 * t80 + t226 * t78 + (-t103 * t9 - t1) * t114 + (-t2 + t220) * t112, t1 * t93 + t226 * t12 - t225 * t13 + t2 * t92 - t227 * t9 + t3 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t229, t76 - t234, -t250, t124 + t217, t165, -t64 * t80 + t138, t64 * t78 + t126, -t32 * t78 + t138 + 0.2e1 * t155 - t230, pkin(4) * t40 + t221 + (t20 - t27) * t80 + (t19 - t200) * t78, -t23 * t78 + t32 * t80 - t126 + t241, -pkin(4) * t7 + qJ(5) * t4 - t19 * t27 + t20 * t200 - t23 * t32, t103 * t16 + t24 * t78 - t139 + t242 + 0.2e1 * t246, -t103 * t15 + t13 * t78 - t24 * t80 - t131 + t241 + t243, -t221 - t233 * t40 + (-t12 + t16) * t80 + (-t9 + t201) * t78, qJ(5) * t1 + t12 * t201 - t13 * t24 - t16 * t9 - t2 * t233; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132, -t250, t240, t230 + t244, t132, t240, t250, t121 - t220 - t242; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t248, t251, -t76 - t234, -t12 * t78 + t80 * t9 + t3;];
tauc_reg  = t17;
