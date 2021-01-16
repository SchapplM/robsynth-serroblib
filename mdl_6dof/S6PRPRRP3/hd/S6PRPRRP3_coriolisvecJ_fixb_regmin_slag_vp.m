% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tauc_reg [6x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:40
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPRRP3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:38:31
% EndTime: 2021-01-16 01:38:41
% DurationCPUTime: 2.68s
% Computational Cost: add. (3744->308), mult. (9684->423), div. (0->0), fcn. (7691->10), ass. (0->164)
t231 = 2 * qJD(4);
t123 = sin(pkin(11));
t125 = cos(pkin(11));
t128 = sin(qJ(4));
t210 = cos(qJ(4));
t140 = -t128 * t123 + t210 * t125;
t201 = qJ(3) + pkin(8);
t108 = t201 * t123;
t109 = t201 * t125;
t217 = -t210 * t108 - t128 * t109;
t47 = t140 * qJD(3) + qJD(4) * t217;
t124 = sin(pkin(6));
t131 = cos(qJ(2));
t181 = t124 * t131;
t137 = t140 * t181;
t77 = qJD(1) * t137;
t230 = -t47 + t77;
t105 = t210 * t123 + t128 * t125;
t100 = t105 * qJD(4);
t129 = sin(qJ(2));
t174 = qJD(1) * t124;
t158 = t129 * t174;
t99 = t140 * qJD(4);
t229 = -t100 * pkin(4) + t99 * pkin(9) + t158;
t107 = qJD(2) * qJ(3) + t158;
t126 = cos(pkin(6));
t173 = qJD(1) * t126;
t114 = t125 * t173;
t192 = pkin(8) * qJD(2);
t71 = t114 + (-t107 - t192) * t123;
t84 = t125 * t107 + t123 * t173;
t72 = t125 * t192 + t84;
t34 = t128 * t71 + t210 * t72;
t228 = t34 * qJD(4);
t227 = t140 * qJD(2);
t127 = sin(qJ(5));
t130 = cos(qJ(5));
t98 = qJD(2) * t105;
t81 = t127 * qJD(4) + t130 * t98;
t190 = t127 * t81;
t136 = qJD(2) * t99;
t225 = (qJD(5) * qJD(4)) + t136;
t30 = qJD(4) * pkin(9) + t34;
t117 = -t125 * pkin(3) - pkin(2);
t157 = t131 * t174;
t149 = qJD(3) - t157;
t90 = t117 * qJD(2) + t149;
t37 = -pkin(4) * t227 - t98 * pkin(9) + t90;
t15 = t127 * t37 + t130 * t30;
t79 = -t130 * qJD(4) + t127 * t98;
t10 = -t79 * qJ(6) + t15;
t91 = qJD(5) - t227;
t209 = t10 * t91;
t103 = (qJD(3) + t157) * qJD(2);
t218 = t140 * t103;
t219 = -t128 * t72 + t210 * t71;
t20 = qJD(4) * t219 + t218;
t112 = qJD(2) * t158;
t89 = qJD(2) * t100;
t43 = t89 * pkin(4) - pkin(9) * t136 + t112;
t39 = t130 * t43;
t134 = -t15 * qJD(5) - t127 * t20 + t39;
t171 = qJD(5) * t127;
t44 = -t225 * t130 + t98 * t171;
t133 = t44 * qJ(6) + t134;
t213 = t89 * pkin(5);
t3 = -t81 * qJD(6) + t133 + t213;
t224 = t3 + t209;
t138 = t105 * t181;
t75 = -t128 * t108 + t210 * t109;
t194 = -qJD(1) * t138 + t105 * qJD(3) + t75 * qJD(4);
t223 = t127 * t77 - t229 * t130;
t222 = t91 * t190;
t221 = t130 * t89 - t91 * t171;
t170 = qJD(5) * t130;
t64 = -pkin(4) * t140 - t105 * pkin(9) + t117;
t220 = t229 * t127 + t230 * t130 - t64 * t170;
t216 = t81 ^ 2;
t14 = -t127 * t30 + t130 * t37;
t9 = -t81 * qJ(6) + t14;
t7 = t91 * pkin(5) + t9;
t215 = t7 - t9;
t214 = t79 * pkin(5);
t147 = -qJ(6) * t99 - qJD(6) * t105;
t70 = t130 * t75;
t212 = t100 * pkin(5) - t127 * t47 + t147 * t130 + (-t70 + (qJ(6) * t105 - t64) * t127) * qJD(5) + t223;
t159 = t105 * t170;
t211 = qJ(6) * t159 - (-qJD(5) * t75 + t147) * t127 + t220;
t208 = t79 * t91;
t207 = t79 * t227;
t21 = t105 * t103 + t228;
t45 = t225 * t127 + t98 * t170;
t8 = t45 * pkin(5) + t21;
t206 = t8 * t127;
t205 = t8 * t130;
t204 = t81 * t91;
t203 = t81 * t98;
t202 = t98 * t79;
t200 = -qJ(6) - pkin(9);
t155 = qJD(5) * t200;
t176 = t130 * qJ(6);
t62 = t98 * pkin(4) - pkin(9) * t227;
t57 = t130 * t62;
t199 = t98 * pkin(5) - t130 * t155 - t227 * t176 + t57 + (qJD(6) - t219) * t127;
t179 = t127 * qJ(6);
t196 = t127 * t62 + t130 * t219;
t198 = -t130 * qJD(6) - t127 * t155 - t179 * t227 + t196;
t186 = t127 * t99;
t197 = (t159 + t186) * pkin(5) + t194;
t195 = -t127 * t45 - t79 * t170;
t193 = t127 * t64 + t70;
t191 = qJD(2) * pkin(2);
t189 = t127 * t89;
t187 = t127 * t227;
t156 = qJD(6) + t214;
t29 = -qJD(4) * pkin(4) - t219;
t23 = t156 + t29;
t185 = t130 * t23;
t184 = t130 * t99;
t183 = t44 * t127;
t182 = t124 * t129;
t132 = qJD(2) ^ 2;
t180 = t124 * t132;
t175 = t123 ^ 2 + t125 ^ 2;
t172 = qJD(2) * t129;
t165 = t129 * t180;
t164 = t131 * t180;
t160 = t124 * t172;
t154 = -t127 * t43 - t130 * t20 - t37 * t170 + t30 * t171;
t153 = t130 * t91;
t152 = t175 * t103;
t151 = qJD(5) * pkin(9) * t91 + t21;
t145 = t45 * qJ(6) + t154;
t4 = -t79 * qJD(6) - t145;
t150 = -t91 * t7 + t4;
t148 = t123 * (-t123 * t107 + t114) - t125 * t84;
t146 = t91 * t187 + t221;
t94 = -t123 * t182 + t125 * t126;
t95 = t126 * t123 + t125 * t182;
t143 = -t128 * t95 + t210 * t94;
t53 = t128 * t94 + t210 * t95;
t40 = -t127 * t53 - t130 * t181;
t142 = t127 * t181 - t130 * t53;
t141 = -pkin(9) * t89 + t91 * t29;
t135 = qJD(2) * t137;
t119 = -t130 * pkin(5) - pkin(4);
t111 = t200 * t130;
t110 = t200 * t127;
t106 = t149 - t191;
t78 = t79 ^ 2;
t60 = t130 * t64;
t49 = t127 * t105 * pkin(5) - t217;
t32 = qJD(2) * t138 + t53 * qJD(4);
t31 = t143 * qJD(4) + t135;
t25 = pkin(5) * t187 + t34;
t24 = -t105 * t179 + t193;
t22 = -pkin(5) * t140 - t105 * t176 - t127 * t75 + t60;
t18 = -t91 ^ 2 * t130 - t189 - t203;
t17 = t146 - t202;
t13 = t142 * qJD(5) - t127 * t31 + t130 * t160;
t12 = t40 * qJD(5) + t127 * t160 + t130 * t31;
t2 = t13 * t91 - t143 * t45 + t32 * t79 + t40 * t89;
t1 = -t12 * t91 + t142 * t89 + t143 * t44 + t32 * t81;
t5 = [0, 0, -t165, -t164, -t125 * t165, t175 * t164, (-t123 * t94 + t125 * t95) * t103 + (t106 * t129 + (-t148 - t158) * t131) * t124 * qJD(2), 0, 0, 0, 0, 0, -t32 * qJD(4) + (-t131 * t89 - t172 * t227) * t124, t98 * t160 + (-t135 - t31) * qJD(4), 0, 0, 0, 0, 0, t2, t1, t2, t1, -t12 * t79 - t13 * t81 + t142 * t45 + t40 * t44, t10 * t12 + t7 * t13 - t142 * t4 - t143 * t8 + t23 * t32 + t3 * t40; 0, 0, 0, 0, 0, qJD(2) * t149 * t175 + t152, -t148 * qJD(3) + qJ(3) * t152 + (t148 * t131 + (-t106 - t191) * t129) * t174, t105 * t136 + t98 * t99, -t98 * t100 - t105 * t89 + t136 * t140 + t227 * t99, t99 * qJD(4), -t100 * qJD(4), 0, -t194 * qJD(4) + t90 * t100 + t117 * t89, t230 * qJD(4) + t105 * t112 + t117 * t136 - t98 * t158 + t90 * t99, t81 * t184 + (-t130 * t44 - t81 * t171) * t105, (-t130 * t79 - t190) * t99 + (t183 - t130 * t45 + (t127 * t79 - t130 * t81) * qJD(5)) * t105, t81 * t100 + t105 * t221 + t140 * t44 + t91 * t184, -t91 * t186 - t79 * t100 + t45 * t140 + (-t91 * t170 - t189) * t105, t91 * t100 - t140 * t89, t60 * t89 - (-t30 * t170 + t39) * t140 + t14 * t100 - t217 * t45 + t29 * t159 + (-t75 * t170 + t223) * t91 + t194 * t79 + ((-qJD(5) * t64 - t47) * t91 - t75 * t89 - (-qJD(5) * t37 - t20) * t140 + t21 * t105 + t29 * t99) * t127, -t193 * t89 - t154 * t140 - t15 * t100 + t217 * t44 + t29 * t184 + (t75 * t171 + t220) * t91 + t194 * t81 + (t21 * t130 - t29 * t171) * t105, t23 * t186 + t7 * t100 - t3 * t140 + t22 * t89 + t49 * t45 + t212 * t91 + t197 * t79 + (t23 * t170 + t206) * t105, t23 * t184 - t10 * t100 + t4 * t140 - t24 * t89 - t49 * t44 + t211 * t91 + t197 * t81 + (-t23 * t171 + t205) * t105, t22 * t44 - t24 * t45 + (-t10 * t127 - t130 * t7) * t99 - t212 * t81 + t211 * t79 + (-t4 * t127 - t3 * t130 + (-t10 * t130 + t127 * t7) * qJD(5)) * t105, -t211 * t10 + t197 * t23 + t212 * t7 + t3 * t22 + t4 * t24 + t8 * t49; 0, 0, 0, 0, 0, -t175 * t132, t148 * qJD(2) + t112, 0, 0, 0, 0, 0, t98 * t231, t227 * t231, 0, 0, 0, 0, 0, t17, t18, t17, t18, (t44 + t207) * t130 + t222 + t195, t150 * t127 + t130 * t224 - t23 * t98; 0, 0, 0, 0, 0, 0, 0, -t98 * t227, -t227 ^ 2 + t98 ^ 2, 0, 0, 0, -t90 * t98 - t21 + t228, -t227 * t90 - t218, t81 * t153 - t183, (-t44 + t207) * t130 - t222 + t195, t91 * t153 + t189 - t203, t146 + t202, -t91 * t98, -pkin(4) * t45 - t14 * t98 - t34 * t79 - t57 * t91 - t151 * t130 + (t219 * t91 + t141) * t127, pkin(4) * t44 + t151 * t127 + t141 * t130 + t15 * t98 + t196 * t91 - t34 * t81, t110 * t89 + t119 * t45 - t205 - t25 * t79 - t7 * t98 - t199 * t91 + (-t23 * t227 + (t23 + t214) * qJD(5)) * t127, -t227 * t185 + t10 * t98 + t111 * t89 - t119 * t44 + t206 - t25 * t81 + t198 * t91 + (pkin(5) * t190 + t185) * qJD(5), t110 * t44 + t111 * t45 - t127 * t224 + t150 * t130 + t198 * t79 + t199 * t81, t3 * t110 - t4 * t111 + t8 * t119 - t199 * t7 + (pkin(5) * t171 - t25) * t23 - t198 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81 * t79, -t78 + t216, -t44 + t208, -t45 + t204, t89, t15 * t91 - t29 * t81 + t134, t14 * t91 + t29 * t79 + t154, 0.2e1 * t213 + t209 + (-t156 - t23) * t81 + t133, -t216 * pkin(5) + t9 * t91 + (qJD(6) + t23) * t79 + t145, t44 * pkin(5) - t215 * t79, t215 * t10 + (-t23 * t81 + t3) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45 + t204, -t44 - t208, -t78 - t216, t10 * t79 + t7 * t81 + t8;];
tauc_reg = t5;
