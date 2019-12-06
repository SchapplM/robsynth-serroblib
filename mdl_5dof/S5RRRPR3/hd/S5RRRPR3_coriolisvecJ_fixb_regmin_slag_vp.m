% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRPR3
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
% tauc_reg [5x22]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:42:58
% EndTime: 2019-12-05 18:43:02
% DurationCPUTime: 1.20s
% Computational Cost: add. (2000->202), mult. (3522->279), div. (0->0), fcn. (2480->8), ass. (0->150)
t142 = sin(qJ(5));
t145 = cos(qJ(5));
t177 = qJD(5) * t142;
t140 = sin(pkin(9));
t141 = cos(pkin(9));
t143 = sin(qJ(3));
t146 = cos(qJ(3));
t113 = -t140 * t143 + t141 * t146;
t137 = qJD(1) + qJD(2);
t94 = t113 * t137;
t81 = t145 * t94;
t114 = t140 * t146 + t141 * t143;
t106 = t114 * qJD(3);
t87 = t137 * t106;
t179 = qJD(3) * t146;
t167 = t137 * t179;
t180 = qJD(3) * t143;
t168 = t137 * t180;
t88 = -t140 * t168 + t141 * t167;
t95 = t114 * t137;
t10 = qJD(5) * t81 - t142 * t87 + t145 * t88 - t95 * t177;
t136 = qJD(3) + qJD(5);
t42 = -t142 * t95 + t81;
t187 = t42 * t136;
t211 = t10 - t187;
t210 = qJ(4) + pkin(7);
t155 = t142 * t94 + t145 * t95;
t209 = t155 * t42;
t11 = t155 * qJD(5) + t142 * t88 + t145 * t87;
t188 = t155 * t136;
t208 = -t11 + t188;
t207 = t155 ^ 2 - t42 ^ 2;
t202 = t94 * pkin(8);
t144 = sin(qJ(2));
t191 = pkin(1) * qJD(1);
t173 = t144 * t191;
t162 = t210 * t137 + t173;
t84 = t162 * t146;
t189 = t141 * t84;
t83 = t162 * t143;
t74 = qJD(3) * pkin(3) - t83;
t34 = t140 * t74 + t189;
t18 = t34 + t202;
t170 = -t146 * pkin(3) - pkin(2);
t147 = cos(qJ(2));
t172 = t147 * t191;
t93 = t170 * t137 + qJD(4) - t172;
t49 = -t94 * pkin(4) + t93;
t206 = t18 * t177 - t49 * t42;
t132 = t146 * qJD(4);
t163 = qJD(3) * t210;
t103 = -t143 * t163 + t132;
t104 = -t143 * qJD(4) - t146 * t163;
t193 = -t140 * t103 + t141 * t104 + t114 * t172;
t192 = t141 * t103 + t140 * t104 - t113 * t172;
t204 = qJD(5) - t136;
t190 = pkin(1) * qJD(2);
t171 = qJD(1) * t190;
t160 = t147 * t171;
t151 = qJD(4) * t137 + t160;
t153 = qJD(3) * t162;
t47 = -t143 * t153 + t151 * t146;
t48 = -t151 * t143 - t146 * t153;
t12 = -t140 * t47 + t141 * t48;
t4 = -t88 * pkin(8) + t12;
t13 = t140 * t48 + t141 * t47;
t5 = -t87 * pkin(8) + t13;
t203 = -t142 * t5 + t145 * t4 - t49 * t155;
t201 = t95 * pkin(8);
t154 = t145 * t113 - t142 * t114;
t107 = t113 * qJD(3);
t58 = t142 * t113 + t145 * t114;
t26 = t58 * qJD(5) + t145 * t106 + t142 * t107;
t126 = t144 * t171;
t105 = pkin(3) * t168 + t126;
t54 = t87 * pkin(4) + t105;
t200 = -t154 * t54 + t49 * t26;
t25 = t154 * qJD(5) - t142 * t106 + t145 * t107;
t199 = t49 * t25 + t54 * t58;
t198 = pkin(3) * t140;
t197 = t107 * pkin(8);
t196 = t114 * pkin(8);
t195 = t147 * pkin(1);
t128 = t144 * pkin(1) + pkin(7);
t183 = -qJ(4) - t128;
t161 = qJD(3) * t183;
t174 = t147 * t190;
t63 = t143 * t161 + t146 * t174 + t132;
t64 = (-qJD(4) - t174) * t143 + t146 * t161;
t30 = t140 * t64 + t141 * t63;
t68 = t140 * t84;
t36 = -t141 * t83 - t68;
t111 = t183 * t143;
t134 = t146 * qJ(4);
t112 = t146 * t128 + t134;
t56 = t140 * t111 + t141 * t112;
t186 = t137 * t143;
t185 = t144 * t146;
t148 = qJD(3) ^ 2;
t184 = t148 * t143;
t133 = t148 * t146;
t120 = -t137 * pkin(2) - t172;
t182 = t120 * t179 + t143 * t126;
t123 = t210 * t143;
t124 = t146 * pkin(7) + t134;
t73 = -t140 * t123 + t141 * t124;
t181 = t143 ^ 2 - t146 ^ 2;
t178 = qJD(3) * t147;
t176 = -qJD(1) - t137;
t175 = -qJD(2) + t137;
t131 = t144 * t190;
t130 = pkin(3) * t180;
t166 = t143 * t178;
t33 = t141 * t74 - t68;
t165 = -t34 * t106 - t33 * t107 + t13 * t113 - t12 * t114;
t80 = t106 * pkin(4) + t130;
t29 = -t140 * t63 + t141 * t64;
t35 = t140 * t83 - t189;
t55 = t141 * t111 - t140 * t112;
t72 = -t141 * t123 - t140 * t124;
t159 = t80 - t173;
t110 = t113 * pkin(8);
t158 = qJD(5) * (t110 + t73) + t197 - t193;
t102 = t106 * pkin(8);
t157 = -qJD(5) * (t72 - t196) + t102 - t192;
t17 = qJD(3) * pkin(4) - t201 + t33;
t156 = -t142 * t17 - t145 * t18;
t89 = -t113 * pkin(4) + t170;
t152 = -t120 * t137 - t160;
t150 = -t144 * t186 + t146 * t178;
t135 = t137 ^ 2;
t129 = -pkin(2) - t195;
t127 = t141 * pkin(3) + pkin(4);
t118 = 0.2e1 * t143 * t167;
t108 = t120 * t180;
t97 = -0.2e1 * t181 * t137 * qJD(3);
t79 = t89 - t195;
t67 = t131 + t80;
t62 = pkin(3) * t186 + t95 * pkin(4);
t38 = t110 + t56;
t37 = t55 - t196;
t22 = t36 - t201;
t21 = t35 - t202;
t20 = t26 * t136;
t19 = t25 * t136;
t15 = -t102 + t30;
t14 = t29 - t197;
t2 = t10 * t58 + t155 * t25;
t1 = t10 * t154 - t58 * t11 - t155 * t26 + t25 * t42;
t3 = [0, 0, 0, 0, -t137 * t131 - t126, t176 * t174, t118, t97, t133, -t184, 0, t129 * t168 - t128 * t133 + t108 + (t176 * t185 - t166) * t190, t128 * t184 + t129 * t167 - t150 * t190 + t182, -t29 * t95 + t30 * t94 - t55 * t88 - t56 * t87 + t165, t13 * t56 + t34 * t30 + t12 * t55 + t33 * t29 + t105 * (t170 - t195) + t93 * (t131 + t130), t2, t1, t19, -t20, 0, -t67 * t42 + t79 * t11 + (t145 * t14 - t142 * t15 + (-t142 * t37 - t145 * t38) * qJD(5)) * t136 + t200, t67 * t155 + t79 * t10 - (t142 * t14 + t145 * t15 + (-t142 * t38 + t145 * t37) * qJD(5)) * t136 + t199; 0, 0, 0, 0, t137 * t173 - t126, t175 * t172, t118, t97, t133, -t184, 0, -pkin(2) * t168 - pkin(7) * t133 + t108 + (t175 * t185 + t166) * t191, -pkin(2) * t167 + pkin(7) * t184 + t150 * t191 + t182, t192 * t94 - t193 * t95 - t72 * t88 - t73 * t87 + t165, t13 * t73 + t12 * t72 + t105 * t170 + (-t173 + t130) * t93 + t192 * t34 + t193 * t33, t2, t1, t19, -t20, 0, t89 * t11 - t159 * t42 + (t157 * t142 - t158 * t145) * t136 + t200, t89 * t10 + t159 * t155 + (t158 * t142 + t157 * t145) * t136 + t199; 0, 0, 0, 0, 0, 0, -t143 * t135 * t146, t181 * t135, 0, 0, 0, t152 * t143, t152 * t146, (t34 + t35) * t95 + (t33 - t36) * t94 + (-t140 * t87 - t141 * t88) * pkin(3), -t33 * t35 - t34 * t36 + (t12 * t141 + t13 * t140 - t93 * t186) * pkin(3), -t209, t207, t211, t208, 0, t62 * t42 - (-t142 * t22 + t145 * t21) * t136 + ((-t127 * t142 - t145 * t198) * t136 + t156) * qJD(5) + t203, -t145 * t5 - t142 * t4 - t62 * t155 + (t142 * t21 + t145 * t22) * t136 + (-(t127 * t145 - t142 * t198) * t136 - t145 * t17) * qJD(5) + t206; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94 ^ 2 - t95 ^ 2, t33 * t95 - t34 * t94 + t105, 0, 0, 0, 0, 0, t11 + t188, t10 + t187; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t209, t207, t211, t208, 0, t204 * t156 + t203, (-t18 * t136 - t4) * t142 + (-t204 * t17 - t5) * t145 + t206;];
tauc_reg = t3;
