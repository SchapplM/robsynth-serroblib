% Calculate minimal parameter regressor of coriolis matrix for
% S5RPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x25]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRPR13_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR13_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:32:55
% EndTime: 2019-12-31 18:32:59
% DurationCPUTime: 1.32s
% Computational Cost: add. (1302->169), mult. (2738->231), div. (0->0), fcn. (3043->6), ass. (0->140)
t179 = pkin(3) + pkin(7);
t96 = sin(pkin(8));
t99 = sin(qJ(3));
t173 = t99 * t96;
t177 = cos(qJ(3));
t97 = cos(pkin(8));
t90 = t177 * t97;
t79 = -t90 + t173;
t183 = t179 * t79;
t81 = t177 * t96 + t97 * t99;
t182 = t179 * t81;
t172 = pkin(6) + qJ(2);
t84 = t172 * t96;
t85 = t172 * t97;
t48 = t177 * t84 + t85 * t99;
t111 = pkin(4) * t81 + t48;
t98 = sin(qJ(5));
t176 = t98 * t111;
t37 = t98 * t81;
t152 = t37 * qJD(1);
t162 = qJD(5) * t98;
t181 = t152 + t162;
t180 = t81 ^ 2;
t75 = t79 ^ 2;
t130 = t180 - t75;
t45 = t75 + t180;
t145 = t180 * qJD(1);
t100 = cos(qJ(5));
t41 = t100 * t79;
t121 = 0.2e1 * t98 * t41;
t94 = t98 ^ 2;
t95 = t100 ^ 2;
t89 = t94 - t95;
t103 = qJD(1) * t121 - qJD(3) * t89;
t178 = t79 * pkin(3);
t49 = t177 * t85 - t84 * t99;
t29 = -pkin(4) * t79 + t49;
t175 = t29 * t98;
t171 = qJ(4) * t79;
t27 = t171 + t182;
t174 = t98 * t27;
t86 = t96 ^ 2 + t97 ^ 2;
t26 = t29 * t100;
t164 = t81 * qJ(4);
t91 = -t97 * pkin(2) - pkin(1);
t112 = t91 - t164;
t19 = t112 + t183;
t9 = -t100 * t111 + t19 * t98;
t1 = (t26 - t174) * t81 + t9 * t79 + (t111 * t79 - t29 * t81) * t100;
t170 = t1 * qJD(1);
t169 = t100 * t27;
t10 = t100 * t19 + t176;
t2 = -t81 * t169 + (t10 - t176) * t79;
t168 = t2 * qJD(1);
t115 = pkin(3) * t81 + t171;
t36 = t112 + t178;
t3 = t36 * t115;
t167 = t3 * qJD(1);
t4 = t26 * t79 + t81 * t9;
t166 = t4 * qJD(1);
t5 = -t10 * t81 + t175 * t79;
t165 = t5 * qJD(1);
t163 = qJD(4) * t98;
t11 = -t115 * t79 - t36 * t81;
t161 = t11 * qJD(1);
t12 = -t115 * t81 + t36 * t79;
t160 = t12 * qJD(1);
t15 = t48 * t81 - t49 * t79;
t159 = t15 * qJD(1);
t158 = t115 * qJD(1);
t18 = t45 * t98;
t157 = t18 * qJD(1);
t22 = t130 * t98;
t156 = t22 * qJD(1);
t23 = t130 * t100;
t155 = t23 * qJD(1);
t24 = t45 * t100;
t154 = t24 * qJD(1);
t153 = t130 * qJD(1);
t151 = t37 * qJD(5);
t150 = t41 * qJD(1);
t149 = t41 * qJD(3);
t148 = t45 * qJD(1);
t147 = t48 * qJD(3);
t44 = t49 * qJD(3);
t74 = t173 / 0.2e1 - t90 / 0.2e1;
t146 = t74 * qJD(1);
t144 = t79 * qJD(1);
t143 = t79 * qJD(2);
t68 = t79 * qJD(3);
t142 = t79 * qJD(4);
t141 = t81 * qJD(1);
t70 = t81 * qJD(2);
t72 = t81 * qJD(3);
t140 = t81 * qJD(4);
t83 = t86 * qJ(2);
t139 = t83 * qJD(1);
t138 = t86 * qJD(1);
t137 = t98 * qJD(3);
t136 = qJD(1) * t100;
t135 = qJD(3) * qJ(4);
t134 = qJD(4) * t100;
t133 = qJD(5) * t100;
t132 = qJD(5) * t179;
t131 = t100 * qJD(3);
t129 = t81 * t162;
t128 = t36 * t141;
t47 = t79 * t141;
t46 = t79 * t72;
t127 = t98 * t144;
t126 = t79 * t137;
t125 = t98 * t145;
t124 = t81 * t133;
t123 = t98 * t133;
t61 = t81 * t136;
t122 = t98 * t131;
t120 = qJD(5) * t74 + t47;
t119 = qJD(1) * t91 + qJD(2);
t118 = qJD(5) + t141;
t116 = qJD(3) * t121;
t114 = t164 - t183;
t113 = t118 * t98;
t110 = t79 * t113;
t109 = t182 / 0.2e1 + t171 / 0.2e1;
t102 = t27 / 0.2e1 + t109;
t14 = t102 * t100;
t108 = -qJD(1) * t14 + t135 * t98;
t39 = (t95 / 0.2e1 - t94 / 0.2e1) * t79;
t107 = -qJD(1) * t39 + t122;
t13 = t102 * t98;
t106 = -qJ(4) * t131 - qJD(1) * t13;
t105 = t136 * t75 * t98 + qJD(3) * t39;
t42 = t89 * t75;
t104 = qJD(1) * t42 + t116;
t62 = t74 * qJD(3);
t50 = -t61 - t133;
t34 = t39 * qJD(5);
t8 = -t175 - t169 / 0.2e1 + t109 * t100;
t7 = t26 - t174 / 0.2e1 + t109 * t98;
t6 = [0, 0, 0, 0, 0, t86 * qJD(2), t83 * qJD(2), -t46, -t130 * qJD(3), 0, 0, 0, t91 * t72, -t91 * t68, t45 * qJD(2), qJD(3) * t11 + t140 * t79, qJD(3) * t12 + qJD(4) * t180, qJD(2) * t15 + qJD(3) * t3 - t140 * t36, t123 * t75 + t46 * t94, -qJD(5) * t42 + t116 * t81, qJD(3) * t22 + t124 * t79, qJD(3) * t23 - t129 * t79, -t46, qJD(2) * t24 + qJD(3) * t1 + qJD(5) * t5 + t163 * t180, -qJD(2) * t18 + qJD(3) * t2 + qJD(5) * t4 + t134 * t180; 0, 0, 0, 0, 0, t138, t139, 0, 0, 0, 0, 0, 0, 0, t148, 0, 0, t159, 0, 0, 0, 0, 0, t154, -t157; 0, 0, 0, 0, 0, 0, 0, -t47, -t153, -t68, -t72, 0, t141 * t91 - t44, -t144 * t91 + t147, (-t164 + t178) * qJD(3) - t142, t44 + t161, -t147 + t160, t167 + (-pkin(3) * t49 - qJ(4) * t48) * qJD(3) + t49 * qJD(4), t34 + (t144 * t94 + t122) * t81, t103 * t81 - 0.2e1 * t79 * t123, -t131 * t79 + t156, t126 + t155, -t120, t170 + (-t100 * t114 - t176) * qJD(3) - t41 * qJD(4) + t7 * qJD(5), -t111 * t131 + t168 + t8 * qJD(5) + (qJD(3) * t114 + t142) * t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, t47, t145, t44 - t128, 0, 0, 0, 0, 0, t125 - t149, t136 * t180 + t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, -t104, t118 * t41, -t110, -t62, qJD(3) * t7 - qJD(5) * t10 + t165, qJD(3) * t8 + qJD(5) * t9 + t166; 0, 0, 0, 0, 0, -t138, -t139, 0, 0, 0, 0, 0, t72, -t68, -t148, -t72, t68, qJD(3) * t115 - t140 - t159, 0, 0, 0, 0, 0, -t124 + t126 - t154, t149 + t151 + t157; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t141, -t144, 0, -t141, t144, t158, 0, 0, 0, 0, 0, t127, t150; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t141, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, t181; 0, 0, 0, 0, 0, 0, 0, t47, t153, 0, 0, 0, -t119 * t81, t119 * t79, 0, t70 - t161, -t143 - t160, -qJD(2) * t115 - t167, -t47 * t94 + t34, -0.2e1 * t100 * t110, -t129 - t156, -t124 - t155, t120, qJD(5) * t13 - t143 * t98 - t170, -qJD(2) * t41 + qJD(5) * t14 - t168; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t141, t144, 0, t141, -t144, -t158, 0, 0, 0, 0, 0, -t127, -t150; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), qJ(4) * qJD(4), -t123, t89 * qJD(5), 0, 0, 0, qJ(4) * t133 + t163, -qJ(4) * t162 + t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t135, 0, 0, 0, 0, 0, t137, t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t107, -t103, -t113, t50, t146, t132 * t98 - t106, t100 * t132 - t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, -t145, t70 + t128, 0, 0, 0, 0, 0, -t125 - t151, (-qJD(5) * t81 - t145) * t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t141, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), -t135, 0, 0, 0, 0, 0, -t137, -t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t181, -t118 * t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t105, t104, (-t136 * t79 + t137) * t81, (t127 + t131) * t81, -t62, -qJD(3) * t13 + qJD(4) * t37 + t100 * t70 - t165, -qJD(2) * t37 - qJD(3) * t14 + t134 * t81 - t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, -t152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, t103, t98 * t141, t61, -t146, t106, t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t152, t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t6;
