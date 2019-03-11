% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RPRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRRPP5_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP5_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:44:42
% EndTime: 2019-03-09 04:44:49
% DurationCPUTime: 2.46s
% Computational Cost: add. (2717->239), mult. (6218->375), div. (0->0), fcn. (5810->6), ass. (0->131)
t97 = sin(qJ(4));
t145 = t97 * qJ(5);
t161 = pkin(4) + pkin(5);
t99 = cos(qJ(4));
t112 = t161 * t99 + t145;
t140 = t99 * qJD(5);
t176 = qJD(4) * t112 - t140;
t158 = cos(qJ(3));
t95 = sin(pkin(9));
t96 = cos(pkin(9));
t98 = sin(qJ(3));
t73 = t158 * t95 + t98 * t96;
t64 = t73 * qJD(3);
t153 = t97 * t64;
t89 = qJD(4) * t99;
t131 = t73 * t89;
t128 = qJD(3) * t158;
t152 = t98 * t95;
t63 = qJD(3) * t152 - t96 * t128;
t40 = -t97 * t63 + t131;
t129 = t158 * t96;
t72 = -t129 + t152;
t12 = t73 * t153 + t40 * t72;
t175 = -0.2e1 * t12;
t150 = t99 * t64;
t151 = t99 * t63;
t88 = qJD(4) * t97;
t37 = t73 * t88 + t151;
t11 = -0.2e1 * t73 * t150 + 0.2e1 * t37 * t72;
t172 = t64 * qJ(5) + t72 * qJD(5);
t146 = qJ(5) * t99;
t166 = t161 * t97 - t146;
t174 = 0.4e1 * t73;
t173 = 0.2e1 * t172;
t93 = t97 ^ 2;
t94 = t99 ^ 2;
t147 = t93 - t94;
t149 = pkin(7) + qJ(2);
t76 = t149 * t95;
t77 = t149 * t96;
t48 = t158 * t77 - t98 * t76;
t171 = t48 * qJD(3);
t170 = t147 * qJD(4);
t139 = t99 * qJD(6);
t59 = t64 * pkin(4);
t123 = t64 * pkin(3) + t63 * pkin(8);
t27 = t76 * t128 - qJD(2) * t129 + (qJD(2) * t95 + qJD(3) * t77) * t98;
t86 = -t96 * pkin(2) - pkin(1);
t41 = t72 * pkin(3) - t73 * pkin(8) + t86;
t7 = t99 * t123 + t97 * t27 - t41 * t88 - t48 * t89;
t5 = -t59 - t7;
t168 = -t37 * qJ(6) + t73 * t139 - t5;
t43 = t97 * t48;
t18 = t99 * t41 - t43;
t19 = t97 * t41 + t99 * t48;
t6 = -t97 * t123 + t99 * t27 - t41 * t89 + t48 * t88;
t165 = t6 * t97 - t7 * t99 + (t18 * t97 - t19 * t99) * qJD(4);
t13 = t72 * qJ(5) + t19;
t14 = -t72 * pkin(4) - t18;
t3 = -t6 + t172;
t164 = t3 * t97 - t5 * t99 + (t13 * t99 + t14 * t97) * qJD(4);
t118 = t99 * pkin(4) + t145;
t163 = t118 * qJD(4) - t140;
t135 = t97 * t151;
t70 = t73 ^ 2;
t16 = t135 * t174 + 0.2e1 * t170 * t70;
t162 = 0.2e1 * t64;
t74 = -0.2e1 * t170;
t101 = 0.2e1 * qJD(5);
t160 = pkin(8) * t64;
t159 = pkin(8) * t72;
t28 = t73 * qJD(2) + t171;
t47 = t158 * t76 + t98 * t77;
t157 = t47 * t28;
t156 = t73 * t63;
t155 = t73 * t97;
t154 = t73 * t99;
t148 = pkin(8) - qJ(6);
t144 = t97 * qJ(6);
t143 = qJD(4) * t73;
t142 = t97 * qJD(5);
t141 = t97 * qJD(6);
t138 = t72 * t162;
t137 = -0.2e1 * pkin(3) * qJD(4);
t133 = pkin(8) * t88;
t132 = pkin(8) * t89;
t130 = t97 * t89;
t80 = t148 * t99;
t127 = -pkin(4) * t88 + t142;
t126 = t70 * t130;
t125 = 0.2e1 * (t95 ^ 2 + t96 ^ 2) * qJD(2);
t4 = -t171 + t166 * t63 + (-qJD(2) - t176) * t73;
t69 = pkin(3) + t112;
t124 = -t69 * t143 + t4;
t122 = pkin(3) * t63 - t160;
t121 = pkin(3) * t73 + t159;
t117 = pkin(4) * t97 - t146;
t115 = t13 * t97 - t14 * t99;
t114 = t18 * t99 + t19 * t97;
t38 = t72 * t89 + t153;
t35 = t72 * t88 - t150;
t15 = -t166 * t73 - t47;
t49 = (-pkin(5) * t97 + t146) * qJD(4) + t127;
t111 = qJD(4) * t15 + t49 * t73 - t63 * t69;
t75 = -pkin(3) - t118;
t8 = -t117 * t63 + t163 * t73 + t28;
t110 = -t8 + (t73 * t75 - t159) * qJD(4);
t108 = qJ(6) * t131 + t73 * t141 - t63 * t144 - t6;
t20 = t117 * t73 + t47;
t61 = qJ(5) * t89 + t127;
t107 = -qJD(4) * t20 + t61 * t73 + t63 * t75 + t160;
t105 = -t166 * qJD(4) + t142;
t1 = -t64 * pkin(5) - t168;
t10 = t73 * t144 + t13;
t2 = t108 + t172;
t9 = t43 + (-qJ(6) * t73 - t41) * t99 - t161 * t72;
t104 = -t1 * t99 + t2 * t97 + (t10 * t99 + t9 * t97) * qJD(4);
t103 = -t115 * qJD(4) + t3 * t99 + t5 * t97;
t102 = -t114 * qJD(4) - t6 * t99 - t7 * t97;
t90 = qJ(5) * t101;
t82 = -0.2e1 * t130;
t81 = 0.2e1 * t130;
t79 = t148 * t97;
t62 = qJD(4) * t80 - t141;
t60 = -t148 * t88 - t139;
t31 = (-t93 - t94) * t63;
t24 = -0.2e1 * t94 * t156 - 0.2e1 * t126;
t23 = -0.2e1 * t93 * t156 + 0.2e1 * t126;
t22 = t143 * t147 + t135;
t21 = t130 * t174 - t147 * t63;
t17 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t125, qJ(2) * t125, -0.2e1 * t156, 0.2e1 * t63 * t72 - 0.2e1 * t73 * t64, 0, t138, 0, 0, t86 * t162, -0.2e1 * t86 * t63, 0.2e1 * t27 * t72 + 0.2e1 * t28 * t73 - 0.2e1 * t47 * t63 - 0.2e1 * t48 * t64, -0.2e1 * t48 * t27 + 0.2e1 * t157, t24, t16, -t11, t23, t175, t138, 0.2e1 * t28 * t155 + 0.2e1 * t18 * t64 + 0.2e1 * t40 * t47 + 0.2e1 * t7 * t72, 0.2e1 * t28 * t154 - 0.2e1 * t19 * t64 - 0.2e1 * t37 * t47 + 0.2e1 * t6 * t72, 0.2e1 * t114 * t63 + 0.2e1 * t165 * t73, 0.2e1 * t18 * t7 - 0.2e1 * t19 * t6 + 0.2e1 * t157, t24, -t11, -t16, t138, 0.2e1 * t12, t23, -0.2e1 * t14 * t64 + 0.2e1 * t8 * t155 + 0.2e1 * t40 * t20 - 0.2e1 * t5 * t72, 0.2e1 * t115 * t63 - 0.2e1 * t164 * t73, 0.2e1 * t13 * t64 - 0.2e1 * t8 * t154 + 0.2e1 * t37 * t20 + 0.2e1 * t3 * t72, 0.2e1 * t13 * t3 + 0.2e1 * t14 * t5 + 0.2e1 * t20 * t8, t24, -t16, t11, t23, t175, t138, -0.2e1 * t1 * t72 - 0.2e1 * t40 * t15 - 0.2e1 * t4 * t155 - 0.2e1 * t9 * t64, 0.2e1 * t10 * t64 - 0.2e1 * t37 * t15 + 0.2e1 * t4 * t154 + 0.2e1 * t2 * t72, 0.2e1 * (-t10 * t97 + t9 * t99) * t63 + 0.2e1 * t104 * t73, 0.2e1 * t9 * t1 + 0.2e1 * t10 * t2 + 0.2e1 * t15 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, -t63, 0, 0, 0, 0, 0, 0, 0, 0, -t35, -t38, -t31, -t165, 0, 0, 0, 0, 0, 0, -t35, -t31, t38, t164, 0, 0, 0, 0, 0, 0, -t35, t38, t31, t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, 0, -t64, 0, -t28, t27, 0, 0, -t22, -t21, t38, t22, -t35, 0, -t28 * t99 + t122 * t97 + (-t121 * t99 + t47 * t97) * qJD(4), t28 * t97 + t122 * t99 + (t121 * t97 + t47 * t99) * qJD(4), t102, -t28 * pkin(3) + t102 * pkin(8), -t22, t38, t21, 0, t35, t22, -t107 * t97 + t110 * t99, t103, t107 * t99 + t110 * t97, t103 * pkin(8) - t20 * t61 + t8 * t75, -t22, t21, -t38, t22, -t35, 0, -t111 * t97 + t124 * t99 - t62 * t72 - t79 * t64, t111 * t99 + t124 * t97 + t60 * t72 + t80 * t64 (-t62 * t73 + t63 * t79 - t2 + (t73 * t80 - t9) * qJD(4)) * t99 + (t60 * t73 - t63 * t80 - t1 + (t73 * t79 + t10) * qJD(4)) * t97, t1 * t79 + t10 * t60 + t15 * t49 + t2 * t80 + t4 * t69 + t9 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97 * t60 - t99 * t62 + (t79 * t97 + t80 * t99) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, t74, 0, t82, 0, 0, t97 * t137, t99 * t137, 0, 0, t81, 0, -t74, 0, 0, t82, 0.2e1 * t61 * t99 + 0.2e1 * t75 * t88, 0, 0.2e1 * t61 * t97 - 0.2e1 * t75 * t89, -0.2e1 * t75 * t61, t81, -t74, 0, t82, 0, 0, 0.2e1 * t49 * t99 - 0.2e1 * t69 * t88, 0.2e1 * t49 * t97 + 0.2e1 * t69 * t89, -0.2e1 * t60 * t99 - 0.2e1 * t62 * t97 + 0.2e1 * (-t79 * t99 + t80 * t97) * qJD(4), 0.2e1 * t69 * t49 + 0.2e1 * t80 * t60 + 0.2e1 * t79 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, 0, -t40, t64, t7, t6, 0, 0, 0, -t37, 0, t64, t40, 0, -t5 + t59, t118 * t63 + (t117 * qJD(4) - t142) * t73, -t6 + t173, -t5 * pkin(4) + t3 * qJ(5) + t13 * qJD(5), 0, 0, t37, 0, -t40, t64 (pkin(5) + t161) * t64 + t168, t108 + t173, t105 * t73 - t112 * t63, t2 * qJ(5) + t10 * qJD(5) - t1 * t161; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t88, -t89, 0, 0, 0, 0, 0, 0, 0, 0, -t88, 0, t89, t61, 0, 0, 0, 0, 0, 0, -t88, t89, 0, t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, 0, -t88, 0, -t132, t133, 0, 0, 0, t89, 0, 0, t88, 0, -t132, -t163, -t133, -t163 * pkin(8), 0, 0, -t89, 0, -t88, 0, -t62, t60, t176, t60 * qJ(5) + t80 * qJD(5) - t161 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, t90, 0, 0, 0, 0, 0, 0, 0, t101, 0, t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, -t37, 0, t5, 0, 0, 0, 0, 0, 0, -t64, 0, t37, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, 0, t132, 0, 0, 0, 0, 0, 0, 0, 0, -t89, t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, -t37, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t88, t89, 0, t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t17;
