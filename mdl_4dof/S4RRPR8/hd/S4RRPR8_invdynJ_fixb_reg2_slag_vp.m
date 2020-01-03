% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRPR8_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR8_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR8_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR8_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:08:26
% EndTime: 2019-12-31 17:08:28
% DurationCPUTime: 1.33s
% Computational Cost: add. (1132->239), mult. (2516->307), div. (0->0), fcn. (1498->6), ass. (0->146)
t103 = sin(qJ(2));
t106 = cos(qJ(2));
t105 = cos(qJ(4));
t149 = qJD(4) * t105;
t102 = sin(qJ(4));
t150 = qJD(4) * t102;
t157 = t106 * t102;
t190 = qJD(2) * t157 + t103 * t149 - t106 * t150;
t147 = qJD(1) * qJD(2);
t139 = t106 * t147;
t85 = t103 * qJDD(1);
t81 = pkin(5) * t85;
t142 = pkin(5) * t139 + qJDD(3) + t81;
t180 = pkin(2) + pkin(3);
t189 = -t139 - t85;
t14 = pkin(6) * t189 - t180 * qJDD(2) + t142;
t140 = t103 * t147;
t146 = t106 * qJDD(1);
t186 = t140 - t146;
t82 = pkin(5) * t146;
t97 = qJDD(2) * qJ(3);
t98 = qJD(2) * qJD(3);
t25 = -pkin(5) * t140 + t82 + t97 + t98;
t15 = pkin(6) * t186 + t25;
t153 = qJD(1) * t103;
t83 = pkin(5) * t153;
t188 = -pkin(6) * t153 + qJD(3) + t83;
t27 = -qJD(2) * t180 + t188;
t152 = qJD(1) * t106;
t84 = pkin(5) * t152;
t50 = -pkin(6) * t152 + t84;
t99 = qJD(2) * qJ(3);
t40 = t50 + t99;
t9 = t102 * t27 + t105 * t40;
t2 = -t9 * qJD(4) - t102 * t15 + t105 * t14;
t45 = t102 * t103 + t105 * t106;
t37 = t45 * qJD(1);
t119 = t45 * qJD(4);
t5 = qJD(1) * t119 - t102 * t186 + t105 * t189;
t96 = qJD(2) - qJD(4);
t185 = t37 * t96 + t5;
t86 = t103 * qJ(3);
t138 = pkin(1) + t86;
t118 = t106 * t180 + t138;
t24 = t118 * qJD(1);
t104 = sin(qJ(1));
t125 = -t103 * t105 + t157;
t32 = t125 * t104;
t107 = cos(qJ(1));
t156 = t106 * t107;
t160 = t103 * t107;
t34 = t102 * t156 - t105 * t160;
t39 = -t102 * t152 + t105 * t153;
t187 = g(1) * t34 + g(2) * t32 + g(3) * t45 - t24 * t39 + t2;
t6 = qJD(1) * t190 + t45 * qJDD(1) - t105 * t140;
t184 = t39 * t96 + t6;
t183 = g(1) * t107 + g(2) * t104;
t182 = t180 * t103;
t165 = pkin(5) * qJDD(2);
t89 = t106 * pkin(2);
t124 = -t138 - t89;
t41 = t124 * qJD(1);
t168 = t89 + t86;
t55 = -pkin(1) - t168;
t181 = (qJD(1) * t55 + t41) * qJD(2) - t165;
t179 = pkin(5) - pkin(6);
t8 = -t102 * t40 + t105 * t27;
t178 = t8 * t96;
t177 = t9 * t96;
t176 = g(1) * t104;
t175 = g(2) * t107;
t88 = t106 * pkin(3);
t173 = t39 * t37;
t52 = -qJ(3) * t102 - t105 * t180;
t171 = t52 * qJD(4) - t102 * t50 + t105 * t188;
t53 = qJ(3) * t105 - t102 * t180;
t170 = -t53 * qJD(4) - t102 * t188 - t105 * t50;
t167 = pkin(1) * t107 + pkin(5) * t104;
t164 = qJ(3) * t106;
t127 = pkin(2) * t103 - t164;
t148 = t103 * qJD(3);
t36 = qJD(2) * t127 - t148;
t163 = qJD(1) * t36;
t162 = qJDD(2) * pkin(2);
t161 = t103 * t104;
t110 = qJD(1) ^ 2;
t159 = t103 * t110;
t158 = t104 * t106;
t100 = t103 ^ 2;
t101 = t106 ^ 2;
t155 = -t100 + t101;
t154 = t100 + t101;
t151 = qJD(2) * t103;
t145 = t88 + t168;
t144 = -g(1) * t160 - g(2) * t161 + g(3) * t106;
t143 = t37 ^ 2 - t39 ^ 2;
t61 = t179 * t106;
t136 = pkin(2) * t156 + qJ(3) * t160 + t167;
t135 = -t81 - t144;
t133 = -qJD(2) * pkin(2) + qJD(3);
t132 = t96 ^ 2;
t131 = t103 * t139;
t130 = t154 * qJDD(1) * pkin(5);
t109 = qJD(2) ^ 2;
t129 = pkin(5) * t109 + t175;
t60 = t179 * t103;
t20 = -t102 * t61 + t105 * t60;
t21 = t102 * t60 + t105 * t61;
t54 = t133 + t83;
t59 = t84 + t99;
t126 = t103 * t59 - t106 * t54;
t30 = t142 - t162;
t123 = t164 - t182;
t122 = -0.2e1 * pkin(1) * t147 - t165;
t1 = t102 * t14 + t105 * t15 + t149 * t27 - t150 * t40;
t117 = 0.2e1 * qJDD(1) * pkin(1) - t129;
t13 = qJDD(1) * t124 + t163;
t115 = -qJDD(1) * t55 - t129 - t13 - t163;
t23 = qJD(2) * t123 + t148;
t114 = -qJD(2) * t126 + t30 * t103 + t25 * t106;
t33 = t45 * t104;
t35 = t45 * t107;
t112 = -g(1) * t35 - g(2) * t33 + g(3) * t125 - t24 * t37 + t1;
t95 = qJDD(2) - qJDD(4);
t90 = t107 * pkin(5);
t78 = g(1) * t158;
t74 = qJ(3) * t156;
t72 = qJ(3) * t158;
t69 = t106 * t159;
t58 = t155 * t110;
t57 = qJDD(2) * t106 - t103 * t109;
t56 = qJDD(2) * t103 + t106 * t109;
t51 = qJD(2) * t61;
t49 = t179 * t151;
t47 = t127 * qJD(1);
t44 = qJDD(1) * t101 - 0.2e1 * t131;
t43 = qJDD(1) * t100 + 0.2e1 * t131;
t42 = pkin(1) + t145;
t31 = t123 * qJD(1);
t26 = t103 * t146 + t147 * t155;
t19 = qJD(2) * t45 - t119;
t18 = -t105 * t151 + t190;
t7 = qJD(1) * t23 + qJDD(1) * t118;
t4 = -qJD(4) * t21 + t102 * t49 + t105 * t51;
t3 = qJD(4) * t20 + t102 * t51 - t105 * t49;
t10 = [0, 0, 0, 0, 0, qJDD(1), -t175 + t176, t183, 0, 0, t43, 0.2e1 * t26, t56, t44, t57, 0, t103 * t122 + t106 * t117 + t78, t122 * t106 + (-t117 - t176) * t103, 0.2e1 * t130 - t183, -g(1) * (-pkin(1) * t104 + t90) - g(2) * t167 + (pkin(5) ^ 2 * t154 + pkin(1) ^ 2) * qJDD(1), t43, t56, -0.2e1 * t26, 0, -t57, t44, t103 * t181 + t106 * t115 + t78, t130 + t114 - t183, -t181 * t106 + (t115 + t176) * t103, pkin(5) * t114 - g(1) * t90 - g(2) * t136 - t124 * t176 + t13 * t55 + t41 * t36, t125 * t5 + t19 * t39, t125 * t6 - t18 * t39 - t19 * t37 + t45 * t5, t125 * t95 - t19 * t96, t18 * t37 + t45 * t6, t18 * t96 + t45 * t95, 0, g(1) * t33 - g(2) * t35 + t18 * t24 - t20 * t95 + t23 * t37 - t4 * t96 + t42 * t6 + t45 * t7, -g(1) * t32 + g(2) * t34 - t125 * t7 + t19 * t24 + t21 * t95 + t23 * t39 + t3 * t96 - t42 * t5, -t1 * t45 + t125 * t2 - t18 * t9 - t19 * t8 + t20 * t5 - t21 * t6 - t3 * t37 - t39 * t4 + t183, t1 * t21 + t9 * t3 + t2 * t20 + t8 * t4 + t7 * t42 + t24 * t23 - g(1) * (-t107 * pkin(6) + t90) - g(2) * (pkin(3) * t156 + t136) + (-g(1) * (t124 - t88) + g(2) * pkin(6)) * t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, -t58, t85, t69, t146, qJDD(2), pkin(1) * t159 + t135, g(3) * t103 - t82 + (pkin(1) * t110 + t183) * t106, 0, 0, -t69, t85, t58, qJDD(2), -t146, t69, 0.2e1 * t162 - qJDD(3) + (-t103 * t41 + t106 * t47) * qJD(1) + t135, -t127 * qJDD(1) + ((t59 - t99) * t103 + (t133 - t54) * t106) * qJD(1), t82 + 0.2e1 * t97 + 0.2e1 * t98 + (qJD(1) * t47 - g(3)) * t103 + (qJD(1) * t41 - t183) * t106, t25 * qJ(3) + t59 * qJD(3) - t30 * pkin(2) - t41 * t47 - g(1) * (-pkin(2) * t160 + t74) - g(2) * (-pkin(2) * t161 + t72) - g(3) * t168 + t126 * qJD(1) * pkin(5), -t173, t143, t185, t173, t184, t95, -t170 * t96 - t31 * t37 - t52 * t95 - t187, t171 * t96 - t31 * t39 + t53 * t95 + t112, t52 * t5 - t53 * t6 + (-t9 - t170) * t39 + (t8 - t171) * t37, -g(1) * t74 - g(2) * t72 - g(3) * t145 + t1 * t53 + t170 * t8 + t171 * t9 + t182 * t183 + t2 * t52 - t24 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(2) - t69, t85, -t100 * t110 - t109, -qJD(2) * t59 + t153 * t41 + t144 + t30, 0, 0, 0, 0, 0, 0, -t102 * t132 - t105 * t95 - t153 * t37, t102 * t95 - t105 * t132 - t153 * t39, -t102 * t184 + t105 * t185, -t24 * t153 + (t2 - t177) * t105 + (t1 + t178) * t102 + t144; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t173, -t143, -t185, -t173, -t184, -t95, -t177 + t187, -t112 - t178, 0, 0;];
tau_reg = t10;
