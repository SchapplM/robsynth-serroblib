% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPRR2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR2_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:45:19
% EndTime: 2019-12-05 15:45:21
% DurationCPUTime: 1.01s
% Computational Cost: add. (1871->188), mult. (3540->253), div. (0->0), fcn. (2713->14), ass. (0->129)
t88 = qJ(2) + pkin(9);
t82 = qJ(4) + t88;
t76 = cos(t82);
t163 = g(3) * t76;
t86 = qJDD(2) + qJDD(4);
t162 = t86 * pkin(4);
t97 = sin(qJ(2));
t143 = qJD(1) * t97;
t100 = cos(qJ(2));
t137 = qJD(1) * t100;
t69 = qJD(2) * pkin(2) + t137;
t91 = sin(pkin(9));
t93 = cos(pkin(9));
t33 = -t91 * t143 + t93 * t69;
t32 = qJD(2) * pkin(3) + t33;
t34 = t93 * t143 + t69 * t91;
t96 = sin(qJ(4));
t99 = cos(qJ(4));
t18 = t32 * t96 + t34 * t99;
t135 = qJD(1) * qJD(2);
t136 = t97 * qJDD(1);
t115 = t100 * t135 + t136;
t141 = qJDD(2) * pkin(2);
t81 = t100 * qJDD(1);
t51 = -t97 * t135 + t141 + t81;
t41 = t93 * t51;
t25 = -t115 * t91 + t41;
t24 = qJDD(2) * pkin(3) + t25;
t26 = t115 * t93 + t91 * t51;
t8 = -t18 * qJD(4) + t99 * t24 - t96 * t26;
t6 = -t162 - t8;
t172 = t6 + t163;
t75 = sin(t82);
t94 = cos(pkin(8));
t153 = t75 * t94;
t92 = sin(pkin(8));
t154 = t75 * t92;
t171 = g(1) * t153 + g(2) * t154;
t166 = pkin(2) * t91;
t77 = pkin(2) * t93 + pkin(3);
t44 = -t96 * t166 + t77 * t99;
t53 = t100 * t91 + t93 * t97;
t48 = t53 * qJD(1);
t52 = t100 * t93 - t91 * t97;
t50 = t52 * qJD(1);
t149 = t44 * qJD(4) + t48 * t96 - t50 * t99;
t87 = qJD(2) + qJD(4);
t170 = t149 * t87;
t101 = qJD(5) ^ 2;
t45 = t99 * t166 + t96 * t77;
t148 = t45 * qJD(4) - t99 * t48 - t50 * t96;
t127 = t148 * t87;
t42 = -pkin(4) - t44;
t43 = pkin(7) + t45;
t169 = t101 * t43 + t42 * t86 + t127;
t168 = -t32 * t99 + t34 * t96;
t16 = pkin(7) * t87 + t18;
t95 = sin(qJ(5));
t98 = cos(qJ(5));
t11 = qJD(3) * t98 - t16 * t95;
t12 = qJD(3) * t95 + t16 * t98;
t140 = t11 * qJD(5);
t124 = t168 * qJD(4) - t96 * t24 - t99 * t26;
t5 = pkin(7) * t86 - t124;
t2 = t95 * qJDD(3) + t98 * t5 + t140;
t139 = t12 * qJD(5);
t80 = t98 * qJDD(3);
t3 = -t95 * t5 - t139 + t80;
t106 = t2 * t98 + (-t11 * t98 - t12 * t95) * qJD(5) - t3 * t95;
t128 = -g(1) * t92 + g(2) * t94;
t167 = -t11 * t95 + t12 * t98;
t165 = pkin(4) * t87;
t72 = g(3) * t75;
t119 = t99 * t52 - t53 * t96;
t47 = t53 * qJD(2);
t49 = t52 * qJD(2);
t9 = t119 * qJD(4) - t96 * t47 + t99 * t49;
t161 = t87 * t9;
t158 = t168 * t87;
t157 = t18 * t87;
t152 = t76 * t92;
t151 = t76 * t94;
t150 = t98 * t86;
t147 = t76 * pkin(4) + t75 * pkin(7);
t79 = cos(t88);
t146 = t100 * pkin(2) + pkin(3) * t79;
t89 = t95 ^ 2;
t90 = t98 ^ 2;
t145 = t89 - t90;
t144 = t89 + t90;
t142 = qJD(5) * t98;
t138 = qJDD(1) - g(3);
t15 = t168 - t165;
t134 = t15 * t142 + t172 * t95;
t85 = t87 ^ 2;
t133 = t95 * t85 * t98;
t132 = t15 * qJD(5) * t95 + t171 * t98;
t131 = -g(1) * t151 - g(2) * t152 - t72;
t78 = sin(t88);
t58 = -pkin(2) * t97 - pkin(3) * t78;
t129 = -pkin(4) * t75 + t58;
t126 = t144 * t86;
t123 = t87 * t95 * t142;
t122 = g(1) * t94 + g(2) * t92;
t28 = t52 * t96 + t53 * t99;
t10 = t28 * qJD(4) + t99 * t47 + t96 * t49;
t121 = -t10 * t87 + t119 * t86;
t117 = t122 * t75;
t116 = t124 - t131;
t114 = pkin(7) * t101 - t157 - t162;
t113 = -qJD(5) * t16 + t128;
t112 = t101 * t28 - t121;
t111 = -pkin(7) * qJDD(5) + (-t168 - t165) * qJD(5);
t110 = -qJDD(5) * t28 + (-t119 * t87 - t9) * qJD(5);
t109 = -g(3) * t100 + t122 * t97;
t108 = -qJDD(5) * t43 + (t42 * t87 - t149) * qJD(5);
t105 = t131 + t106;
t104 = t8 - t163 + t171;
t103 = -qJD(5) * qJD(3) + t122 * t76 - t15 * t87 - t5 + t72;
t102 = qJD(2) ^ 2;
t66 = qJDD(5) * t98 - t101 * t95;
t65 = qJDD(5) * t95 + t101 * t98;
t60 = pkin(7) * t151;
t59 = pkin(7) * t152;
t57 = qJDD(3) + t128;
t40 = t86 * t90 - 0.2e1 * t123;
t39 = t86 * t89 + 0.2e1 * t123;
t29 = -0.2e1 * t145 * t87 * qJD(5) + 0.2e1 * t95 * t150;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t138, 0, 0, 0, 0, 0, 0, qJDD(2) * t100 - t102 * t97, -qJDD(2) * t97 - t100 * t102, 0, -g(3) + (t100 ^ 2 + t97 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, -qJD(2) * t47 + qJDD(2) * t52, -qJD(2) * t49 - qJDD(2) * t53, 0, t25 * t52 + t26 * t53 - t33 * t47 + t34 * t49 - g(3), 0, 0, 0, 0, 0, 0, t121, -t28 * t86 - t161, 0, t10 * t168 + t119 * t8 - t124 * t28 + t18 * t9 - g(3), 0, 0, 0, 0, 0, 0, t110 * t95 - t112 * t98, t110 * t98 + t112 * t95, t28 * t126 + t144 * t161, t15 * t10 + t106 * t28 - t119 * t6 + t167 * t9 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t81 + t109, t122 * t100 - t138 * t97, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t93 * t141 - t91 * t136 - g(3) * t79 + t41 + t122 * t78 + (-t91 * t137 + t48) * qJD(2), -t93 * t136 + g(3) * t78 + (-t51 - t141) * t91 + t122 * t79 + (-t93 * t137 + t50) * qJD(2), 0, t33 * t48 - t34 * t50 + (t25 * t93 + t26 * t91 + t109) * pkin(2), 0, 0, 0, 0, 0, t86, t44 * t86 + t104 - t127, -t45 * t86 + t116 - t170, 0, -g(3) * t146 - t122 * t58 - t124 * t45 + t148 * t168 + t149 * t18 + t8 * t44, t39, t29, t65, t40, t66, 0, t108 * t95 + (-t172 - t169) * t98 + t132, t108 * t98 + (-t117 + t169) * t95 + t134, t43 * t126 + t144 * t170 + t105, t6 * t42 - g(1) * (t129 * t94 + t60) - g(2) * (t129 * t92 + t59) - g(3) * (t146 + t147) + t148 * t15 + t106 * t43 + t167 * t149; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, 0, 0, 0, 0, 0, t66, -t65, 0, qJD(5) * t167 + t2 * t95 + t3 * t98 + t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, t104 + t157, t116 - t158, 0, 0, t39, t29, t65, t40, t66, 0, t111 * t95 + (-t114 - t172) * t98 + t132, t111 * t98 + (-t117 + t114) * t95 + t134, pkin(7) * t126 + t144 * t158 + t105, -t6 * pkin(4) - t15 * t18 - g(1) * (-pkin(4) * t153 + t60) - g(2) * (-pkin(4) * t154 + t59) - g(3) * t147 + t167 * t168 + t106 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t133, t145 * t85, t95 * t86, t133, t150, qJDD(5), t103 * t95 + t113 * t98 + t139 + t80, t140 + (-qJDD(3) - t113) * t95 + t103 * t98, 0, 0;];
tau_reg = t1;
