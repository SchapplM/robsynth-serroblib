% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPPRP2
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRP2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP2_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP2_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP2_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:49:28
% EndTime: 2019-12-31 17:49:30
% DurationCPUTime: 1.03s
% Computational Cost: add. (1530->219), mult. (3152->253), div. (0->0), fcn. (2228->12), ass. (0->125)
t96 = qJ(1) + pkin(7);
t89 = sin(t96);
t167 = g(1) * t89;
t91 = cos(t96);
t135 = g(2) * t91 - t167;
t100 = cos(pkin(7));
t83 = -t100 * pkin(1) - pkin(2);
t144 = qJDD(1) * t83;
t67 = qJDD(3) + t144;
t177 = -t135 - t67;
t129 = g(1) * t91 + g(2) * t89;
t162 = cos(qJ(4));
t132 = qJDD(1) * t162;
t102 = sin(qJ(4));
t141 = qJDD(1) * t102;
t97 = sin(pkin(8));
t99 = cos(pkin(8));
t124 = -t99 * t132 + t97 * t141;
t65 = t102 * t99 + t162 * t97;
t59 = t65 * qJD(4);
t26 = qJD(1) * t59 + t124;
t137 = t162 * t99;
t130 = qJD(1) * t137;
t151 = t102 * t97;
t136 = qJD(1) * t151;
t54 = -t130 + t136;
t134 = qJD(4) * t162;
t145 = qJD(4) * t102;
t58 = -t99 * t134 + t97 * t145;
t156 = -t65 * t26 + t58 * t54;
t115 = t137 - t151;
t171 = t65 * qJD(1);
t138 = qJD(4) * t130 + t97 * t132 + t99 * t141;
t25 = qJD(4) * t136 - t138;
t174 = t115 * t25 + t171 * t59;
t176 = t174 - t156;
t175 = t174 + t156;
t168 = t171 ^ 2;
t51 = t54 ^ 2;
t173 = -t51 - t168;
t172 = -t51 + t168;
t149 = pkin(1) * qJDD(1);
t98 = sin(pkin(7));
t79 = t98 * pkin(1) + qJ(3);
t163 = pkin(6) + t79;
t60 = t163 * t97;
t61 = t163 * t99;
t116 = -t102 * t61 - t162 * t60;
t13 = t115 * qJD(3) + t116 * qJD(4);
t23 = -t102 * t60 + t162 * t61;
t95 = pkin(8) + qJ(4);
t88 = sin(t95);
t170 = -t13 * qJD(4) - t23 * qJDD(4) + t135 * t88;
t62 = qJD(1) * qJD(3) + t79 * qJDD(1);
t85 = t99 * qJDD(2);
t39 = t85 + (-pkin(6) * qJDD(1) - t62) * t97;
t143 = t99 * qJDD(1);
t44 = t97 * qJDD(2) + t99 * t62;
t40 = pkin(6) * t143 + t44;
t153 = pkin(6) * qJD(1);
t72 = t79 * qJD(1);
t87 = t99 * qJD(2);
t41 = t87 + (-t72 - t153) * t97;
t139 = t102 * t39 + t41 * t134 + t162 * t40;
t90 = cos(t95);
t169 = -g(3) * t88 - t129 * t90 + t139;
t164 = t99 * pkin(3);
t103 = sin(qJ(1));
t161 = t103 * pkin(1);
t160 = t171 * t54;
t158 = t90 * t91;
t46 = t97 * qJD(2) + t99 * t72;
t82 = pkin(2) + t164;
t104 = cos(qJ(1));
t92 = t104 * pkin(1);
t155 = t91 * t82 + t92;
t93 = t97 ^ 2;
t94 = t99 ^ 2;
t154 = t93 + t94;
t42 = t99 * t153 + t46;
t152 = t102 * t42;
t150 = t88 * qJ(5);
t148 = qJDD(4) * pkin(4);
t12 = t102 * t41 + t162 * t42;
t147 = t12 * qJD(4);
t11 = t162 * t41 - t152;
t146 = qJD(5) - t11;
t140 = qJDD(4) * qJ(5);
t133 = t11 + t152;
t131 = t102 * t40 + t42 * t134 + t41 * t145 - t162 * t39;
t127 = t90 * pkin(4) + t150;
t125 = g(1) * t103 - g(2) * t104;
t122 = -t115 * t26 + t54 * t59;
t43 = -t97 * t62 + t85;
t121 = -t43 * t97 + t44 * t99;
t120 = (-t97 * t72 + t87) * t97 - t46 * t99;
t101 = -pkin(6) - qJ(3);
t119 = -t91 * t101 - t161;
t71 = t83 - t164;
t29 = qJD(4) * t59 - qJDD(4) * t115;
t114 = -g(3) * t90 + t129 * t88 - t131;
t112 = -t144 + t177;
t50 = t71 * qJD(1) + qJD(3);
t48 = t71 * qJDD(1) + qJDD(3);
t14 = t65 * qJD(3) + t23 * qJD(4);
t111 = -g(2) * t158 - t14 * qJD(4) + qJDD(4) * t116 + t90 * t167;
t110 = t116 * t25 - t13 * t54 + t14 * t171 - t23 * t26 - t129;
t17 = t54 * pkin(4) - qJ(5) * t171 + t50;
t109 = t17 * t171 + qJDD(5) - t114;
t108 = t26 * pkin(4) + t25 * qJ(5) + t48;
t107 = 0.2e1 * t171 * qJD(4) + t124;
t27 = -t58 * qJD(4) + t65 * qJDD(4);
t24 = pkin(4) * t171 + t54 * qJ(5);
t21 = -pkin(4) * t115 - t65 * qJ(5) + t71;
t18 = t59 * pkin(4) + t58 * qJ(5) - t65 * qJD(5);
t16 = (t54 - t136) * qJD(4) + t138;
t15 = (t54 + t136) * qJD(4) - t138;
t10 = qJD(4) * qJ(5) + t12;
t9 = -qJD(4) * pkin(4) + t146;
t6 = -t171 * t58 - t25 * t65;
t5 = -qJD(5) * t171 + t108;
t3 = -t42 * t145 + t139;
t2 = qJDD(5) + t131 - t148;
t1 = t140 + (qJD(5) - t152) * qJD(4) + t139;
t4 = [0, 0, 0, 0, 0, qJDD(1), t125, g(1) * t104 + g(2) * t103, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t100 * t149 - t135, -0.2e1 * t98 * t149 + t129, 0, (t125 + (t100 ^ 2 + t98 ^ 2) * t149) * pkin(1), t93 * qJDD(1), 0.2e1 * t97 * t143, 0, t94 * qJDD(1), 0, 0, t112 * t99, -t112 * t97, t62 * t154 + t121 - t129, t67 * t83 - g(1) * (-t89 * pkin(2) + t91 * qJ(3) - t161) - g(2) * (t91 * pkin(2) + t89 * qJ(3) + t92) + t121 * t79 - t120 * qJD(3), t6, -t176, t27, t122, -t29, 0, -t115 * t48 + t71 * t26 + t50 * t59 + t111, -t71 * t25 + t48 * t65 - t50 * t58 + t170, t11 * t58 + t115 * t3 - t12 * t59 + t131 * t65 + t110, t3 * t23 + t12 * t13 - t131 * t116 - t11 * t14 + t48 * t71 - g(1) * (-t89 * t82 + t119) - g(2) * (-t89 * t101 + t155), t6, t27, t176, 0, t29, t122, -t115 * t5 + t17 * t59 + t18 * t54 + t21 * t26 + t111, t1 * t115 - t10 * t59 + t2 * t65 - t9 * t58 + t110, t17 * t58 - t171 * t18 + t21 * t25 - t5 * t65 - t170, t1 * t23 + t10 * t13 + t5 * t21 + t17 * t18 - t2 * t116 + t9 * t14 - g(1) * t119 - g(2) * (pkin(4) * t158 + t91 * t150 + t155) + (-g(1) * (-t127 - t82) + g(2) * t101) * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t43 * t99 + t44 * t97 - g(3), 0, 0, 0, 0, 0, 0, -t29, -t27, t175, -t11 * t59 - t115 * t131 - t12 * t58 + t3 * t65 - g(3), 0, 0, 0, 0, 0, 0, -t29, t175, t27, t1 * t65 - t10 * t58 - t115 * t2 + t9 * t59 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t143, t97 * qJDD(1), -t154 * qJD(1) ^ 2, t120 * qJD(1) - t177, 0, 0, 0, 0, 0, 0, t107, -t15, t173, t11 * t171 + t12 * t54 + t135 + t48, 0, 0, 0, 0, 0, 0, t107, t173, t15, t10 * t54 + (-qJD(5) - t9) * t171 + t108 + t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t160, t172, t16, -t160, -t124, qJDD(4), -t171 * t50 + t114 + t147, t133 * qJD(4) + t50 * t54 - t169, 0, 0, t160, t16, -t172, qJDD(4), t124, -t160, -t24 * t54 - t109 + t147 + 0.2e1 * t148, pkin(4) * t25 - t26 * qJ(5) + (t10 - t12) * t171 + (t9 - t146) * t54, 0.2e1 * t140 - t17 * t54 + t24 * t171 + (0.2e1 * qJD(5) - t133) * qJD(4) + t169, -t2 * pkin(4) - g(3) * t127 + t1 * qJ(5) + t146 * t10 - t9 * t12 - t17 * t24 + t129 * (pkin(4) * t88 - qJ(5) * t90); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) + t160, t16, -qJD(4) ^ 2 - t168, -t10 * qJD(4) + t109 - t148;];
tau_reg = t4;
