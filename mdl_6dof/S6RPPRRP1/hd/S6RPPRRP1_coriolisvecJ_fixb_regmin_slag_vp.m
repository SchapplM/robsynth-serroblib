% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
% 
% Output:
% tauc_reg [6x24]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRRP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:58:43
% EndTime: 2019-03-09 01:58:47
% DurationCPUTime: 1.34s
% Computational Cost: add. (2678->233), mult. (6446->315), div. (0->0), fcn. (4672->8), ass. (0->126)
t108 = sin(qJ(4));
t110 = cos(qJ(4));
t103 = sin(pkin(10));
t95 = sin(pkin(9)) * pkin(1) + qJ(3);
t88 = t95 * qJD(1);
t105 = cos(pkin(10));
t99 = t105 * qJD(2);
t61 = t99 + (-pkin(7) * qJD(1) - t88) * t103;
t144 = qJD(1) * t105;
t70 = t103 * qJD(2) + t105 * t88;
t62 = pkin(7) * t144 + t70;
t26 = t108 * t61 + t110 * t62;
t176 = qJD(4) * t26;
t86 = t103 * t110 + t105 * t108;
t107 = sin(qJ(5));
t149 = t103 * t108;
t136 = qJD(1) * t149;
t94 = t110 * t144;
t78 = t94 - t136;
t74 = qJD(5) - t78;
t133 = t107 * t74;
t109 = cos(qJ(5));
t79 = t86 * qJD(1);
t65 = qJD(4) * t107 + t109 * t79;
t175 = t65 * t133;
t174 = -t108 * t62 + t110 * t61;
t169 = pkin(7) + t95;
t82 = t169 * t103;
t83 = t169 * t105;
t49 = t108 * t83 + t110 * t82;
t23 = qJD(4) * pkin(8) + t26;
t87 = -cos(pkin(9)) * pkin(1) - pkin(3) * t105 - pkin(2);
t75 = t87 * qJD(1) + qJD(3);
t32 = -pkin(4) * t78 - pkin(8) * t79 + t75;
t12 = t107 * t32 + t109 * t23;
t85 = -t110 * t105 + t149;
t115 = t85 * qJD(3);
t17 = -qJD(1) * t115 + qJD(4) * t174;
t93 = qJD(4) * t94;
t118 = qJD(4) * t136 - t93;
t81 = t86 * qJD(4);
t71 = qJD(1) * t81;
t41 = t71 * pkin(4) + t118 * pkin(8);
t37 = t109 * t41;
t113 = -t12 * qJD(5) - t107 * t17 + t37;
t141 = t109 * qJD(4);
t143 = qJD(5) * t107;
t33 = -qJD(5) * t141 + t109 * t118 + t79 * t143;
t1 = pkin(5) * t71 + qJ(6) * t33 - qJD(6) * t65 + t113;
t63 = t107 * t79 - t141;
t7 = -qJ(6) * t63 + t12;
t173 = t74 * t7 + t1;
t163 = t86 * t71;
t80 = t85 * qJD(4);
t123 = -t74 * t80 + t163;
t138 = t86 * t143;
t172 = -t123 * t109 + t74 * t138;
t171 = t65 ^ 2;
t11 = -t107 * t23 + t109 * t32;
t6 = -qJ(6) * t65 + t11;
t5 = pkin(5) * t74 + t6;
t170 = t5 - t6;
t162 = -qJ(6) - pkin(8);
t135 = qJD(5) * t162;
t151 = qJ(6) * t109;
t52 = pkin(4) * t79 - pkin(8) * t78;
t47 = t109 * t52;
t168 = -pkin(5) * t79 + t109 * t135 + t78 * t151 - t47 + (-qJD(6) + t174) * t107;
t167 = t63 * t78;
t166 = t63 * t79;
t165 = t65 * t79;
t164 = t65 * t80;
t114 = t107 * t118;
t34 = t65 * qJD(5) - t114;
t161 = (-t34 * t86 + t63 * t80) * t109;
t160 = t107 * t52 + t109 * t174;
t142 = qJD(5) * t109;
t159 = -t107 * t34 - t63 * t142;
t158 = -t33 * t85 + t65 * t81;
t50 = -t108 * t82 + t110 * t83;
t43 = t109 * t50;
t44 = pkin(4) * t85 - pkin(8) * t86 + t87;
t157 = t107 * t44 + t43;
t152 = qJ(6) * t107;
t156 = t109 * qJD(6) + t107 * t135 + t78 * t152 - t160;
t155 = t107 * t71;
t150 = qJD(5) * t63;
t146 = t80 * qJD(4);
t145 = t103 ^ 2 + t105 ^ 2;
t140 = qJD(1) * qJD(3);
t28 = -t49 * qJD(4) - t115;
t53 = pkin(4) * t81 + pkin(8) * t80;
t139 = t107 * t53 + t109 * t28 + t44 * t142;
t137 = t86 * t142;
t18 = t86 * t140 + t176;
t134 = t109 * t74;
t132 = qJD(1) * t145;
t130 = t65 * t137;
t129 = pkin(8) * qJD(5) * t74 + t18;
t116 = t107 * t41 + t109 * t17 + t32 * t142 - t23 * t143;
t2 = -qJ(6) * t34 - qJD(6) * t63 + t116;
t128 = -t74 * t5 + t2;
t127 = -t107 * t7 - t109 * t5;
t126 = t107 * t5 - t109 * t7;
t22 = -qJD(4) * pkin(4) - t174;
t125 = t18 * t86 - t22 * t80;
t124 = -t85 * t34 - t81 * t63;
t122 = t103 * (-t103 * t88 + t99) - t105 * t70;
t120 = qJ(6) * t80 - qJD(6) * t86;
t119 = t109 * t71 + t78 * t133 - t74 * t143;
t9 = pkin(5) * t34 + t18;
t117 = -pkin(8) * t71 + t74 * t22;
t112 = -t123 * t107 - t74 * t137;
t29 = t86 * qJD(3) + t50 * qJD(4);
t90 = t162 * t109;
t89 = t162 * t107;
t73 = t81 * qJD(4);
t60 = t63 ^ 2;
t48 = t109 * t53;
t40 = t109 * t44;
t15 = pkin(5) * t63 + qJD(6) + t22;
t14 = -t86 * t152 + t157;
t13 = pkin(5) * t85 - t107 * t50 - t86 * t151 + t40;
t4 = -qJ(6) * t137 + (-qJD(5) * t50 + t120) * t107 + t139;
t3 = pkin(5) * t81 - t107 * t28 + t48 + t120 * t109 + (-t43 + (qJ(6) * t86 - t44) * t107) * qJD(5);
t8 = [0, 0, 0, 0, 0, 0, 0.2e1 * qJD(3) * t132 (t95 * t132 - t122) * qJD(3), -t118 * t86 - t79 * t80, t118 * t85 - t80 * t78 - t79 * t81 - t163, -t146, -t73, 0, -qJD(4) * t29 + t71 * t87 + t75 * t81, -t28 * qJD(4) - t87 * t118 - t75 * t80, -t65 * t138 + (-t33 * t86 - t164) * t109, -t130 + (t164 + (t33 + t150) * t86) * t107 + t161, t158 - t172, t112 + t124, t71 * t85 + t74 * t81 (-t50 * t142 + t48) * t74 + t40 * t71 + (-t23 * t142 + t37) * t85 + t11 * t81 + t29 * t63 + t49 * t34 + t22 * t137 + ((-qJD(5) * t44 - t28) * t74 - t50 * t71 + (-qJD(5) * t32 - t17) * t85 + t125) * t107 -(-t50 * t143 + t139) * t74 - t157 * t71 - t116 * t85 - t12 * t81 + t29 * t65 - t49 * t33 - t22 * t138 + t125 * t109, t13 * t33 - t14 * t34 - t3 * t65 - t4 * t63 - t127 * t80 + (t126 * qJD(5) - t1 * t109 - t107 * t2) * t86, t2 * t14 + t7 * t4 + t1 * t13 + t5 * t3 + t9 * (pkin(5) * t107 * t86 + t49) + t15 * ((-t107 * t80 + t137) * pkin(5) + t29); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73, t146, 0, 0, 0, 0, 0, t112 - t124, t158 + t172, t130 + (-t164 + (-t33 + t150) * t86) * t107 + t161, t15 * t81 + t85 * t9 + t126 * t80 + (t127 * qJD(5) - t1 * t107 + t109 * t2) * t86; 0, 0, 0, 0, 0, 0, -t145 * qJD(1) ^ 2, t122 * qJD(1), 0, 0, 0, 0, 0, 0.2e1 * t79 * qJD(4), t93 + (t78 - t136) * qJD(4), 0, 0, 0, 0, 0, t119 - t166, -t74 ^ 2 * t109 - t155 - t165 (t33 + t167) * t109 + t175 + t159, t128 * t107 + t173 * t109 - t15 * t79; 0, 0, 0, 0, 0, 0, 0, 0, -t79 * t78, -t78 ^ 2 + t79 ^ 2, t93 + (-t78 - t136) * qJD(4), 0, 0, -t75 * t79 + t176 - t18, t85 * t140 - t75 * t78, -t33 * t107 + t65 * t134 (-t33 + t167) * t109 - t175 + t159, t74 * t134 + t155 - t165, t119 + t166, -t74 * t79, -pkin(4) * t34 - t11 * t79 - t26 * t63 - t47 * t74 - t129 * t109 + (t174 * t74 + t117) * t107, pkin(4) * t33 + t129 * t107 + t117 * t109 + t12 * t79 + t160 * t74 - t26 * t65, -t173 * t107 + t128 * t109 - t156 * t63 - t168 * t65 + t33 * t89 + t34 * t90, -t2 * t90 + t1 * t89 + t9 * (-pkin(5) * t109 - pkin(4)) + t156 * t7 + t168 * t5 + (pkin(5) * t133 - t26) * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65 * t63, -t60 + t171, t63 * t74 - t33, t114 + (-qJD(5) + t74) * t65, t71, t12 * t74 - t22 * t65 + t113, t11 * t74 + t22 * t63 - t116, pkin(5) * t33 - t170 * t63, t170 * t7 + (-t15 * t65 + t1) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60 - t171, t5 * t65 + t63 * t7 + t9;];
tauc_reg  = t8;
