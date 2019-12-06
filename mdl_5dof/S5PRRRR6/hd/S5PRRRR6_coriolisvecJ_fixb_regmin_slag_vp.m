% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% tauc_reg [5x21]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRR6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:10:12
% EndTime: 2019-12-05 17:10:16
% DurationCPUTime: 0.80s
% Computational Cost: add. (928->138), mult. (1769->205), div. (0->0), fcn. (1280->8), ass. (0->115)
t80 = sin(qJ(5));
t81 = sin(qJ(4));
t84 = cos(qJ(5));
t85 = cos(qJ(4));
t49 = t80 * t85 + t84 * t81;
t76 = qJD(4) + qJD(5);
t150 = t49 * t76;
t77 = qJD(2) + qJD(3);
t15 = t150 * t77;
t153 = pkin(7) + pkin(8);
t148 = qJD(5) - t76;
t121 = qJD(1) * qJD(2);
t83 = sin(qJ(2));
t126 = qJD(1) * t83;
t152 = qJD(3) * t126 + t83 * t121;
t132 = t84 * t85;
t136 = t80 * t81;
t47 = -t132 + t136;
t94 = t47 * t76;
t82 = sin(qJ(3));
t86 = cos(qJ(3));
t87 = cos(qJ(2));
t48 = t82 * t83 - t86 * t87;
t151 = t48 * t77;
t109 = t87 * t121;
t125 = qJD(3) * t82;
t122 = t87 * qJD(1);
t68 = qJD(2) * pkin(2) + t122;
t19 = t82 * t109 + t68 * t125 + t152 * t86;
t50 = t82 * t87 + t86 * t83;
t105 = pkin(2) * t125 - t50 * qJD(1);
t96 = t105 * t77;
t149 = -t96 - t19;
t128 = t152 * t82;
t18 = (qJD(3) * t68 + t109) * t86 - t128;
t146 = t77 * pkin(3);
t145 = t86 * pkin(2);
t70 = t82 * pkin(2) + pkin(7);
t144 = -pkin(8) - t70;
t124 = qJD(4) * t81;
t115 = pkin(4) * t124;
t13 = t77 * t115 + t19;
t69 = t82 * t126;
t36 = t86 * t68 - t69;
t72 = -t85 * pkin(4) - pkin(3);
t25 = t72 * t77 - t36;
t143 = t13 * t47 + t150 * t25;
t142 = t13 * t49 - t25 * t94;
t29 = t77 * t50;
t141 = t29 * t77;
t37 = t86 * t126 + t82 * t68;
t140 = t37 * t77;
t118 = t77 * t132;
t119 = t77 * t136;
t39 = -t118 + t119;
t41 = t49 * t77;
t139 = t41 * t39;
t138 = t50 * t76;
t88 = qJD(4) ^ 2;
t137 = t70 * t88;
t112 = t153 * t77 + t37;
t22 = t112 * t85;
t134 = t84 * t22;
t131 = t88 * t81;
t123 = qJD(4) * t85;
t34 = -t36 - t146;
t130 = t34 * t123 + t19 * t81;
t129 = t115 + t105;
t127 = t81 ^ 2 - t85 ^ 2;
t120 = pkin(4) * t77 * t81;
t117 = qJD(3) * t145;
t114 = t77 * t123;
t21 = t112 * t81;
t20 = qJD(4) * pkin(4) - t21;
t113 = -pkin(4) * t76 - t20;
t111 = qJD(4) * t153;
t107 = -t34 * t77 - t18;
t106 = qJD(4) * t144;
t104 = -t37 + t115;
t103 = pkin(7) * t88 - t140;
t101 = t50 * t88 + t141;
t100 = qJD(4) * (t36 - t146);
t99 = qJD(4) * t112;
t98 = 0.2e1 * qJD(4) * t151;
t4 = t85 * t18 - t81 * t99;
t5 = -t81 * t18 - t85 * t99;
t97 = -t25 * t41 - t80 * t4 + t84 * t5;
t14 = qJD(5) * t118 + t84 * t114 - t76 * t119;
t43 = t86 * t122 - t69;
t92 = qJD(4) * ((-pkin(3) - t145) * t77 - t117 + t43);
t91 = t25 * t39 + (t148 * t22 - t5) * t80;
t89 = qJD(2) ^ 2;
t75 = t77 ^ 2;
t74 = t85 * pkin(8);
t73 = t88 * t85;
t62 = t85 * pkin(7) + t74;
t61 = t153 * t81;
t57 = t72 - t145;
t54 = 0.2e1 * t81 * t114;
t52 = t85 * t111;
t51 = t81 * t111;
t45 = t85 * t70 + t74;
t44 = t144 * t81;
t38 = -0.2e1 * t127 * t77 * qJD(4);
t33 = t85 * t106 - t81 * t117;
t32 = t81 * t106 + t85 * t117;
t30 = t34 * t124;
t24 = t150 * t76;
t23 = t94 * t76;
t12 = -t39 ^ 2 + t41 ^ 2;
t7 = t41 * t76 - t15;
t6 = t39 * t76 + t14;
t2 = t14 * t49 - t41 * t94;
t1 = -t14 * t47 - t49 * t15 - t150 * t41 + t39 * t94;
t3 = [0, 0, -t89 * t83, -t89 * t87, 0, -t141, t151 * t77, 0, 0, 0, 0, 0, -t101 * t85 + t81 * t98, t101 * t81 + t85 * t98, 0, 0, 0, 0, 0, t138 * t94 + t48 * t15 + t150 * t151 + t29 * t39, t138 * t150 + t48 * t14 - t151 * t94 + t29 * t41; 0, 0, 0, 0, 0, t149, t43 * t77 + (-t109 + (-pkin(2) * t77 - t68) * qJD(3)) * t86 + t128, t54, t38, t73, -t131, 0, t30 + t81 * t92 + (-t137 + t149) * t85, (t137 + t96) * t81 + t85 * t92 + t130, t2, t1, -t23, -t24, 0, (-t80 * t32 + t84 * t33 + (-t44 * t80 - t45 * t84) * qJD(5)) * t76 + t57 * t15 + t43 * t150 + t129 * t39 + t143, -(t84 * t32 + t80 * t33 + (t44 * t84 - t45 * t80) * qJD(5)) * t76 + t57 * t14 - t43 * t94 + t129 * t41 + t142; 0, 0, 0, 0, 0, -t19 + t140, t36 * t77 - t18, t54, t38, t73, -t131, 0, t30 + t81 * t100 + (-t103 - t19) * t85, t100 * t85 + t103 * t81 + t130, t2, t1, -t23, -t24, 0, (t80 * t51 - t84 * t52 + (t61 * t80 - t62 * t84) * qJD(5)) * t76 + t72 * t15 + t104 * t39 + t36 * t150 + t143, -(-t84 * t51 - t80 * t52 + (-t61 * t84 - t62 * t80) * qJD(5)) * t76 + t72 * t14 + t104 * t41 - t36 * t94 + t142; 0, 0, 0, 0, 0, 0, 0, -t81 * t75 * t85, t127 * t75, 0, 0, 0, t107 * t81, t107 * t85, t139, t12, t6, t7, 0, -(t80 * t21 - t134) * t76 - t39 * t120 + (t113 * t80 - t134) * qJD(5) + t97, -t41 * t120 + (qJD(5) * t113 - t21 * t76 - t4) * t84 + t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, t12, t6, t7, 0, t97 + t148 * (-t80 * t20 - t134), (-t148 * t20 - t4) * t84 + t91;];
tauc_reg = t3;
