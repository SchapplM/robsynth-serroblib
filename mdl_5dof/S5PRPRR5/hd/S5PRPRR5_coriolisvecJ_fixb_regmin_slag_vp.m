% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tauc_reg [5x22]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPRR5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:54:41
% EndTime: 2019-12-05 15:54:45
% DurationCPUTime: 0.78s
% Computational Cost: add. (832->144), mult. (2255->212), div. (0->0), fcn. (1764->8), ass. (0->94)
t68 = cos(pkin(9));
t73 = cos(qJ(4));
t110 = t73 * t68;
t98 = qJD(2) * t110;
t67 = sin(pkin(9));
t70 = sin(qJ(4));
t113 = t70 * t67;
t99 = qJD(2) * t113;
t44 = -t98 + t99;
t72 = cos(qJ(5));
t35 = t72 * t44;
t52 = t73 * t67 + t70 * t68;
t46 = qJD(2) * t52;
t69 = sin(qJ(5));
t19 = t69 * t46 + t35;
t66 = qJD(4) + qJD(5);
t115 = t19 * t66;
t102 = qJD(5) * t69;
t57 = qJD(4) * t98;
t40 = -qJD(4) * t99 + t57;
t48 = t52 * qJD(4);
t41 = qJD(2) * t48;
t4 = -qJD(5) * t35 - t46 * t102 + t72 * t40 - t69 * t41;
t129 = t4 + t115;
t85 = t69 * t44 - t72 * t46;
t116 = t85 * t66;
t77 = t85 * qJD(5) - t69 * t40 - t72 * t41;
t128 = t77 - t116;
t74 = cos(qJ(2));
t100 = t74 * qJD(1);
t90 = qJD(3) - t100;
t127 = t85 * t19;
t126 = -t19 ^ 2 + t85 ^ 2;
t71 = sin(qJ(2));
t101 = t71 * qJD(1);
t58 = qJD(2) * qJ(3) + t101;
t92 = pkin(6) * qJD(2) + t58;
t38 = t92 * t67;
t39 = t92 * t68;
t86 = t70 * t38 - t73 * t39;
t13 = -t44 * pkin(7) - t86;
t62 = -t68 * pkin(3) - pkin(2);
t49 = t62 * qJD(2) + t90;
t25 = t44 * pkin(4) + t49;
t53 = (qJD(3) + t100) * qJD(2);
t80 = t52 * t53;
t3 = -t40 * pkin(7) + t86 * qJD(4) - t80;
t125 = t25 * t19 + t13 * t102 + (-t13 * t66 - t3) * t69;
t120 = -t73 * t38 - t70 * t39;
t51 = -t110 + t113;
t123 = t51 * t53;
t2 = -t41 * pkin(7) + t120 * qJD(4) - t123;
t124 = -t69 * t2 + t25 * t85 + t72 * t3;
t107 = pkin(6) + qJ(3);
t54 = t107 * t67;
t55 = t107 * t68;
t84 = t70 * t54 - t73 * t55;
t122 = t84 * qJD(4) - t52 * t90;
t105 = t67 ^ 2 + t68 ^ 2;
t121 = t105 * t58;
t47 = t51 * qJD(4);
t119 = qJD(5) - t66;
t111 = t73 * t54;
t79 = t51 * t74;
t118 = -qJD(1) * t79 + (qJD(3) * t67 + qJD(4) * t55) * t70 - qJD(3) * t110 + qJD(4) * t111;
t117 = pkin(4) * t46;
t112 = t72 * t13;
t75 = qJD(2) ^ 2;
t109 = t75 * t71;
t108 = t75 * t74;
t104 = qJD(2) * pkin(2);
t103 = qJD(2) * t71;
t12 = -t46 * pkin(7) + t120;
t9 = qJD(4) * pkin(4) + t12;
t97 = -pkin(4) * t66 - t9;
t95 = t105 * t74;
t94 = t105 * t53;
t91 = t105 * qJD(3);
t89 = t48 * pkin(7) - qJD(5) * (-t52 * pkin(7) - t70 * t55 - t111) + t118;
t88 = -t47 * pkin(7) + qJD(5) * (-t51 * pkin(7) - t84) - t122;
t23 = t72 * t51 + t69 * t52;
t24 = -t69 * t51 + t72 * t52;
t83 = t48 * pkin(4) - t101;
t42 = t52 * t71;
t63 = qJD(2) * t101;
t56 = t90 - t104;
t43 = t51 * t71;
t33 = t51 * pkin(4) + t62;
t26 = t41 * pkin(4) + t63;
t17 = -t46 * t74 + t71 * t47;
t16 = -qJD(2) * t79 - qJD(4) * t42;
t7 = t24 * qJD(5) - t69 * t47 + t72 * t48;
t6 = -t23 * qJD(5) - t72 * t47 - t69 * t48;
t1 = [0, 0, -t109, -t108, -t68 * t109, t67 * t109, t105 * t108, t71 * t94 + (t56 * t71 + (-t101 + t121) * t74) * qJD(2), 0, 0, 0, 0, 0, t17 * qJD(4) + t44 * t103 - t74 * t41, -t16 * qJD(4) + t46 * t103 - t74 * t40, 0, 0, 0, 0, 0, (-t69 * t16 + t72 * t17 + (t42 * t69 + t43 * t72) * qJD(5)) * t66 + t19 * t103 + t74 * t77, -(t72 * t16 + t69 * t17 + (-t42 * t72 + t43 * t69) * qJD(5)) * t66 - t85 * t103 - t74 * t4; 0, 0, 0, 0, 0, 0, t94 + (-qJD(1) * t95 + t91) * qJD(2), t58 * t91 + qJ(3) * t94 + ((-t56 - t104) * t71 - t58 * t95) * qJD(1), t40 * t52 - t46 * t47, -t40 * t51 - t52 * t41 + t47 * t44 - t46 * t48, -t47 * qJD(4), -t48 * qJD(4), 0, t62 * t41 + t49 * t48 + (qJD(2) * t51 - t44) * t101 + t122 * qJD(4), t118 * qJD(4) + t62 * t40 - t49 * t47, t4 * t24 - t6 * t85, -t6 * t19 - t4 * t23 + t24 * t77 + t7 * t85, t6 * t66, -t7 * t66, 0, t26 * t23 + t25 * t7 - t33 * t77 + (t89 * t69 - t88 * t72) * t66 + t83 * t19, t26 * t24 + t25 * t6 + t33 * t4 + (t88 * t69 + t89 * t72) * t66 - t83 * t85; 0, 0, 0, 0, 0, 0, -t105 * t75, -qJD(2) * t121 + t63, 0, 0, 0, 0, 0, 0.2e1 * t46 * qJD(4), t57 + (-t44 - t99) * qJD(4), 0, 0, 0, 0, 0, -t77 - t116, t4 - t115; 0, 0, 0, 0, 0, 0, 0, 0, t46 * t44, -t44 ^ 2 + t46 ^ 2, t57 + (t44 - t99) * qJD(4), 0, 0, -t49 * t46 - t80, t49 * t44 + t123, -t127, t126, t129, t128, 0, -(-t69 * t12 - t112) * t66 - t19 * t117 + (t97 * t69 - t112) * qJD(5) + t124, t85 * t117 + (t97 * qJD(5) + t12 * t66 - t2) * t72 + t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t127, t126, t129, t128, 0, t119 * (-t69 * t9 - t112) + t124, (-t119 * t9 - t2) * t72 + t125;];
tauc_reg = t1;
