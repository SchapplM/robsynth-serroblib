% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tauc_reg [5x20]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPR7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:19:51
% EndTime: 2019-12-31 18:19:54
% DurationCPUTime: 0.76s
% Computational Cost: add. (1181->162), mult. (2908->239), div. (0->0), fcn. (2016->8), ass. (0->99)
t61 = sin(pkin(8)) * pkin(1) + pkin(6);
t111 = qJ(4) + t61;
t74 = sin(qJ(3));
t76 = cos(qJ(3));
t92 = t111 * qJD(1);
t36 = t76 * qJD(2) - t92 * t74;
t103 = qJD(1) * qJD(4);
t37 = t74 * qJD(2) + t92 * t76;
t130 = -t37 * qJD(3) - t74 * t103;
t107 = t74 * qJD(1);
t110 = qJD(1) * t76;
t71 = cos(pkin(9));
t57 = t71 * t110;
t69 = sin(pkin(9));
t46 = -t69 * t107 + t57;
t44 = qJD(5) - t46;
t129 = -qJD(5) + t44;
t27 = t36 * qJD(3) + t76 * t103;
t3 = -t71 * t130 + t69 * t27;
t53 = t69 * t76 + t71 * t74;
t48 = t53 * qJD(1);
t60 = t69 * pkin(3) + pkin(7);
t128 = (pkin(3) * t107 + t48 * pkin(4) - t46 * pkin(7) + qJD(5) * t60) * t44 + t3;
t47 = t53 * qJD(3);
t41 = qJD(1) * t47;
t75 = cos(qJ(5));
t39 = t75 * t41;
t73 = sin(qJ(5));
t108 = qJD(5) * t73;
t52 = t69 * t74 - t71 * t76;
t50 = t52 * qJD(3);
t85 = t53 * t108 + t75 * t50;
t127 = t53 * t39 - t44 * t85;
t35 = t73 * qJD(3) + t75 * t48;
t104 = qJD(1) * qJD(3);
t99 = t74 * t104;
t42 = qJD(3) * t57 - t69 * t99;
t17 = t35 * qJD(5) + t73 * t42;
t93 = qJD(3) * t111;
t38 = t76 * qJD(4) - t74 * t93;
t82 = -t74 * qJD(4) - t76 * t93;
t13 = t71 * t38 + t69 * t82;
t63 = -cos(pkin(8)) * pkin(1) - pkin(2);
t88 = -t76 * pkin(3) + t63;
t83 = t88 * qJD(1);
t45 = qJD(4) + t83;
t14 = -t46 * pkin(4) - t48 * pkin(7) + t45;
t22 = t52 * pkin(4) - t53 * pkin(7) + t88;
t4 = t130 * t69 + t71 * t27;
t118 = t69 * t37;
t112 = qJD(3) * pkin(3);
t31 = t36 + t112;
t8 = t71 * t31 - t118;
t6 = -qJD(3) * pkin(4) - t8;
t51 = t111 * t76;
t98 = t111 * t74;
t25 = t71 * t51 - t69 * t98;
t91 = -t25 * t41 + t3 * t53;
t126 = -(qJD(5) * t22 + t13) * t44 - t6 * t50 - (qJD(5) * t14 + t4) * t52 + t91;
t105 = t75 * qJD(3);
t16 = qJD(5) * t105 - t48 * t108 + t75 * t42;
t125 = t16 * t52 + t35 * t47;
t124 = t16 * t73;
t123 = t22 * t41;
t33 = t73 * t48 - t105;
t122 = t33 * t44;
t121 = t35 * t44;
t120 = t35 * t48;
t119 = t48 * t33;
t29 = t71 * t37;
t117 = t73 * t41;
t77 = qJD(3) ^ 2;
t115 = t77 * t74;
t114 = t77 * t76;
t9 = t69 * t31 + t29;
t113 = t74 ^ 2 - t76 ^ 2;
t55 = qJD(1) * t63;
t109 = qJD(5) * t53;
t101 = t74 * t112;
t96 = t44 * t75;
t7 = qJD(3) * pkin(7) + t9;
t1 = t75 * t14 - t73 * t7;
t2 = t73 * t14 + t75 * t7;
t90 = -t52 * t17 - t47 * t33;
t87 = 0.2e1 * qJD(3) * t55;
t86 = t39 + (t46 * t73 - t108) * t44;
t11 = t71 * t36 - t118;
t81 = -t60 * t41 + (t11 + t6) * t44;
t79 = (-t75 * t109 + t73 * t50) * t44 - t53 * t117;
t78 = qJD(1) ^ 2;
t62 = -t71 * pkin(3) - pkin(4);
t58 = pkin(3) * t99;
t24 = t69 * t51 + t71 * t98;
t20 = t47 * pkin(4) + t50 * pkin(7) + t101;
t18 = t41 * pkin(4) - t42 * pkin(7) + t58;
t15 = t75 * t18;
t12 = t69 * t38 - t71 * t82;
t10 = t69 * t36 + t29;
t5 = [0, 0, 0, 0, 0.2e1 * t76 * t99, -0.2e1 * t113 * t104, t114, -t115, 0, -t61 * t114 + t74 * t87, t61 * t115 + t76 * t87, t12 * t48 + t13 * t46 + t24 * t42 - t4 * t52 - t9 * t47 + t8 * t50 + t91, -t8 * t12 + t9 * t13 + t3 * t24 + t4 * t25 + (t45 + t83) * t101, t16 * t75 * t53 - t35 * t85, -(-t33 * t75 - t35 * t73) * t50 + (-t124 - t17 * t75 + (t33 * t73 - t35 * t75) * qJD(5)) * t53, t125 + t127, t79 + t90, t41 * t52 + t44 * t47, t1 * t47 + t12 * t33 + t15 * t52 + t24 * t17 + (t20 * t44 + t123 + (-t25 * t44 - t7 * t52 + t6 * t53) * qJD(5)) * t75 + t126 * t73, t12 * t35 + t24 * t16 - t2 * t47 + (-(-qJD(5) * t25 + t20) * t44 - t123 - (-qJD(5) * t7 + t18) * t52 - t6 * t109) * t73 + t126 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t115, -t114, -t53 * t41 + t52 * t42 - t50 * t46 + t47 * t48, t3 * t52 + t4 * t53 - t8 * t47 - t9 * t50, 0, 0, 0, 0, 0, t79 - t90, t125 - t127; 0, 0, 0, 0, -t74 * t78 * t76, t113 * t78, 0, 0, 0, -t55 * t107, -t55 * t110, (-t10 + t9) * t48 + (-t11 + t8) * t46 + (-t41 * t69 - t42 * t71) * pkin(3), t8 * t10 - t9 * t11 + (-t45 * t107 - t3 * t71 + t4 * t69) * pkin(3), t35 * t96 + t124, (t16 - t122) * t75 + (-t17 - t121) * t73, t44 * t96 + t117 - t120, t86 + t119, -t44 * t48, -t1 * t48 - t10 * t33 - t128 * t75 + t62 * t17 + t81 * t73, -t10 * t35 + t128 * t73 + t62 * t16 + t2 * t48 + t81 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46 ^ 2 - t48 ^ 2, -t9 * t46 + t8 * t48 + t58, 0, 0, 0, 0, 0, t86 - t119, -t44 ^ 2 * t75 - t117 - t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35 * t33, -t33 ^ 2 + t35 ^ 2, t16 + t122, t121 - t17, t41, t129 * t2 - t6 * t35 - t73 * t4 + t15, t129 * t1 - t73 * t18 + t6 * t33 - t75 * t4;];
tauc_reg = t5;
