% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPR10_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR10_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:26:05
% EndTime: 2019-12-31 18:26:07
% DurationCPUTime: 0.56s
% Computational Cost: add. (1538->103), mult. (2550->157), div. (0->0), fcn. (1291->6), ass. (0->87)
t55 = sin(pkin(8));
t56 = cos(pkin(8));
t58 = sin(qJ(3));
t60 = cos(qJ(3));
t38 = t55 * t58 - t56 * t60;
t84 = qJD(1) - qJD(3);
t111 = t38 * t84;
t116 = t111 * t84;
t115 = t84 ^ 2;
t39 = t55 * t60 + t56 * t58;
t62 = qJD(5) ^ 2;
t95 = t84 * t39;
t78 = t95 * t84;
t114 = t39 * t62 + t78;
t61 = -pkin(1) - pkin(2);
t46 = t61 * qJD(1) + qJD(2);
t86 = qJD(1) * qJ(2);
t79 = t58 * t86;
t85 = qJD(1) * qJD(2);
t89 = qJD(3) * t60;
t20 = -qJD(3) * t79 + t46 * t89 + t60 * t85;
t29 = t60 * t46 - t79;
t113 = t29 * t84 + t20;
t87 = t58 * qJD(2);
t98 = t58 * t46;
t21 = -(qJ(2) * t89 + t87) * qJD(1) - qJD(3) * t98;
t30 = t60 * t86 + t98;
t112 = -t30 * t84 + t21;
t25 = t56 * t30;
t13 = t55 * t29 + t25;
t5 = t55 * t20 - t56 * t21;
t110 = -t13 * t84 - t5;
t77 = -t58 * qJ(2) + t60 * t61;
t27 = t60 * qJD(2) + t77 * qJD(3);
t45 = t60 * qJ(2) + t58 * t61;
t28 = -t45 * qJD(3) - t87;
t11 = t55 * t27 - t56 * t28;
t109 = t11 * t84 + t5;
t57 = sin(qJ(5));
t59 = cos(qJ(5));
t23 = -pkin(3) * t84 + t29;
t10 = t55 * t23 + t25;
t8 = -pkin(7) * t84 + t10;
t3 = t59 * qJD(4) - t57 * t8;
t4 = t57 * qJD(4) + t59 * t8;
t72 = t3 * t57 - t4 * t59;
t6 = t56 * t20 + t55 * t21;
t1 = t3 * qJD(5) + t59 * t6;
t2 = -t4 * qJD(5) - t57 * t6;
t64 = -(t3 * t59 + t4 * t57) * qJD(5) + t1 * t59 - t2 * t57;
t106 = t5 * t38;
t12 = t56 * t27 + t55 * t28;
t104 = t12 * t84;
t100 = t55 * t30;
t14 = t56 * t29 - t100;
t102 = t14 * t84;
t99 = t57 * t59;
t97 = t62 * t57;
t96 = t62 * t59;
t41 = -pkin(3) + t77;
t93 = t55 * t41 + t56 * t45;
t53 = t57 ^ 2;
t54 = t59 ^ 2;
t92 = t53 - t54;
t91 = t53 + t54;
t88 = qJD(5) * t84;
t83 = t115 * t99;
t81 = 0.2e1 * t85;
t9 = t56 * t23 - t100;
t7 = pkin(4) * t84 - t9;
t80 = t7 * t84 - t6;
t75 = t88 * t99;
t71 = t56 * t41 - t55 * t45;
t16 = -pkin(7) + t93;
t70 = t16 * t62 - t109;
t47 = t55 * pkin(3) + pkin(7);
t69 = -t47 * t62 + t110;
t15 = pkin(4) - t71;
t68 = qJD(5) * (-t15 * t84 - t12 - t7);
t48 = -t56 * pkin(3) - pkin(4);
t67 = qJD(5) * (-t48 * t84 + t14 + t7);
t66 = -0.2e1 * qJD(5) * t111;
t63 = qJD(1) ^ 2;
t43 = 0.2e1 * t75;
t42 = -0.2e1 * t75;
t31 = t92 * t88;
t17 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, qJ(2) * t81, 0, 0, 0, 0, 0, 0, -t28 * t84 - t21, t27 * t84 + t20, 0, t20 * t45 + t21 * t77 + t30 * t27 + t29 * t28, 0, 0, 0, 0, 0, 0, t109, t6 + t104, 0, t10 * t12 - t9 * t11 - t5 * t71 + t6 * t93, t43, -0.2e1 * t31, -t96, t42, t97, 0, t57 * t68 - t59 * t70, t70 * t57 + t59 * t68, -t91 * t104 - t64, t7 * t11 - t12 * t72 + t5 * t15 + t16 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, -t63 * qJ(2), 0, 0, 0, 0, 0, 0, -t58 * t115, -t60 * t115, 0, t112 * t60 + t113 * t58, 0, 0, 0, 0, 0, 0, -t78, t116, 0, t10 * t111 + t6 * t39 + t9 * t95 + t106, 0, 0, 0, 0, 0, 0, -t114 * t59 + t57 * t66, t114 * t57 + t59 * t66, -t91 * t116, -t111 * t72 + t39 * t64 - t7 * t95 + t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, -t113, 0, 0, 0, 0, 0, 0, 0, 0, t110, -t6 - t102, 0, -t10 * t14 + t9 * t13 + (-t5 * t56 + t55 * t6) * pkin(3), t42, 0.2e1 * t31, t96, t43, -t97, 0, t57 * t67 + t59 * t69, -t69 * t57 + t59 * t67, t91 * t102 + t64, -t7 * t13 + t14 * t72 + t47 * t64 + t5 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97, -t96, 0, -qJD(5) * t72 + t1 * t57 + t2 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83, t92 * t115, 0, t83, 0, 0, t80 * t57, t80 * t59, 0, 0;];
tauc_reg = t17;
