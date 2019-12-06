% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRP1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:40:13
% EndTime: 2019-12-05 16:40:18
% DurationCPUTime: 0.71s
% Computational Cost: add. (717->146), mult. (1369->198), div. (0->0), fcn. (661->4), ass. (0->97)
t58 = cos(qJ(4));
t53 = qJD(2) + qJD(3);
t57 = sin(qJ(3));
t96 = pkin(2) * qJD(2);
t82 = t57 * t96;
t32 = t53 * pkin(7) + t82;
t75 = qJ(5) * t53 + t32;
t65 = t75 * t58;
t56 = sin(qJ(4));
t89 = t56 * qJD(1);
t11 = t65 + t89;
t59 = cos(qJ(3));
t95 = pkin(2) * qJD(3);
t77 = qJD(2) * t95;
t69 = t59 * t77;
t113 = -qJD(1) * qJD(4) - t69;
t49 = t58 * qJD(1);
t15 = -t56 * t32 + t49;
t16 = t58 * t32 + t89;
t6 = -t16 * qJD(4) - t56 * t69;
t92 = qJD(4) * t56;
t85 = t113 * t58 + t32 * t92;
t112 = -t85 * t58 + (-t15 * t58 - t16 * t56) * qJD(4) - t6 * t56;
t66 = t15 * t56 - t16 * t58;
t111 = t66 * t59;
t54 = t56 ^ 2;
t110 = pkin(4) * t54;
t109 = t58 * pkin(4);
t108 = t59 * pkin(2);
t10 = -t75 * t56 + t49;
t94 = qJD(4) * pkin(4);
t7 = t10 + t94;
t107 = -t10 + t7;
t52 = t53 ^ 2;
t106 = t52 * t58;
t105 = t53 * t56;
t104 = t53 * t58;
t60 = qJD(4) ^ 2;
t103 = t60 * t56;
t50 = t60 * t58;
t102 = -qJ(5) - pkin(7);
t47 = -pkin(3) - t109;
t81 = t59 * t96;
t14 = t47 * t53 + qJD(5) - t81;
t43 = t57 * t77;
t80 = t53 * t92;
t21 = pkin(4) * t80 + t43;
t91 = qJD(4) * t58;
t101 = t14 * t91 + t21 * t56;
t33 = -t53 * pkin(3) - t81;
t100 = t33 * t91 + t56 * t43;
t70 = qJD(4) * t81;
t72 = t53 * t82;
t99 = t56 * t70 + t58 * t72;
t55 = t58 ^ 2;
t98 = t54 - t55;
t97 = t54 + t55;
t45 = t57 * pkin(2) + pkin(7);
t93 = -qJ(5) - t45;
t90 = qJD(4) * t59;
t48 = t58 * qJD(5);
t88 = -qJD(2) - t53;
t87 = -qJD(5) - t14;
t84 = t59 * t95;
t83 = t57 * t95;
t12 = t14 * t92;
t79 = t53 * t91;
t76 = qJ(5) * t92;
t2 = (-t76 + t48) * t53 - t85;
t3 = (-qJD(5) * t53 - t69) * t56 - t11 * qJD(4);
t78 = t2 * t58 - t3 * t56;
t74 = qJD(4) * t102;
t73 = qJD(4) * t93;
t71 = t56 * t79;
t68 = t11 * t58 - t56 * t7;
t67 = -t11 * t56 - t58 * t7;
t64 = t97 * t81;
t51 = t58 * qJ(5);
t46 = -pkin(3) - t108;
t41 = t56 * t106;
t40 = t58 * pkin(7) + t51;
t39 = t102 * t56;
t36 = t58 * t70;
t34 = t47 - t108;
t31 = -0.2e1 * t71;
t30 = 0.2e1 * t71;
t28 = pkin(4) * t92 + t83;
t27 = t58 * t45 + t51;
t26 = t98 * t52;
t25 = t93 * t56;
t22 = t33 * t92;
t20 = -t56 * qJD(5) + t58 * t74;
t19 = t56 * t74 + t48;
t18 = -0.2e1 * t98 * t53 * qJD(4);
t9 = (-qJD(5) - t84) * t56 + t58 * t73;
t8 = t56 * t73 + t58 * t84 + t48;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t103, -t50, 0, -t66 * qJD(4) - t56 * t85 + t6 * t58, 0, 0, 0, 0, 0, 0, -t103, -t50, 0, t68 * qJD(4) + t2 * t56 + t3 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53 * t83 - t43, t88 * t84, 0, 0, t30, t18, t50, t31, -t103, 0, t46 * t80 - t45 * t50 + t22 + (t88 * t58 * t57 - t56 * t90) * t95, t46 * t79 + t45 * t103 + (t57 * t105 - t58 * t90) * t95 + t100, t97 * t53 * t84 + t112, t112 * t45 + (-t111 + (qJD(2) * t46 + t33) * t57) * t95, t30, t18, t50, t31, -t103, 0, t12 + (-t28 * t53 - t21) * t58 + (t34 * t105 + t9) * qJD(4), t28 * t105 + (t34 * t104 - t8) * qJD(4) + t101, (-t56 * t9 + t58 * t8) * t53 + ((-t25 * t58 - t27 * t56) * t53 + t67) * qJD(4) + t78, t11 * t8 + t14 * t28 + t2 * t27 + t21 * t34 + t3 * t25 + t7 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43 + t72, (-qJD(3) + t53) * t81, 0, 0, t30, t18, t50, t31, -t103, 0, -pkin(3) * t80 + t22 + (-pkin(7) * t60 - t43) * t58 + t99, pkin(7) * t103 + t36 + (-pkin(3) * t91 - t56 * t82) * t53 + t100, -t53 * t64 + t112, t112 * pkin(7) + (t111 + (-pkin(3) * qJD(3) - t33) * t57) * t96, t30, t18, t50, t31, -t103, 0, -t21 * t58 + t12 + (t20 + (t47 - t109) * t105) * qJD(4) + t99, -t56 * t72 + t36 + (-t19 + (t47 * t58 + t110) * t53) * qJD(4) + t101, t67 * qJD(4) + (t19 * t58 - t20 * t56 + (-t39 * t58 - t40 * t56) * qJD(4) - t64) * t53 + t78, pkin(4) * t12 + t11 * t19 + t2 * t40 + t7 * t20 + t21 * t47 + t3 * t39 + (-t14 * t57 - t68 * t59) * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, t26, 0, t41, 0, 0, (-t33 * t53 - t69) * t56, t15 * qJD(4) - t33 * t104 + t85, 0, 0, -t41, t26, 0, t41, 0, 0, (t11 - t65) * qJD(4) + (pkin(4) * t106 + t87 * t53 + t113) * t56, -t52 * t110 + t10 * qJD(4) + (t87 * t58 + t76) * t53 + t85, (-t94 + t107) * t104, t107 * t11 + (-t14 * t105 + t3) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t80, 0.2e1 * t79, -t97 * t52, -t68 * t53 + t21;];
tauc_reg = t1;
