% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRRP2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:13:09
% EndTime: 2019-12-31 17:13:11
% DurationCPUTime: 0.52s
% Computational Cost: add. (525->119), mult. (1042->178), div. (0->0), fcn. (472->4), ass. (0->90)
t49 = sin(qJ(3));
t47 = t49 ^ 2;
t51 = cos(qJ(3));
t48 = t51 ^ 2;
t85 = t47 + t48;
t46 = qJD(1) + qJD(2);
t50 = sin(qJ(2));
t84 = pkin(1) * qJD(1);
t74 = t50 * t84;
t26 = t46 * pkin(6) + t74;
t68 = qJ(4) * t46 + t26;
t9 = t68 * t49;
t10 = t68 * t51;
t82 = qJD(3) * pkin(3);
t6 = -t9 + t82;
t99 = t6 + t9;
t98 = pkin(3) * t47;
t97 = t51 * pkin(3);
t52 = cos(qJ(2));
t96 = t52 * pkin(1);
t41 = -pkin(2) - t97;
t73 = t52 * t84;
t11 = t41 * t46 + qJD(4) - t73;
t83 = pkin(1) * qJD(2);
t69 = qJD(1) * t83;
t38 = t50 * t69;
t80 = qJD(3) * t49;
t72 = t46 * t80;
t16 = pkin(3) * t72 + t38;
t79 = qJD(3) * t51;
t95 = t11 * t79 + t16 * t49;
t45 = t46 ^ 2;
t94 = t45 * t51;
t93 = t46 * t49;
t92 = t46 * t51;
t53 = qJD(3) ^ 2;
t91 = t53 * t49;
t43 = t53 * t51;
t90 = -qJ(4) - pkin(6);
t27 = -t46 * pkin(2) - t73;
t89 = t27 * t79 + t49 * t38;
t62 = qJD(3) * t73;
t64 = t46 * t74;
t88 = t49 * t62 + t51 * t64;
t61 = t52 * t69;
t87 = t85 * t61;
t86 = t47 - t48;
t39 = t50 * pkin(1) + pkin(6);
t81 = -qJ(4) - t39;
t78 = qJD(3) * t52;
t77 = -qJD(1) - t46;
t76 = t52 * t83;
t75 = t50 * t83;
t7 = t11 * t80;
t71 = t46 * t79;
t55 = qJD(4) * t46 + t61;
t58 = qJD(3) * t68;
t2 = -t49 * t58 + t51 * t55;
t3 = -t49 * t55 - t51 * t58;
t70 = t2 * t51 - t3 * t49;
t67 = qJD(3) * t90;
t66 = t85 * qJD(2);
t65 = qJD(3) * t81;
t63 = t49 * t71;
t60 = -t10 * t51 + t49 * t6;
t59 = -t10 * t49 - t51 * t6;
t57 = t85 * t73;
t56 = -t27 * t46 - t61;
t54 = -t61 + (-qJD(4) - t11) * t46;
t44 = t51 * qJ(4);
t42 = t51 * qJD(4);
t40 = -pkin(2) - t96;
t36 = t49 * t94;
t35 = t51 * pkin(6) + t44;
t34 = t90 * t49;
t32 = t51 * t62;
t28 = t41 - t96;
t25 = -0.2e1 * t63;
t24 = 0.2e1 * t63;
t22 = pkin(3) * t80 + t75;
t21 = t51 * t39 + t44;
t20 = t86 * t45;
t19 = t81 * t49;
t17 = t27 * t80;
t15 = -t49 * qJD(4) + t51 * t67;
t14 = t49 * t67 + t42;
t13 = -0.2e1 * t86 * t46 * qJD(3);
t5 = (-qJD(4) - t76) * t49 + t51 * t65;
t4 = t49 * t65 + t51 * t76 + t42;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46 * t75 - t38, t77 * t76, 0, 0, t24, t13, t43, t25, -t91, 0, t40 * t72 - t39 * t43 + t17 + (t50 * t51 * t77 - t49 * t78) * t83, t40 * t71 + t39 * t91 + (t50 * t93 - t51 * t78) * t83 + t89, t46 * t66 * t96 + t87, ((qJD(1) * t40 + t27) * t50 + (qJD(1) * t39 + t26) * t52 * t85) * t83, t24, t13, t43, t25, -t91, 0, t7 + (-t22 * t46 - t16) * t51 + (t28 * t93 + t5) * qJD(3), t22 * t93 + (t28 * t92 - t4) * qJD(3) + t95, (t4 * t51 - t49 * t5) * t46 + ((-t19 * t51 - t21 * t49) * t46 + t59) * qJD(3) + t70, t10 * t4 + t11 * t22 + t16 * t28 + t3 * t19 + t2 * t21 + t6 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38 + t64, (-qJD(2) + t46) * t73, 0, 0, t24, t13, t43, t25, -t91, 0, -pkin(2) * t72 + t17 + (-pkin(6) * t53 - t38) * t51 + t88, pkin(6) * t91 + t32 + (-pkin(2) * t79 - t49 * t74) * t46 + t89, -t46 * t57 + t87, ((-pkin(2) * qJD(2) - t27) * t50 + (pkin(6) * t66 - t26 * t85) * t52) * t84, t24, t13, t43, t25, -t91, 0, -t16 * t51 + t7 + (t15 + (t41 - t97) * t93) * qJD(3) + t88, -t49 * t64 + t32 + (-t14 + (t41 * t51 + t98) * t46) * qJD(3) + t95, t59 * qJD(3) + (t14 * t51 - t15 * t49 + (-t34 * t51 - t35 * t49) * qJD(3) - t57) * t46 + t70, pkin(3) * t7 + t10 * t14 + t6 * t15 + t16 * t41 + t2 * t35 + t3 * t34 + (-t11 * t50 + t52 * t60) * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, t20, 0, t36, 0, 0, t56 * t49, t56 * t51, 0, 0, -t36, t20, 0, t36, 0, 0, (pkin(3) * t94 + t54) * t49, -t45 * t98 + t54 * t51, (-t82 + t99) * t92, t99 * t10 + (-t11 * t93 + t3) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t72, 0.2e1 * t71, -t85 * t45, t46 * t60 + t16;];
tauc_reg = t1;
