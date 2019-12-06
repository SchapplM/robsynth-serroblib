% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPPR2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:24:42
% EndTime: 2019-12-05 15:24:45
% DurationCPUTime: 0.51s
% Computational Cost: add. (883->119), mult. (2195->170), div. (0->0), fcn. (1707->8), ass. (0->84)
t66 = sin(pkin(8));
t70 = sin(qJ(2));
t89 = qJD(1) * t70;
t56 = t66 * t89;
t68 = cos(pkin(8));
t72 = cos(qJ(2));
t84 = t72 * qJD(1);
t83 = t68 * t84;
t40 = -t56 + t83;
t100 = t40 - qJD(4);
t53 = qJD(2) * t83;
t26 = t53 + (qJD(4) - t56) * qJD(2);
t67 = cos(pkin(9));
t71 = cos(qJ(5));
t65 = sin(pkin(9));
t69 = sin(qJ(5));
t93 = t69 * t65;
t48 = -t71 * t67 + t93;
t101 = t26 * t48;
t49 = t66 * t72 + t68 * t70;
t12 = t48 * t49;
t36 = t49 * qJD(1);
t50 = t71 * t65 + t69 * t67;
t41 = t50 * qJD(2);
t99 = t41 ^ 2;
t59 = t66 * pkin(2) + qJ(4);
t96 = pkin(6) + t59;
t45 = t96 * t65;
t46 = t96 * t67;
t16 = -t69 * t45 + t71 * t46;
t98 = -t16 * qJD(5) + t100 * t50;
t97 = t68 * pkin(2);
t15 = -t71 * t45 - t69 * t46;
t95 = -t15 * qJD(5) - t100 * t48;
t35 = t49 * qJD(2);
t29 = qJD(1) * t35;
t47 = t66 * t70 - t68 * t72;
t19 = t29 * t47;
t88 = qJD(2) * t67;
t81 = t71 * t88;
t82 = qJD(2) * t93;
t37 = -t81 + t82;
t94 = t41 * t37;
t44 = t50 * qJD(5);
t32 = qJD(2) * t44;
t43 = t48 * qJD(5);
t91 = -t50 * t32 + t43 * t37;
t55 = qJD(2) * pkin(2) + t84;
t28 = t66 * t55 + t68 * t89;
t25 = qJD(2) * qJ(4) + t28;
t18 = t65 * qJD(3) + t67 * t25;
t90 = t65 ^ 2 + t67 ^ 2;
t87 = t35 * qJD(2);
t39 = t47 * qJD(2);
t86 = t39 * qJD(2);
t85 = t43 * qJD(5);
t80 = -t67 * pkin(4) - pkin(3);
t79 = t90 * t26;
t27 = t68 * t55 - t56;
t78 = t36 * qJD(2) - t29;
t77 = qJD(4) - t27;
t62 = t67 * qJD(3);
t13 = t62 + (-pkin(6) * qJD(2) - t25) * t65;
t14 = pkin(6) * t88 + t18;
t5 = t71 * t13 - t69 * t14;
t6 = t69 * t13 + t71 * t14;
t76 = (-t65 * t25 + t62) * t65 - t18 * t67;
t54 = qJD(5) * t81;
t31 = qJD(5) * t82 - t54;
t75 = -t48 * t31 + t41 * t44;
t74 = t50 * t26;
t11 = t50 * t49;
t73 = qJD(2) ^ 2;
t52 = t80 - t97;
t34 = t37 ^ 2;
t33 = t44 * qJD(5);
t30 = -qJD(2) * t56 + t53;
t24 = -qJD(2) * pkin(3) + t77;
t22 = qJD(2) * t80 + t77;
t4 = qJD(5) * t12 + t39 * t50;
t3 = -qJD(5) * t11 + t39 * t48;
t2 = -t6 * qJD(5) - t74;
t1 = t5 * qJD(5) - t101;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73 * t70, -t73 * t72, 0, 0, 0, 0, 0, 0, 0, 0, -t87, t86, 0, -t27 * t35 - t28 * t39 + t30 * t49 + t19, 0, 0, 0, 0, 0, 0, -t67 * t87, t65 * t87, -t90 * t86, t24 * t35 + t39 * t76 + t49 * t79 + t19, 0, 0, 0, 0, 0, 0, t4 * qJD(5) + t47 * t32 + t35 * t37, -t3 * qJD(5) - t47 * t31 + t35 * t41, -t11 * t31 + t12 * t32 - t3 * t37 - t4 * t41, -t1 * t12 - t2 * t11 + t22 * t35 + t6 * t3 + t5 * t4 + t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53 + (t40 + t56) * qJD(2), 0, t27 * t36 - t28 * t40 + (-t29 * t68 + t30 * t66) * pkin(2), 0, 0, 0, 0, 0, 0, t78 * t67, -t78 * t65, -t100 * qJD(2) * t90 + t79, t29 * (-pkin(3) - t97) - t24 * t36 + t59 * t79 + t100 * t76, -t31 * t50 - t41 * t43, -t75 + t91, -t85, t32 * t48 + t37 * t44, -t33, 0, qJD(5) * t98 + t22 * t44 + t29 * t48 + t52 * t32 - t36 * t37, qJD(5) * t95 - t22 * t43 + t29 * t50 - t52 * t31 - t36 * t41, -t1 * t48 + t15 * t31 - t16 * t32 - t2 * t50 + t37 * t95 - t41 * t98 + t5 * t43 - t6 * t44, t1 * t16 + t2 * t15 - t22 * t36 + t29 * t52 + t5 * t98 - t6 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t85, t75 + t91, t1 * t50 - t2 * t48 - t6 * t43 - t5 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t90 * t73, (t36 + t76) * qJD(2), 0, 0, 0, 0, 0, 0, 0.2e1 * t41 * qJD(5), t54 + (-t37 - t82) * qJD(5), -t34 - t99, t6 * t37 + t5 * t41 + t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, -t34 + t99, t54 + (t37 - t82) * qJD(5), -t94, 0, 0, -t22 * t41 - t74, t22 * t37 + t101, 0, 0;];
tauc_reg = t7;
