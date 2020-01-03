% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPPR3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:26:36
% EndTime: 2019-12-31 19:26:37
% DurationCPUTime: 0.37s
% Computational Cost: add. (751->101), mult. (1445->131), div. (0->0), fcn. (769->6), ass. (0->84)
t43 = sin(pkin(8));
t44 = cos(pkin(8));
t48 = cos(qJ(2));
t80 = pkin(1) * qJD(2);
t70 = qJD(1) * t80;
t64 = t48 * t70;
t46 = sin(qJ(2));
t65 = t46 * t70;
t16 = t43 * t64 + t44 * t65;
t40 = qJD(1) + qJD(2);
t81 = pkin(1) * qJD(1);
t73 = t48 * t81;
t28 = t40 * pkin(2) + t73;
t74 = t46 * t81;
t14 = t43 * t28 + t44 * t74;
t9 = t40 * qJ(4) + t14;
t91 = t9 * t40;
t93 = t16 - t91;
t47 = cos(qJ(5));
t15 = t47 * t16;
t45 = sin(qJ(5));
t33 = t43 * t74;
t13 = t44 * t28 - t33;
t59 = qJD(4) - t13;
t6 = (-pkin(3) - pkin(7)) * t40 + t59;
t5 = t47 * qJD(3) + t45 * t6;
t79 = t5 * qJD(5);
t2 = t15 - t79;
t4 = -t45 * qJD(3) + t47 * t6;
t78 = qJD(5) * t45;
t1 = t4 * qJD(5) + t45 * t16;
t92 = t1 * t45;
t51 = -(t2 + t79) * t47 + t4 * t78 - t92;
t39 = t40 ^ 2;
t17 = -t43 * t65 + t44 * t64;
t12 = t40 * qJD(4) + t17;
t77 = qJD(5) * t47;
t90 = t12 * t45 + t9 * t77;
t87 = t44 * t46;
t53 = pkin(1) * (t43 * t48 + t87);
t22 = qJD(1) * t53;
t89 = t22 * t40;
t23 = qJD(2) * t53;
t88 = t23 * t40;
t49 = qJD(5) ^ 2;
t86 = t49 * t45;
t85 = t49 * t47;
t84 = -t39 - t49;
t41 = t45 ^ 2;
t42 = t47 ^ 2;
t83 = t41 - t42;
t82 = t41 + t42;
t24 = t44 * t73 - t33;
t76 = qJD(4) - t24;
t35 = t43 * t46 * pkin(1);
t75 = t47 * t39 * t45;
t72 = -t44 * pkin(2) - pkin(3);
t38 = t48 * pkin(1) + pkin(2);
t54 = pkin(1) * t87 + t43 * t38;
t21 = qJ(4) + t54;
t68 = t21 * t40 + t23;
t67 = t44 * t38 - t35;
t66 = t40 * t45 * t77;
t63 = -pkin(3) - t67;
t62 = -t16 + t89;
t61 = t16 + t88;
t60 = t4 * t47 + t45 * t5;
t58 = (-qJD(2) + t40) * t81;
t57 = (-qJD(1) - t40) * t80;
t25 = t44 * t48 * t80 - qJD(2) * t35;
t19 = qJD(4) + t25;
t56 = t12 * t21 + t9 * t19;
t20 = -pkin(7) + t63;
t55 = t19 * t40 - t20 * t49;
t37 = t43 * pkin(2) + qJ(4);
t52 = t12 * t37 + t76 * t9;
t50 = t92 + t2 * t47 + (-t4 * t45 + t47 * t5) * qJD(5);
t36 = -pkin(7) + t72;
t27 = -0.2e1 * t66;
t26 = 0.2e1 * t66;
t18 = 0.2e1 * t83 * t40 * qJD(5);
t11 = t12 * t47;
t8 = -t40 * pkin(3) + t59;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46 * t57, t48 * t57, 0, 0, 0, 0, 0, 0, 0, 0, -t61, -t25 * t40 - t17, 0, -t13 * t23 + t14 * t25 - t16 * t67 + t17 * t54, 0, 0, 0, 0, 0, 0, 0, t61, (qJD(4) + t19) * t40 + t17, t16 * t63 + t8 * t23 + t56, t27, t18, -t86, t26, -t85, 0, t55 * t45 + t68 * t77 + t90, t11 + t55 * t47 + (-t68 - t9) * t78, -t82 * t88 + t51, t50 * t20 + t60 * t23 + t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46 * t58, t48 * t58, 0, 0, 0, 0, 0, 0, 0, 0, t62, t24 * t40 - t17, 0, t13 * t22 - t14 * t24 + (-t16 * t44 + t17 * t43) * pkin(2), 0, 0, 0, 0, 0, 0, 0, -t62, (0.2e1 * qJD(4) - t24) * t40 + t17, t16 * t72 - t8 * t22 + t52, t27, t18, -t86, t26, -t85, 0, -t22 * t77 - t36 * t86 + (t37 * t77 + t76 * t45) * t40 + t90, t11 + (-t36 * t49 + t76 * t40) * t47 + (-t37 * t40 + t22 - t9) * t78, t82 * t89 + t51, -t60 * t22 + t50 * t36 + t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85, t86, 0, -t60 * qJD(5) + t1 * t47 - t2 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, t93, 0, 0, 0, 0, 0, 0, t84 * t45, t84 * t47, 0, -t51 - t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, -t83 * t39, 0, -t75, 0, 0, -t47 * t91 + t15, -t93 * t45, 0, 0;];
tauc_reg = t3;
