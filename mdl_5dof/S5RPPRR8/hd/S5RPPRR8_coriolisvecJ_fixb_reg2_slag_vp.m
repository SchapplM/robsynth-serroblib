% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPRR8_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR8_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:01:14
% EndTime: 2019-12-31 18:01:16
% DurationCPUTime: 0.55s
% Computational Cost: add. (1419->89), mult. (2403->140), div. (0->0), fcn. (1296->6), ass. (0->79)
t40 = sin(pkin(8));
t41 = cos(pkin(8));
t43 = sin(qJ(4));
t45 = cos(qJ(4));
t26 = t43 * t40 - t45 * t41;
t94 = qJD(1) - qJD(4);
t92 = t94 * t26;
t96 = t92 * t94;
t27 = t45 * t40 + t43 * t41;
t47 = qJD(5) ^ 2;
t78 = t94 * t27;
t66 = t78 * t94;
t95 = -t27 * t47 - t66;
t51 = t27 * qJD(2);
t46 = -pkin(1) - pkin(2);
t73 = t40 * qJ(2);
t63 = t41 * t46 - t73;
t28 = -pkin(3) + t63;
t29 = t41 * qJ(2) + t40 * t46;
t76 = t43 * t28 + t45 * t29;
t10 = t76 * qJD(4) + t51;
t34 = t46 * qJD(1) + qJD(2);
t30 = t41 * t34;
t16 = t30 + (-pkin(3) - t73) * qJD(1);
t71 = qJD(1) * qJ(2);
t18 = t40 * t34 + t41 * t71;
t12 = t43 * t16 + t45 * t18;
t91 = t12 * qJD(4);
t4 = qJD(1) * t51 + t91;
t93 = t10 * t94 + t4;
t70 = qJD(1) * qJD(2);
t64 = t41 * t70;
t65 = t40 * t70;
t82 = t45 * t16;
t3 = -(qJD(4) * t18 + t65) * t43 + qJD(4) * t82 + t45 * t64;
t44 = cos(qJ(5));
t42 = sin(qJ(5));
t8 = -pkin(7) * t94 + t12;
t5 = t44 * qJD(3) - t42 * t8;
t1 = t5 * qJD(5) + t44 * t3;
t6 = t42 * qJD(3) + t44 * t8;
t2 = -t6 * qJD(5) - t42 * t3;
t49 = -(t42 * t6 + t44 * t5) * qJD(5) + t1 * t44 - t2 * t42;
t90 = t94 * pkin(4);
t89 = t4 * t26;
t57 = t45 * t28 - t43 * t29;
t9 = -t26 * qJD(2) + t57 * qJD(4);
t88 = t9 * t94;
t11 = -t43 * t18 + t82;
t86 = t11 * t94;
t85 = t12 * t94;
t83 = t42 * t44;
t81 = t47 * t42;
t80 = t47 * t44;
t38 = t42 ^ 2;
t39 = t44 ^ 2;
t75 = t38 - t39;
t74 = t38 + t39;
t72 = qJD(5) * t94;
t36 = t94 ^ 2;
t69 = t36 * t83;
t68 = 0.2e1 * t70;
t7 = -t11 + t90;
t67 = t7 * t94 - t3;
t62 = t72 * t83;
t59 = t42 * t5 - t44 * t6;
t58 = (-t40 * t71 + t30) * t40 - t18 * t41;
t56 = pkin(7) * t47 + t4 + t85;
t14 = -pkin(7) + t76;
t55 = t14 * t47 - t93;
t54 = qJD(5) * (t11 + t7 + t90);
t13 = pkin(4) - t57;
t53 = qJD(5) * (-t13 * t94 - t7 - t9);
t52 = -0.2e1 * qJD(5) * t92;
t48 = qJD(1) ^ 2;
t32 = 0.2e1 * t62;
t31 = -0.2e1 * t62;
t19 = t75 * t72;
t15 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, qJ(2) * t68, 0, 0, 0, 0, 0, 0, 0.2e1 * t65, 0.2e1 * t64, 0, ((t41 * t29 - t40 * t63) * qJD(1) - t58) * qJD(2), 0, 0, 0, 0, 0, 0, t93, t3 + t88, 0, -t11 * t10 + t12 * t9 + t3 * t76 - t4 * t57, t32, -0.2e1 * t19, -t80, t31, t81, 0, t42 * t53 - t55 * t44, t55 * t42 + t44 * t53, -t74 * t88 - t49, t7 * t10 + t4 * t13 + t49 * t14 - t59 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, -t48 * qJ(2), 0, 0, 0, 0, 0, 0, -t40 * t48, -t41 * t48, 0, t58 * qJD(1), 0, 0, 0, 0, 0, 0, -t66, t96, 0, t78 * t11 + t12 * t92 + t3 * t27 + t89, 0, 0, 0, 0, 0, 0, t42 * t52 + t95 * t44, -t95 * t42 + t44 * t52, -t74 * t96, t49 * t27 - t92 * t59 - t78 * t7 + t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, -t80, 0, -t59 * qJD(5) + t1 * t42 + t2 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27 * t70 - t85 - t91, -t3 - t86, 0, 0, t31, 0.2e1 * t19, t80, t32, -t81, 0, t42 * t54 - t56 * t44, t56 * t42 + t44 * t54, t74 * t86 + t49, -t4 * pkin(4) + t49 * pkin(7) + t59 * t11 - t7 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, t75 * t36, 0, t69, 0, 0, t67 * t42, t67 * t44, 0, 0;];
tauc_reg = t15;
