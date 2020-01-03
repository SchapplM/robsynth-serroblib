% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S4RRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRRR3_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR3_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR3_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR3_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:24:37
% EndTime: 2019-12-31 17:24:39
% DurationCPUTime: 0.75s
% Computational Cost: add. (1095->98), mult. (2572->196), div. (0->0), fcn. (2238->6), ass. (0->64)
t81 = qJD(2) + qJD(3);
t45 = cos(qJ(3));
t44 = sin(qJ(2));
t78 = -pkin(6) - pkin(5);
t61 = t78 * t44;
t35 = t45 * t61;
t43 = sin(qJ(3));
t46 = cos(qJ(2));
t79 = t78 * t46;
t74 = t43 * t79;
t22 = t35 + t74;
t80 = -t43 * t61 + t45 * t79;
t33 = t43 * t46 + t45 * t44;
t77 = t33 * pkin(7);
t76 = t44 * pkin(2);
t75 = cos(qJ(4));
t73 = qJD(3) * t43;
t72 = qJD(3) * t45;
t42 = sin(qJ(4));
t71 = qJD(4) * t42;
t70 = t44 * qJD(2);
t69 = t46 * qJD(2);
t68 = t42 * t43 * pkin(2);
t67 = -0.2e1 * pkin(1) * qJD(2);
t66 = t22 * qJD(2) + qJD(3) * t35;
t65 = pkin(2) * t70;
t64 = pkin(2) * t73;
t63 = pkin(2) * t72;
t62 = pkin(3) * t71;
t59 = t44 * t69;
t41 = -t46 * pkin(2) - pkin(1);
t58 = t75 * t43;
t57 = t75 * qJD(4);
t56 = pkin(3) * t57;
t32 = t43 * t44 - t45 * t46;
t55 = -t22 + t77;
t20 = -t42 * t32 + t75 * t33;
t54 = t41 * t33;
t53 = t33 * qJD(2);
t52 = qJD(3) * t33;
t40 = t45 * pkin(2) + pkin(3);
t14 = -t40 * t57 - t75 * t63 + (qJD(3) + qJD(4)) * t68;
t51 = t75 * t55;
t12 = -t32 * pkin(7) - t80;
t8 = t75 * t12 - t42 * t55;
t50 = -t52 - t53;
t49 = (-t43 * t57 + (-t42 * t45 - t58) * qJD(3)) * pkin(2);
t10 = t81 * t80;
t48 = -pkin(7) * t53 + (t74 - t77) * qJD(3) + t66;
t21 = t81 * t32;
t47 = t21 * pkin(7) + t10;
t27 = pkin(2) * t58 + t42 * t40;
t26 = t75 * t40 - t68;
t24 = t32 * pkin(3) + t41;
t19 = t75 * t32 + t42 * t33;
t15 = -t40 * t71 + t49;
t13 = pkin(3) * t52 + (t33 * pkin(3) + t76) * qJD(2);
t9 = -t73 * t79 - t66;
t7 = -t42 * t12 - t51;
t4 = qJD(4) * t20 - t42 * t21 - t75 * t50;
t3 = t75 * t21 + t32 * t57 + t33 * t71 - t42 * t50;
t2 = -t8 * qJD(4) - t42 * t48 + t75 * t47;
t1 = qJD(4) * t51 + t12 * t71 - t42 * t47 - t75 * t48;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t59, 0.2e1 * (-t44 ^ 2 + t46 ^ 2) * qJD(2), 0, -0.2e1 * t59, 0, 0, t44 * t67, t46 * t67, 0, 0, -0.2e1 * t33 * t21, 0.2e1 * t21 * t32 + 0.2e1 * t33 * t50, 0, -0.2e1 * t32 * t50, 0, 0, 0.2e1 * qJD(3) * t54 + 0.2e1 * (t32 * t76 + t54) * qJD(2), -0.2e1 * t41 * t21 + 0.2e1 * t33 * t65, -0.2e1 * t10 * t33 + 0.2e1 * t22 * t21 + 0.2e1 * t9 * t32 - 0.2e1 * t50 * t80, 0.2e1 * t22 * t10 + 0.2e1 * t41 * t65 + 0.2e1 * t80 * t9, -0.2e1 * t20 * t3, 0.2e1 * t3 * t19 - 0.2e1 * t20 * t4, 0, 0.2e1 * t19 * t4, 0, 0, 0.2e1 * t13 * t19 + 0.2e1 * t24 * t4, 0.2e1 * t13 * t20 - 0.2e1 * t24 * t3, 0.2e1 * t1 * t19 - 0.2e1 * t2 * t20 + 0.2e1 * t7 * t3 - 0.2e1 * t8 * t4, -0.2e1 * t8 * t1 + 0.2e1 * t24 * t13 + 0.2e1 * t7 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, 0, -t70, 0, -pkin(5) * t69, pkin(5) * t70, 0, 0, 0, 0, -t21, 0, t50, 0, t10, t9, (t45 * t21 - t32 * t72 - t43 * t53) * pkin(2), (t10 * t45 - t43 * t9 + (-t22 * t43 - t45 * t80) * qJD(3)) * pkin(2), 0, 0, -t3, 0, -t4, 0, t2, t1, t14 * t19 - t15 * t20 + t26 * t3 - t27 * t4, -t1 * t27 - t8 * t14 + t7 * t15 + t2 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t64, -0.2e1 * t63, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t15, 0.2e1 * t14, 0, -0.2e1 * t27 * t14 + 0.2e1 * t26 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0, t50, 0, t10, t9, 0, 0, 0, 0, -t3, 0, -t4, 0, t2, t1, (t75 * t3 - t4 * t42 + (-t75 * t19 + t20 * t42) * qJD(4)) * pkin(3), (t75 * t2 - t1 * t42 + (-t42 * t7 + t75 * t8) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, -t63, 0, 0, 0, 0, 0, 0, 0, 0, (-pkin(3) - t40) * t71 + t49, -t56 + t14, 0, (t75 * t15 - t14 * t42 + (-t26 * t42 + t75 * t27) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t62, -0.2e1 * t56, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, -t4, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, t14, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, -t56, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;
