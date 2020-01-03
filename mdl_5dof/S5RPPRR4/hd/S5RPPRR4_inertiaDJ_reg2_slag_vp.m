% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPRR4_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:31:23
% EndTime: 2020-01-03 11:31:28
% DurationCPUTime: 0.90s
% Computational Cost: add. (1679->121), mult. (3922->245), div. (0->0), fcn. (3870->8), ass. (0->76)
t62 = sin(pkin(8));
t63 = cos(pkin(9));
t96 = cos(qJ(4));
t80 = t96 * t63;
t76 = qJD(4) * t80;
t51 = t62 * t76;
t61 = sin(pkin(9));
t66 = sin(qJ(4));
t87 = qJD(4) * t66;
t81 = t61 * t87;
t74 = t62 * t81 - t51;
t99 = -0.2e1 * t74;
t47 = -t66 * t61 + t80;
t48 = t96 * t61 + t66 * t63;
t45 = t48 * qJD(4);
t59 = t62 ^ 2;
t85 = t59 * qJD(2);
t64 = cos(pkin(8));
t98 = 0.2e1 * t64;
t97 = t64 * pkin(4);
t95 = cos(qJ(5));
t94 = t61 * t62;
t50 = -t64 * pkin(2) - t62 * qJ(3) - pkin(1);
t90 = qJ(2) * t64;
t91 = t61 * t50 + t63 * t90;
t28 = -pkin(6) * t94 + t91;
t70 = (-qJ(2) * t61 - pkin(3)) * t64 + (-pkin(6) * t62 + t50) * t63;
t18 = t96 * t28 + t66 * t70;
t73 = t48 * t62;
t14 = -pkin(7) * t73 + t18;
t65 = sin(qJ(5));
t93 = t65 * t14;
t49 = pkin(3) * t94 + t62 * qJ(2);
t89 = qJD(2) * t64;
t88 = qJD(3) * t62;
t86 = qJD(5) * t65;
t57 = t62 * qJD(2);
t84 = qJ(2) * qJD(2);
t23 = t96 * t70;
t17 = -t66 * t28 + t23;
t40 = t47 * t62;
t13 = -t40 * pkin(7) + t17 - t97;
t35 = t62 * t45;
t42 = -t61 * t89 - t63 * t88;
t43 = -t61 * t88 + t63 * t89;
t9 = -t18 * qJD(4) + t96 * t42 - t66 * t43;
t67 = t35 * pkin(7) + t9;
t8 = -qJD(4) * t23 + t28 * t87 - t66 * t42 - t96 * t43;
t69 = t74 * pkin(7) - t8;
t78 = qJD(5) * t95;
t83 = -t13 * t78 - t65 * t67 - t95 * t69;
t82 = pkin(4) * t86;
t79 = t95 * t14;
t77 = pkin(4) * t78;
t75 = t42 * t63 + t43 * t61;
t4 = t65 * t13 + t79;
t27 = t65 * t47 + t95 * t48;
t72 = t65 * t73;
t71 = t95 * t73;
t68 = -t65 * t69 + t95 * t67;
t60 = t64 ^ 2;
t56 = t59 * t84;
t44 = -t76 + t81;
t29 = -t74 * pkin(4) + t57;
t26 = t95 * t47 - t65 * t48;
t25 = pkin(4) * t73 + t49;
t20 = t95 * t40 - t72;
t19 = t65 * t40 + t71;
t16 = -t27 * qJD(5) + t65 * t44 - t95 * t45;
t15 = t95 * t44 + t65 * t45 - t47 * t78 + t48 * t86;
t12 = -qJD(5) * t72 - t65 * t35 + t40 * t78 - t95 * t74;
t11 = qJD(5) * t71 + t95 * t35 + t40 * t86 - t65 * t74;
t3 = t95 * t13 - t93;
t2 = -t4 * qJD(5) + t68;
t1 = t14 * t86 + t83;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * (t59 + t60) * qJD(2), 0.2e1 * t60 * t84 + 0.2e1 * t56, 0, 0, 0, 0, 0, 0, -0.2e1 * t42 * t64 + 0.2e1 * t61 * t85, 0.2e1 * t43 * t64 + 0.2e1 * t63 * t85, -0.2e1 * t75 * t62, 0.2e1 * t91 * t43 + 0.2e1 * (t63 * t50 - t61 * t90) * t42 + 0.2e1 * t56, -0.2e1 * t40 * t35, -0.2e1 * t40 * t51 + 0.2e1 * (t35 * t48 + t40 * t81) * t62, t35 * t98, t73 * t99, t64 * t99, 0, 0.2e1 * t48 * t85 - 0.2e1 * t49 * t74 - 0.2e1 * t9 * t64, -0.2e1 * t49 * t35 + 0.2e1 * t40 * t57 - 0.2e1 * t8 * t64, 0.2e1 * t17 * t35 + 0.2e1 * t18 * t74 - 0.2e1 * t9 * t40 + 0.2e1 * t8 * t73, 0.2e1 * t17 * t9 - 0.2e1 * t18 * t8 + 0.2e1 * t49 * t57, -0.2e1 * t20 * t11, 0.2e1 * t11 * t19 - 0.2e1 * t20 * t12, t11 * t98, 0.2e1 * t19 * t12, t12 * t98, 0, 0.2e1 * t25 * t12 + 0.2e1 * t29 * t19 - 0.2e1 * t2 * t64, -0.2e1 * t1 * t64 - 0.2e1 * t25 * t11 + 0.2e1 * t29 * t20, 0.2e1 * t1 * t19 + 0.2e1 * t3 * t11 - 0.2e1 * t4 * t12 - 0.2e1 * t2 * t20, -0.2e1 * t4 * t1 + 0.2e1 * t3 * t2 + 0.2e1 * t25 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, 0, 0, 0, 0, 0, 0, t45 * t64, -t44 * t64, t47 * t35 + t45 * t40 + t44 * t73 + t48 * t74, -t17 * t45 - t18 * t44 + t9 * t47 - t8 * t48, 0, 0, 0, 0, 0, 0, -t16 * t64, -t15 * t64, t26 * t11 - t27 * t12 + t15 * t19 - t16 * t20, -t1 * t27 - t4 * t15 + t3 * t16 + t2 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t48 * t44 - 0.2e1 * t47 * t45, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t27 * t15 + 0.2e1 * t26 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, 0, 0, 0, 0, 0, -t74, -t35, 0, t57, 0, 0, 0, 0, 0, 0, t12, -t11, 0, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, 0, t74, 0, t9, t8, 0, 0, 0, 0, -t11, 0, -t12, 0, (-t79 + (-t13 + t97) * t65) * qJD(5) + t68, (t95 * t97 + t93) * qJD(5) + t83, (t95 * t11 - t12 * t65 + (-t95 * t19 + t20 * t65) * qJD(5)) * pkin(4), (t95 * t2 - t1 * t65 + (-t3 * t65 + t95 * t4) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, t44, 0, 0, 0, 0, 0, 0, 0, 0, t16, t15, 0, (t95 * t16 - t15 * t65 + (-t26 * t65 + t95 * t27) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t82, -0.2e1 * t77, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, -t12, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t15, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82, -t77, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;
