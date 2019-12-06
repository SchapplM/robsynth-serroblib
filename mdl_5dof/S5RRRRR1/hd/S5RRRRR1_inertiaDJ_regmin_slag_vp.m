% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x31]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRRR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR1_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_inertiaDJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:51:00
% EndTime: 2019-12-05 18:51:04
% DurationCPUTime: 0.78s
% Computational Cost: add. (1197->104), mult. (2896->204), div. (0->0), fcn. (2953->8), ass. (0->79)
t41 = sin(qJ(4));
t42 = sin(qJ(3));
t46 = cos(qJ(3));
t45 = cos(qJ(4));
t80 = qJD(4) * t45;
t85 = t42 * t45;
t93 = ((t41 * t46 + t85) * qJD(3) + t42 * t80) * pkin(2);
t44 = cos(qJ(5));
t39 = t44 ^ 2;
t40 = sin(qJ(5));
t83 = t40 ^ 2 - t39;
t64 = t83 * qJD(5);
t43 = sin(qJ(2));
t47 = cos(qJ(2));
t59 = t42 * t43 - t46 * t47;
t60 = t42 * t47 + t46 * t43;
t16 = t41 * t59 - t45 * t60;
t15 = -t41 * t60 - t45 * t59;
t17 = (qJD(2) + qJD(3)) * t59;
t51 = t60 * qJD(3);
t48 = t60 * qJD(2) + t51;
t6 = -t15 * qJD(4) + t45 * t17 + t41 * t48;
t91 = t16 * t6;
t7 = t16 * qJD(4) + t41 * t17 - t45 * t48;
t90 = t40 * t7;
t89 = t43 * pkin(2);
t88 = t44 * t7;
t35 = t47 * pkin(2) + pkin(1);
t34 = t46 * pkin(2) + pkin(3);
t81 = qJD(4) * t41;
t12 = t34 * t81 + t93;
t76 = t41 * t42 * pkin(2);
t21 = -t45 * t34 - pkin(4) + t76;
t36 = qJD(5) * t44;
t87 = t12 * t40 + t21 * t36;
t86 = t40 * t44;
t33 = -t45 * pkin(3) - pkin(4);
t70 = pkin(3) * t81;
t84 = t33 * t36 + t40 * t70;
t82 = pkin(2) * qJD(3);
t79 = qJD(5) * t40;
t78 = t43 * qJD(2);
t77 = t47 * qJD(2);
t75 = -0.2e1 * pkin(1) * qJD(2);
t74 = pkin(4) * t79;
t73 = pkin(4) * t36;
t72 = t42 * t82;
t71 = t46 * t82;
t69 = pkin(3) * t80;
t67 = t40 * t36;
t66 = -0.4e1 * t16 * t86;
t18 = t21 * t79;
t65 = -t12 * t44 + t18;
t22 = pkin(2) * t85 + t41 * t34 + pkin(6);
t63 = t15 * t22 - t16 * t21;
t32 = t41 * pkin(3) + pkin(6);
t62 = t15 * t32 - t16 * t33;
t26 = t33 * t79;
t58 = -t44 * t70 + t26;
t10 = -pkin(3) * t51 + (-t60 * pkin(3) - t89) * qJD(2);
t2 = t7 * pkin(4) - t6 * pkin(6) + t10;
t20 = -t59 * pkin(3) + t35;
t8 = t15 * pkin(4) - t16 * pkin(6) + t20;
t57 = -t40 * t2 - t8 * t36;
t56 = t44 * t2 - t8 * t79;
t5 = t15 * t36 + t90;
t55 = t15 * t79 - t88;
t54 = -t16 * t36 - t40 * t6;
t53 = t16 * t79 - t44 * t6;
t52 = t35 * t60;
t11 = -t34 * t80 - t45 * t71 + (qJD(3) + qJD(4)) * t76;
t50 = t11 * t15 + t12 * t16 + t21 * t6 - t22 * t7;
t49 = -t32 * t7 + t33 * t6 + (-t15 * t45 + t16 * t41) * qJD(4) * pkin(3);
t31 = 0.2e1 * t67;
t24 = -0.2e1 * t64;
t14 = t16 ^ 2;
t3 = -t16 * t64 + t6 * t86;
t1 = qJD(5) * t66 - t83 * t6;
t4 = [0, 0, 0, 0.2e1 * t43 * t77, 0.2e1 * (-t43 ^ 2 + t47 ^ 2) * qJD(2), 0, 0, 0, t43 * t75, t47 * t75, -0.2e1 * t60 * t17, 0.2e1 * t59 * t17 - 0.2e1 * t48 * t60, 0, 0, 0, -0.2e1 * qJD(3) * t52 + 0.2e1 * (t59 * t89 - t52) * qJD(2), 0.2e1 * pkin(2) * t60 * t78 + 0.2e1 * t35 * t17, 0.2e1 * t91, -0.2e1 * t6 * t15 - 0.2e1 * t16 * t7, 0, 0, 0, 0.2e1 * t10 * t15 + 0.2e1 * t20 * t7, 0.2e1 * t10 * t16 + 0.2e1 * t20 * t6, -0.2e1 * t14 * t67 + 0.2e1 * t39 * t91, 0.2e1 * t14 * t64 + t6 * t66, -0.2e1 * t53 * t15 + 0.2e1 * t16 * t88, 0.2e1 * t54 * t15 - 0.2e1 * t16 * t90, 0.2e1 * t15 * t7, 0.2e1 * t56 * t15 + 0.2e1 * t8 * t88, 0.2e1 * t57 * t15 - 0.2e1 * t8 * t90; 0, 0, 0, 0, 0, -t77, t78, 0, 0, 0, 0, 0, t17, t48, 0, 0, 0, 0, 0, t6, -t7, 0, 0, 0, t3, t1, t5, -t55, 0, -t63 * t36 + t50 * t40, t50 * t44 + t63 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t72, -0.2e1 * t71, 0, 0, 0, 0, 0, -0.2e1 * t12, 0.2e1 * t11, t31, t24, 0, 0, 0, 0.2e1 * t65, 0.2e1 * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t48, 0, 0, 0, 0, 0, t6, -t7, 0, 0, 0, t3, t1, t5, -t55, 0, -t62 * t36 + t49 * t40, t49 * t44 + t62 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, -t71, 0, 0, 0, 0, 0, (-pkin(3) - t34) * t81 - t93, t11 - t69, t31, t24, 0, 0, 0, t18 + t26 + (-t12 - t70) * t44, t84 + t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t70, -0.2e1 * t69, t31, t24, 0, 0, 0, 0.2e1 * t58, 0.2e1 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t7, 0, 0, 0, t3, t1, t5, -t55, 0, t54 * pkin(4) - t5 * pkin(6), t53 * pkin(4) + t55 * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, t11, t31, t24, 0, 0, 0, t65 - t74, -t73 + t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, -t69, t31, t24, 0, 0, 0, t58 - t74, -t73 + t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t24, 0, 0, 0, -0.2e1 * t74, -0.2e1 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, t54, t7, t56, t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t79, 0, t40 * t11 - t22 * t36, t44 * t11 + t22 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t79, 0, -t32 * t36 - t40 * t69, t32 * t79 - t44 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t79, 0, -pkin(6) * t36, pkin(6) * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t4;
