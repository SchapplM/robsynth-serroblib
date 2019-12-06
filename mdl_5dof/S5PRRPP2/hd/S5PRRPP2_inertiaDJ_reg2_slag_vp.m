% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRPP2_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP2_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:10:09
% EndTime: 2019-12-05 16:10:12
% DurationCPUTime: 0.58s
% Computational Cost: add. (406->76), mult. (1161->150), div. (0->0), fcn. (1012->6), ass. (0->60)
t44 = sin(qJ(2));
t45 = cos(qJ(3));
t46 = cos(qJ(2));
t70 = t46 * qJD(2);
t43 = sin(qJ(3));
t72 = t43 * qJD(3);
t78 = t44 * t72 - t45 * t70;
t77 = 2 * qJD(5);
t42 = sin(pkin(8));
t76 = t42 * t43;
t75 = -qJ(4) - pkin(6);
t40 = t43 ^ 2;
t41 = t45 ^ 2;
t74 = t40 + t41;
t73 = cos(pkin(8));
t39 = t44 * qJD(2);
t71 = t45 * qJD(3);
t59 = t73 * t43;
t26 = t42 * t45 + t59;
t23 = t26 * qJD(3);
t58 = t73 * t45;
t50 = t58 - t76;
t69 = -0.2e1 * t50 * t23;
t68 = -0.2e1 * pkin(2) * qJD(3);
t38 = pkin(3) * t72;
t67 = t43 * t71;
t65 = t44 * t70;
t64 = t46 * t72;
t63 = t43 * t70;
t37 = -t45 * pkin(3) - pkin(2);
t57 = qJD(3) * t75;
t22 = t45 * qJD(4) + t43 * t57;
t49 = -t43 * qJD(4) + t45 * t57;
t11 = t42 * t22 - t73 * t49;
t12 = t73 * t22 + t42 * t49;
t32 = t75 * t45;
t17 = -t42 * t32 - t75 * t59;
t18 = -t73 * t32 + t75 * t76;
t61 = t17 * t11 + t18 * t12;
t60 = t74 * t46;
t56 = t73 * t70;
t55 = qJD(3) * t58;
t24 = -t42 * t72 + t55;
t54 = t26 * t23 - t24 * t50;
t53 = -t46 * t24 + t26 * t39;
t10 = t44 * t23 + t42 * t63 - t45 * t56;
t20 = t26 * t44;
t21 = t50 * t44;
t9 = t78 * t42 - t43 * t56 - t44 * t55;
t52 = -t10 * t18 + t20 * t11 + t21 * t12 - t9 * t17;
t51 = -t10 * t50 + t20 * t24 - t21 * t23 - t9 * t26;
t48 = 0.2e1 * t11 * t26 + 0.2e1 * t12 * t50 + 0.2e1 * t17 * t24 - 0.2e1 * t18 * t23;
t47 = -0.2e1 * t21 * t10 - 0.2e1 * t20 * t9 - 0.2e1 * t65;
t36 = -t73 * pkin(3) - pkin(4);
t34 = t42 * pkin(3) + qJ(5);
t16 = 0.2e1 * t26 * t24;
t14 = -pkin(4) * t50 - t26 * qJ(5) + t37;
t13 = -t46 * t23 - t39 * t50;
t5 = t23 * pkin(4) - t24 * qJ(5) - t26 * qJD(5) + t38;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * (-0.1e1 + t74) * t65, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, -t70, 0, 0, 0, 0, 0, 0, 0, 0, -t45 * t39 - t64, t43 * t39 - t46 * t71, qJD(2) * t60, (-pkin(2) * t44 + pkin(6) * t60) * qJD(2), 0, 0, 0, 0, 0, 0, t13, t53, t51, -pkin(3) * t64 + t37 * t39 + t52, 0, 0, 0, 0, 0, 0, t13, t51, -t53, t14 * t39 - t46 * t5 + t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t67, 0.2e1 * (-t40 + t41) * qJD(3), 0, -0.2e1 * t67, 0, 0, t43 * t68, t45 * t68, 0, 0, t16, -0.2e1 * t54, 0, t69, 0, 0, 0.2e1 * t37 * t23 - 0.2e1 * t38 * t50, 0.2e1 * t37 * t24 + 0.2e1 * t26 * t38, t48, 0.2e1 * t37 * t38 + 0.2e1 * t61, t16, 0, 0.2e1 * t54, 0, 0, t69, 0.2e1 * t14 * t23 - 0.2e1 * t5 * t50, t48, -0.2e1 * t14 * t24 - 0.2e1 * t5 * t26, 0.2e1 * t14 * t5 + 0.2e1 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44 * t71 - t63, t78, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, (-t10 * t42 + t73 * t9) * pkin(3), 0, 0, 0, 0, 0, 0, t9, 0, -t10, t21 * qJD(5) - t10 * t34 - t9 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, 0, -t72, 0, -pkin(6) * t71, pkin(6) * t72, 0, 0, 0, 0, t24, 0, -t23, 0, -t11, -t12, (-t23 * t42 - t73 * t24) * pkin(3), (-t73 * t11 + t12 * t42) * pkin(3), 0, t24, 0, 0, t23, 0, -t11, qJD(5) * t50 - t34 * t23 + t36 * t24, t12, t18 * qJD(5) + t11 * t36 + t12 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, t34 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t24, 0, t38, 0, 0, 0, 0, 0, 0, t23, 0, -t24, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
