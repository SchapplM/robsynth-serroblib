% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRP7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP7_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP7_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:57:44
% EndTime: 2019-12-31 21:57:46
% DurationCPUTime: 0.54s
% Computational Cost: add. (537->79), mult. (1028->148), div. (0->0), fcn. (1121->6), ass. (0->67)
t41 = sin(qJ(3));
t73 = t41 * pkin(2);
t33 = pkin(8) + t73;
t40 = sin(qJ(4));
t38 = t40 ^ 2;
t43 = cos(qJ(4));
t39 = t43 ^ 2;
t57 = t38 + t39;
t59 = t57 * t33;
t68 = cos(qJ(2));
t35 = -t68 * pkin(2) - pkin(1);
t79 = 0.2e1 * t35;
t78 = -0.2e1 * t40;
t77 = -0.2e1 * t43;
t42 = sin(qJ(2));
t67 = cos(qJ(3));
t21 = t41 * t42 - t67 * t68;
t76 = pkin(8) * t21;
t75 = t21 * pkin(4);
t74 = t40 * pkin(8);
t72 = t43 * pkin(8);
t26 = (-pkin(7) - pkin(6)) * t42;
t53 = t68 * pkin(6);
t27 = t68 * pkin(7) + t53;
t13 = -t67 * t26 + t41 * t27;
t22 = t41 * t68 + t67 * t42;
t47 = pkin(4) * t40 - t43 * qJ(5);
t6 = t47 * t22 + t13;
t71 = t6 * t40;
t70 = t6 * t43;
t52 = t67 * pkin(2);
t34 = -t52 - pkin(3);
t69 = pkin(3) - t34;
t14 = t41 * t26 + t67 * t27;
t9 = t21 * pkin(3) - t22 * pkin(8) + t35;
t5 = t43 * t14 + t40 * t9;
t66 = t13 * t43;
t65 = t21 * t33;
t64 = t40 * t22;
t63 = t40 * t33;
t62 = t40 * t43;
t18 = t43 * t22;
t61 = t43 * t33;
t48 = -t43 * pkin(4) - t40 * qJ(5);
t24 = -pkin(3) + t48;
t19 = -t52 + t24;
t60 = -t19 - t24;
t58 = t57 * pkin(8);
t56 = t21 * qJ(5);
t55 = 0.2e1 * t68;
t54 = -0.2e1 * t22 * t21;
t51 = t40 * t14 - t43 * t9;
t50 = -pkin(3) * t22 - t76;
t2 = t56 + t5;
t3 = t51 - t75;
t1 = t2 * t43 + t3 * t40;
t49 = -t22 * t24 + t76;
t46 = t19 * t22 - t65;
t45 = t22 * t34 - t65;
t30 = 0.2e1 * t62;
t20 = t22 ^ 2;
t17 = t43 * t21;
t16 = t40 * t21;
t15 = t40 * t18;
t12 = t13 * t40;
t10 = (-t38 + t39) * t22;
t4 = [1, 0, 0, t42 ^ 2, t42 * t55, 0, 0, 0, pkin(1) * t55, -0.2e1 * pkin(1) * t42, t20, t54, 0, 0, 0, t21 * t79, t22 * t79, t39 * t20, -0.2e1 * t20 * t62, 0.2e1 * t21 * t18, t40 * t54, t21 ^ 2, 0.2e1 * t13 * t64 - 0.2e1 * t21 * t51, 0.2e1 * t13 * t18 - 0.2e1 * t5 * t21, -0.2e1 * t3 * t21 + 0.2e1 * t6 * t64, 0.2e1 * (-t2 * t40 + t3 * t43) * t22, -0.2e1 * t6 * t18 + 0.2e1 * t2 * t21, t2 ^ 2 + t3 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, t42, t68, 0, -t42 * pkin(6), -t53, 0, 0, t22, -t21, 0, -t13, -t14, t15, t10, t16, t17, 0, t45 * t40 - t66, t45 * t43 + t12, t46 * t40 - t70, t1, -t46 * t43 - t71, t1 * t33 + t6 * t19; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t52, -0.2e1 * t73, t38, t30, 0, 0, 0, t34 * t77, 0.2e1 * t34 * t40, t19 * t77, 0.2e1 * t59, t19 * t78, t57 * t33 ^ 2 + t19 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t21, 0, -t13, -t14, t15, t10, t16, t17, 0, t50 * t40 - t66, t50 * t43 + t12, -t49 * t40 - t70, t1, t49 * t43 - t71, t1 * pkin(8) + t6 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t52, -t73, t38, t30, 0, 0, 0, t69 * t43, -t69 * t40, t60 * t43, t58 + t59, t60 * t40, pkin(8) * t59 + t19 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t38, t30, 0, 0, 0, 0.2e1 * pkin(3) * t43, pkin(3) * t78, t24 * t77, 0.2e1 * t58, t24 * t78, t57 * pkin(8) ^ 2 + t24 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t64, t21, -t51, -t5, -t51 + 0.2e1 * t75, t48 * t22, 0.2e1 * t56 + t5, -t3 * pkin(4) + t2 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t43, 0, -t63, -t61, -t63, -t47, t61, -t47 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t43, 0, -t74, -t72, -t74, -t47, t72, -t47 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, t18, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t4;
