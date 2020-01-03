% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPR4_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t60 = sin(qJ(3));
t57 = t60 ^ 2;
t63 = cos(qJ(3));
t58 = t63 ^ 2;
t90 = t57 + t58;
t61 = sin(qJ(2));
t81 = t61 * pkin(1);
t45 = pkin(7) + t81;
t72 = t90 * t45;
t59 = sin(qJ(5));
t62 = cos(qJ(5));
t23 = t59 * t60 + t62 * t63;
t89 = 0.2e1 * t23;
t25 = -t59 * t63 + t60 * t62;
t88 = 0.2e1 * t25;
t87 = -0.2e1 * t60;
t86 = -0.2e1 * t63;
t85 = -pkin(3) - pkin(4);
t40 = t63 * t45;
t80 = t63 * pkin(8);
t20 = t40 - t80;
t67 = (-pkin(8) + t45) * t60;
t6 = t20 * t59 - t62 * t67;
t8 = t62 * t20 + t59 * t67;
t84 = -t23 * t8 + t25 * t6;
t53 = t63 * pkin(7);
t35 = t53 - t80;
t68 = (pkin(7) - pkin(8)) * t60;
t12 = t35 * t59 - t62 * t68;
t14 = t62 * t35 + t59 * t68;
t83 = t12 * t25 - t14 * t23;
t82 = t60 * pkin(7);
t64 = cos(qJ(2));
t56 = t64 * pkin(1);
t46 = -t56 - pkin(2);
t79 = pkin(2) - t46;
t78 = t60 * t45;
t77 = t60 * t63;
t69 = t63 * pkin(3) + t60 * qJ(4) + pkin(2);
t18 = -t56 - t69;
t54 = t63 * pkin(4);
t16 = t54 - t18;
t19 = t54 + t69;
t76 = t16 + t19;
t75 = -t18 + t69;
t74 = t72 * pkin(7);
t73 = t90 * t45 ^ 2;
t71 = t90 * pkin(7) ^ 2;
t70 = t90 * pkin(7);
t34 = -pkin(3) * t60 + qJ(4) * t63;
t42 = -0.2e1 * t77;
t41 = 0.2e1 * t77;
t30 = qJ(4) * t62 + t59 * t85;
t28 = qJ(4) * t59 - t62 * t85;
t27 = 0.2e1 * t70;
t22 = t25 ^ 2;
t21 = t23 ^ 2;
t15 = 0.2e1 * t72;
t11 = -0.2e1 * t25 * t23;
t10 = t70 + t72;
t9 = -t23 * t59 - t25 * t62;
t1 = -t23 * t30 + t25 * t28;
t2 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t56, -0.2e1 * t81, 0, (t61 ^ 2 + t64 ^ 2) * pkin(1) ^ 2, t57, t41, 0, t58, 0, 0, t46 * t86, 0.2e1 * t46 * t60, t15, t46 ^ 2 + t73, t57, 0, t42, 0, 0, t58, t18 * t86, t15, t18 * t87, t18 ^ 2 + t73, t22, t11, 0, t21, 0, 0, t16 * t89, t16 * t88, 0.2e1 * t84, t16 ^ 2 + t6 ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t56, -t81, 0, 0, t57, t41, 0, t58, 0, 0, t79 * t63, -t79 * t60, t10, -pkin(2) * t46 + t74, t57, 0, t42, 0, 0, t58, t75 * t63, t10, t75 * t60, -t18 * t69 + t74, t22, t11, 0, t21, 0, 0, t76 * t23, t76 * t25, t83 + t84, t12 * t6 + t14 * t8 + t16 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t57, t41, 0, t58, 0, 0, 0.2e1 * pkin(2) * t63, pkin(2) * t87, t27, pkin(2) ^ 2 + t71, t57, 0, t42, 0, 0, t58, -t69 * t86, t27, -t69 * t87, t69 ^ 2 + t71, t22, t11, 0, t21, 0, 0, t19 * t89, t19 * t88, 0.2e1 * t83, t12 ^ 2 + t14 ^ 2 + t19 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, t63, 0, -t78, -t40, 0, 0, 0, t60, 0, 0, -t63, 0, -t78, t34, t40, t34 * t45, 0, 0, -t25, 0, t23, 0, t6, t8, t1, t28 * t6 + t30 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, t63, 0, -t82, -t53, 0, 0, 0, t60, 0, 0, -t63, 0, -t82, t34, t53, t34 * pkin(7), 0, 0, -t25, 0, t23, 0, t12, t14, t1, t12 * t28 + t14 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, 0.2e1 * qJ(4), pkin(3) ^ 2 + qJ(4) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t28, 0.2e1 * t30, 0, t28 ^ 2 + t30 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, t78, 0, 0, 0, 0, 0, 0, 0, 0, t9, t59 * t8 - t6 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, t82, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t12 * t62 + t14 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3), 0, 0, 0, 0, 0, 0, -t62, t59, 0, -t28 * t62 + t30 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59 ^ 2 + t62 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, -t23, 0, -t6, -t8, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, -t23, 0, -t12, -t14, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t28, -t30, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, -t59, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t2;
