% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRPPR7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:16:28
% EndTime: 2019-05-05 17:16:29
% DurationCPUTime: 0.41s
% Computational Cost: add. (344->68), mult. (569->99), div. (0->0), fcn. (664->6), ass. (0->48)
t37 = sin(pkin(9));
t38 = cos(pkin(9));
t41 = cos(qJ(3));
t58 = sin(qJ(3));
t20 = -t37 * t58 + t38 * t41;
t21 = -t37 * t41 - t38 * t58;
t64 = (t20 * t38 - t21 * t37) * pkin(3);
t18 = t21 ^ 2;
t62 = t20 ^ 2;
t63 = t18 + t62;
t61 = -0.2e1 * t20;
t26 = t37 * pkin(3) + qJ(5);
t60 = 0.2e1 * t26;
t59 = 2 * qJ(2);
t57 = t21 * t26;
t39 = sin(qJ(6));
t56 = t39 * t20;
t12 = t39 * t21;
t40 = cos(qJ(6));
t55 = t40 * t20;
t14 = t40 * t21;
t54 = t40 * t39;
t32 = t58 * pkin(3) + qJ(2);
t53 = t21 * t61;
t42 = -pkin(1) - pkin(7);
t51 = t58 * t42;
t24 = -t58 * qJ(4) + t51;
t33 = t41 * t42;
t50 = -t41 * qJ(4) + t33;
t10 = t38 * t24 + t37 * t50;
t8 = t37 * t24 - t38 * t50;
t52 = t10 ^ 2 + t8 ^ 2;
t31 = -t38 * pkin(3) - pkin(4);
t49 = -t20 * qJ(5) + t32;
t48 = t10 * t21 + t8 * t20;
t25 = -pkin(8) + t31;
t47 = t20 * t25 + t57;
t46 = t31 * t20 + t57;
t44 = 0.2e1 * t48;
t36 = t40 ^ 2;
t35 = t39 ^ 2;
t7 = -t21 * pkin(4) + t49;
t6 = t21 * pkin(5) + t10;
t5 = t20 * pkin(5) + t8;
t4 = (-pkin(4) - pkin(8)) * t21 + t49;
t2 = t39 * t5 + t40 * t4;
t1 = -t39 * t4 + t40 * t5;
t3 = [1, 0, 0, -2 * pkin(1), t59, pkin(1) ^ 2 + qJ(2) ^ 2, t41 ^ 2, -0.2e1 * t41 * t58, 0, 0, 0, t58 * t59, t41 * t59, t44, t32 ^ 2 + t52, t44, 0.2e1 * t7 * t21, t7 * t61, t7 ^ 2 + t52, t35 * t18, 0.2e1 * t18 * t54, t39 * t53, t40 * t53, t62, 0.2e1 * t1 * t20 + 0.2e1 * t6 * t14, -0.2e1 * t6 * t12 - 0.2e1 * t2 * t20; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, -t63, -t48, -t63, 0, 0, -t48, 0, 0, 0, 0, 0, -t63 * t40, t63 * t39; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t63, 0, 0, 0, t63, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t41, -t58, 0, t33, -t51, -t64 (t10 * t37 - t38 * t8) * pkin(3), t46, t8, t10, t10 * t26 + t8 * t31, -t21 * t54 (t35 - t36) * t21, t55, -t56, 0, t6 * t39 + t47 * t40, -t47 * t39 + t6 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t58, 0, t64, 0, -t20, -t21, -t46, 0, 0, 0, 0, 0, -t12, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t37 ^ 2 + t38 ^ 2) * pkin(3) ^ 2, 0, 0.2e1 * t31, t60, t26 ^ 2 + t31 ^ 2, t36, -0.2e1 * t54, 0, 0, 0, t39 * t60, t40 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, t21, -t20, t7, 0, 0, 0, 0, 0, -t56, -t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, 0, t8, 0, 0, 0, 0, 0, t55, -t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, t31, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, -t14, t20, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t39, 0, t40 * t25, -t39 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, -t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
