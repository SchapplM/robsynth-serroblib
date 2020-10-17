% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
% 
% Output:
% MM_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRPRP4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:44:27
% EndTime: 2019-05-05 17:44:28
% DurationCPUTime: 0.46s
% Computational Cost: add. (346->80), mult. (556->130), div. (0->0), fcn. (523->6), ass. (0->56)
t37 = sin(qJ(5));
t39 = cos(qJ(5));
t17 = pkin(5) * t39 + t37 * qJ(6);
t40 = cos(qJ(3));
t63 = t17 * t40;
t41 = -pkin(3) - pkin(8);
t62 = t41 * t40;
t61 = 0.2e1 * t40;
t60 = 2 * qJ(4);
t35 = sin(pkin(9));
t23 = t35 * pkin(1) + pkin(7);
t38 = sin(qJ(3));
t18 = t38 * t23;
t12 = t38 * pkin(4) + t18;
t36 = cos(pkin(9));
t24 = -t36 * pkin(1) - pkin(2);
t50 = t38 * qJ(4);
t43 = t24 - t50;
t8 = t43 + t62;
t5 = t37 * t12 + t39 * t8;
t59 = t38 * pkin(5);
t58 = t40 * pkin(3);
t57 = t37 * t40;
t56 = t37 * t41;
t55 = t38 * t40;
t54 = t38 * t41;
t53 = t39 * t37;
t52 = t39 * t40;
t28 = t39 * t41;
t19 = t40 * t23;
t13 = t40 * pkin(4) + t19;
t30 = t37 ^ 2;
t32 = t39 ^ 2;
t21 = t30 + t32;
t31 = t38 ^ 2;
t33 = t40 ^ 2;
t51 = t31 + t33;
t49 = t38 * qJ(6);
t48 = t40 * qJ(4);
t47 = -0.2e1 * t55;
t46 = -t39 * t12 + t37 * t8;
t2 = t49 + t5;
t3 = t46 - t59;
t1 = t2 * t37 - t3 * t39;
t45 = -t38 * pkin(3) + t48;
t44 = t37 * pkin(5) - t39 * qJ(6);
t27 = t39 * t38;
t26 = t37 * t38;
t25 = t30 * t33;
t20 = t38 * t28;
t16 = qJ(4) + t44;
t15 = t21 * t41;
t14 = t21 * t40;
t11 = t43 - t58;
t6 = t13 + t63;
t4 = [1, 0, 0 (t35 ^ 2 + t36 ^ 2) * pkin(1) ^ 2, t31, 0.2e1 * t55, 0, 0, 0, -0.2e1 * t24 * t40, 0.2e1 * t24 * t38, 0.2e1 * t51 * t23, t11 * t61, -0.2e1 * t11 * t38, t51 * t23 ^ 2 + t11 ^ 2, t25, 0.2e1 * t33 * t53, t37 * t47, t39 * t47, t31, 0.2e1 * t13 * t52 - 0.2e1 * t38 * t46, -0.2e1 * t13 * t57 - 0.2e1 * t5 * t38, -0.2e1 * t3 * t38 + 0.2e1 * t6 * t52 (-t2 * t39 - t3 * t37) * t61, 0.2e1 * t2 * t38 + 0.2e1 * t6 * t57, t2 ^ 2 + t3 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t40 + t6 * t38; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32 * t33 + t25 + t31; 0, 0, 0, 0, 0, 0, t38, t40, 0, -t18, -t19, t45, t18, t19, t45 * t23, -t37 * t52 (t30 - t32) * t40, t27, -t26, 0, t13 * t37 + t39 * t48 + t20, t13 * t39 + (-t48 - t54) * t37, t16 * t52 + t6 * t37 + t20, -t1, -t6 * t39 + (t16 * t40 + t54) * t37, t1 * t41 + t6 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t38, 0, -t40, t38, t50 + t58, 0, 0, 0, 0, 0, t26, t27, t26, t14, -t27, t38 * t16 - t21 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(3), t60, pkin(3) ^ 2 + (qJ(4) ^ 2) t32, -0.2e1 * t53, 0, 0, 0, t37 * t60, t39 * t60, 0.2e1 * t16 * t37, -0.2e1 * t15, -0.2e1 * t16 * t39, t21 * t41 ^ 2 + t16 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, 0, t18, 0, 0, 0, 0, 0, t27, -t26, t27, 0, t26, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(3), 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, -t52, t38, -t46, -t5, -t46 + 0.2e1 * t59, t44 * t40, 0.2e1 * t49 + t5, -t3 * pkin(5) + t2 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, t57, -t52, 0, -t57, -t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t37, 0, t28, -t56, t28, -t17, t56, t17 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t37, t39, 0, t37, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, -t57, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t4;
