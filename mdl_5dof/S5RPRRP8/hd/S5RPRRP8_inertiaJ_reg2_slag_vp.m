% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRP8_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_inertiaJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:30
% EndTime: 2019-12-31 18:47:32
% DurationCPUTime: 0.63s
% Computational Cost: add. (245->61), mult. (371->77), div. (0->0), fcn. (306->4), ass. (0->44)
t35 = sin(qJ(3));
t34 = sin(qJ(4));
t30 = t34 ^ 2;
t36 = cos(qJ(4));
t32 = t36 ^ 2;
t40 = t30 + t32;
t62 = t40 * t35;
t37 = cos(qJ(3));
t38 = -pkin(1) - pkin(2);
t18 = t37 * qJ(2) + t35 * t38;
t14 = -pkin(7) + t18;
t61 = t40 * t14;
t59 = -0.2e1 * t34;
t58 = 0.2e1 * t36;
t57 = t62 * t14;
t56 = pkin(7) * t61;
t54 = t34 * pkin(7);
t53 = t36 * pkin(7);
t16 = t35 * qJ(2) - t37 * t38;
t13 = pkin(3) + t16;
t52 = pkin(3) + t13;
t51 = t40 * t14 ^ 2;
t19 = -t36 * pkin(4) - t34 * qJ(5) - pkin(3);
t3 = t16 - t19;
t50 = t19 - t3;
t47 = t34 * t14;
t46 = t34 * t35;
t45 = t34 * t36;
t44 = t36 * t14;
t43 = t36 * t35;
t42 = t62 * pkin(7);
t41 = t40 * pkin(7) ^ 2;
t20 = -t34 * pkin(4) + t36 * qJ(5);
t33 = t37 ^ 2;
t31 = t35 ^ 2;
t26 = t37 * t36;
t25 = t37 * t34;
t24 = -0.2e1 * t45;
t23 = 0.2e1 * t45;
t15 = 0.2e1 * t40 * pkin(7);
t8 = t40 * t31 + t33;
t2 = -0.2e1 * t61;
t1 = (-pkin(7) + t14) * t40;
t4 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2 * pkin(1), 0, 0.2e1 * qJ(2), (pkin(1) ^ 2) + qJ(2) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t16, 0.2e1 * t18, 0, t16 ^ 2 + t18 ^ 2, t30, t23, 0, t32, 0, 0, t13 * t58, t13 * t59, t2, t13 ^ 2 + t51, t30, 0, t24, 0, 0, t32, t3 * t58, t2, 0.2e1 * t3 * t34, t3 ^ 2 + t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(1), 0, 0, 0, 0, 0, 0, -t37, t35, 0, -t16 * t37 + t18 * t35, 0, 0, 0, 0, 0, 0, -t26, t25, -t62, -t13 * t37 + t57, 0, 0, 0, 0, 0, 0, -t26, -t62, -t25, -t3 * t37 + t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31 + t33, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t16, -t18, 0, 0, -t30, t24, 0, -t32, 0, 0, -t52 * t36, t52 * t34, t1, -t13 * pkin(3) + t56, -t30, 0, t23, 0, 0, -t32, t50 * t36, t1, t50 * t34, t3 * t19 + t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t35, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t25, t62, t37 * pkin(3) + t42, 0, 0, 0, 0, 0, 0, t26, t62, t25, -t37 * t19 + t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t30, t23, 0, t32, 0, 0, pkin(3) * t58, pkin(3) * t59, t15, pkin(3) ^ 2 + t41, t30, 0, t24, 0, 0, t32, -0.2e1 * t19 * t36, t15, t19 * t59, t19 ^ 2 + t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, 0, -t36, 0, -t47, -t44, 0, 0, 0, -t34, 0, 0, t36, 0, -t47, -t20, t44, t20 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, -t43, 0, 0, 0, 0, 0, 0, 0, 0, -t46, 0, t43, t20 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, t36, 0, -t54, -t53, 0, 0, 0, t34, 0, 0, -t36, 0, -t54, t20, t53, t20 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, 0, t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t4;
