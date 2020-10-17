% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% MM_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRRP8_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:00:36
% EndTime: 2019-12-05 17:00:37
% DurationCPUTime: 0.32s
% Computational Cost: add. (218->72), mult. (518->134), div. (0->0), fcn. (568->8), ass. (0->50)
t24 = sin(qJ(4));
t27 = cos(qJ(4));
t32 = -t27 * pkin(4) - t24 * qJ(5);
t14 = -pkin(3) + t32;
t55 = -0.2e1 * t14;
t25 = sin(qJ(3));
t54 = -0.2e1 * t25;
t28 = cos(qJ(3));
t53 = 0.2e1 * t28;
t52 = pkin(3) * t24;
t51 = pkin(3) * t27;
t50 = pkin(7) * t24;
t49 = pkin(7) * t27;
t48 = t24 * pkin(8);
t47 = t27 * pkin(8);
t23 = cos(pkin(5));
t22 = sin(pkin(5));
t44 = t22 * sin(qJ(2));
t10 = -t23 * t28 + t25 * t44;
t46 = t10 * t24;
t45 = t10 * t27;
t43 = t22 * cos(qJ(2));
t42 = t24 * t25;
t41 = t24 * t27;
t40 = t24 * t28;
t18 = t27 * t25;
t39 = t27 * t28;
t15 = -t28 * pkin(3) - t25 * pkin(8) - pkin(2);
t8 = pkin(7) * t39 + t24 * t15;
t19 = t24 ^ 2;
t21 = t27 ^ 2;
t38 = t19 + t21;
t37 = t28 * qJ(5);
t36 = t25 * t53;
t11 = t23 * t25 + t28 * t44;
t2 = t11 * t24 + t27 * t43;
t35 = t10 * t42 + t2 * t28;
t3 = t11 * t27 - t24 * t43;
t34 = t2 * t24 + t3 * t27;
t4 = -t37 + t8;
t13 = t27 * t15;
t5 = -t13 + (pkin(4) + t50) * t28;
t33 = t5 * t24 + t4 * t27;
t31 = -pkin(4) * t24 + t27 * qJ(5);
t20 = t25 ^ 2;
t16 = pkin(8) * t40;
t9 = (pkin(7) - t31) * t25;
t7 = -pkin(7) * t40 + t13;
t1 = t10 * t18 + t3 * t28;
t6 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, t43, -t44, 0, 0, 0, 0, 0, t28 * t43, -t25 * t43, 0, 0, 0, 0, 0, t35, t1, t35, (t2 * t27 - t24 * t3) * t25, -t1, t10 * t9 + t2 * t5 + t3 * t4; 0, 1, 0, 0, t20, t36, 0, 0, 0, pkin(2) * t53, pkin(2) * t54, t21 * t20, -0.2e1 * t20 * t41, t39 * t54, t24 * t36, t28 ^ 2, 0.2e1 * t20 * t50 - 0.2e1 * t7 * t28, 0.2e1 * t20 * t49 + 0.2e1 * t8 * t28, 0.2e1 * t5 * t28 + 0.2e1 * t9 * t42, 0.2e1 * (-t24 * t4 + t27 * t5) * t25, -0.2e1 * t9 * t18 - 0.2e1 * t4 * t28, t4 ^ 2 + t5 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, 0, 0, 0, 0, 0, -t45, t46, -t45, t34, -t46, t34 * pkin(8) + t10 * t14; 0, 0, 0, 0, 0, 0, t25, t28, 0, -t25 * pkin(7), -t28 * pkin(7), t24 * t18, (-t19 + t21) * t25, -t40, -t39, 0, t16 + (-t49 - t52) * t25, pkin(8) * t39 + (t50 - t51) * t25, t14 * t42 - t9 * t27 + t16, t33, -t9 * t24 + (-pkin(8) * t28 - t14 * t25) * t27, t33 * pkin(8) + t9 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t19, 0.2e1 * t41, 0, 0, 0, 0.2e1 * t51, -0.2e1 * t52, t27 * t55, 0.2e1 * t38 * pkin(8), t24 * t55, t38 * pkin(8) ^ 2 + t14 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, -t3, -t2, 0, t3, -t2 * pkin(4) + t3 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t42, -t28, t7, -t8, t13 + (-0.2e1 * pkin(4) - t50) * t28, t32 * t25, -0.2e1 * t37 + t8, -t5 * pkin(4) + t4 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t27, 0, -t48, -t47, -t48, t31, t47, t31 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t18, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t6;
