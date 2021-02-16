% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRRP1
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
% Datum: 2021-01-15 23:53
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:52:44
% EndTime: 2021-01-15 23:52:46
% DurationCPUTime: 0.36s
% Computational Cost: add. (490->55), mult. (936->94), div. (0->0), fcn. (1082->6), ass. (0->44)
t36 = cos(qJ(2));
t27 = -t36 * pkin(2) - pkin(1);
t32 = sin(qJ(3));
t33 = sin(qJ(2));
t35 = cos(qJ(3));
t39 = t32 * t33 - t35 * t36;
t14 = t39 * pkin(3) + t27;
t20 = t32 * t36 + t35 * t33;
t31 = sin(qJ(4));
t34 = cos(qJ(4));
t9 = t31 * t20 + t34 * t39;
t5 = t9 * pkin(4) + t14;
t47 = 0.2e1 * t5;
t46 = 0.2e1 * t14;
t45 = 0.2e1 * t27;
t44 = 0.2e1 * t36;
t43 = pkin(6) + pkin(7);
t42 = t31 * pkin(3);
t41 = t32 * pkin(2);
t40 = t34 * t41;
t22 = t43 * t33;
t23 = t43 * t36;
t11 = -t35 * t22 - t32 * t23;
t7 = -t20 * pkin(8) + t11;
t12 = t32 * t22 - t35 * t23;
t8 = -t39 * pkin(8) - t12;
t3 = -t31 * t8 + t34 * t7;
t30 = t35 * pkin(2);
t26 = t30 + pkin(3);
t17 = t34 * t26 - t31 * t41;
t4 = -t31 * t7 - t34 * t8;
t37 = 0.2e1 * pkin(4);
t38 = t17 + t37;
t29 = t34 * pkin(3);
t28 = -0.2e1 * t42;
t25 = t29 + pkin(4);
t18 = t31 * t26 + t40;
t16 = pkin(4) + t17;
t15 = 0.2e1 * t18;
t13 = -t40 + (-pkin(3) - t26) * t31;
t10 = t34 * t20 - t31 * t39;
t2 = -t9 * qJ(5) - t4;
t1 = -t10 * qJ(5) + t3;
t6 = [1, 0, 0, t33 ^ 2, t33 * t44, 0, 0, 0, pkin(1) * t44, -0.2e1 * pkin(1) * t33, t20 ^ 2, -0.2e1 * t20 * t39, 0, 0, 0, t39 * t45, t20 * t45, t10 ^ 2, -0.2e1 * t10 * t9, 0, 0, 0, t9 * t46, t10 * t46, t9 * t47, t10 * t47, -0.2e1 * t1 * t10 - 0.2e1 * t2 * t9, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t33, t36, 0, -t33 * pkin(6), -t36 * pkin(6), 0, 0, t20, -t39, 0, t11, t12, 0, 0, t10, -t9, 0, t3, t4, t1, -t2, -t16 * t10 - t18 * t9, t1 * t16 + t2 * t18; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t30, -0.2e1 * t41, 0, 0, 0, 0, 1, 0.2e1 * t17, -t15, 0.2e1 * t16, -t15, 0, t16 ^ 2 + t18 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t39, 0, t11, t12, 0, 0, t10, -t9, 0, t3, t4, t1, -t2, -t25 * t10 - t9 * t42, t1 * t25 + t2 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t30, -t41, 0, 0, 0, 0, 1, t17 + t29, t13, t29 + t38, t13, 0, t16 * t25 + t18 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t29, t28, 0.2e1 * t25, t28, 0, t31 ^ 2 * pkin(3) ^ 2 + t25 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, 0, t3, t4, t1, -t2, -t10 * pkin(4), t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t17, -t18, t38, -t18, 0, t16 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t29, -t42, t37 + t29, -t42, 0, t25 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t37, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t6;
