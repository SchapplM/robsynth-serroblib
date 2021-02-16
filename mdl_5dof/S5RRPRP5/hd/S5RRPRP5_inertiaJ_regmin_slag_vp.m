% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:19
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRP5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:18:18
% EndTime: 2021-01-15 20:18:21
% DurationCPUTime: 0.30s
% Computational Cost: add. (405->52), mult. (772->96), div. (0->0), fcn. (873->6), ass. (0->36)
t25 = sin(pkin(8));
t26 = cos(pkin(8));
t28 = sin(qJ(2));
t29 = cos(qJ(2));
t16 = t25 * t28 - t26 * t29;
t17 = t25 * t29 + t26 * t28;
t27 = sin(qJ(4));
t36 = cos(qJ(4));
t7 = -t27 * t16 + t36 * t17;
t42 = -0.2e1 * t7;
t24 = -t29 * pkin(2) - pkin(1);
t10 = t16 * pkin(3) + t24;
t41 = 0.2e1 * t10;
t40 = 0.2e1 * t24;
t39 = 0.2e1 * t29;
t38 = t25 * pkin(2);
t37 = t26 * pkin(2);
t35 = -qJ(3) - pkin(6);
t23 = pkin(3) + t37;
t14 = -t27 * t23 - t36 * t38;
t34 = -t36 * t23 + t27 * t38;
t18 = t35 * t28;
t19 = t35 * t29;
t8 = t26 * t18 + t25 * t19;
t9 = t25 * t18 - t26 * t19;
t33 = -t17 * pkin(7) + t8;
t31 = 2 * pkin(4);
t30 = 2 * qJ(5);
t12 = -pkin(4) + t34;
t11 = qJ(5) - t14;
t6 = t36 * t16 + t27 * t17;
t5 = -t16 * pkin(7) + t9;
t3 = t27 * t33 + t36 * t5;
t2 = t27 * t5 - t36 * t33;
t1 = t6 * pkin(4) - t7 * qJ(5) + t10;
t4 = [1, 0, 0, t28 ^ 2, t28 * t39, 0, 0, 0, pkin(1) * t39, -0.2e1 * pkin(1) * t28, t16 * t40, t17 * t40, -0.2e1 * t9 * t16 - 0.2e1 * t8 * t17, t24 ^ 2 + t8 ^ 2 + t9 ^ 2, t7 ^ 2, t6 * t42, 0, 0, 0, t6 * t41, t7 * t41, 0.2e1 * t1 * t6, 0.2e1 * t2 * t7 - 0.2e1 * t3 * t6, t1 * t42, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, t28, t29, 0, -t28 * pkin(6), -t29 * pkin(6), t8, -t9, (-t16 * t25 - t17 * t26) * pkin(2), (t25 * t9 + t26 * t8) * pkin(2), 0, 0, t7, -t6, 0, -t2, -t3, -t2, -t11 * t6 + t12 * t7, t3, t3 * t11 + t2 * t12; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * t37, -0.2e1 * t38, 0, (t25 ^ 2 + t26 ^ 2) * pkin(2) ^ 2, 0, 0, 0, 0, 1, -0.2e1 * t34, 0.2e1 * t14, -0.2e1 * t12, 0, 0.2e1 * t11, t11 ^ 2 + t12 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t17, 0, t24, 0, 0, 0, 0, 0, t6, t7, t6, 0, -t7, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t6, 0, -t2, -t3, -t2, -pkin(4) * t7 - t6 * qJ(5), t3, -t2 * pkin(4) + t3 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -t34, t14, t31 - t34, 0, t30 - t14, -t12 * pkin(4) + t11 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t31, 0, t30, pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t4;
