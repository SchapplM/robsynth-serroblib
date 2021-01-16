% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:52
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRP8_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:52:19
% EndTime: 2021-01-15 20:52:21
% DurationCPUTime: 0.26s
% Computational Cost: add. (194->52), mult. (338->83), div. (0->0), fcn. (348->4), ass. (0->30)
t28 = sin(qJ(4));
t29 = sin(qJ(2));
t30 = cos(qJ(4));
t31 = cos(qJ(2));
t9 = -t29 * t28 - t31 * t30;
t41 = -0.2e1 * t9;
t10 = -t31 * t28 + t29 * t30;
t40 = 0.2e1 * t10;
t39 = -0.2e1 * t29;
t38 = 0.2e1 * t31;
t37 = -pkin(2) - pkin(3);
t36 = t29 * pkin(6);
t26 = t29 ^ 2;
t35 = t31 ^ 2 + t26;
t16 = -t31 * pkin(2) - t29 * qJ(3) - pkin(1);
t34 = (pkin(6) - pkin(7)) * t29;
t23 = t31 * pkin(6);
t17 = -t31 * pkin(7) + t23;
t5 = t28 * t17 - t30 * t34;
t7 = t31 * pkin(3) - t16;
t14 = t28 * qJ(3) - t30 * t37;
t33 = -t29 * pkin(2) + t31 * qJ(3);
t6 = t30 * t17 + t28 * t34;
t15 = t30 * qJ(3) + t28 * t37;
t12 = pkin(4) + t14;
t8 = 0.2e1 * t15;
t4 = -t9 * pkin(4) + t7;
t3 = t9 * qJ(5) + t6;
t1 = t10 * qJ(5) + t5;
t2 = [1, 0, 0, t26, t29 * t38, 0, 0, 0, pkin(1) * t38, pkin(1) * t39, -0.2e1 * t16 * t31, 0.2e1 * t35 * pkin(6), t16 * t39, t35 * pkin(6) ^ 2 + t16 ^ 2, t10 ^ 2, t9 * t40, 0, 0, 0, t7 * t41, t7 * t40, t4 * t41, t4 * t40, 0.2e1 * t1 * t10 + 0.2e1 * t3 * t9, t1 ^ 2 + t3 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, t29, t31, 0, -t36, -t23, -t36, t33, t23, t33 * pkin(6), 0, 0, -t10, -t9, 0, t5, t6, t1, t3, t12 * t10 + t15 * t9, t1 * t12 + t3 * t15; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2), 0, 0.2e1 * qJ(3), pkin(2) ^ 2 + qJ(3) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t14, t8, 0.2e1 * t12, t8, 0, t12 ^ 2 + t15 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, t36, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30 * t10 + t28 * t9, -t1 * t30 + t3 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(2), 0, 0, 0, 0, 0, -t30, t28, -t30, t28, 0, -t12 * t30 + t15 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28 ^ 2 + t30 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t9, 0, -t5, -t6, -t1, -t3, -t10 * pkin(4), -t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t14, -t15, -0.2e1 * pkin(4) - t14, -t15, 0, -t12 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t28, t30, -t28, 0, t30 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, t10, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t2;
