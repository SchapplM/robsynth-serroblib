% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% MM_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:34
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRRP5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:33:42
% EndTime: 2021-01-15 16:33:44
% DurationCPUTime: 0.24s
% Computational Cost: add. (167->51), mult. (360->76), div. (0->0), fcn. (422->6), ass. (0->32)
t23 = cos(qJ(3));
t16 = -t23 * pkin(3) - pkin(2);
t19 = sin(qJ(4));
t20 = sin(qJ(3));
t22 = cos(qJ(4));
t8 = t19 * t20 - t22 * t23;
t5 = t8 * pkin(4) + t16;
t34 = 0.2e1 * t5;
t33 = 0.2e1 * t16;
t32 = 0.2e1 * t23;
t31 = pkin(6) + pkin(7);
t30 = t19 * pkin(3);
t24 = cos(qJ(2));
t29 = t24 * t8;
t9 = t19 * t23 + t22 * t20;
t28 = t24 * t9;
t21 = sin(qJ(2));
t27 = t20 * t21;
t26 = t23 * t21;
t11 = t31 * t20;
t12 = t31 * t23;
t3 = -t22 * t11 - t19 * t12;
t4 = t19 * t11 - t22 * t12;
t25 = 0.2e1 * pkin(4);
t18 = t22 * pkin(3);
t17 = -0.2e1 * t30;
t15 = t18 + pkin(4);
t7 = -t19 * t27 + t22 * t26;
t6 = t9 * t21;
t2 = -t8 * qJ(5) - t4;
t1 = -t9 * qJ(5) + t3;
t10 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24 ^ 2 + t6 ^ 2 + t7 ^ 2; 0, 0, t24, -t21, 0, 0, 0, 0, 0, t24 * t23, -t24 * t20, 0, 0, 0, 0, 0, -t29, -t28, -t29, -t28, t6 * t9 - t7 * t8, -t6 * t1 + t7 * t2 - t24 * t5; 0, 1, 0, 0, t20 ^ 2, t20 * t32, 0, 0, 0, pkin(2) * t32, -0.2e1 * pkin(2) * t20, t9 ^ 2, -0.2e1 * t9 * t8, 0, 0, 0, t8 * t33, t9 * t33, t8 * t34, t9 * t34, -0.2e1 * t1 * t9 - 0.2e1 * t2 * t8, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, -t26, 0, 0, 0, 0, 0, -t6, -t7, -t6, -t7, 0, -t6 * t15 + t7 * t30; 0, 0, 0, 0, 0, 0, t20, t23, 0, -t20 * pkin(6), -t23 * pkin(6), 0, 0, t9, -t8, 0, t3, t4, t1, -t2, -t15 * t9 - t8 * t30, t1 * t15 + t2 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t18, t17, 0.2e1 * t15, t17, 0, t19 ^ 2 * pkin(3) ^ 2 + t15 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, -t7, -t6, -t7, 0, -t6 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t8, 0, t3, t4, t1, -t2, -t9 * pkin(4), t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t18, -t30, t25 + t18, -t30, 0, t15 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t25, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t9, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t10;
