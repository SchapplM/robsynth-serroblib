% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPP5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t24 = sin(qJ(3));
t26 = cos(qJ(2));
t25 = sin(qJ(2));
t38 = cos(qJ(3));
t34 = t38 * t25;
t9 = t24 * t26 + t34;
t43 = -0.2e1 * t9;
t22 = -t26 * pkin(2) - pkin(1);
t42 = 0.2e1 * t22;
t41 = 0.2e1 * t26;
t27 = pkin(3) + pkin(4);
t40 = -pkin(7) - pkin(6);
t23 = t24 * pkin(2);
t17 = t23 + qJ(4);
t36 = t24 * t25;
t8 = -t38 * t26 + t36;
t39 = t17 * t8;
t37 = qJ(4) * t8;
t35 = t38 * pkin(2);
t11 = t40 * t26;
t5 = -t24 * t11 - t40 * t34;
t20 = t35 + pkin(3);
t31 = 0.2e1 * pkin(3);
t33 = t31 + t35;
t32 = t9 * qJ(4) - t22;
t6 = -t38 * t11 + t40 * t36;
t30 = qJ(4) ^ 2;
t29 = 0.2e1 * qJ(4);
t18 = t29 + t23;
t15 = pkin(4) + t20;
t14 = t17 ^ 2;
t13 = t17 * qJ(4);
t12 = 0.2e1 * t17;
t4 = t8 * pkin(3) - t32;
t3 = t8 * qJ(5) + t6;
t2 = -t9 * qJ(5) + t5;
t1 = -t27 * t8 + t32;
t7 = [1, 0, 0, t25 ^ 2, t25 * t41, 0, 0, 0, pkin(1) * t41, -0.2e1 * pkin(1) * t25, t9 ^ 2, t8 * t43, 0, 0, 0, t8 * t42, t9 * t42, 0.2e1 * t4 * t8, 0.2e1 * t5 * t9 - 0.2e1 * t6 * t8, t4 * t43, t4 ^ 2 + t5 ^ 2 + t6 ^ 2, -0.2e1 * t1 * t8, 0.2e1 * t1 * t9, -0.2e1 * t2 * t9 + 0.2e1 * t3 * t8, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, t25, t26, 0, -t25 * pkin(6), -t26 * pkin(6), 0, 0, t9, -t8, 0, -t5, -t6, -t5, -t20 * t9 - t39, t6, t6 * t17 - t5 * t20, -t2, t3, t15 * t9 + t39, -t2 * t15 + t3 * t17; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t35, -0.2e1 * t23, 0.2e1 * t20, 0, t12, t20 ^ 2 + t14, 0.2e1 * t15, t12, 0, t15 ^ 2 + t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t8, 0, -t5, -t6, -t5, -pkin(3) * t9 - t37, t6, -t5 * pkin(3) + t6 * qJ(4), -t2, t3, t27 * t9 + t37, t3 * qJ(4) - t2 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t35, -t23, t33, 0, t18, t20 * pkin(3) + t13, 0.2e1 * pkin(4) + t33, t18, 0, t15 * t27 + t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t31, 0, t29, pkin(3) ^ 2 + t30, 0.2e1 * t27, t29, 0, t27 ^ 2 + t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, t5, 0, 0, -t9, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t20, -1, 0, 0, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3), -1, 0, 0, -t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t9, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t7;
