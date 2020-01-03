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
% MM_reg [((5+1)*5/2)x23]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
t24 = sin(pkin(8));
t25 = cos(pkin(8));
t27 = sin(qJ(2));
t28 = cos(qJ(2));
t16 = t24 * t28 + t25 * t27;
t26 = sin(qJ(4));
t33 = t24 * t27 - t25 * t28;
t38 = cos(qJ(4));
t7 = t38 * t16 - t26 * t33;
t42 = -0.2e1 * t7;
t23 = -t28 * pkin(2) - pkin(1);
t10 = t33 * pkin(3) + t23;
t41 = 0.2e1 * t10;
t40 = 0.2e1 * t28;
t39 = t24 * pkin(2);
t37 = -qJ(3) - pkin(6);
t35 = t37 * t27;
t36 = t37 * t28;
t9 = t24 * t35 - t25 * t36;
t22 = t25 * pkin(2) + pkin(3);
t14 = -t26 * t22 - t38 * t39;
t34 = -t38 * t22 + t26 * t39;
t8 = t24 * t36 + t25 * t35;
t32 = -t16 * pkin(7) + t8;
t30 = 2 * pkin(4);
t29 = 2 * qJ(5);
t12 = -pkin(4) + t34;
t11 = qJ(5) - t14;
t6 = t26 * t16 + t38 * t33;
t5 = -t33 * pkin(7) + t9;
t3 = t26 * t32 + t38 * t5;
t2 = t26 * t5 - t38 * t32;
t1 = t6 * pkin(4) - t7 * qJ(5) + t10;
t4 = [1, 0, 0, t27 ^ 2, t27 * t40, 0, 0, 0, pkin(1) * t40, -0.2e1 * pkin(1) * t27, -0.2e1 * t8 * t16 - 0.2e1 * t9 * t33, t23 ^ 2 + t8 ^ 2 + t9 ^ 2, t7 ^ 2, t6 * t42, 0, 0, 0, t6 * t41, t7 * t41, 0.2e1 * t1 * t6, 0.2e1 * t2 * t7 - 0.2e1 * t3 * t6, t1 * t42, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, t27, t28, 0, -t27 * pkin(6), -t28 * pkin(6), (-t25 * t16 - t24 * t33) * pkin(2), (t24 * t9 + t25 * t8) * pkin(2), 0, 0, t7, -t6, 0, -t2, -t3, -t2, -t11 * t6 + t12 * t7, t3, t3 * t11 + t2 * t12; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, (t24 ^ 2 + t25 ^ 2) * pkin(2) ^ 2, 0, 0, 0, 0, 1, -0.2e1 * t34, 0.2e1 * t14, -0.2e1 * t12, 0, 0.2e1 * t11, t11 ^ 2 + t12 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, 0, 0, 0, 0, t6, t7, t6, 0, -t7, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t6, 0, -t2, -t3, -t2, -pkin(4) * t7 - t6 * qJ(5), t3, -t2 * pkin(4) + t3 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -t34, t14, t30 - t34, 0, t29 - t14, -t12 * pkin(4) + t11 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t30, 0, t29, pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t4;
