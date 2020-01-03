% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x23]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPP4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t37 = cos(qJ(2));
t29 = -t37 * pkin(2) - pkin(1);
t49 = 0.2e1 * t29;
t48 = 0.2e1 * t37;
t47 = pkin(6) + pkin(7);
t33 = cos(pkin(8));
t46 = t33 * pkin(3);
t34 = sin(qJ(3));
t45 = t34 * pkin(2);
t36 = cos(qJ(3));
t31 = t36 * pkin(2);
t28 = t31 + pkin(3);
t32 = sin(pkin(8));
t17 = t32 * t28 + t33 * t45;
t42 = t47 * t37;
t35 = sin(qJ(2));
t43 = t47 * t35;
t11 = -t34 * t42 - t36 * t43;
t20 = t34 * t37 + t36 * t35;
t39 = -t20 * qJ(4) + t11;
t12 = t34 * t43 - t36 * t42;
t19 = t34 * t35 - t36 * t37;
t8 = -t19 * qJ(4) - t12;
t4 = t32 * t8 - t33 * t39;
t6 = t32 * t39 + t33 * t8;
t44 = t4 ^ 2 + t6 ^ 2;
t41 = -t33 * t28 + t32 * t45;
t10 = -t32 * t19 + t33 * t20;
t9 = t33 * t19 + t32 * t20;
t40 = 0.2e1 * t4 * t10 - 0.2e1 * t6 * t9;
t13 = t19 * pkin(3) + t29;
t30 = t32 * pkin(3);
t25 = pkin(4) + t46;
t24 = t30 + qJ(5);
t15 = -pkin(4) + t41;
t14 = qJ(5) + t17;
t2 = t9 * pkin(4) - t10 * qJ(5) + t13;
t1 = [1, 0, 0, t35 ^ 2, t35 * t48, 0, 0, 0, pkin(1) * t48, -0.2e1 * pkin(1) * t35, t20 ^ 2, -0.2e1 * t20 * t19, 0, 0, 0, t19 * t49, t20 * t49, t40, t13 ^ 2 + t44, 0.2e1 * t2 * t9, t40, -0.2e1 * t2 * t10, t2 ^ 2 + t44; 0, 0, 0, 0, 0, t35, t37, 0, -t35 * pkin(6), -t37 * pkin(6), 0, 0, t20, -t19, 0, t11, t12, t10 * t41 - t17 * t9, t6 * t17 + t4 * t41, -t4, t15 * t10 - t14 * t9, t6, t6 * t14 + t4 * t15; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t31, -0.2e1 * t45, 0, t17 ^ 2 + t41 ^ 2, -0.2e1 * t15, 0, 0.2e1 * t14, t14 ^ 2 + t15 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, 0, t11, t12, (-t10 * t33 - t32 * t9) * pkin(3), (t32 * t6 - t33 * t4) * pkin(3), -t4, -t25 * t10 - t24 * t9, t6, t6 * t24 - t4 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t31, -t45, 0, (t17 * t32 - t33 * t41) * pkin(3), 0.2e1 * pkin(4) - t41 + t46, 0, t30 + 0.2e1 * qJ(5) + t17, t14 * t24 - t15 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, (t32 ^ 2 + t33 ^ 2) * pkin(3) ^ 2, 0.2e1 * t25, 0, 0.2e1 * t24, t24 ^ 2 + t25 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t9, 0, -t10, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t1;
