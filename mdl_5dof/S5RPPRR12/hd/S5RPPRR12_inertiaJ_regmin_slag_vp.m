% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x24]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPPRR12_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:07:21
% EndTime: 2019-12-31 18:07:22
% DurationCPUTime: 0.26s
% Computational Cost: add. (153->40), mult. (295->73), div. (0->0), fcn. (350->6), ass. (0->37)
t24 = sin(pkin(8));
t25 = cos(pkin(8));
t28 = sin(qJ(4));
t37 = cos(qJ(4));
t11 = t37 * t24 + t28 * t25;
t7 = t11 ^ 2;
t10 = t28 * t24 - t37 * t25;
t8 = t10 ^ 2;
t42 = -t7 - t8;
t41 = -0.2e1 * t10;
t40 = 0.2e1 * t11;
t39 = 2 * qJ(2);
t26 = -pkin(1) - qJ(3);
t38 = -pkin(6) + t26;
t27 = sin(qJ(5));
t36 = t10 * t27;
t29 = cos(qJ(5));
t35 = t10 * t29;
t34 = t27 * t11;
t33 = t27 * t29;
t6 = t29 * t11;
t15 = t24 ^ 2 + t25 ^ 2;
t16 = t24 * pkin(3) + qJ(2);
t32 = t10 * t40;
t31 = pkin(4) * t10 - pkin(7) * t11;
t30 = qJ(2) ^ 2;
t23 = t29 ^ 2;
t22 = t27 ^ 2;
t14 = t38 * t25;
t13 = t38 * t24;
t9 = t15 * t26;
t5 = t37 * t13 + t28 * t14;
t4 = t28 * t13 - t37 * t14;
t3 = t11 * pkin(4) + t10 * pkin(7) + t16;
t2 = t27 * t3 + t29 * t5;
t1 = -t27 * t5 + t29 * t3;
t12 = [1, 0, 0, -2 * pkin(1), t39, pkin(1) ^ 2 + t30, t24 * t39, t25 * t39, -0.2e1 * t9, t15 * t26 ^ 2 + t30, t8, t32, 0, 0, 0, t16 * t40, t16 * t41, t23 * t8, -0.2e1 * t8 * t33, t6 * t41, t27 * t32, t7, 0.2e1 * t1 * t11 - 0.2e1 * t4 * t36, -0.2e1 * t2 * t11 - 0.2e1 * t4 * t35; 0, 0, 0, 1, 0, -pkin(1), 0, 0, -t15, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42 * t27, t42 * t29; 0, 0, 0, 0, 0, 1, 0, 0, 0, t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t24, t25, 0, qJ(2), 0, 0, 0, 0, 0, t11, -t10, 0, 0, 0, 0, 0, t6, -t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, 0, -t4, -t5, -t10 * t33, (t22 - t23) * t10, t34, t6, 0, t31 * t27 - t4 * t29, t4 * t27 + t31 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, 0, 0, 0, 0, 0, -t35, t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t22, 0.2e1 * t33, 0, 0, 0, 0.2e1 * pkin(4) * t29, -0.2e1 * pkin(4) * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, t36, t11, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t29, 0, -t27 * pkin(7), -t29 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t12;
