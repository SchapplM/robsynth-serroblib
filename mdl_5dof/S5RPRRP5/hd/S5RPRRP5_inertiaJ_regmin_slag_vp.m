% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x18]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRP5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:40:57
% EndTime: 2019-12-31 18:40:58
% DurationCPUTime: 0.19s
% Computational Cost: add. (139->28), mult. (254->57), div. (0->0), fcn. (214->6), ass. (0->31)
t20 = sin(qJ(4));
t16 = t20 ^ 2;
t22 = cos(qJ(4));
t27 = t22 ^ 2 + t16;
t19 = cos(pkin(8));
t13 = t19 * pkin(1) + pkin(2);
t21 = sin(qJ(3));
t23 = cos(qJ(3));
t18 = sin(pkin(8));
t33 = pkin(1) * t18;
t8 = -t21 * t13 - t23 * t33;
t6 = pkin(7) - t8;
t34 = t27 * t6;
t39 = -0.2e1 * t20;
t38 = -0.2e1 * t22;
t37 = 0.2e1 * t22;
t7 = t23 * t13 - t21 * t33;
t5 = -pkin(3) - t7;
t36 = pkin(3) - t5;
t26 = t22 * pkin(4) + t20 * qJ(5);
t9 = -pkin(3) - t26;
t1 = -t7 + t9;
t35 = -t1 - t9;
t32 = t20 * pkin(7);
t31 = t20 * t6;
t30 = t22 * pkin(7);
t29 = t22 * t6;
t28 = t27 * pkin(7);
t10 = -t20 * pkin(4) + t22 * qJ(5);
t12 = t20 * t37;
t2 = [1, 0, 0, (t18 ^ 2 + t19 ^ 2) * pkin(1) ^ 2, 1, 0.2e1 * t7, 0.2e1 * t8, t16, t12, 0, 0, 0, t5 * t38, 0.2e1 * t5 * t20, t1 * t38, 0.2e1 * t34, t1 * t39, t27 * t6 ^ 2 + t1 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27; 0, 0, 0, 0, 1, t7, t8, t16, t12, 0, 0, 0, t36 * t22, -t36 * t20, t35 * t22, t28 + t34, t35 * t20, pkin(7) * t34 + t1 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0, 0, t16, t12, 0, 0, 0, pkin(3) * t37, pkin(3) * t39, t9 * t38, 0.2e1 * t28, t9 * t39, t27 * pkin(7) ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, t22, 0, -t31, -t29, -t31, t10, t29, t10 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t20, t22, 0, t20, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, t22, 0, -t32, -t30, -t32, t10, t30, t10 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t2;
