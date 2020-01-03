% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x21]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPPRR7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t14 = sin(pkin(8));
t5 = t14 * pkin(1) + qJ(3);
t30 = 0.2e1 * t5;
t18 = cos(qJ(5));
t29 = 0.2e1 * t18;
t19 = cos(qJ(4));
t13 = t19 ^ 2;
t15 = cos(pkin(8));
t7 = -t15 * pkin(1) - pkin(2);
t4 = -pkin(6) + t7;
t28 = t13 * t4;
t16 = sin(qJ(5));
t17 = sin(qJ(4));
t8 = t16 * t17;
t27 = t16 * t18;
t26 = t16 * t19;
t25 = t18 * t17;
t9 = t18 * t19;
t24 = t19 * t17;
t11 = t17 ^ 2;
t23 = -t11 - t13;
t22 = -0.2e1 * t24;
t21 = -pkin(4) * t19 - pkin(7) * t17;
t12 = t18 ^ 2;
t10 = t16 ^ 2;
t3 = t17 * pkin(4) - t19 * pkin(7) + t5;
t2 = t16 * t3 + t4 * t25;
t1 = t18 * t3 - t4 * t8;
t6 = [1, 0, 0, (t14 ^ 2 + t15 ^ 2) * pkin(1) ^ 2, 0.2e1 * t7, t30, t5 ^ 2 + t7 ^ 2, t13, t22, 0, 0, 0, t17 * t30, t19 * t30, t12 * t13, -0.2e1 * t13 * t27, t24 * t29, t16 * t22, t11, 0.2e1 * t1 * t17 - 0.2e1 * t16 * t28, -0.2e1 * t2 * t17 - 0.2e1 * t18 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23 * t16, t23 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t17, 0, t19 * t4, -t17 * t4, t16 * t9, (-t10 + t12) * t19, t8, t25, 0, t21 * t16 + t4 * t9, t21 * t18 - t4 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, -t19, 0, 0, 0, 0, 0, -t25, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t17, 0, 0, 0, 0, 0, t9, -t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t10, 0.2e1 * t27, 0, 0, 0, pkin(4) * t29, -0.2e1 * pkin(4) * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t26, t17, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, -t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t18, 0, -t16 * pkin(7), -t18 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t6;
