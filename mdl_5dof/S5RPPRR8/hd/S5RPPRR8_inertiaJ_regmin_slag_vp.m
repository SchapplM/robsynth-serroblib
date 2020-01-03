% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPPRR8
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
% MM_reg [((5+1)*5/2)x19]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPPRR8_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t16 = sin(qJ(5));
t26 = -0.2e1 * t16;
t18 = cos(qJ(5));
t25 = 0.2e1 * t18;
t24 = -pkin(1) - pkin(2);
t14 = sin(pkin(8));
t15 = cos(pkin(8));
t10 = t15 * qJ(2) + t14 * t24;
t17 = sin(qJ(4));
t19 = cos(qJ(4));
t8 = t14 * qJ(2) - t15 * t24;
t7 = -pkin(3) - t8;
t3 = t17 * t10 - t19 * t7;
t1 = pkin(4) + t3;
t23 = pkin(4) + t1;
t5 = t17 * t14 - t19 * t15;
t22 = t5 * t16;
t21 = t5 * t18;
t20 = t16 * t18;
t4 = t19 * t10 + t17 * t7;
t13 = t16 ^ 2;
t11 = 0.2e1 * t20;
t6 = t19 * t14 + t17 * t15;
t2 = -pkin(7) + t4;
t9 = [1, 0, 0, 2 * pkin(1), 0.2e1 * qJ(2), (pkin(1) ^ 2) + qJ(2) ^ 2, 0.2e1 * t8, 0.2e1 * t10, t10 ^ 2 + t8 ^ 2, 1, 0.2e1 * t3, 0.2e1 * t4, t13, t11, 0, 0, 0, t1 * t25, t1 * t26; 0, 0, 0, -1, 0, -pkin(1), -t15, t14, t10 * t14 - t8 * t15, 0, t5, t6, 0, 0, 0, 0, 0, t21, -t22; 0, 0, 0, 0, 0, 1, 0, 0, t14 ^ 2 + t15 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t3, -t4, -t13, -0.2e1 * t20, 0, 0, 0, -t23 * t18, t23 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6, 0, 0, 0, 0, 0, -t21, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t13, t11, 0, 0, 0, pkin(4) * t25, pkin(4) * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t18, 0, -t16 * t2, -t18 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16 * t6, -t18 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t18, 0, -t16 * pkin(7), -t18 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t9;
