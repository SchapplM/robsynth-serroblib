% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x23]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPPRR11_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t11 = (pkin(1) + qJ(3));
t30 = 2 * t11;
t15 = cos(qJ(5));
t29 = 0.2e1 * t15;
t14 = sin(qJ(4));
t7 = t14 ^ 2;
t16 = cos(qJ(4));
t9 = t16 ^ 2;
t28 = -t7 - t9;
t10 = -pkin(6) + qJ(2);
t27 = t10 * t9;
t13 = sin(qJ(5));
t26 = t13 * t15;
t25 = t13 * t16;
t24 = t14 * t10;
t23 = t15 * t14;
t5 = t15 * t16;
t22 = t16 * t10;
t21 = t16 * t14;
t20 = -0.2e1 * t21;
t19 = -pkin(4) * t16 - pkin(7) * t14;
t18 = (qJ(2) ^ 2);
t17 = 2 * qJ(2);
t8 = t15 ^ 2;
t6 = t13 ^ 2;
t4 = t13 * t14;
t3 = t14 * pkin(4) - t16 * pkin(7) + t11;
t2 = t10 * t23 + t13 * t3;
t1 = -t13 * t24 + t15 * t3;
t12 = [1, 0, 0, -2 * pkin(1), t17, pkin(1) ^ 2 + t18, t17, t30, t11 ^ 2 + t18, t9, t20, 0, 0, 0, t14 * t30, t16 * t30, t8 * t9, -0.2e1 * t9 * t26, t21 * t29, t13 * t20, t7, 0.2e1 * t1 * t14 - 0.2e1 * t13 * t27, -0.2e1 * t2 * t14 - 0.2e1 * t15 * t27; 0, 0, 0, 1, 0, -pkin(1), 0, -1, -t11, 0, 0, 0, 0, 0, -t14, -t16, 0, 0, 0, 0, 0, -t23, t4; 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 1, 0, qJ(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28 * t13, t28 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t14, 0, t22, -t24, t13 * t5, (-t6 + t8) * t16, t4, t23, 0, t19 * t13 + t15 * t22, -t13 * t22 + t19 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t14, 0, 0, 0, 0, 0, t5, -t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t6, 0.2e1 * t26, 0, 0, 0, pkin(4) * t29, -0.2e1 * pkin(4) * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t25, t14, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, -t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t15, 0, -t13 * pkin(7), -t15 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t12;
