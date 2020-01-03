% Calculate minimal parameter regressor of joint inertia matrix for
% S4PRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% MM_reg [((4+1)*4/2)x18]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4PRRR7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t14 = sin(qJ(3));
t30 = -0.2e1 * t14;
t17 = cos(qJ(3));
t29 = 0.2e1 * t17;
t16 = cos(qJ(4));
t28 = pkin(3) * t16;
t13 = sin(qJ(4));
t27 = pkin(6) * t13;
t11 = sin(pkin(4));
t26 = t11 * sin(qJ(2));
t25 = t11 * cos(qJ(2));
t24 = t13 * t14;
t23 = t13 * t16;
t22 = t13 * t17;
t21 = t16 * t14;
t20 = t16 * t17;
t19 = t14 * t29;
t12 = cos(pkin(4));
t10 = t16 ^ 2;
t9 = t14 ^ 2;
t8 = t13 ^ 2;
t7 = -t17 * pkin(3) - t14 * pkin(7) - pkin(2);
t6 = t12 * t14 + t17 * t26;
t5 = -t12 * t17 + t14 * t26;
t4 = pkin(6) * t20 + t13 * t7;
t3 = -pkin(6) * t22 + t16 * t7;
t2 = -t13 * t25 + t6 * t16;
t1 = -t6 * t13 - t16 * t25;
t15 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t25, -t26, 0, 0, 0, 0, 0, t17 * t25, -t14 * t25, 0, 0, 0, 0, 0, -t1 * t17 + t5 * t24, t2 * t17 + t5 * t21; 0, 1, 0, 0, t9, t19, 0, 0, 0, pkin(2) * t29, pkin(2) * t30, t10 * t9, -0.2e1 * t9 * t23, t20 * t30, t13 * t19, t17 ^ 2, -0.2e1 * t3 * t17 + 0.2e1 * t9 * t27, 0.2e1 * t9 * pkin(6) * t16 + 0.2e1 * t4 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6, 0, 0, 0, 0, 0, -t5 * t16, t5 * t13; 0, 0, 0, 0, 0, 0, t14, t17, 0, -t14 * pkin(6), -t17 * pkin(6), t13 * t21, (t10 - t8) * t14, -t22, -t20, 0, -pkin(6) * t21 + (-pkin(3) * t14 + pkin(7) * t17) * t13, pkin(7) * t20 + (t27 - t28) * t14; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t8, 0.2e1 * t23, 0, 0, 0, 0.2e1 * t28, -0.2e1 * pkin(3) * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t24, -t17, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t16, 0, -t13 * pkin(7), -t16 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t15;
