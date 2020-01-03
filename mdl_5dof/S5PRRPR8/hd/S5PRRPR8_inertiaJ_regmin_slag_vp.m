% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x15]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRPR8_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t22 = sin(qJ(5));
t33 = 0.2e1 * t22;
t25 = cos(qJ(5));
t32 = -0.2e1 * t25;
t23 = sin(qJ(3));
t24 = sin(qJ(2));
t26 = cos(qJ(3));
t27 = cos(qJ(2));
t10 = -t23 * t24 + t26 * t27;
t11 = t23 * t27 + t26 * t24;
t20 = sin(pkin(9));
t21 = cos(pkin(9));
t2 = -t21 * t10 + t20 * t11;
t31 = t2 * t25;
t30 = t23 * pkin(2);
t16 = -t21 * pkin(3) - pkin(4);
t18 = t26 * pkin(2);
t17 = t18 + pkin(3);
t7 = t21 * t17 - t20 * t30;
t5 = -pkin(4) - t7;
t29 = t16 + t5;
t8 = t20 * t17 + t21 * t30;
t19 = t22 ^ 2;
t15 = t20 * pkin(3) + pkin(7);
t14 = t25 * t33;
t6 = pkin(7) + t8;
t4 = t20 * t10 + t21 * t11;
t1 = t2 * t22;
t3 = [1, 0, 0, 0, 0, 0, 0, t2 ^ 2 + t4 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, 0, t27, -t24, 0, t10, -t11, -t2 * t7 + t4 * t8, 0, 0, 0, 0, 0, -t31, t1; 0, 1, 0, 0, 1, 0.2e1 * t18, -0.2e1 * t30, t7 ^ 2 + t8 ^ 2, t19, t14, 0, 0, 0, t5 * t32, t5 * t33; 0, 0, 0, 0, 0, t10, -t11, (-t2 * t21 + t20 * t4) * pkin(3), 0, 0, 0, 0, 0, -t31, t1; 0, 0, 0, 0, 1, t18, -t30, (t20 * t8 + t21 * t7) * pkin(3), t19, t14, 0, 0, 0, -t29 * t25, t29 * t22; 0, 0, 0, 0, 1, 0, 0, (t20 ^ 2 + t21 ^ 2) * pkin(3) ^ 2, t19, t14, 0, 0, 0, t16 * t32, t16 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22 * t4, -t25 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t25, 0, -t22 * t6, -t25 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t25, 0, -t22 * t15, -t25 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t3;
