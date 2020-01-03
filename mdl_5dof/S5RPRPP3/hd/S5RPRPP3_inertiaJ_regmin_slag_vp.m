% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPP3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t17 = sin(pkin(7));
t18 = cos(pkin(7));
t21 = sin(qJ(3));
t29 = cos(qJ(3));
t10 = t29 * t17 + t21 * t18;
t31 = -0.2e1 * t10;
t14 = -t18 * pkin(2) - pkin(1);
t30 = 0.2e1 * t14;
t9 = t21 * t17 - t29 * t18;
t28 = qJ(4) * t9;
t19 = pkin(3) + qJ(5);
t27 = pkin(6) + qJ(2);
t26 = t17 ^ 2 + t18 ^ 2;
t11 = t27 * t17;
t12 = t27 * t18;
t5 = t29 * t11 + t21 * t12;
t25 = -t10 * qJ(4) + t14;
t6 = -t21 * t11 + t29 * t12;
t23 = qJ(4) ^ 2;
t22 = 0.2e1 * qJ(4);
t4 = t9 * pkin(3) + t25;
t3 = -t9 * pkin(4) + t6;
t2 = t10 * pkin(4) + t5;
t1 = t19 * t9 + t25;
t7 = [1, 0, 0, 0.2e1 * pkin(1) * t18, -0.2e1 * pkin(1) * t17, 0.2e1 * t26 * qJ(2), t26 * qJ(2) ^ 2 + pkin(1) ^ 2, t10 ^ 2, t9 * t31, 0, 0, 0, t9 * t30, t10 * t30, 0.2e1 * t5 * t10 - 0.2e1 * t6 * t9, -0.2e1 * t4 * t9, t4 * t31, t4 ^ 2 + t5 ^ 2 + t6 ^ 2, 0.2e1 * t2 * t10 - 0.2e1 * t3 * t9, t1 * t31, 0.2e1 * t1 * t9, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, -t18, t17, 0, -pkin(1), 0, 0, 0, 0, 0, t9, t10, 0, -t9, -t10, t4, 0, -t10, t9, t1; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, 0, -t5, -t6, -pkin(3) * t10 - t28, t5, t6, -t5 * pkin(3) + t6 * qJ(4), -t19 * t10 - t28, t3, -t2, t3 * qJ(4) - t2 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(3), t22, pkin(3) ^ 2 + t23, 0, t22, 0.2e1 * t19, t19 ^ 2 + t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, 0, t5, t10, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(3), 0, 0, -1, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, qJ(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t7;
