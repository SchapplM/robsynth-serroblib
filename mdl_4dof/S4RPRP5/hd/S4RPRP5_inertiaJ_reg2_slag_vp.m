% Calculate inertial parameters regressor of joint inertia matrix for
% S4RPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% MM_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RPRP5_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_inertiaJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t17 = sin(pkin(6));
t18 = cos(pkin(6));
t19 = sin(qJ(3));
t26 = cos(qJ(3));
t9 = t19 * t17 - t26 * t18;
t30 = t9 ^ 2;
t14 = -t18 * pkin(2) - pkin(1);
t29 = 0.2e1 * t14;
t28 = 0.2e1 * t18;
t11 = t26 * t17 + t19 * t18;
t27 = t11 * t9;
t25 = pkin(5) + qJ(2);
t15 = t17 ^ 2;
t16 = t18 ^ 2;
t24 = t15 + t16;
t12 = t25 * t18;
t22 = t25 * t17;
t4 = t19 * t12 + t26 * t22;
t6 = t26 * t12 - t19 * t22;
t23 = t4 ^ 2 + t6 ^ 2;
t21 = 0.2e1 * t4 * t11 - 0.2e1 * t6 * t9;
t7 = t11 ^ 2;
t2 = t9 * pkin(3) - t11 * qJ(4) + t14;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t15, t17 * t28, 0, t16, 0, 0, pkin(1) * t28, -0.2e1 * pkin(1) * t17, 0.2e1 * t24 * qJ(2), t24 * qJ(2) ^ 2 + pkin(1) ^ 2, t7, -0.2e1 * t27, 0, t30, 0, 0, t9 * t29, t11 * t29, t21, t14 ^ 2 + t23, t7, 0, 0.2e1 * t27, 0, 0, t30, 0.2e1 * t2 * t9, t21, -0.2e1 * t2 * t11, t2 ^ 2 + t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, t17, 0, -pkin(1), 0, 0, 0, 0, 0, 0, t9, t11, 0, t14, 0, 0, 0, 0, 0, 0, t9, 0, -t11, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, -t9, 0, -t4, -t6, 0, 0, 0, t11, 0, 0, t9, 0, -t4, -pkin(3) * t11 - t9 * qJ(4), t6, -t4 * pkin(3) + t6 * qJ(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, 0.2e1 * qJ(4), pkin(3) ^ 2 + qJ(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t1;
