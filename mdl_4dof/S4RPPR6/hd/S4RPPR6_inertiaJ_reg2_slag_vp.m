% Calculate inertial parameters regressor of joint inertia matrix for
% S4RPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% MM_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RPPR6_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_inertiaJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t17 = sin(pkin(6));
t15 = t17 ^ 2;
t18 = cos(pkin(6));
t16 = t18 ^ 2;
t27 = t15 + t16;
t22 = t17 * qJ(3) + pkin(1);
t3 = (pkin(2) + pkin(3)) * t18 + t22;
t26 = 0.2e1 * t3;
t25 = -0.2e1 * t17;
t24 = t17 * t18;
t23 = t27 * qJ(2) ^ 2;
t20 = cos(qJ(4));
t19 = sin(qJ(4));
t13 = t17 * qJ(2);
t10 = (-pkin(5) + qJ(2)) * t18;
t9 = -t17 * pkin(5) + t13;
t8 = -t18 * pkin(2) - t22;
t7 = 0.2e1 * t27 * qJ(2);
t6 = t17 * t20 - t18 * t19;
t4 = t17 * t19 + t18 * t20;
t2 = t20 * t10 + t19 * t9;
t1 = -t19 * t10 + t20 * t9;
t5 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t15, 0.2e1 * t24, 0, t16, 0, 0, 0.2e1 * pkin(1) * t18, pkin(1) * t25, t7, pkin(1) ^ 2 + t23, t15, 0, -0.2e1 * t24, 0, 0, t16, -0.2e1 * t8 * t18, t7, t8 * t25, t8 ^ 2 + t23, t6 ^ 2, -0.2e1 * t6 * t4, 0, t4 ^ 2, 0, 0, t4 * t26, t6 * t26, -0.2e1 * t1 * t6 - 0.2e1 * t2 * t4, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, t17, 0, -pkin(1), 0, 0, 0, 0, 0, 0, -t18, 0, -t17, t8, 0, 0, 0, 0, 0, 0, -t4, -t6, 0, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, t13, 0, 0, 0, 0, 0, 0, 0, 0, -t19 * t4 - t20 * t6, t1 * t20 + t2 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19 ^ 2 + t20 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, -t4, 0, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t5;
