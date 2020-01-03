% Calculate minimal parameter regressor of joint inertia matrix for
% S4RRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% MM_reg [((4+1)*4/2)x21]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRRP5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_inertiaJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t14 = sin(qJ(3));
t16 = cos(qJ(2));
t15 = sin(qJ(2));
t22 = cos(qJ(3));
t19 = t22 * t15;
t5 = t14 * t16 + t19;
t26 = -0.2e1 * t5;
t12 = -t16 * pkin(2) - pkin(1);
t25 = 0.2e1 * t12;
t24 = 0.2e1 * t16;
t23 = -pkin(6) - pkin(5);
t21 = t14 * t15;
t20 = t22 * pkin(2);
t18 = 2 * pkin(3);
t17 = 2 * qJ(4);
t13 = t14 * pkin(2);
t10 = t20 + pkin(3);
t8 = t13 + qJ(4);
t7 = t23 * t16;
t4 = -t22 * t16 + t21;
t3 = t23 * t21 - t22 * t7;
t2 = -t14 * t7 - t23 * t19;
t1 = t4 * pkin(3) - t5 * qJ(4) + t12;
t6 = [1, 0, 0, t15 ^ 2, t15 * t24, 0, 0, 0, pkin(1) * t24, -0.2e1 * pkin(1) * t15, t5 ^ 2, t4 * t26, 0, 0, 0, t4 * t25, t5 * t25, 0.2e1 * t1 * t4, 0.2e1 * t2 * t5 - 0.2e1 * t3 * t4, t1 * t26, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, t15, t16, 0, -t15 * pkin(5), -t16 * pkin(5), 0, 0, t5, -t4, 0, -t2, -t3, -t2, -t10 * t5 - t8 * t4, t3, -t2 * t10 + t3 * t8; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t20, -0.2e1 * t13, 0.2e1 * t10, 0, 0.2e1 * t8, t10 ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t4, 0, -t2, -t3, -t2, -pkin(3) * t5 - t4 * qJ(4), t3, -t2 * pkin(3) + t3 * qJ(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t20, -t13, t18 + t20, 0, t17 + t13, t10 * pkin(3) + t8 * qJ(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t18, 0, t17, pkin(3) ^ 2 + qJ(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t6;
