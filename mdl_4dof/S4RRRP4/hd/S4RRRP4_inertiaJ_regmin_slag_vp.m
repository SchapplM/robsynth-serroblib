% Calculate minimal parameter regressor of joint inertia matrix for
% S4RRRP4
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
% MM_reg [((4+1)*4/2)x19]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRRP4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_inertiaJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t17 = cos(qJ(2));
t13 = -t17 * pkin(2) - pkin(1);
t22 = 0.2e1 * t13;
t21 = 0.2e1 * t17;
t20 = pkin(5) + pkin(6);
t15 = sin(qJ(3));
t19 = t15 * pkin(2);
t18 = cos(qJ(3));
t10 = t20 * t17;
t16 = sin(qJ(2));
t9 = t20 * t16;
t3 = -t15 * t10 - t18 * t9;
t4 = -t18 * t10 + t15 * t9;
t14 = t18 * pkin(2);
t12 = t14 + pkin(3);
t7 = t15 * t17 + t18 * t16;
t6 = t15 * t16 - t18 * t17;
t5 = t6 * pkin(3) + t13;
t2 = -t6 * qJ(4) - t4;
t1 = -t7 * qJ(4) + t3;
t8 = [1, 0, 0, t16 ^ 2, t16 * t21, 0, 0, 0, pkin(1) * t21, -0.2e1 * pkin(1) * t16, t7 ^ 2, -0.2e1 * t7 * t6, 0, 0, 0, t6 * t22, t7 * t22, -0.2e1 * t1 * t7 - 0.2e1 * t2 * t6, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t16, t17, 0, -t16 * pkin(5), -t17 * pkin(5), 0, 0, t7, -t6, 0, t3, t4, -t12 * t7 - t6 * t19, t1 * t12 + t2 * t19; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t14, -0.2e1 * t19, 0, t15 ^ 2 * pkin(2) ^ 2 + t12 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t6, 0, t3, t4, -pkin(3) * t7, t1 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t14, -t19, 0, t12 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(3) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t8;
