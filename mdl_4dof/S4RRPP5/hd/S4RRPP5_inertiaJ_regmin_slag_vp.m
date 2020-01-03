% Calculate minimal parameter regressor of joint inertia matrix for
% S4RRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% MM_reg [((4+1)*4/2)x18]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRPP5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_inertiaJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t12 = sin(qJ(2));
t13 = cos(qJ(2));
t5 = t13 * qJ(3);
t22 = -t12 * pkin(2) + t5;
t21 = -0.2e1 * t12;
t20 = 0.2e1 * t13;
t8 = t12 ^ 2;
t19 = t13 ^ 2 + t8;
t10 = pkin(2) + qJ(4);
t17 = -t12 * qJ(3) - pkin(1);
t15 = qJ(3) ^ 2;
t14 = 0.2e1 * qJ(3);
t7 = t13 * pkin(5);
t6 = t12 * pkin(5);
t4 = t13 * pkin(3) + t7;
t3 = t12 * pkin(3) + t6;
t2 = -t13 * pkin(2) + t17;
t1 = -t10 * t13 + t17;
t9 = [1, 0, 0, t8, t12 * t20, 0, 0, 0, pkin(1) * t20, pkin(1) * t21, 0.2e1 * t19 * pkin(5), t2 * t20, t2 * t21, t19 * pkin(5) ^ 2 + t2 ^ 2, 0.2e1 * t3 * t12 + 0.2e1 * t4 * t13, t1 * t21, -0.2e1 * t1 * t13, t1 ^ 2 + t3 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, t12, t13, 0, -t6, -t7, t22, t6, t7, t22 * pkin(5), -t10 * t12 + t5, t4, -t3, t4 * qJ(3) - t3 * t10; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(2), t14, pkin(2) ^ 2 + t15, 0, t14, 0.2e1 * t10, t10 ^ 2 + t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, 0, t6, t12, 0, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, -1, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, qJ(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t9;
