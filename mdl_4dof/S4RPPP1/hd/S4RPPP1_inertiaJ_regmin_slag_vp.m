% Calculate minimal parameter regressor of joint inertia matrix for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
% 
% Output:
% MM_reg [((4+1)*4/2)x15]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RPPP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_inertiaJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t16 = sin(pkin(4));
t24 = 0.2e1 * t16;
t17 = cos(pkin(6));
t23 = pkin(1) * t17;
t18 = cos(pkin(4));
t22 = pkin(1) * t18;
t15 = sin(pkin(6));
t12 = t16 * t15;
t13 = t16 * t17;
t21 = qJ(2) * t16;
t8 = t15 * t22 + t17 * t21;
t20 = -pkin(2) - t23;
t19 = -qJ(3) * t15 - pkin(1);
t4 = -t18 * qJ(3) - t8;
t14 = t16 ^ 2;
t9 = t15 * t21;
t7 = t17 * t22 - t9;
t6 = (-pkin(2) * t17 + t19) * t16;
t5 = t20 * t18 + t9;
t3 = ((-pkin(2) - qJ(4)) * t17 + t19) * t16;
t2 = pkin(3) * t13 - t4;
t1 = pkin(3) * t12 + t9 + (-qJ(4) + t20) * t18;
t10 = [1, 0, 0, 0.2e1 * t14 * t23 + 0.2e1 * t7 * t18, -0.2e1 * t14 * pkin(1) * t15 - 0.2e1 * t8 * t18 (-t15 * t7 + t17 * t8) * t24, t14 * pkin(1) ^ 2 + t7 ^ 2 + t8 ^ 2 (t15 * t5 - t17 * t4) * t24, 0.2e1 * t6 * t13 + 0.2e1 * t5 * t18, -0.2e1 * t6 * t12 - 0.2e1 * t4 * t18, t4 ^ 2 + t5 ^ 2 + t6 ^ 2 (t1 * t15 + t17 * t2) * t24, -0.2e1 * t3 * t12 + 0.2e1 * t2 * t18, -0.2e1 * t1 * t18 - 0.2e1 * t3 * t13, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, -t13, t12, 0, -t16 * pkin(1), 0, t13, -t12, t6, 0, -t12, -t13, t3; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, t12, t18, 0, t5, t12, 0, -t18, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t18, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t10;
