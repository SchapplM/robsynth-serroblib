% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x17]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t25 = sin(qJ(3));
t37 = 0.2e1 * t25;
t26 = cos(qJ(3));
t22 = sin(pkin(7));
t34 = t22 * pkin(1) + pkin(6);
t31 = qJ(4) + t34;
t10 = t31 * t26;
t21 = sin(pkin(8));
t23 = cos(pkin(8));
t29 = t31 * t25;
t5 = t21 * t10 + t23 * t29;
t7 = t23 * t10 - t21 * t29;
t36 = t5 ^ 2 + t7 ^ 2;
t12 = t21 * t25 - t23 * t26;
t14 = t21 * t26 + t23 * t25;
t35 = t12 ^ 2 + t14 ^ 2;
t24 = cos(pkin(7));
t20 = -t24 * pkin(1) - pkin(2);
t33 = t5 * t12 + t7 * t14;
t30 = -0.2e1 * t7 * t12 + 0.2e1 * t5 * t14;
t15 = -t26 * pkin(3) + t20;
t18 = t23 * pkin(3) + pkin(4);
t16 = t21 * pkin(3) + qJ(5);
t3 = t12 * pkin(4) - t14 * qJ(5) + t15;
t1 = [1, 0, 0, (t22 ^ 2 + t24 ^ 2) * pkin(1) ^ 2, t25 ^ 2, t26 * t37, 0, 0, 0, -0.2e1 * t20 * t26, t20 * t37, t30, t15 ^ 2 + t36, 0.2e1 * t3 * t12, t30, -0.2e1 * t3 * t14, t3 ^ 2 + t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, 0, 0, t33; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, 0, 0, t35; 0, 0, 0, 0, 0, 0, t25, t26, 0, -t25 * t34, -t26 * t34, (-t12 * t21 - t14 * t23) * pkin(3), (t21 * t7 - t23 * t5) * pkin(3), -t5, -t16 * t12 - t18 * t14, t7, t7 * t16 - t5 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t25, 0, (-t12 * t23 + t14 * t21) * pkin(3), -t12, 0, t14, -t12 * t18 + t14 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, (t21 ^ 2 + t23 ^ 2) * pkin(3) ^ 2, 0.2e1 * t18, 0, 0.2e1 * t16, t16 ^ 2 + t18 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, t12, 0, -t14, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t1;
