% Calculate minimal parameter regressor of joint inertia matrix for
% S4RRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% MM_reg [((4+1)*4/2)x20]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRRR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR2_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR2_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t25 = cos(qJ(3));
t18 = -t25 * pkin(3) - pkin(2);
t30 = cos(qJ(2)) * pkin(1);
t11 = t18 - t30;
t37 = 0.2e1 * t11;
t36 = 0.2e1 * t18;
t35 = 0.2e1 * t25;
t21 = sin(qJ(4));
t34 = t21 * pkin(3);
t33 = sin(qJ(2)) * pkin(1);
t24 = cos(qJ(4));
t32 = t24 * pkin(3);
t31 = t25 * pkin(6);
t17 = -pkin(2) - t30;
t29 = pkin(2) - t17;
t16 = pkin(6) + t33;
t28 = t25 * t16;
t27 = t11 + t18;
t22 = sin(qJ(3));
t20 = t22 ^ 2;
t19 = t25 * pkin(7);
t14 = t22 * t35;
t13 = t19 + t31;
t12 = (-pkin(6) - pkin(7)) * t22;
t10 = t21 * t25 + t24 * t22;
t9 = t21 * t22 - t24 * t25;
t8 = t10 ^ 2;
t7 = t19 + t28;
t6 = (-pkin(7) - t16) * t22;
t5 = -t21 * t12 - t24 * t13;
t4 = t24 * t12 - t21 * t13;
t3 = -0.2e1 * t10 * t9;
t2 = -t21 * t6 - t24 * t7;
t1 = -t21 * t7 + t24 * t6;
t15 = [1, 0, 0, 1, 0.2e1 * t30, -0.2e1 * t33, t20, t14, 0, 0, 0, -0.2e1 * t17 * t25, 0.2e1 * t17 * t22, t8, t3, 0, 0, 0, t9 * t37, t10 * t37; 0, 0, 0, 1, t30, -t33, t20, t14, 0, 0, 0, t29 * t25, -t29 * t22, t8, t3, 0, 0, 0, t27 * t9, t27 * t10; 0, 0, 0, 1, 0, 0, t20, t14, 0, 0, 0, pkin(2) * t35, -0.2e1 * pkin(2) * t22, t8, t3, 0, 0, 0, t9 * t36, t10 * t36; 0, 0, 0, 0, 0, 0, 0, 0, t22, t25, 0, -t22 * t16, -t28, 0, 0, t10, -t9, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, t22, t25, 0, -t22 * pkin(6), -t31, 0, 0, t10, -t9, 0, t4, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t32, -0.2e1 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, 0, t4, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t32, -t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t15;
