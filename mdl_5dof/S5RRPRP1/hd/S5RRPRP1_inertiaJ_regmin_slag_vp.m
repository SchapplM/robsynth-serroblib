% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x16]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t25 = sin(qJ(4));
t37 = 0.2e1 * t25;
t27 = cos(qJ(4));
t36 = -0.2e1 * t27;
t35 = pkin(4) * t25;
t34 = sin(qJ(2)) * pkin(1);
t33 = t27 * pkin(4);
t21 = cos(qJ(2)) * pkin(1);
t19 = t21 + pkin(2);
t23 = sin(pkin(8));
t24 = cos(pkin(8));
t9 = t23 * t19 + t24 * t34;
t6 = pkin(7) + t9;
t32 = t27 * t6;
t18 = -t24 * pkin(2) - pkin(3);
t8 = t24 * t19 - t23 * t34;
t5 = -pkin(3) - t8;
t31 = t18 + t5;
t17 = t23 * pkin(2) + pkin(7);
t30 = t27 * t17;
t22 = t25 ^ 2;
t20 = t27 * qJ(5);
t16 = t27 * t37;
t12 = t18 - t33;
t11 = t20 + t30;
t10 = (-qJ(5) - t17) * t25;
t7 = t11 * t27;
t4 = t5 - t33;
t3 = t20 + t32;
t2 = (-qJ(5) - t6) * t25;
t1 = t3 * t27;
t13 = [1, 0, 0, 1, 0.2e1 * t21, -0.2e1 * t34, t8 ^ 2 + t9 ^ 2, t22, t16, 0, 0, 0, t5 * t36, t5 * t37, -0.2e1 * t2 * t25 + 0.2e1 * t1, t2 ^ 2 + t3 ^ 2 + t4 ^ 2; 0, 0, 0, 1, t21, -t34, (t23 * t9 + t24 * t8) * pkin(2), t22, t16, 0, 0, 0, -t31 * t27, t31 * t25, t1 + t7 + (-t10 - t2) * t25, t2 * t10 + t3 * t11 + t4 * t12; 0, 0, 0, 1, 0, 0, (t23 ^ 2 + t24 ^ 2) * pkin(2) ^ 2, t22, t16, 0, 0, 0, t18 * t36, t18 * t37, -0.2e1 * t10 * t25 + 0.2e1 * t7, t10 ^ 2 + t11 ^ 2 + t12 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2 * t27 + t3 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 * t27 + t11 * t25; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t27 ^ 2 + t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t27, 0, -t25 * t6, -t32, -t35, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t27, 0, -t25 * t17, -t30, -t35, t10 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t25, 0, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t13;
