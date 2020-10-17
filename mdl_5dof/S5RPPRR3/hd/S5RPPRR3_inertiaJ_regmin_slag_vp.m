% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPPRR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:28:15
% EndTime: 2020-01-03 11:28:17
% DurationCPUTime: 0.21s
% Computational Cost: add. (153->29), mult. (291->56), div. (0->0), fcn. (354->8), ass. (0->30)
t21 = sin(pkin(9));
t23 = cos(pkin(9));
t26 = sin(qJ(4));
t28 = cos(qJ(4));
t12 = t26 * t21 - t28 * t23;
t24 = cos(pkin(8));
t18 = -t24 * pkin(1) - pkin(2);
t14 = -t23 * pkin(3) + t18;
t37 = 0.2e1 * t12 * pkin(4) + 0.2e1 * t14;
t36 = 0.2e1 * t14;
t25 = sin(qJ(5));
t35 = t25 * pkin(4);
t27 = cos(qJ(5));
t34 = t27 * pkin(4);
t22 = sin(pkin(8));
t16 = t22 * pkin(1) + qJ(3);
t33 = pkin(6) + t16;
t32 = t21 ^ 2 + t23 ^ 2;
t10 = t33 * t21;
t11 = t33 * t23;
t31 = -t28 * t10 - t26 * t11;
t30 = t26 * t10 - t28 * t11;
t13 = t28 * t21 + t26 * t23;
t6 = -t25 * t12 + t27 * t13;
t5 = t27 * t12 + t25 * t13;
t4 = -t12 * pkin(7) - t30;
t3 = -t13 * pkin(7) + t31;
t2 = -t25 * t3 - t27 * t4;
t1 = -t25 * t4 + t27 * t3;
t7 = [1, 0, 0, (t22 ^ 2 + t24 ^ 2) * pkin(1) ^ 2, -0.2e1 * t18 * t23, 0.2e1 * t18 * t21, 0.2e1 * t32 * t16, t32 * t16 ^ 2 + t18 ^ 2, t13 ^ 2, -0.2e1 * t13 * t12, 0, 0, 0, t12 * t36, t13 * t36, t6 ^ 2, -0.2e1 * t6 * t5, 0, 0, 0, t5 * t37, t6 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, t32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t23, t21, 0, t18, 0, 0, 0, 0, 0, t12, t13, 0, 0, 0, 0, 0, t5, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t12, 0, t31, t30, 0, 0, t6, -t5, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, -t13, 0, 0, 0, 0, 0, -t5, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t34, -0.2e1 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t34, -t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t7;
