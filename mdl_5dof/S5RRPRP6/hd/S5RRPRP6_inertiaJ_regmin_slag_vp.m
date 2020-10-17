% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPRP6
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
% MM_reg [((5+1)*5/2)x21]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRP6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:58:27
% EndTime: 2019-12-31 19:58:28
% DurationCPUTime: 0.28s
% Computational Cost: add. (324->55), mult. (625->119), div. (0->0), fcn. (688->6), ass. (0->43)
t31 = cos(qJ(2));
t47 = 0.2e1 * t31;
t30 = cos(qJ(4));
t46 = t30 * pkin(4);
t26 = sin(pkin(8));
t27 = cos(pkin(8));
t29 = sin(qJ(2));
t15 = t26 * t29 - t27 * t31;
t28 = sin(qJ(4));
t45 = t28 * t15;
t16 = t26 * t31 + t27 * t29;
t44 = t28 * t16;
t43 = t28 * t30;
t40 = -qJ(3) - pkin(6);
t19 = t40 * t31;
t36 = t40 * t29;
t10 = -t27 * t19 + t26 * t36;
t42 = t30 * t10;
t41 = t30 * t16;
t24 = t28 ^ 2;
t25 = t30 ^ 2;
t39 = t24 + t25;
t38 = qJ(5) * t16;
t21 = t26 * pkin(2) + pkin(7);
t37 = qJ(5) + t21;
t23 = -t31 * pkin(2) - pkin(1);
t22 = -t27 * pkin(2) - pkin(3);
t7 = t15 * pkin(3) - t16 * pkin(7) + t23;
t3 = -t28 * t10 + t30 * t7;
t8 = -t26 * t19 - t27 * t36;
t1 = t15 * pkin(4) - t30 * t38 + t3;
t2 = t42 + (t7 - t38) * t28;
t35 = t1 * t30 + t2 * t28;
t12 = t37 * t28;
t13 = t37 * t30;
t34 = -t12 * t30 + t13 * t28;
t33 = -t15 * t21 + t16 * t22;
t18 = t22 - t46;
t14 = t16 ^ 2;
t11 = t30 * t15;
t5 = pkin(4) * t44 + t8;
t4 = t28 * t7 + t42;
t6 = [1, 0, 0, t29 ^ 2, t29 * t47, 0, 0, 0, pkin(1) * t47, -0.2e1 * pkin(1) * t29, -0.2e1 * t10 * t15 + 0.2e1 * t8 * t16, t10 ^ 2 + t23 ^ 2 + t8 ^ 2, t25 * t14, -0.2e1 * t14 * t43, 0.2e1 * t15 * t41, -0.2e1 * t15 * t44, t15 ^ 2, 0.2e1 * t3 * t15 + 0.2e1 * t8 * t44, -0.2e1 * t4 * t15 + 0.2e1 * t8 * t41, -0.2e1 * t35 * t16, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t29, t31, 0, -t29 * pkin(6), -t31 * pkin(6), (-t15 * t26 - t16 * t27) * pkin(2), (t10 * t26 - t27 * t8) * pkin(2), t28 * t41, (-t24 + t25) * t16, t45, t11, 0, t33 * t28 - t8 * t30, t8 * t28 + t33 * t30, -t1 * t28 - t34 * t16 + t2 * t30, -t1 * t12 + t2 * t13 + t5 * t18; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, (t26 ^ 2 + t27 ^ 2) * pkin(2) ^ 2, t24, 0.2e1 * t43, 0, 0, 0, -0.2e1 * t22 * t30, 0.2e1 * t22 * t28, 0.2e1 * t12 * t28 + 0.2e1 * t13 * t30, t12 ^ 2 + t13 ^ 2 + t18 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, 0, 0, 0, 0, t11, -t45, -t39 * t16, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t44, t15, t3, -t4, -pkin(4) * t41, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t30, 0, -t28 * t21, -t30 * t21, -t28 * pkin(4), -t12 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t28, 0, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t6;
