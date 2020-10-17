% Calculate minimal parameter regressor of joint inertia matrix for
% S4RRRR5
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
% MM_reg [((4+1)*4/2)x24]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRRR5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR5_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR5_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:28:13
% EndTime: 2019-12-31 17:28:14
% DurationCPUTime: 0.24s
% Computational Cost: add. (159->58), mult. (380->112), div. (0->0), fcn. (417->6), ass. (0->46)
t25 = sin(qJ(4));
t26 = sin(qJ(3));
t28 = cos(qJ(4));
t29 = cos(qJ(3));
t13 = t25 * t26 - t28 * t29;
t27 = sin(qJ(2));
t11 = t13 * t27;
t47 = 0.2e1 * t11;
t20 = -t29 * pkin(3) - pkin(2);
t46 = 0.2e1 * t20;
t45 = -0.2e1 * t27;
t30 = cos(qJ(2));
t44 = 0.2e1 * t30;
t43 = pkin(6) + pkin(7);
t42 = pkin(2) * t29;
t41 = pkin(5) * t26;
t40 = t25 * pkin(3);
t39 = t28 * pkin(3);
t16 = -t30 * pkin(2) - t27 * pkin(6) - pkin(1);
t33 = t29 * t30;
t31 = pkin(5) * t33;
t5 = t31 + (-pkin(7) * t27 + t16) * t26;
t38 = t28 * t5;
t37 = t30 * pkin(3);
t36 = t26 * t29;
t35 = t26 * t30;
t34 = t29 * t27;
t32 = t27 * t44;
t12 = t29 * t16;
t4 = -pkin(7) * t34 + t12 + (-pkin(3) - t41) * t30;
t1 = -t25 * t5 + t28 * t4;
t14 = t25 * t29 + t28 * t26;
t24 = t30 ^ 2;
t23 = t29 ^ 2;
t22 = t27 ^ 2;
t21 = t26 ^ 2;
t18 = t43 * t29;
t17 = t43 * t26;
t15 = (pkin(3) * t26 + pkin(5)) * t27;
t10 = t14 * t27;
t9 = t26 * t16 + t31;
t8 = -pkin(5) * t35 + t12;
t7 = -t25 * t17 + t28 * t18;
t6 = -t28 * t17 - t25 * t18;
t2 = t25 * t4 + t38;
t3 = [1, 0, 0, t22, t32, 0, 0, 0, pkin(1) * t44, pkin(1) * t45, t23 * t22, -0.2e1 * t22 * t36, t33 * t45, t26 * t32, t24, 0.2e1 * t22 * t41 - 0.2e1 * t8 * t30, 0.2e1 * t22 * pkin(5) * t29 + 0.2e1 * t9 * t30, t11 ^ 2, t10 * t47, t30 * t47, t10 * t44, t24, -0.2e1 * t1 * t30 + 0.2e1 * t15 * t10, -0.2e1 * t15 * t11 + 0.2e1 * t2 * t30; 0, 0, 0, 0, 0, t27, t30, 0, -t27 * pkin(5), -t30 * pkin(5), t26 * t34, (-t21 + t23) * t27, -t35, -t33, 0, -pkin(5) * t34 + (-pkin(2) * t27 + pkin(6) * t30) * t26, pkin(6) * t33 + (t41 - t42) * t27, -t11 * t14, -t14 * t10 + t11 * t13, -t14 * t30, t13 * t30, 0, t20 * t10 + t15 * t13 - t6 * t30, -t20 * t11 + t15 * t14 + t7 * t30; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t21, 0.2e1 * t36, 0, 0, 0, 0.2e1 * t42, -0.2e1 * pkin(2) * t26, t14 ^ 2, -0.2e1 * t14 * t13, 0, 0, 0, t13 * t46, t14 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t26 * t27, -t30, t8, -t9, 0, 0, -t11, -t10, -t30, -t28 * t37 + t1, -t38 + (-t4 + t37) * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, t29, 0, -t26 * pkin(6), -t29 * pkin(6), 0, 0, t14, -t13, 0, t6, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t39, -0.2e1 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t10, -t30, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13, 0, t6, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t39, -t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t3;
