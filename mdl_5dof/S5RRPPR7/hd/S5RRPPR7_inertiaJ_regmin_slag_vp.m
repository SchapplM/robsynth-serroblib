% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:00
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPPR7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:59:01
% EndTime: 2021-01-15 19:59:02
% DurationCPUTime: 0.28s
% Computational Cost: add. (240->49), mult. (475->99), div. (0->0), fcn. (530->6), ass. (0->44)
t25 = sin(pkin(8));
t26 = cos(pkin(8));
t28 = sin(qJ(2));
t30 = cos(qJ(2));
t13 = t25 * t28 - t26 * t30;
t14 = t25 * t30 + t26 * t28;
t22 = -t30 * pkin(2) - pkin(1);
t33 = -t14 * qJ(4) + t22;
t6 = t13 * pkin(3) + t33;
t49 = -0.2e1 * t6;
t45 = t25 * pkin(2);
t18 = qJ(4) + t45;
t48 = 0.2e1 * t18;
t47 = 0.2e1 * t22;
t46 = 0.2e1 * t30;
t44 = t26 * pkin(2);
t43 = t18 * t13;
t27 = sin(qJ(5));
t42 = t27 * t13;
t41 = t27 * t14;
t29 = cos(qJ(5));
t40 = t29 * t13;
t11 = t29 * t14;
t39 = t29 * t27;
t38 = -qJ(3) - pkin(6);
t37 = 0.2e1 * t13 * t14;
t16 = t38 * t30;
t35 = t38 * t28;
t7 = -t25 * t16 - t26 * t35;
t9 = -t26 * t16 + t25 * t35;
t36 = t7 ^ 2 + t9 ^ 2;
t21 = -pkin(3) - t44;
t17 = -pkin(7) + t21;
t34 = -t14 * t17 + t43;
t32 = -0.2e1 * t9 * t13 + 0.2e1 * t7 * t14;
t24 = t29 ^ 2;
t23 = t27 ^ 2;
t12 = t13 ^ 2;
t5 = -t13 * pkin(4) + t9;
t4 = t14 * pkin(4) + t7;
t3 = (pkin(3) + pkin(7)) * t13 + t33;
t2 = t27 * t4 + t29 * t3;
t1 = -t27 * t3 + t29 * t4;
t8 = [1, 0, 0, t28 ^ 2, t28 * t46, 0, 0, 0, pkin(1) * t46, -0.2e1 * pkin(1) * t28, t13 * t47, t14 * t47, t32, t22 ^ 2 + t36, t32, t13 * t49, t14 * t49, t6 ^ 2 + t36, t23 * t12, 0.2e1 * t12 * t39, t27 * t37, t29 * t37, t14 ^ 2, 0.2e1 * t1 * t14 - 0.2e1 * t5 * t40, -0.2e1 * t2 * t14 + 0.2e1 * t5 * t42; 0, 0, 0, 0, 0, t28, t30, 0, -t28 * pkin(6), -t30 * pkin(6), -t7, -t9, (-t13 * t25 - t14 * t26) * pkin(2), (t25 * t9 - t26 * t7) * pkin(2), t21 * t14 - t43, t7, t9, t9 * t18 + t7 * t21, t13 * t39, (-t23 + t24) * t13, t11, -t41, 0, t5 * t27 - t34 * t29, t34 * t27 + t5 * t29; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * t44, -0.2e1 * t45, 0, (t25 ^ 2 + t26 ^ 2) * pkin(2) ^ 2, 0, 0.2e1 * t21, t48, t18 ^ 2 + t21 ^ 2, t24, -0.2e1 * t39, 0, 0, 0, t27 * t48, t29 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t14, 0, t22, 0, -t13, -t14, t6, 0, 0, 0, 0, 0, -t41, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, 0, t7, 0, 0, 0, 0, 0, t11, -t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, t21, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, t40, t14, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t27, 0, t29 * t17, -t27 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, -t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t8;
