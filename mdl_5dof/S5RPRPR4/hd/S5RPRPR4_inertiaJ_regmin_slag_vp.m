% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:23
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPR4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:23:10
% EndTime: 2022-01-23 09:23:11
% DurationCPUTime: 0.23s
% Computational Cost: add. (217->40), mult. (409->80), div. (0->0), fcn. (477->8), ass. (0->33)
t24 = sin(pkin(9));
t26 = cos(pkin(9));
t29 = sin(qJ(3));
t31 = cos(qJ(3));
t16 = t24 * t29 - t26 * t31;
t27 = cos(pkin(8));
t23 = -t27 * pkin(1) - pkin(2);
t19 = -t31 * pkin(3) + t23;
t39 = 0.2e1 * t16 * pkin(4) + 0.2e1 * t19;
t38 = 0.2e1 * t19;
t37 = 0.2e1 * t29;
t36 = t24 * pkin(3);
t35 = t26 * pkin(3);
t25 = sin(pkin(8));
t21 = t25 * pkin(1) + pkin(6);
t34 = qJ(4) + t21;
t14 = t34 * t29;
t15 = t34 * t31;
t5 = -t26 * t14 - t24 * t15;
t6 = -t24 * t14 + t26 * t15;
t30 = cos(qJ(5));
t28 = sin(qJ(5));
t22 = pkin(4) + t35;
t18 = t24 * t31 + t26 * t29;
t12 = -t28 * t22 - t30 * t36;
t11 = t30 * t22 - t28 * t36;
t8 = -t28 * t16 + t30 * t18;
t7 = t30 * t16 + t28 * t18;
t4 = -t16 * pkin(7) + t6;
t3 = -t18 * pkin(7) + t5;
t2 = -t28 * t3 - t30 * t4;
t1 = -t28 * t4 + t30 * t3;
t9 = [1, 0, 0, (t25 ^ 2 + t27 ^ 2) * pkin(1) ^ 2, t29 ^ 2, t31 * t37, 0, 0, 0, -0.2e1 * t23 * t31, t23 * t37, t16 * t38, t18 * t38, -0.2e1 * t6 * t16 - 0.2e1 * t5 * t18, t19 ^ 2 + t5 ^ 2 + t6 ^ 2, t8 ^ 2, -0.2e1 * t8 * t7, 0, 0, 0, t7 * t39, t8 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5 * t16 + t6 * t18, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16 ^ 2 + t18 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t29, t31, 0, -t29 * t21, -t31 * t21, t5, -t6, (-t16 * t24 - t18 * t26) * pkin(3), (t24 * t6 + t26 * t5) * pkin(3), 0, 0, t8, -t7, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t29, -t16, -t18, 0, (-t16 * t26 + t18 * t24) * pkin(3), 0, 0, 0, 0, 0, -t7, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * t35, -0.2e1 * t36, 0, (t24 ^ 2 + t26 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t11, 0.2e1 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t18, 0, t19, 0, 0, 0, 0, 0, t7, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t11, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t9;
