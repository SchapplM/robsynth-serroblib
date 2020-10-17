% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPR14_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:35:16
% EndTime: 2019-12-31 18:35:17
% DurationCPUTime: 0.26s
% Computational Cost: add. (197->47), mult. (353->87), div. (0->0), fcn. (411->6), ass. (0->35)
t24 = sin(pkin(8));
t25 = cos(pkin(8));
t28 = cos(qJ(3));
t41 = sin(qJ(3));
t11 = -t24 * t41 + t25 * t28;
t12 = -t24 * t28 - t25 * t41;
t44 = (t11 * t25 - t12 * t24) * pkin(3);
t10 = t12 ^ 2;
t9 = t11 ^ 2;
t43 = t10 + t9;
t42 = 2 * qJ(2);
t26 = sin(qJ(5));
t40 = t26 * t11;
t39 = t26 * t12;
t27 = cos(qJ(5));
t38 = t26 * t27;
t37 = t27 * t11;
t7 = t27 * t12;
t20 = t41 * pkin(3) + qJ(2);
t29 = -pkin(1) - pkin(6);
t35 = t41 * t29;
t34 = (-qJ(4) + t29) * t28;
t14 = -t41 * qJ(4) + t35;
t4 = t24 * t14 - t25 * t34;
t6 = t25 * t14 + t24 * t34;
t33 = t4 * t11 + t6 * t12;
t18 = t24 * pkin(3) + pkin(7);
t19 = -t25 * pkin(3) - pkin(4);
t32 = t11 * t19 + t12 * t18;
t23 = t27 ^ 2;
t22 = t26 ^ 2;
t3 = -t12 * pkin(4) - t11 * pkin(7) + t20;
t2 = t26 * t3 + t27 * t6;
t1 = -t26 * t6 + t27 * t3;
t5 = [1, 0, 0, -2 * pkin(1), t42, pkin(1) ^ 2 + qJ(2) ^ 2, t28 ^ 2, -0.2e1 * t28 * t41, 0, 0, 0, t41 * t42, t28 * t42, 0.2e1 * t33, t20 ^ 2 + t4 ^ 2 + t6 ^ 2, t23 * t9, -0.2e1 * t9 * t38, -0.2e1 * t11 * t7, 0.2e1 * t11 * t39, t10, -0.2e1 * t1 * t12 + 0.2e1 * t4 * t40, 0.2e1 * t2 * t12 + 0.2e1 * t4 * t37; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, -t43, -t33, 0, 0, 0, 0, 0, -t43 * t26, -t43 * t27; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t43, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t28, -t41, 0, t28 * t29, -t35, -t44, (t24 * t6 - t25 * t4) * pkin(3), t26 * t37, (-t22 + t23) * t11, -t39, -t7, 0, t32 * t26 - t4 * t27, t4 * t26 + t32 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t41, 0, t44, 0, 0, 0, 0, 0, t37, -t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, (t24 ^ 2 + t25 ^ 2) * pkin(3) ^ 2, t22, 0.2e1 * t38, 0, 0, 0, -0.2e1 * t19 * t27, 0.2e1 * t19 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, 0, 0, 0, 0, -t7, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t40, -t12, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, t27, 0, -t26 * t18, -t27 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t5;
