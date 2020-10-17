% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPRP2
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
% MM_reg [((5+1)*5/2)x18]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRP2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:49:43
% EndTime: 2019-12-31 19:49:43
% DurationCPUTime: 0.19s
% Computational Cost: add. (167->34), mult. (287->63), div. (0->0), fcn. (245->6), ass. (0->34)
t24 = sin(pkin(8));
t18 = t24 * pkin(2) + pkin(7);
t26 = sin(qJ(4));
t22 = t26 ^ 2;
t28 = cos(qJ(4));
t33 = t28 ^ 2 + t22;
t34 = t33 * t18;
t46 = -0.2e1 * t26;
t45 = 0.2e1 * t26;
t44 = -0.2e1 * t28;
t32 = t28 * pkin(4) + t26 * qJ(5);
t31 = -pkin(3) - t32;
t21 = cos(qJ(2)) * pkin(1);
t20 = t21 + pkin(2);
t25 = cos(pkin(8));
t39 = sin(qJ(2)) * pkin(1);
t7 = t25 * t20 - t24 * t39;
t1 = t31 - t7;
t41 = t25 * pkin(2);
t9 = t31 - t41;
t43 = -t1 - t9;
t8 = t24 * t20 + t25 * t39;
t6 = pkin(7) + t8;
t42 = t33 * t6;
t40 = t26 * t6;
t38 = t28 * t6;
t19 = -pkin(3) - t41;
t5 = -pkin(3) - t7;
t37 = t19 + t5;
t36 = t26 * t18;
t35 = t28 * t18;
t10 = -t26 * pkin(4) + t28 * qJ(5);
t17 = t28 * t45;
t2 = [1, 0, 0, 1, 0.2e1 * t21, -0.2e1 * t39, t7 ^ 2 + t8 ^ 2, t22, t17, 0, 0, 0, t5 * t44, t5 * t45, t1 * t44, 0.2e1 * t42, t1 * t46, t33 * t6 ^ 2 + t1 ^ 2; 0, 0, 0, 1, t21, -t39, (t24 * t8 + t25 * t7) * pkin(2), t22, t17, 0, 0, 0, -t37 * t28, t37 * t26, t43 * t28, t34 + t42, t43 * t26, t1 * t9 + t6 * t34; 0, 0, 0, 1, 0, 0, (t24 ^ 2 + t25 ^ 2) * pkin(2) ^ 2, t22, t17, 0, 0, 0, t19 * t44, t19 * t45, t9 * t44, 0.2e1 * t34, t9 * t46, t33 * t18 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, t28, 0, -t40, -t38, -t40, t10, t38, t10 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, t28, 0, -t36, -t35, -t36, t10, t35, t10 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t26, t28, 0, t26, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t2;
