% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x23]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:03
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRR5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:02:55
% EndTime: 2022-01-20 11:02:57
% DurationCPUTime: 0.33s
% Computational Cost: add. (306->45), mult. (590->80), div. (0->0), fcn. (707->8), ass. (0->53)
t40 = sin(pkin(9));
t41 = cos(pkin(9));
t53 = t40 ^ 2 + t41 ^ 2;
t54 = t53 * qJ(3);
t43 = sin(qJ(4));
t46 = cos(qJ(4));
t23 = t43 * t40 - t46 * t41;
t33 = -t41 * pkin(3) - pkin(2);
t15 = t23 * pkin(4) + t33;
t58 = cos(qJ(2)) * pkin(1);
t14 = t15 - t58;
t66 = 0.2e1 * t14;
t65 = 0.2e1 * t15;
t25 = t33 - t58;
t64 = 0.2e1 * t25;
t63 = 0.2e1 * t33;
t24 = t46 * t40 + t43 * t41;
t62 = t24 * pkin(8);
t42 = sin(qJ(5));
t61 = t42 * pkin(4);
t60 = sin(qJ(2)) * pkin(1);
t45 = cos(qJ(5));
t59 = t45 * pkin(4);
t57 = t14 + t15;
t56 = t25 + t33;
t32 = qJ(3) + t60;
t55 = t53 * t32;
t18 = (-pkin(7) - t32) * t40;
t37 = t41 * pkin(7);
t19 = t41 * t32 + t37;
t52 = t46 * t18 - t43 * t19;
t26 = (-pkin(7) - qJ(3)) * t40;
t27 = t41 * qJ(3) + t37;
t51 = t46 * t26 - t43 * t27;
t50 = -t43 * t18 - t46 * t19;
t49 = -t43 * t26 - t46 * t27;
t34 = -pkin(2) - t58;
t21 = t24 ^ 2;
t20 = t23 * pkin(8);
t13 = -0.2e1 * t24 * t23;
t12 = -t42 * t23 + t45 * t24;
t11 = t45 * t23 + t42 * t24;
t10 = t12 ^ 2;
t9 = -t20 - t49;
t8 = t51 - t62;
t7 = -t20 - t50;
t6 = t52 - t62;
t5 = -0.2e1 * t12 * t11;
t4 = -t42 * t8 - t45 * t9;
t3 = -t42 * t9 + t45 * t8;
t2 = -t42 * t6 - t45 * t7;
t1 = -t42 * t7 + t45 * t6;
t16 = [1, 0, 0, 1, 0.2e1 * t58, -0.2e1 * t60, -0.2e1 * t34 * t41, 0.2e1 * t55, t53 * t32 ^ 2 + t34 ^ 2, t21, t13, 0, 0, 0, t23 * t64, t24 * t64, t10, t5, 0, 0, 0, t11 * t66, t12 * t66; 0, 0, 0, 1, t58, -t60, (pkin(2) - t34) * t41, t54 + t55, -t34 * pkin(2) + t32 * t54, t21, t13, 0, 0, 0, t56 * t23, t56 * t24, t10, t5, 0, 0, 0, t57 * t11, t57 * t12; 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2) * t41, 0.2e1 * t54, t53 * qJ(3) ^ 2 + pkin(2) ^ 2, t21, t13, 0, 0, 0, t23 * t63, t24 * t63, t10, t5, 0, 0, 0, t11 * t65, t12 * t65; 0, 0, 0, 0, 0, 0, -t41, 0, t34, 0, 0, 0, 0, 0, t23, t24, 0, 0, 0, 0, 0, t11, t12; 0, 0, 0, 0, 0, 0, -t41, 0, -pkin(2), 0, 0, 0, 0, 0, t23, t24, 0, 0, 0, 0, 0, t11, t12; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t23, 0, t52, t50, 0, 0, t12, -t11, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t23, 0, t51, t49, 0, 0, t12, -t11, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t59, -0.2e1 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t11, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t11, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t59, -t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t16;
