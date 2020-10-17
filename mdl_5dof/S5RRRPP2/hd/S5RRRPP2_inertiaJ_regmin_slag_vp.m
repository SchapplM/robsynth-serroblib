% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x21]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPP2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP2_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP2_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:51:58
% EndTime: 2019-12-31 20:51:59
% DurationCPUTime: 0.29s
% Computational Cost: add. (192->53), mult. (315->84), div. (0->0), fcn. (256->4), ass. (0->40)
t45 = sin(qJ(2)) * pkin(1);
t16 = pkin(7) + t45;
t28 = sin(qJ(3));
t26 = t28 ^ 2;
t30 = cos(qJ(3));
t40 = t30 ^ 2 + t26;
t42 = t40 * t16;
t52 = -0.2e1 * t28;
t51 = 0.2e1 * t28;
t50 = -0.2e1 * t30;
t49 = 0.2e1 * t30;
t25 = cos(qJ(2)) * pkin(1);
t37 = t30 * pkin(3) + t28 * qJ(4) + pkin(2);
t2 = -t25 - t37;
t23 = t30 * pkin(4);
t1 = t23 - t2;
t3 = t23 + t37;
t48 = t1 + t3;
t47 = -t2 + t37;
t46 = t28 * pkin(7);
t17 = -t25 - pkin(2);
t44 = pkin(2) - t17;
t43 = t28 * t16;
t41 = t40 * pkin(7);
t39 = qJ(4) * t30;
t38 = qJ(5) * t30;
t9 = -pkin(3) * t28 + t39;
t35 = qJ(4) ^ 2;
t34 = 0.2e1 * qJ(4);
t32 = pkin(3) + pkin(4);
t22 = t30 * pkin(7);
t18 = t28 * qJ(5);
t14 = t28 * t49;
t13 = t30 * t16;
t10 = t22 - t38;
t8 = -t18 + t46;
t6 = t28 * t32 - t39;
t5 = t13 - t38;
t4 = -t18 + t43;
t7 = [1, 0, 0, 1, 0.2e1 * t25, -0.2e1 * t45, t26, t14, 0, 0, 0, t17 * t50, t17 * t51, t2 * t50, 0.2e1 * t42, t2 * t52, t40 * t16 ^ 2 + t2 ^ 2, t1 * t49, t1 * t51, -0.2e1 * t4 * t28 - 0.2e1 * t5 * t30, t1 ^ 2 + t4 ^ 2 + t5 ^ 2; 0, 0, 0, 1, t25, -t45, t26, t14, 0, 0, 0, t44 * t30, -t44 * t28, t47 * t30, t41 + t42, t47 * t28, pkin(7) * t42 - t2 * t37, t48 * t30, t48 * t28, (-t10 - t5) * t30 + (-t4 - t8) * t28, t1 * t3 + t10 * t5 + t4 * t8; 0, 0, 0, 1, 0, 0, t26, t14, 0, 0, 0, pkin(2) * t49, pkin(2) * t52, -t37 * t50, 0.2e1 * t41, -t37 * t52, t40 * pkin(7) ^ 2 + t37 ^ 2, t3 * t49, t3 * t51, -0.2e1 * t10 * t30 - 0.2e1 * t8 * t28, t10 ^ 2 + t3 ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, t28, t30, 0, -t43, -t13, -t43, t9, t13, t9 * t16, -t4, t5, t6, qJ(4) * t5 - t32 * t4; 0, 0, 0, 0, 0, 0, 0, 0, t28, t30, 0, -t46, -t22, -t46, t9, t22, t9 * pkin(7), -t8, t10, t6, qJ(4) * t10 - t32 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, t34, pkin(3) ^ 2 + t35, 0.2e1 * t32, t34, 0, t32 ^ 2 + t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, t43, 0, 0, -t28, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, t46, 0, 0, -t28, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3), -1, 0, 0, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t28, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t28, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t7;
