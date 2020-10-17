% Calculate inertial parameters regressor of joint inertia matrix for
% S4RRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% MM_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRRP3_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_inertiaJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:14:15
% EndTime: 2019-12-31 17:14:16
% DurationCPUTime: 0.28s
% Computational Cost: add. (101->29), mult. (227->56), div. (0->0), fcn. (161->4), ass. (0->34)
t24 = sin(qJ(3));
t22 = t24 ^ 2;
t26 = cos(qJ(3));
t23 = t26 ^ 2;
t46 = t22 + t23;
t25 = sin(qJ(2));
t40 = t25 * pkin(1);
t16 = pkin(6) + t40;
t32 = t46 * t16;
t45 = -0.2e1 * t24;
t44 = -0.2e1 * t26;
t27 = cos(qJ(2));
t38 = t27 * pkin(1);
t5 = -t26 * pkin(3) - t24 * qJ(4) - pkin(2);
t3 = t5 - t38;
t43 = -t3 - t5;
t42 = t32 * pkin(6);
t41 = t24 * pkin(6);
t39 = t26 * pkin(6);
t17 = -pkin(2) - t38;
t37 = pkin(2) - t17;
t36 = t46 * t16 ^ 2;
t35 = t24 * t16;
t34 = t24 * t26;
t33 = t26 * t16;
t31 = t46 * pkin(6) ^ 2;
t30 = t46 * pkin(6);
t8 = -t24 * pkin(3) + t26 * qJ(4);
t14 = -0.2e1 * t34;
t13 = 0.2e1 * t34;
t4 = 0.2e1 * t30;
t2 = 0.2e1 * t32;
t1 = t30 + t32;
t6 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t38, -0.2e1 * t40, 0, (t25 ^ 2 + t27 ^ 2) * pkin(1) ^ 2, t22, t13, 0, t23, 0, 0, t17 * t44, 0.2e1 * t17 * t24, t2, t17 ^ 2 + t36, t22, 0, t14, 0, 0, t23, t3 * t44, t2, t3 * t45, t3 ^ 2 + t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t38, -t40, 0, 0, t22, t13, 0, t23, 0, 0, t37 * t26, -t37 * t24, t1, -t17 * pkin(2) + t42, t22, 0, t14, 0, 0, t23, t43 * t26, t1, t43 * t24, t3 * t5 + t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t22, t13, 0, t23, 0, 0, 0.2e1 * pkin(2) * t26, pkin(2) * t45, t4, pkin(2) ^ 2 + t31, t22, 0, t14, 0, 0, t23, t5 * t44, t4, t5 * t45, t5 ^ 2 + t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, t26, 0, -t35, -t33, 0, 0, 0, t24, 0, 0, -t26, 0, -t35, t8, t33, t8 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, t26, 0, -t41, -t39, 0, 0, 0, t24, 0, 0, -t26, 0, -t41, t8, t39, t8 * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, 0.2e1 * qJ(4), pkin(3) ^ 2 + qJ(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t6;
