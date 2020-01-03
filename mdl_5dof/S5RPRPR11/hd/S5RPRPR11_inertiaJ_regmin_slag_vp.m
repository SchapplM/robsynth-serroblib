% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPR11_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t23 = sin(pkin(8));
t24 = cos(pkin(8));
t26 = sin(qJ(3));
t33 = cos(qJ(3));
t13 = t26 * t23 - t33 * t24;
t28 = -pkin(3) - pkin(4);
t14 = t33 * t23 + t26 * t24;
t20 = -t24 * pkin(2) - pkin(1);
t30 = t14 * qJ(4) - t20;
t36 = 0.2e1 * t28 * t13 + 0.2e1 * t30;
t35 = -0.2e1 * t14;
t34 = 0.2e1 * t20;
t32 = pkin(6) + qJ(2);
t31 = t23 ^ 2 + t24 ^ 2;
t17 = t32 * t23;
t18 = t32 * t24;
t9 = t33 * t17 + t26 * t18;
t10 = -t26 * t17 + t33 * t18;
t27 = cos(qJ(5));
t25 = sin(qJ(5));
t16 = t27 * qJ(4) + t25 * t28;
t15 = t25 * qJ(4) - t27 * t28;
t8 = t25 * t13 + t27 * t14;
t7 = -t27 * t13 + t25 * t14;
t6 = t13 * pkin(3) - t30;
t5 = t13 * pkin(7) + t10;
t4 = -t14 * pkin(7) + t9;
t2 = t25 * t4 + t27 * t5;
t1 = t25 * t5 - t27 * t4;
t3 = [1, 0, 0, 0.2e1 * pkin(1) * t24, -0.2e1 * pkin(1) * t23, 0.2e1 * t31 * qJ(2), t31 * qJ(2) ^ 2 + pkin(1) ^ 2, t14 ^ 2, t13 * t35, 0, 0, 0, t13 * t34, t14 * t34, 0.2e1 * t6 * t13, -0.2e1 * t10 * t13 + 0.2e1 * t9 * t14, t6 * t35, t10 ^ 2 + t6 ^ 2 + t9 ^ 2, t8 ^ 2, -0.2e1 * t8 * t7, 0, 0, 0, t7 * t36, t8 * t36; 0, 0, 0, -t24, t23, 0, -pkin(1), 0, 0, 0, 0, 0, t13, t14, t13, 0, -t14, t6, 0, 0, 0, 0, 0, -t7, -t8; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13, 0, -t9, -t10, -t9, -pkin(3) * t14 - t13 * qJ(4), t10, -t9 * pkin(3) + t10 * qJ(4), 0, 0, -t8, t7, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, 0.2e1 * qJ(4), pkin(3) ^ 2 + qJ(4) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t15, 0.2e1 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, t9, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3), 0, 0, 0, 0, 0, -t27, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, 0, -t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t15, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t3;
