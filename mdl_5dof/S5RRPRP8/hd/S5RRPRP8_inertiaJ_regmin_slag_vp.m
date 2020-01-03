% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x23]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRP8_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t26 = sin(qJ(2));
t28 = cos(qJ(2));
t13 = -t28 * pkin(2) - t26 * qJ(3) - pkin(1);
t6 = t28 * pkin(3) - t13;
t35 = 0.2e1 * t6;
t34 = -0.2e1 * t26;
t33 = 0.2e1 * t28;
t23 = t26 ^ 2;
t32 = t28 ^ 2 + t23;
t19 = t26 * pkin(6);
t14 = -t26 * pkin(7) + t19;
t20 = t28 * pkin(6);
t15 = -t28 * pkin(7) + t20;
t25 = sin(qJ(4));
t27 = cos(qJ(4));
t4 = -t27 * t14 + t25 * t15;
t29 = -pkin(2) - pkin(3);
t11 = t25 * qJ(3) - t27 * t29;
t31 = -t26 * pkin(2) + t28 * qJ(3);
t5 = t25 * t14 + t27 * t15;
t12 = t27 * qJ(3) + t25 * t29;
t10 = -pkin(4) - t11;
t8 = -t28 * t25 + t26 * t27;
t7 = t26 * t25 + t28 * t27;
t3 = t7 * pkin(4) + t6;
t2 = -t7 * qJ(5) + t5;
t1 = -t8 * qJ(5) - t4;
t9 = [1, 0, 0, t23, t26 * t33, 0, 0, 0, pkin(1) * t33, pkin(1) * t34, -0.2e1 * t13 * t28, 0.2e1 * t32 * pkin(6), t13 * t34, t32 * pkin(6) ^ 2 + t13 ^ 2, t8 ^ 2, -0.2e1 * t8 * t7, 0, 0, 0, t7 * t35, t8 * t35, -0.2e1 * t1 * t8 - 0.2e1 * t2 * t7, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, t26, t28, 0, -t19, -t20, -t19, t31, t20, t31 * pkin(6), 0, 0, -t8, t7, 0, t4, t5, -t10 * t8 - t12 * t7, t1 * t10 + t2 * t12; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2), 0, 0.2e1 * qJ(3), pkin(2) ^ 2 + qJ(3) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t11, 0.2e1 * t12, 0, t10 ^ 2 + t12 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, t19, 0, 0, 0, 0, 0, 0, 0, -t25 * t7 - t27 * t8, t1 * t27 + t2 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(2), 0, 0, 0, 0, 0, -t27, t25, 0, t10 * t27 + t12 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t25 ^ 2 + t27 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, 0, -t4, -t5, -pkin(4) * t8, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t11, -t12, 0, t10 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t25, 0, t27 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t9;
