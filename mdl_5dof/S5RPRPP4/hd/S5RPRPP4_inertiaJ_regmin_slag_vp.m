% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x19]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPP4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t24 = sin(pkin(7));
t25 = cos(pkin(7));
t26 = cos(qJ(3));
t36 = sin(qJ(3));
t11 = -t24 * t36 + t25 * t26;
t12 = -t24 * t26 - t25 * t36;
t44 = (t11 * t25 - t12 * t24) * pkin(3);
t34 = t11 ^ 2 + t12 ^ 2;
t39 = -pkin(1) - pkin(6);
t31 = t36 * t39;
t15 = -t36 * qJ(4) + t31;
t32 = (-qJ(4) + t39) * t26;
t6 = t24 * t15 - t25 * t32;
t8 = t25 * t15 + t24 * t32;
t43 = t6 * t11 + t8 * t12;
t22 = t36 * pkin(3) + qJ(2);
t4 = -t12 * pkin(4) - t11 * qJ(5) + t22;
t42 = -0.2e1 * t4;
t40 = 0.2e1 * qJ(2);
t35 = t6 ^ 2 + t8 ^ 2;
t16 = t24 * pkin(3) + qJ(5);
t20 = t25 * pkin(3) + pkin(4);
t30 = t20 * t11 - t16 * t12;
t28 = 0.2e1 * t43;
t1 = [1, 0, 0, -2 * pkin(1), t40, (pkin(1) ^ 2) + qJ(2) ^ 2, t26 ^ 2, -0.2e1 * t26 * t36, 0, 0, 0, t36 * t40, t26 * t40, t28, t22 ^ 2 + t35, t12 * t42, t28, t11 * t42, t4 ^ 2 + t35; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, -t34, -t43, 0, -t34, 0, -t43; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, 0, 0, t34; 0, 0, 0, 0, 0, 0, 0, 0, t26, -t36, 0, t26 * t39, -t31, -t44, (t24 * t8 - t25 * t6) * pkin(3), -t6, -t30, t8, t8 * t16 - t6 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t36, 0, t44, t11, 0, -t12, t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, (t24 ^ 2 + t25 ^ 2) * pkin(3) ^ 2, 0.2e1 * t20, 0, 0.2e1 * t16, t16 ^ 2 + t20 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t12, 0, -t11, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t1;
