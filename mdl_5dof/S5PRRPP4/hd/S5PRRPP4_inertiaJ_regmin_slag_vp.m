% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
% 
% Output:
% MM_reg [((5+1)*5/2)x19]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRPP4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t11 = sin(qJ(3));
t23 = -0.2e1 * t11;
t12 = cos(qJ(3));
t22 = 0.2e1 * t12;
t6 = t11 * qJ(4);
t21 = t12 * pkin(3) + t6;
t20 = t11 * pkin(6);
t9 = t11 ^ 2;
t5 = t12 ^ 2 + t9;
t19 = t12 * qJ(4);
t2 = -pkin(2) - t21;
t18 = -t11 * pkin(3) + t19;
t16 = qJ(4) ^ 2;
t15 = 0.2e1 * qJ(4);
t13 = pkin(3) + pkin(4);
t7 = t12 * pkin(6);
t4 = -t12 * qJ(5) + t7;
t3 = (pkin(6) - qJ(5)) * t11;
t1 = t12 * pkin(4) - t2;
t8 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11 * t4 - t12 * t3; 0, 1, 0, 0, t9, t11 * t22, 0, 0, 0, pkin(2) * t22, pkin(2) * t23, -0.2e1 * t2 * t12, 0.2e1 * t5 * pkin(6), t2 * t23, t5 * pkin(6) ^ 2 + t2 ^ 2, t1 * t22, 0.2e1 * t1 * t11, -0.2e1 * t3 * t11 - 0.2e1 * t4 * t12, t1 ^ 2 + t3 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t11, t12, 0, t11, t21, t12, t11, 0, t12 * t13 + t6; 0, 0, 0, 0, 0, 0, t11, t12, 0, -t20, -t7, -t20, t18, t7, t18 * pkin(6), -t3, t4, t13 * t11 - t19, t4 * qJ(4) - t3 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, t15, pkin(3) ^ 2 + t16, 0.2e1 * t13, t15, 0, t13 ^ 2 + t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, 0, 0, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, t20, 0, 0, -t11, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3), -1, 0, 0, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t11, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t8;
