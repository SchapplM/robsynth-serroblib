% Calculate minimal parameter regressor of joint inertia matrix for
% S4RPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% MM_reg [((4+1)*4/2)x18]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RPRP5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_inertiaJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t12 = sin(pkin(6));
t13 = cos(pkin(6));
t14 = sin(qJ(3));
t19 = cos(qJ(3));
t6 = t19 * t12 + t14 * t13;
t21 = -0.2e1 * t6;
t9 = -t13 * pkin(2) - pkin(1);
t20 = 0.2e1 * t9;
t18 = pkin(5) + qJ(2);
t17 = t12 ^ 2 + t13 ^ 2;
t16 = t18 * t12;
t7 = t18 * t13;
t5 = t14 * t12 - t19 * t13;
t3 = -t14 * t16 + t19 * t7;
t2 = t14 * t7 + t19 * t16;
t1 = t5 * pkin(3) - t6 * qJ(4) + t9;
t4 = [1, 0, 0, 0.2e1 * pkin(1) * t13, -0.2e1 * pkin(1) * t12, 0.2e1 * t17 * qJ(2), t17 * qJ(2) ^ 2 + pkin(1) ^ 2, t6 ^ 2, t5 * t21, 0, 0, 0, t5 * t20, t6 * t20, 0.2e1 * t1 * t5, 0.2e1 * t2 * t6 - 0.2e1 * t3 * t5, t1 * t21, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, -t13, t12, 0, -pkin(1), 0, 0, 0, 0, 0, t5, t6, t5, 0, -t6, t1; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, 0, -t2, -t3, -t2, -pkin(3) * t6 - t5 * qJ(4), t3, -t2 * pkin(3) + t3 * qJ(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, 0.2e1 * qJ(4), pkin(3) ^ 2 + qJ(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t4;
