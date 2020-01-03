% Calculate minimal parameter regressor of joint inertia matrix for
% S4RPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% 
% Output:
% MM_reg [((4+1)*4/2)x20]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RPRR8_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_inertiaJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t10 = sin(qJ(3));
t17 = 0.2e1 * t10 * pkin(3) + (2 * qJ(2));
t16 = 2 * qJ(2);
t9 = sin(qJ(4));
t15 = t9 * pkin(3);
t11 = cos(qJ(4));
t14 = t11 * pkin(3);
t13 = -pkin(1) - pkin(5);
t12 = cos(qJ(3));
t8 = t12 * t13;
t6 = -t12 * pkin(6) + t8;
t5 = (-pkin(6) + t13) * t10;
t4 = -t9 * t10 + t11 * t12;
t3 = t11 * t10 + t9 * t12;
t2 = -t11 * t5 - t9 * t6;
t1 = t11 * t6 - t9 * t5;
t7 = [1, 0, 0, -2 * pkin(1), t16, pkin(1) ^ 2 + qJ(2) ^ 2, t12 ^ 2, -0.2e1 * t12 * t10, 0, 0, 0, t10 * t16, t12 * t16, t4 ^ 2, -0.2e1 * t4 * t3, 0, 0, 0, t3 * t17, t4 * t17; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t12, -t10, 0, t8, -t10 * t13, 0, 0, t4, -t3, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t10, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t14, -0.2e1 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t14, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t7;
