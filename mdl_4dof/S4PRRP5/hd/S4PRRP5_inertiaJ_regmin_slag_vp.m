% Calculate minimal parameter regressor of joint inertia matrix for
% S4PRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% MM_reg [((4+1)*4/2)x13]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4PRRP5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_inertiaJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:29:22
% EndTime: 2019-12-31 16:29:22
% DurationCPUTime: 0.09s
% Computational Cost: add. (27->17), mult. (58->32), div. (0->0), fcn. (56->4), ass. (0->14)
t9 = cos(qJ(3));
t15 = 0.2e1 * t9;
t7 = sin(qJ(3));
t8 = sin(qJ(2));
t14 = t7 * t8;
t4 = t7 ^ 2;
t13 = t9 ^ 2 + t4;
t12 = -qJ(4) - pkin(5);
t1 = t12 * t7;
t2 = t12 * t9;
t11 = -t1 * t7 - t2 * t9;
t10 = cos(qJ(2));
t3 = -t9 * pkin(3) - pkin(2);
t5 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13 * t8 ^ 2 + t10 ^ 2; 0, 0, t10, -t8, 0, 0, 0, 0, 0, t10 * t9, -t10 * t7, t13 * t8, -t10 * t3 + t11 * t8; 0, 1, 0, 0, t4, t7 * t15, 0, 0, 0, pkin(2) * t15, -0.2e1 * pkin(2) * t7, 0.2e1 * t11, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, -t9 * t8, 0, -pkin(3) * t14; 0, 0, 0, 0, 0, 0, t7, t9, 0, -t7 * pkin(5), -t9 * pkin(5), -t7 * pkin(3), t1 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(3) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t5;
