% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x17]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPPRR5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR5_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR5_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:56:38
% EndTime: 2019-12-31 17:56:39
% DurationCPUTime: 0.13s
% Computational Cost: add. (69->28), mult. (92->36), div. (0->0), fcn. (90->6), ass. (0->22)
t13 = sin(qJ(5));
t23 = -0.2e1 * t13;
t15 = cos(qJ(5));
t22 = 0.2e1 * t15;
t14 = sin(qJ(4));
t16 = cos(qJ(4));
t12 = cos(pkin(8));
t8 = t12 * pkin(1) + pkin(2);
t5 = -pkin(3) - t8;
t11 = sin(pkin(8));
t7 = t11 * pkin(1) + qJ(3);
t3 = t14 * t7 - t16 * t5;
t1 = pkin(4) + t3;
t21 = pkin(4) + t1;
t20 = t13 * t15;
t19 = t16 * t13;
t18 = t16 * t15;
t4 = t14 * t5 + t16 * t7;
t10 = t13 ^ 2;
t6 = 0.2e1 * t20;
t2 = -pkin(7) + t4;
t9 = [1, 0, 0, (t11 ^ 2 + t12 ^ 2) * pkin(1) ^ 2, 0.2e1 * t8, 0.2e1 * t7, t7 ^ 2 + t8 ^ 2, 1, 0.2e1 * t3, 0.2e1 * t4, t10, t6, 0, 0, 0, t1 * t22, t1 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -1, 0, -t8, 0, -t16, t14, 0, 0, 0, 0, 0, -t18, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, -1, -t3, -t4, -t10, -0.2e1 * t20, 0, 0, 0, -t21 * t15, t21 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t16, -t14, 0, 0, 0, 0, 0, t18, -t19; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t10, t6, 0, 0, 0, pkin(4) * t22, pkin(4) * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t15, 0, -t13 * t2, -t15 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13 * t14, -t15 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t15, 0, -t13 * pkin(7), -t15 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t9;
