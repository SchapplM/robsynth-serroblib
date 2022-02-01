% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x18]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:28
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRP2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:28:12
% EndTime: 2022-01-23 09:28:13
% DurationCPUTime: 0.22s
% Computational Cost: add. (155->41), mult. (267->67), div. (0->0), fcn. (246->6), ass. (0->32)
t19 = sin(qJ(4));
t33 = 0.2e1 * t19;
t21 = cos(qJ(4));
t32 = -0.2e1 * t21;
t31 = 0.2e1 * t21;
t18 = cos(pkin(8));
t14 = t18 * pkin(1) + pkin(2);
t20 = sin(qJ(3));
t22 = cos(qJ(3));
t17 = sin(pkin(8));
t29 = pkin(1) * t17;
t7 = t22 * t14 - t20 * t29;
t5 = -pkin(3) - t7;
t30 = pkin(3) - t5;
t28 = t19 * pkin(4);
t27 = t21 * pkin(4);
t15 = -pkin(3) - t27;
t4 = t5 - t27;
t26 = t15 + t4;
t25 = -qJ(5) - pkin(7);
t8 = -t20 * t14 - t22 * t29;
t6 = pkin(7) - t8;
t24 = qJ(5) + t6;
t16 = t19 ^ 2;
t13 = t19 * t31;
t11 = t25 * t21;
t10 = t25 * t19;
t9 = t11 * t21;
t3 = t24 * t21;
t2 = t24 * t19;
t1 = t3 * t21;
t12 = [1, 0, 0, (t17 ^ 2 + t18 ^ 2) * pkin(1) ^ 2, 1, 0.2e1 * t7, 0.2e1 * t8, t16, t13, 0, 0, 0, t5 * t32, t5 * t33, t4 * t32, t4 * t33, 0.2e1 * t2 * t19 + 0.2e1 * t1, t2 ^ 2 + t3 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3 * t19 - t2 * t21; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21 ^ 2 + t16; 0, 0, 0, 0, 1, t7, t8, t16, t13, 0, 0, 0, t30 * t21, -t30 * t19, -t26 * t21, t26 * t19, t1 - t9 + (-t10 + t2) * t19, -t2 * t10 - t3 * t11 + t4 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21 * t10 - t19 * t11; 0, 0, 0, 0, 1, 0, 0, t16, t13, 0, 0, 0, pkin(3) * t31, -0.2e1 * pkin(3) * t19, t15 * t32, t15 * t33, -0.2e1 * t10 * t19 - 0.2e1 * t9, t10 ^ 2 + t11 ^ 2 + t15 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, t21, 0, -t19 * t6, -t21 * t6, -t2, -t3, -t28, -t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t19, t21, -t19, 0, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, t21, 0, -t19 * pkin(7), -t21 * pkin(7), t10, t11, -t28, t10 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, t19, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, t19, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t12;
