% Calculate minimal parameter regressor of joint inertia matrix for
% S4RRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% MM_reg [((4+1)*4/2)x21]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:46
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRPR6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:46:22
% EndTime: 2021-01-15 10:46:23
% DurationCPUTime: 0.14s
% Computational Cost: add. (135->27), mult. (287->67), div. (0->0), fcn. (325->6), ass. (0->29)
t21 = sin(pkin(7));
t22 = cos(pkin(7));
t24 = sin(qJ(2));
t26 = cos(qJ(2));
t13 = t21 * t24 - t22 * t26;
t20 = -t26 * pkin(2) - pkin(1);
t33 = 0.2e1 * t13 * pkin(3) + 0.2e1 * t20;
t32 = 0.2e1 * t20;
t31 = 0.2e1 * t26;
t30 = t21 * pkin(2);
t29 = t22 * pkin(2);
t28 = -qJ(3) - pkin(5);
t16 = t28 * t24;
t17 = t28 * t26;
t7 = t22 * t16 + t21 * t17;
t8 = t21 * t16 - t22 * t17;
t25 = cos(qJ(4));
t23 = sin(qJ(4));
t19 = pkin(3) + t29;
t14 = t21 * t26 + t22 * t24;
t11 = -t23 * t19 - t25 * t30;
t10 = t25 * t19 - t23 * t30;
t6 = -t23 * t13 + t25 * t14;
t5 = t25 * t13 + t23 * t14;
t4 = -t13 * pkin(6) + t8;
t3 = -t14 * pkin(6) + t7;
t2 = -t23 * t3 - t25 * t4;
t1 = -t23 * t4 + t25 * t3;
t9 = [1, 0, 0, t24 ^ 2, t24 * t31, 0, 0, 0, pkin(1) * t31, -0.2e1 * pkin(1) * t24, t13 * t32, t14 * t32, -0.2e1 * t8 * t13 - 0.2e1 * t7 * t14, t20 ^ 2 + t7 ^ 2 + t8 ^ 2, t6 ^ 2, -0.2e1 * t6 * t5, 0, 0, 0, t5 * t33, t6 * t33; 0, 0, 0, 0, 0, t24, t26, 0, -t24 * pkin(5), -t26 * pkin(5), t7, -t8, (-t13 * t21 - t14 * t22) * pkin(2), (t21 * t8 + t22 * t7) * pkin(2), 0, 0, t6, -t5, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * t29, -0.2e1 * t30, 0, (t21 ^ 2 + t22 ^ 2) * pkin(2) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t10, 0.2e1 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t14, 0, t20, 0, 0, 0, 0, 0, t5, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t10, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t9;
