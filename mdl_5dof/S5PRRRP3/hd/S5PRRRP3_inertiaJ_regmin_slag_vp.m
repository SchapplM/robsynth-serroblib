% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% MM_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:24
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRRP3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:23:27
% EndTime: 2021-01-15 16:23:29
% DurationCPUTime: 0.21s
% Computational Cost: add. (150->40), mult. (304->61), div. (0->0), fcn. (352->4), ass. (0->26)
t24 = cos(qJ(3));
t15 = -t24 * pkin(3) - pkin(2);
t18 = sin(qJ(4));
t19 = sin(qJ(3));
t23 = cos(qJ(4));
t7 = t18 * t19 - t23 * t24;
t26 = t7 * pkin(4);
t5 = t15 + t26;
t28 = 0.2e1 * t5;
t27 = 0.2e1 * t15;
t25 = t18 * pkin(3);
t22 = 0.2e1 * t24;
t21 = t24 * pkin(6);
t11 = (-pkin(6) - pkin(7)) * t19;
t12 = t24 * pkin(7) + t21;
t3 = t23 * t11 - t18 * t12;
t4 = -t18 * t11 - t23 * t12;
t20 = 0.2e1 * pkin(4);
t17 = t23 * pkin(3);
t16 = -0.2e1 * t25;
t14 = t17 + pkin(4);
t9 = t18 * t24 + t23 * t19;
t6 = t9 ^ 2;
t2 = -t7 * qJ(5) - t4;
t1 = -t9 * qJ(5) + t3;
t8 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 ^ 2 + t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7 * t1 + t9 * t2; 0, 1, 0, 0, t19 ^ 2, t19 * t22, 0, 0, 0, pkin(2) * t22, -0.2e1 * pkin(2) * t19, t6, -0.2e1 * t9 * t7, 0, 0, 0, t7 * t27, t9 * t27, t7 * t28, t9 * t28, -0.2e1 * t1 * t9 - 0.2e1 * t2 * t7, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t19, 0, 0, 0, 0, 0, -t7, -t9, -t7, -t9, 0, -t7 * t14 + t9 * t25; 0, 0, 0, 0, 0, 0, t19, t24, 0, -t19 * pkin(6), -t21, 0, 0, t9, -t7, 0, t3, t4, t1, -t2, -t14 * t9 - t7 * t25, t1 * t14 + t2 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t17, t16, 0.2e1 * t14, t16, 0, t18 ^ 2 * pkin(3) ^ 2 + t14 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t9, -t7, -t9, 0, -t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t7, 0, t3, t4, t1, -t2, -t9 * pkin(4), t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t17, -t25, t20 + t17, -t25, 0, t14 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t20, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t9, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t8;
