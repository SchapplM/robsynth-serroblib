% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x18]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 17:13
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPPRP4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 17:13:15
% EndTime: 2021-01-15 17:13:16
% DurationCPUTime: 0.15s
% Computational Cost: add. (114->35), mult. (140->49), div. (0->0), fcn. (137->4), ass. (0->26)
t19 = sin(qJ(4));
t29 = -0.2e1 * t19;
t20 = cos(qJ(4));
t28 = 0.2e1 * t20;
t17 = sin(pkin(7));
t18 = cos(pkin(7));
t21 = -pkin(1) - pkin(2);
t7 = t18 * qJ(2) + t17 * t21;
t27 = t18 * t20;
t26 = t19 * t17;
t25 = t20 * t17;
t5 = -pkin(6) + t7;
t24 = qJ(5) - t5;
t15 = t19 ^ 2;
t23 = t20 ^ 2 + t15;
t6 = -t17 * qJ(2) + t18 * t21;
t4 = pkin(3) - t6;
t1 = t24 * t19;
t2 = t24 * t20;
t22 = t1 * t19 + t2 * t20;
t14 = t18 ^ 2;
t13 = t17 ^ 2;
t12 = t20 * pkin(4);
t8 = t19 * t18;
t3 = t12 + t4;
t9 = [1, 0, 0, 2 * pkin(1), 0.2e1 * qJ(2), (pkin(1) ^ 2) + qJ(2) ^ 2, t6 ^ 2 + t7 ^ 2, t15, t19 * t28, 0, 0, 0, t4 * t28, t4 * t29, t3 * t28, t3 * t29, 0.2e1 * t22, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, -1, 0, -pkin(1), t7 * t17 + t6 * t18, 0, 0, 0, 0, 0, -t27, t8, -t27, t8, -t23 * t17, -t22 * t17 - t3 * t18; 0, 0, 0, 0, 0, 1, t13 + t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23 * t13 + t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 * t20 - t2 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t20, 0, -t19 * t5, -t20 * t5, t1, t2, t19 * pkin(4), t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, -t25, -t26, -t25, 0, -pkin(4) * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, t20, -t19, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t9;
