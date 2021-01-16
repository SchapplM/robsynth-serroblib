% Calculate minimal parameter regressor of joint inertia matrix for
% S4RRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% MM_reg [((4+1)*4/2)x21]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:30
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRRP4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_inertiaJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:30:22
% EndTime: 2021-01-15 14:30:23
% DurationCPUTime: 0.13s
% Computational Cost: add. (115->26), mult. (237->56), div. (0->0), fcn. (252->4), ass. (0->24)
t18 = cos(qJ(2));
t13 = -t18 * pkin(2) - pkin(1);
t16 = sin(qJ(3));
t17 = sin(qJ(2));
t20 = cos(qJ(3));
t6 = t16 * t17 - t20 * t18;
t5 = t6 * pkin(3) + t13;
t25 = 0.2e1 * t5;
t24 = 0.2e1 * t13;
t23 = 0.2e1 * t18;
t22 = pkin(5) + pkin(6);
t21 = t16 * pkin(2);
t10 = t22 * t18;
t9 = t22 * t17;
t3 = -t16 * t10 - t20 * t9;
t4 = -t20 * t10 + t16 * t9;
t19 = 0.2e1 * pkin(3);
t15 = t20 * pkin(2);
t14 = -0.2e1 * t21;
t12 = t15 + pkin(3);
t7 = t16 * t18 + t20 * t17;
t2 = -t6 * qJ(4) - t4;
t1 = -t7 * qJ(4) + t3;
t8 = [1, 0, 0, t17 ^ 2, t17 * t23, 0, 0, 0, pkin(1) * t23, -0.2e1 * pkin(1) * t17, t7 ^ 2, -0.2e1 * t7 * t6, 0, 0, 0, t6 * t24, t7 * t24, t6 * t25, t7 * t25, -0.2e1 * t1 * t7 - 0.2e1 * t2 * t6, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t17, t18, 0, -t17 * pkin(5), -t18 * pkin(5), 0, 0, t7, -t6, 0, t3, t4, t1, -t2, -t12 * t7 - t6 * t21, t1 * t12 + t2 * t21; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t15, t14, 0.2e1 * t12, t14, 0, t16 ^ 2 * pkin(2) ^ 2 + t12 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t6, 0, t3, t4, t1, -t2, -t7 * pkin(3), t1 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t15, -t21, t19 + t15, -t21, 0, t12 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t19, 0, 0, pkin(3) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t7, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t8;
