% Calculate inertial parameters regressor of joint inertia matrix for
% S4RRRP1
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
% MM_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRRP1_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP1_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP1_inertiaJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:24:40
% EndTime: 2019-05-04 19:24:40
% DurationCPUTime: 0.15s
% Computational Cost: add. (67->25), mult. (153->36), div. (0->0), fcn. (111->4), ass. (0->24)
t15 = sin(qJ(3));
t25 = t15 * pkin(2);
t16 = sin(qJ(2));
t24 = t16 * pkin(1);
t17 = cos(qJ(3));
t13 = t17 * pkin(2);
t23 = t17 * t24;
t18 = cos(qJ(2));
t14 = t18 * pkin(1);
t10 = t14 + pkin(2);
t6 = t17 * t10 - t15 * t24;
t19 = 2 * pkin(3);
t22 = t19 + t6;
t20 = pkin(2) ^ 2;
t12 = t15 ^ 2 * t20;
t11 = -0.2e1 * t25;
t9 = t13 + pkin(3);
t7 = t15 * t10 + t23;
t5 = t7 ^ 2;
t4 = pkin(3) + t6;
t3 = 0.2e1 * t7;
t2 = t7 * t25;
t1 = -t23 + (-pkin(2) - t10) * t15;
t8 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t14, -0.2e1 * t24, 0 (t16 ^ 2 + t18 ^ 2) * pkin(1) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t6, -t3, 0, t6 ^ 2 + t5, 0, 0, 0, 0, 0, 1, 0.2e1 * t4, -t3, 0, t4 ^ 2 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t14, -t24, 0, 0, 0, 0, 0, 0, 0, 1, t13 + t6, t1, 0, t6 * t13 + t2, 0, 0, 0, 0, 0, 1, t13 + t22, t1, 0, t4 * t9 + t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t13, t11, 0, t17 ^ 2 * t20 + t12, 0, 0, 0, 0, 0, 1, 0.2e1 * t9, t11, 0, t9 ^ 2 + t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t6, -t7, 0, 0, 0, 0, 0, 0, 0, 1, t22, -t7, 0, t4 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t13, -t25, 0, 0, 0, 0, 0, 0, 0, 1, t19 + t13, -t25, 0, t9 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t19, 0, 0, pkin(3) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t8;
