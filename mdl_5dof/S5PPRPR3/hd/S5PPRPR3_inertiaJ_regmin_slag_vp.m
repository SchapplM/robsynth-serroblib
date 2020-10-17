% Calculate minimal parameter regressor of joint inertia matrix for
% S5PPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x13]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PPRPR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:05:22
% EndTime: 2019-12-05 15:05:22
% DurationCPUTime: 0.12s
% Computational Cost: add. (45->25), mult. (108->46), div. (0->0), fcn. (143->8), ass. (0->17)
t15 = sin(qJ(5));
t20 = 0.2e1 * t15;
t11 = sin(pkin(9));
t13 = cos(pkin(9));
t16 = sin(qJ(3));
t18 = cos(qJ(3));
t6 = t11 * t18 + t13 * t16;
t4 = t11 * t16 - t13 * t18;
t17 = cos(qJ(5));
t14 = cos(pkin(8));
t12 = sin(pkin(8));
t10 = t14 ^ 2;
t9 = -t13 * pkin(3) - pkin(4);
t8 = t11 * pkin(3) + pkin(6);
t3 = t4 * t12;
t1 = t6 * t12;
t2 = [1, t12 ^ 2 + t10, 0, 0, 0, t1 ^ 2 + t3 ^ 2 + t10, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t1 * t4 - t3 * t6, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, t4 ^ 2 + t6 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t16 * t12, -t18 * t12, (-t1 * t13 - t11 * t3) * pkin(3), 0, 0, 0, 0, 0, -t1 * t17, t1 * t15; 0, 0, 0, t18, -t16, (t11 * t6 - t13 * t4) * pkin(3), 0, 0, 0, 0, 0, -t4 * t17, t4 * t15; 0, 0, 1, 0, 0, (t11 ^ 2 + t13 ^ 2) * pkin(3) ^ 2, t15 ^ 2, t17 * t20, 0, 0, 0, -0.2e1 * t9 * t17, t9 * t20; 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14 * t17 + t15 * t3, t14 * t15 + t17 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15 * t6, -t17 * t6; 0, 0, 0, 0, 0, 0, 0, 0, t15, t17, 0, -t15 * t8, -t17 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t2;
