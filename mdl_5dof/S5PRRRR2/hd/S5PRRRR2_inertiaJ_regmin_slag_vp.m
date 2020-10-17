% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x17]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRRR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_inertiaJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:04:51
% EndTime: 2019-12-05 17:04:51
% DurationCPUTime: 0.14s
% Computational Cost: add. (47->22), mult. (122->31), div. (0->0), fcn. (118->6), ass. (0->21)
t14 = sin(qJ(4));
t23 = t14 * pkin(3);
t22 = sin(qJ(3)) * pkin(2);
t17 = cos(qJ(4));
t10 = t17 * pkin(3);
t16 = cos(qJ(5));
t11 = cos(qJ(3)) * pkin(2);
t9 = t11 + pkin(3);
t3 = t14 * t22 - t17 * t9;
t21 = t3 * t16;
t13 = sin(qJ(5));
t20 = t13 * t10;
t19 = t17 * t22;
t4 = -t14 * t9 - t19;
t12 = t13 ^ 2;
t8 = pkin(6) + t23;
t7 = 0.2e1 * t13 * t16;
t6 = t16 * t10;
t2 = pkin(6) - t4;
t1 = t3 * t13;
t5 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 1, 0.2e1 * t11, -0.2e1 * t22, 1, -0.2e1 * t3, 0.2e1 * t4, t12, t7, 0, 0, 0, -0.2e1 * t21, 0.2e1 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, t11, -t22, 1, t10 - t3, -t19 + (-pkin(3) - t9) * t14, t12, t7, 0, 0, 0, t6 - t21, t1 - t20; 0, 0, 0, 0, 1, 0, 0, 1, 0.2e1 * t10, -0.2e1 * t23, t12, t7, 0, 0, 0, 0.2e1 * t6, -0.2e1 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 1, -t3, t4, t12, t7, 0, 0, 0, -t21, t1; 0, 0, 0, 0, 0, 0, 0, 1, t10, -t23, t12, t7, 0, 0, 0, t6, -t20; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t12, t7, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t16, 0, -t13 * t2, -t16 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t16, 0, -t13 * t8, -t16 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t16, 0, -t13 * pkin(6), -t16 * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t5;
