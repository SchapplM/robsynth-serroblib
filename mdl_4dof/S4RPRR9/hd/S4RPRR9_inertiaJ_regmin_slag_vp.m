% Calculate minimal parameter regressor of joint inertia matrix for
% S4RPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% 
% Output:
% MM_reg [((4+1)*4/2)x20]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RPRR9_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_inertiaJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:28
% EndTime: 2019-12-31 16:56:29
% DurationCPUTime: 0.12s
% Computational Cost: add. (42->25), mult. (102->57), div. (0->0), fcn. (102->4), ass. (0->27)
t11 = cos(qJ(4));
t26 = 0.2e1 * t11;
t25 = 2 * qJ(2);
t9 = sin(qJ(4));
t24 = pkin(3) * t9;
t10 = sin(qJ(3));
t6 = t10 ^ 2;
t12 = cos(qJ(3));
t8 = t12 ^ 2;
t23 = -t6 - t8;
t22 = t9 * t10;
t21 = t9 * t11;
t20 = t9 * t12;
t13 = -pkin(1) - pkin(5);
t19 = t10 * t13;
t18 = t11 * t10;
t4 = t11 * t12;
t17 = t11 * t13;
t16 = t12 * t10;
t15 = t12 * t13;
t14 = -0.2e1 * t16;
t7 = t11 ^ 2;
t5 = t9 ^ 2;
t3 = t10 * pkin(3) - t12 * pkin(6) + qJ(2);
t2 = t10 * t17 + t9 * t3;
t1 = t11 * t3 - t9 * t19;
t27 = [1, 0, 0, -2 * pkin(1), t25, pkin(1) ^ 2 + qJ(2) ^ 2, t8, t14, 0, 0, 0, t10 * t25, t12 * t25, t7 * t8, -0.2e1 * t8 * t21, t16 * t26, t9 * t14, t6, -0.2e1 * t8 * t13 * t9 + 0.2e1 * t1 * t10, -0.2e1 * t2 * t10 - 0.2e1 * t8 * t17; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23 * t9, t23 * t11; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t12, -t10, 0, t15, -t19, t9 * t4, (-t5 + t7) * t12, t22, t18, 0, -pkin(6) * t22 + (t17 - t24) * t12, -t9 * t15 + (-pkin(3) * t12 - pkin(6) * t10) * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t10, 0, 0, 0, 0, 0, t4, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t5, 0.2e1 * t21, 0, 0, 0, pkin(3) * t26, -0.2e1 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t20, t10, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t11, 0, -t9 * pkin(6), -t11 * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t27;
