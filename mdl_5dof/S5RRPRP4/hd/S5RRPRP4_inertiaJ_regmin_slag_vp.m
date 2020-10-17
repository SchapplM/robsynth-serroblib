% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x20]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRP4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:56
% EndTime: 2019-12-31 19:52:56
% DurationCPUTime: 0.20s
% Computational Cost: add. (134->37), mult. (196->58), div. (0->0), fcn. (144->4), ass. (0->28)
t21 = cos(qJ(4));
t17 = t21 ^ 2;
t19 = sin(qJ(4));
t7 = t19 ^ 2 + t17;
t30 = sin(qJ(2)) * pkin(1);
t11 = qJ(3) + t30;
t35 = 0.2e1 * t11;
t34 = 0.2e1 * t19;
t33 = -0.2e1 * t21;
t24 = 0.2e1 * qJ(3);
t4 = t19 * pkin(4) - t21 * qJ(5) + qJ(3);
t2 = t4 + t30;
t32 = t2 + t4;
t29 = cos(qJ(2)) * pkin(1);
t13 = -pkin(2) - t29;
t9 = -pkin(7) + t13;
t31 = t19 * t9;
t23 = -pkin(2) - pkin(7);
t27 = t19 * t23;
t26 = qJ(3) + t11;
t3 = t7 * t23;
t6 = t21 * pkin(4) + t19 * qJ(5);
t25 = -0.2e1 * pkin(2);
t14 = t21 * t23;
t10 = t19 * t33;
t5 = t21 * t9;
t1 = t7 * t9;
t8 = [1, 0, 0, 1, 0.2e1 * t29, -0.2e1 * t30, 0.2e1 * t13, t35, t11 ^ 2 + t13 ^ 2, t17, t10, 0, 0, 0, t11 * t34, t21 * t35, t2 * t34, -0.2e1 * t1, t2 * t33, t7 * t9 ^ 2 + t2 ^ 2; 0, 0, 0, 1, t29, -t30, t25 - t29, t24 + t30, -t13 * pkin(2) + t11 * qJ(3), t17, t10, 0, 0, 0, t26 * t19, t26 * t21, t32 * t19, t7 * (-t23 - t9), -t32 * t21, t2 * t4 + t9 * t3; 0, 0, 0, 1, 0, 0, t25, t24, pkin(2) ^ 2 + qJ(3) ^ 2, t17, t10, 0, 0, 0, t19 * t24, t21 * t24, t4 * t34, -0.2e1 * t3, t4 * t33, t7 * t23 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, 0, 1, 0, t13, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, t1; 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t19, 0, t5, -t31, t5, -t6, t31, t6 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t19, 0, t14, -t27, t14, -t6, t27, t6 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t19, t21, 0, t19, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t8;
