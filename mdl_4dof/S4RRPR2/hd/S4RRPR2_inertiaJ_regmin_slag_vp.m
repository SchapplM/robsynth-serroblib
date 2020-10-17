% Calculate minimal parameter regressor of joint inertia matrix for
% S4RRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% MM_reg [((4+1)*4/2)x12]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRPR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_inertiaJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-19 15:06:30
% EndTime: 2019-07-19 15:06:31
% DurationCPUTime: 0.11s
% Computational Cost: add. (55->27), mult. (73->30), div. (0->0), fcn. (55->4), ass. (0->17)
t19 = cos(qJ(2)) * pkin(1);
t14 = -pkin(2) - pkin(3);
t7 = pkin(2) + t19;
t5 = -pkin(3) - t7;
t18 = t14 + t5;
t9 = sin(qJ(2)) * pkin(1);
t6 = t9 + qJ(3);
t17 = qJ(3) + t6;
t16 = 0.2e1 * pkin(2);
t15 = 0.2e1 * qJ(3);
t12 = cos(qJ(4));
t10 = sin(qJ(4));
t4 = t12 * qJ(3) + t10 * t14;
t3 = t10 * qJ(3) - t12 * t14;
t2 = t10 * t5 + t12 * t6;
t1 = t10 * t6 - t12 * t5;
t8 = [1, 0, 0, 1, 0.2e1 * t19, -0.2e1 * t9, 0.2e1 * t7, 0.2e1 * t6, t6 ^ 2 + t7 ^ 2, 1, 0.2e1 * t1, 0.2e1 * t2; 0, 0, 0, 1, t19, -t9, t16 + t19, t15 + t9, t7 * pkin(2) + t6 * qJ(3), 1, t17 * t10 - t18 * t12, t18 * t10 + t17 * t12; 0, 0, 0, 1, 0, 0, t16, t15, pkin(2) ^ 2 + qJ(3) ^ 2, 1, 0.2e1 * t3, 0.2e1 * t4; 0, 0, 0, 0, 0, 0, -1, 0, -t7, 0, -t12, t10; 0, 0, 0, 0, 0, 0, -1, 0, -pkin(2), 0, -t12, t10; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t8;
