% Calculate inertial parameters regressor of joint inertia matrix for
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
% MM_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRPR2_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_inertiaJ_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-19 15:06:30
% EndTime: 2019-07-19 15:06:34
% DurationCPUTime: 0.20s
% Computational Cost: add. (98->36), mult. (127->47), div. (0->0), fcn. (99->4), ass. (0->20)
t23 = -pkin(2) - pkin(3);
t16 = cos(qJ(2));
t22 = t16 * pkin(1);
t14 = sin(qJ(2));
t12 = t14 * pkin(1);
t8 = t12 + qJ(3);
t21 = qJ(3) + t8;
t9 = pkin(2) + t22;
t20 = -pkin(3) - t9;
t18 = 0.2e1 * pkin(2);
t17 = 0.2e1 * qJ(3);
t15 = cos(qJ(4));
t13 = sin(qJ(4));
t11 = t15 * t23;
t7 = t15 * t20;
t6 = t15 * qJ(3) + t13 * t23;
t4 = t13 * qJ(3) - t11;
t3 = t13 * t20 + t15 * t8;
t1 = t13 * t8 - t7;
t2 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t22, -0.2e1 * t12, 0, (t14 ^ 2 + t16 ^ 2) * pkin(1) ^ 2, 0, 0, 0, 1, 0, 0, 0.2e1 * t9, 0, 0.2e1 * t8, t8 ^ 2 + t9 ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t1, 0.2e1 * t3, 0, t1 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t22, -t12, 0, 0, 0, 0, 0, 1, 0, 0, t18 + t22, 0, t17 + t12, t9 * pkin(2) + t8 * qJ(3), 0, 0, 0, 0, 0, 1, t21 * t13 - t11 - t7, (-0.2e1 * pkin(2) - 0.2e1 * pkin(3) - t22) * t13 + t21 * t15, 0, t1 * t4 + t3 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t18, 0, t17, pkin(2) ^ 2 + qJ(3) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t4, 0.2e1 * t6, 0, t4 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t9, 0, 0, 0, 0, 0, 0, -t15, t13, 0, -t1 * t15 + t3 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(2), 0, 0, 0, 0, 0, 0, -t15, t13, 0, t6 * t13 - t4 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13 ^ 2 + t15 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t1, -t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t4, -t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t13, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t2;
