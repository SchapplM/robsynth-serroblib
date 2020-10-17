% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPPPRR5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:50:44
% EndTime: 2019-05-05 13:50:45
% DurationCPUTime: 0.34s
% Computational Cost: add. (153->59), mult. (217->97), div. (0->0), fcn. (238->6), ass. (0->43)
t23 = sin(qJ(5));
t44 = 0.2e1 * t23;
t24 = cos(qJ(6));
t43 = pkin(5) * t24;
t22 = sin(qJ(6));
t17 = sin(pkin(9));
t18 = cos(pkin(9));
t19 = pkin(3) + qJ(2);
t20 = pkin(1) + qJ(3);
t11 = t17 * t19 - t18 * t20;
t9 = pkin(7) + t11;
t42 = t22 * t9;
t41 = t22 * t17;
t40 = t22 * t18;
t39 = t22 * t23;
t38 = t22 * t24;
t25 = cos(qJ(5));
t37 = t22 * t25;
t36 = t23 * t17;
t35 = t23 * t18;
t34 = t24 * t17;
t33 = t24 * t18;
t32 = t24 * t23;
t31 = t24 * t25;
t30 = t25 * t17;
t29 = t25 * t18;
t28 = t25 * t44;
t10 = t17 * t20 + t18 * t19;
t8 = -pkin(4) - t10;
t27 = (qJ(2) ^ 2);
t26 = 2 * qJ(2);
t16 = t24 ^ 2;
t15 = t23 ^ 2;
t14 = t22 ^ 2;
t12 = t17 ^ 2 + t18 ^ 2;
t7 = t24 * t29 + t41;
t6 = t24 * t30 - t40;
t5 = -t22 * t29 + t34;
t4 = -t22 * t30 - t33;
t3 = -t25 * pkin(5) - t23 * pkin(8) + t8;
t2 = t22 * t3 + t9 * t31;
t1 = t24 * t3 - t9 * t37;
t13 = [1, 0, 0, -2 * pkin(1), t26, pkin(1) ^ 2 + t27, t26, 0.2e1 * t20, t20 ^ 2 + t27, 0.2e1 * t10, -0.2e1 * t11, t10 ^ 2 + t11 ^ 2, t15, t28, 0, 0, 0, -0.2e1 * t8 * t25, t8 * t44, t16 * t15, -0.2e1 * t15 * t38, -0.2e1 * t23 * t31, t22 * t28, t25 ^ 2, -0.2e1 * t1 * t25 + 0.2e1 * t15 * t42, 0.2e1 * t15 * t9 * t24 + 0.2e1 * t2 * t25; 0, 0, 0, 1, 0, -pkin(1), 0, -1, -t20, -t17, -t18, -t10 * t17 + t11 * t18, 0, 0, 0, 0, 0, -t30, t36, 0, 0, 0, 0, 0, t15 * t40 - t5 * t25, t15 * t33 + t7 * t25; 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 1, 0, qJ(2), t18, -t17, t10 * t18 + t11 * t17, 0, 0, 0, 0, 0, t29, -t35, 0, 0, 0, 0, 0, t15 * t41 - t4 * t25, t15 * t34 + t6 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t25, 0, -t23 * t9, -t25 * t9, t22 * t32 (-t14 + t16) * t23, -t37, -t31, 0, -t9 * t32 + (-pkin(5) * t23 + pkin(8) * t25) * t22, pkin(8) * t31 + (t42 - t43) * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, -t29, 0, 0, 0, 0, 0, -t18 * t32, t22 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t30, 0, 0, 0, 0, 0, -t17 * t32, t22 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t23, 0, 0, 0, 0, 0, t31, -t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t14, 0.2e1 * t38, 0, 0, 0, 0.2e1 * t43, -0.2e1 * pkin(5) * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t39, -t25, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t24, 0, -t22 * pkin(8), -t24 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t13;
