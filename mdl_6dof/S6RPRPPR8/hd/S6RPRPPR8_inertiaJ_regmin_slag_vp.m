% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
% 
% Output:
% MM_reg [((6+1)*6/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRPPR8_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:23:32
% EndTime: 2019-05-05 17:23:33
% DurationCPUTime: 0.41s
% Computational Cost: add. (166->59), mult. (234->90), div. (0->0), fcn. (221->4), ass. (0->42)
t26 = sin(qJ(3));
t45 = 0.2e1 * t26;
t27 = cos(qJ(6));
t44 = 0.2e1 * t27;
t43 = 2 * qJ(2);
t29 = -pkin(3) - pkin(4);
t25 = sin(qJ(6));
t42 = t25 * t26;
t41 = t25 * t27;
t28 = cos(qJ(3));
t40 = t25 * t28;
t30 = -pkin(1) - pkin(7);
t39 = t26 * t30;
t14 = t27 * t26;
t38 = t28 * t26;
t20 = t26 ^ 2;
t22 = t28 ^ 2;
t37 = t20 + t22;
t36 = t26 * qJ(4);
t35 = t28 * qJ(4) - qJ(2);
t18 = -pkin(8) + t29;
t34 = -0.2e1 * t38;
t12 = t28 * pkin(3) + t36;
t24 = qJ(4) + pkin(5);
t33 = -t18 * t28 + t24 * t26;
t32 = qJ(4) ^ 2;
t31 = 0.2e1 * qJ(4);
t21 = t27 ^ 2;
t19 = t25 ^ 2;
t16 = t28 * t30;
t15 = t27 * t28;
t11 = t26 * pkin(3) - t35;
t10 = -t29 * t28 + t36;
t9 = -t28 * qJ(5) - t16;
t7 = (-qJ(5) - t30) * t26;
t6 = t37 * t30;
t5 = t29 * t26 + t35;
t4 = t28 * pkin(5) + t18 * t26 + t35;
t3 = -t7 * t26 - t9 * t28;
t2 = t25 * t4 + t27 * t9;
t1 = -t25 * t9 + t27 * t4;
t8 = [1, 0, 0, -2 * pkin(1), t43, pkin(1) ^ 2 + qJ(2) ^ 2, t22, t34, 0, 0, 0, t26 * t43, t28 * t43, t11 * t45, -0.2e1 * t6, -0.2e1 * t11 * t28, t37 * t30 ^ 2 + t11 ^ 2, 0.2e1 * t5 * t28, t5 * t45, 0.2e1 * t3, t5 ^ 2 + t7 ^ 2 + t9 ^ 2, t21 * t20, -0.2e1 * t20 * t41, t38 * t44, t25 * t34, t22, 0.2e1 * t1 * t28 - 0.2e1 * t7 * t42, -0.2e1 * t7 * t14 - 0.2e1 * t2 * t28; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, -t37, 0, t6, 0, 0, t37, t3, 0, 0, 0, 0, 0, t37 * t25, t37 * t27; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, 0, 0, t37, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t28, -t26, 0, t16, -t39, t16, -t12, t39, t12 * t30, -t7, t9, t10, -t7 * qJ(4) + t9 * t29, -t25 * t14 (t19 - t21) * t26, -t40, -t15, 0, t33 * t25 - t7 * t27, t7 * t25 + t33 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t26, t28, 0, t26, t12, t26, -t28, 0, t10, 0, 0, 0, 0, 0, t14, -t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, t31, pkin(3) ^ 2 + t32, t31, 0.2e1 * t29, 0, t29 ^ 2 + t32, t19, 0.2e1 * t41, 0, 0, 0, t24 * t44, -0.2e1 * t24 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, -t16, 0, 0, -t28, t9, 0, 0, 0, 0, 0, -t40, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, 0, 0, 0, -t28, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3), 0, 1, 0, t29, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t26, 0, t5, 0, 0, 0, 0, 0, t15, -t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t42, t28, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, -t27, 0, -t25 * t18, -t27 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, -t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t8;
