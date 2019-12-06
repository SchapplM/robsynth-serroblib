% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x17]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRRR3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR3_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:06:18
% EndTime: 2019-12-05 17:06:19
% DurationCPUTime: 0.24s
% Computational Cost: add. (128->44), mult. (392->64), div. (0->0), fcn. (242->6), ass. (0->36)
t20 = sin(qJ(4));
t21 = sin(qJ(3));
t24 = cos(qJ(3));
t23 = cos(qJ(4));
t37 = qJD(4) * t23;
t41 = t21 * t23;
t44 = ((t20 * t24 + t41) * qJD(3) + t21 * t37) * pkin(2);
t22 = cos(qJ(5));
t18 = qJD(5) * t22;
t19 = sin(qJ(5));
t17 = t24 * pkin(2) + pkin(3);
t38 = qJD(4) * t20;
t3 = t17 * t38 + t44;
t35 = t20 * t21 * pkin(2);
t6 = -t23 * t17 - pkin(4) + t35;
t42 = t6 * t18 + t3 * t19;
t16 = -t23 * pkin(3) - pkin(4);
t30 = pkin(3) * t38;
t40 = t16 * t18 + t19 * t30;
t39 = pkin(2) * qJD(3);
t36 = qJD(5) * t19;
t34 = pkin(4) * t36;
t33 = pkin(4) * t18;
t32 = t21 * t39;
t31 = t24 * t39;
t29 = pkin(3) * t37;
t4 = t6 * t36;
t27 = -t3 * t22 + t4;
t9 = t16 * t36;
t25 = -t22 * t30 + t9;
t2 = -t17 * t37 - t23 * t31 + (qJD(3) + qJD(4)) * t35;
t15 = t20 * pkin(3) + pkin(8);
t14 = 0.2e1 * t19 * t18;
t8 = 0.2e1 * (-t19 ^ 2 + t22 ^ 2) * qJD(5);
t7 = pkin(2) * t41 + t20 * t17 + pkin(8);
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -0.2e1 * t32, -0.2e1 * t31, 0, -0.2e1 * t3, 0.2e1 * t2, t14, t8, 0, 0, 0, 0.2e1 * t27, 0.2e1 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t32, -t31, 0, (-pkin(3) - t17) * t38 - t44, t2 - t29, t14, t8, 0, 0, 0, t4 + t9 + (-t3 - t30) * t22, t40 + t42; 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t30, -0.2e1 * t29, t14, t8, 0, 0, 0, 0.2e1 * t25, 0.2e1 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t3, t2, t14, t8, 0, 0, 0, t27 - t34, -t33 + t42; 0, 0, 0, 0, 0, 0, 0, 0, -t30, -t29, t14, t8, 0, 0, 0, t25 - t34, -t33 + t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t8, 0, 0, 0, -0.2e1 * t34, -0.2e1 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t36, 0, -t7 * t18 + t19 * t2, t22 * t2 + t7 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t36, 0, -t15 * t18 - t19 * t29, t15 * t36 - t22 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t36, 0, -pkin(8) * t18, pkin(8) * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
