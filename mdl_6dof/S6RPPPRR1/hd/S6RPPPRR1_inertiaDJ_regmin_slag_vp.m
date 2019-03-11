% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x24]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPPPRR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:30:14
% EndTime: 2019-03-09 01:30:14
% DurationCPUTime: 0.34s
% Computational Cost: add. (161->55), mult. (392->112), div. (0->0), fcn. (276->6), ass. (0->49)
t21 = sin(qJ(5));
t15 = t21 ^ 2;
t23 = cos(qJ(5));
t17 = t23 ^ 2;
t34 = (t15 - t17) * qJD(5);
t12 = t21 * qJD(5);
t11 = sin(pkin(9)) * pkin(1) + qJ(3);
t9 = -pkin(7) + t11;
t28 = t23 * qJD(3) - t9 * t12;
t31 = pkin(5) * t23 + pkin(8) * t21;
t53 = t31 * qJD(6) - t28;
t22 = cos(qJ(6));
t16 = t22 ^ 2;
t20 = sin(qJ(6));
t52 = t20 ^ 2 - t16;
t50 = t15 + t17;
t49 = qJD(6) * t17;
t48 = qJD(6) * t20;
t47 = qJD(6) * t21;
t13 = qJD(6) * t22;
t46 = qJD(6) * t23;
t45 = t11 * qJD(3);
t44 = t17 * qJD(3);
t43 = t23 * qJD(5);
t42 = -0.2e1 * pkin(5) * qJD(6);
t41 = t20 * t46;
t40 = t22 * t46;
t39 = t20 * t13;
t38 = t22 * t12;
t37 = t22 * t43;
t36 = t21 * t43;
t35 = t52 * qJD(6);
t33 = t20 * t38;
t10 = cos(pkin(9)) * pkin(1) + pkin(2) + qJ(4);
t32 = -t31 * qJD(5) + t9 * t47 - qJD(4);
t30 = t21 * pkin(5) - t23 * pkin(8);
t29 = -t21 * qJD(3) - t9 * t43;
t3 = t10 + t30;
t27 = t3 * t43 - t9 * t49;
t26 = -qJD(6) * t3 + t29;
t25 = t30 * qJD(5) - t9 * t46;
t24 = 0.2e1 * qJD(3);
t7 = t20 * t12 - t40;
t6 = t21 * t13 + t20 * t43;
t5 = t38 + t41;
t4 = t20 * t47 - t37;
t2 = t26 * t20 - t32 * t22;
t1 = t32 * t20 + t26 * t22;
t8 = [0, 0, 0, 0, 0, t24, 0.2e1 * t45, t24, 0.2e1 * qJD(4), 0.2e1 * t10 * qJD(4) + 0.2e1 * t45, -0.2e1 * t36, 0.2e1 * t34, 0, 0, 0, 0.2e1 * qJD(4) * t21 + 0.2e1 * t10 * t43, 0.2e1 * qJD(4) * t23 - 0.2e1 * t10 * t12, -0.2e1 * t16 * t36 - 0.2e1 * t17 * t39, 0.4e1 * t23 * t33 + 0.2e1 * t52 * t49, -0.2e1 * t21 * t41 - 0.2e1 * t22 * t34, 0.2e1 * t20 * t34 - 0.2e1 * t21 * t40, 0.2e1 * t36, 0.2e1 * t2 * t21 + 0.2e1 * t27 * t22 + 0.2e1 * (t9 * t36 - t44) * t20, -0.2e1 * t22 * t44 + 0.2e1 * (t9 * t37 + t1) * t21 - 0.2e1 * t27 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4), 0, 0, 0, 0, 0, -t43, t12, 0, 0, 0, 0, 0, t4, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50 * t13, t50 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, -t43, 0, t28, t29, -t23 * t35 - t33, t52 * t12 - 0.4e1 * t23 * t39, t6, -t4, 0, t25 * t20 - t53 * t22, t53 * t20 + t25 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, t12, 0, 0, 0, 0, 0, t4, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, -t43, 0, 0, 0, 0, 0, -t5, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t39, -0.2e1 * t35, 0, 0, 0, t20 * t42, t22 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, t7, t43, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t48, 0, -pkin(8) * t13, pkin(8) * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t8;
