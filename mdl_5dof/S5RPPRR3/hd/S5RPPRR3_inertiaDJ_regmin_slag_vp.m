% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x21]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:15
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPRR3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:14:36
% EndTime: 2022-01-23 09:14:37
% DurationCPUTime: 0.25s
% Computational Cost: add. (312->48), mult. (716->84), div. (0->0), fcn. (704->8), ass. (0->40)
t28 = sin(pkin(9));
t29 = cos(pkin(9));
t32 = sin(qJ(4));
t33 = cos(qJ(4));
t35 = t32 * t28 - t33 * t29;
t25 = sin(pkin(8)) * pkin(1) + qJ(3);
t48 = pkin(6) + t25;
t17 = t48 * t28;
t18 = t48 * t29;
t37 = -t33 * t17 - t32 * t18;
t50 = t35 * qJD(3) - t37 * qJD(4);
t44 = qJD(4) * t33;
t45 = qJD(4) * t32;
t16 = t28 * t45 - t29 * t44;
t49 = -0.2e1 * t16;
t47 = cos(qJ(5));
t46 = t28 * t44 + t29 * t45;
t31 = sin(qJ(5));
t43 = qJD(5) * t31;
t42 = pkin(4) * t43;
t41 = pkin(4) * t46;
t40 = qJD(5) * t47;
t39 = pkin(4) * t40;
t38 = 0.2e1 * (t28 ^ 2 + t29 ^ 2) * qJD(3);
t36 = t32 * t17 - t33 * t18;
t20 = t33 * t28 + t32 * t29;
t21 = -cos(pkin(8)) * pkin(1) - t29 * pkin(3) - pkin(2);
t10 = t47 * t20 - t31 * t35;
t34 = -t20 * qJD(3) + t36 * qJD(4);
t11 = pkin(4) * t35 + t21;
t9 = t31 * t20 + t35 * t47;
t8 = -pkin(7) * t35 - t36;
t7 = -t20 * pkin(7) + t37;
t6 = t16 * pkin(7) + t34;
t5 = -t46 * pkin(7) - t50;
t4 = t10 * qJD(5) - t31 * t16 + t47 * t46;
t3 = t47 * t16 + t20 * t43 + t31 * t46 + t35 * t40;
t2 = t47 * t6 - t31 * t5 + (-t31 * t7 - t47 * t8) * qJD(5);
t1 = -t47 * t5 - t31 * t6 + (t31 * t8 - t47 * t7) * qJD(5);
t12 = [0, 0, 0, 0, 0, t38, t25 * t38, t20 * t49, 0.2e1 * t16 * t35 - 0.2e1 * t20 * t46, 0, 0, 0, 0.2e1 * t21 * t46, t21 * t49, -0.2e1 * t10 * t3, -0.2e1 * t10 * t4 + 0.2e1 * t3 * t9, 0, 0, 0, 0.2e1 * t11 * t4 + 0.2e1 * t41 * t9, 0.2e1 * t10 * t41 - 0.2e1 * t11 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t16, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t46, 0, t34, t50, 0, 0, -t3, -t4, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, t16, 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t42, -0.2e1 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, -t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t12;
