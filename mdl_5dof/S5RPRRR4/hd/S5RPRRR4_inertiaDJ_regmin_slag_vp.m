% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x17]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:35
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRR4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:34:46
% EndTime: 2022-01-23 09:34:46
% DurationCPUTime: 0.23s
% Computational Cost: add. (278->47), mult. (654->68), div. (0->0), fcn. (470->8), ass. (0->42)
t22 = cos(pkin(9)) * pkin(1) + pkin(2);
t29 = sin(qJ(3));
t32 = cos(qJ(3));
t51 = pkin(1) * sin(pkin(9));
t53 = -t32 * t22 + t29 * t51;
t30 = cos(qJ(5));
t25 = t30 * qJD(5);
t27 = sin(qJ(5));
t14 = pkin(3) - t53;
t28 = sin(qJ(4));
t15 = t29 * t22 + t32 * t51;
t31 = cos(qJ(4));
t48 = t31 * t15;
t34 = t28 * t14 + t48;
t12 = t53 * qJD(3);
t13 = t15 * qJD(3);
t35 = -t28 * t12 + t31 * t13;
t3 = t34 * qJD(4) + t35;
t49 = t28 * t15;
t6 = -t31 * t14 - pkin(4) + t49;
t52 = t6 * t25 + t3 * t27;
t50 = t31 * pkin(3);
t24 = -pkin(4) - t50;
t45 = qJD(4) * t28;
t38 = pkin(3) * t45;
t46 = t24 * t25 + t27 * t38;
t44 = qJD(4) * t31;
t43 = t27 * qJD(5);
t42 = t31 * t12 + t28 * t13 - t14 * t44;
t40 = pkin(4) * t43;
t39 = pkin(4) * t25;
t37 = pkin(3) * t44;
t4 = t6 * t43;
t36 = -t3 * t30 + t4;
t18 = t24 * t43;
t33 = -t30 * t38 + t18;
t23 = t28 * pkin(3) + pkin(8);
t21 = 0.2e1 * t27 * t25;
t16 = 0.2e1 * (-t27 ^ 2 + t30 ^ 2) * qJD(5);
t7 = pkin(8) + t34;
t2 = t15 * t45 + t42;
t1 = [0, 0, 0, 0, 0, -0.2e1 * t13, 0.2e1 * t12, 0, -0.2e1 * t3, 0.2e1 * t2, t21, t16, 0, 0, 0, 0.2e1 * t36, 0.2e1 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t13, t12, 0, (-t48 + (-pkin(3) - t14) * t28) * qJD(4) - t35, (t49 - t50) * qJD(4) + t42, t21, t16, 0, 0, 0, t18 + t4 + (-t3 - t38) * t30, t46 + t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t38, -0.2e1 * t37, t21, t16, 0, 0, 0, 0.2e1 * t33, 0.2e1 * t46; 0, 0, 0, 0, 0, 0, 0, 0, -t3, t2, t21, t16, 0, 0, 0, t36 - t40, -t39 + t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t38, -t37, t21, t16, 0, 0, 0, t33 - t40, -t39 + t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t16, 0, 0, 0, -0.2e1 * t40, -0.2e1 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t43, 0, t27 * t2 - t7 * t25, t30 * t2 + t7 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t43, 0, -t23 * t25 - t27 * t37, t23 * t43 - t30 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t43, 0, -pkin(8) * t25, pkin(8) * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
