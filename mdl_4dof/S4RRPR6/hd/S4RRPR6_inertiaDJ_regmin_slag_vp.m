% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x21]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:46
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRPR6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR6_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:46:22
% EndTime: 2021-01-15 10:46:23
% DurationCPUTime: 0.22s
% Computational Cost: add. (300->50), mult. (746->115), div. (0->0), fcn. (655->6), ass. (0->42)
t37 = sin(pkin(7));
t48 = pkin(2) * t37;
t47 = -qJ(3) - pkin(5);
t40 = sin(qJ(2));
t30 = t47 * t40;
t42 = cos(qJ(2));
t31 = t47 * t42;
t38 = cos(pkin(7));
t14 = t37 * t30 - t38 * t31;
t46 = t40 * qJD(2);
t45 = t42 * qJD(2);
t44 = -0.2e1 * pkin(1) * qJD(2);
t36 = pkin(2) * t46;
t35 = -t42 * pkin(2) - pkin(1);
t43 = qJD(2) * t47;
t21 = t42 * qJD(3) + t40 * t43;
t22 = -t40 * qJD(3) + t42 * t43;
t7 = -t37 * t21 + t38 * t22;
t13 = t38 * t30 + t37 * t31;
t8 = t38 * t21 + t37 * t22;
t25 = t37 * t40 - t38 * t42;
t26 = t37 * t42 + t38 * t40;
t39 = sin(qJ(4));
t41 = cos(qJ(4));
t11 = t41 * t25 + t39 * t26;
t12 = -t39 * t25 + t41 * t26;
t34 = t38 * pkin(2) + pkin(3);
t24 = -t37 * t46 + t38 * t45;
t23 = t26 * qJD(2);
t19 = (-t34 * t39 - t41 * t48) * qJD(4);
t18 = (-t34 * t41 + t39 * t48) * qJD(4);
t16 = t25 * pkin(3) + t35;
t15 = t23 * pkin(3) + t36;
t10 = -t25 * pkin(6) + t14;
t9 = -t26 * pkin(6) + t13;
t6 = -t23 * pkin(6) + t8;
t5 = -t24 * pkin(6) + t7;
t4 = t12 * qJD(4) + t41 * t23 + t39 * t24;
t3 = -t11 * qJD(4) - t39 * t23 + t41 * t24;
t2 = -t39 * t6 + t41 * t5 + (-t10 * t41 - t39 * t9) * qJD(4);
t1 = -t39 * t5 - t41 * t6 + (t10 * t39 - t41 * t9) * qJD(4);
t17 = [0, 0, 0, 0.2e1 * t40 * t45, 0.2e1 * (-t40 ^ 2 + t42 ^ 2) * qJD(2), 0, 0, 0, t40 * t44, t42 * t44, 0.2e1 * t35 * t23 + 0.2e1 * t25 * t36, 0.2e1 * t35 * t24 + 0.2e1 * t26 * t36, -0.2e1 * t13 * t24 - 0.2e1 * t14 * t23 - 0.2e1 * t8 * t25 - 0.2e1 * t7 * t26, 0.2e1 * t13 * t7 + 0.2e1 * t14 * t8 + 0.2e1 * t35 * t36, 0.2e1 * t12 * t3, -0.2e1 * t3 * t11 - 0.2e1 * t12 * t4, 0, 0, 0, 0.2e1 * t15 * t11 + 0.2e1 * t16 * t4, 0.2e1 * t15 * t12 + 0.2e1 * t16 * t3; 0, 0, 0, 0, 0, t45, -t46, 0, -pkin(5) * t45, pkin(5) * t46, t7, -t8, (-t23 * t37 - t24 * t38) * pkin(2), (t37 * t8 + t38 * t7) * pkin(2), 0, 0, t3, -t4, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t19, 0.2e1 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t24, 0, t36, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t17;
