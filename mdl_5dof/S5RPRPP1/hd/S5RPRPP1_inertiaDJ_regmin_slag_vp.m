% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x19]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:14
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPP1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP1_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:14:03
% EndTime: 2021-01-15 11:14:05
% DurationCPUTime: 0.25s
% Computational Cost: add. (381->55), mult. (842->106), div. (0->0), fcn. (707->6), ass. (0->38)
t51 = 2 * qJD(5);
t50 = cos(pkin(8));
t30 = sin(pkin(8));
t32 = sin(qJ(3));
t33 = cos(qJ(3));
t22 = t30 * t33 + t50 * t32;
t49 = t22 * qJD(5);
t48 = t32 * qJD(3);
t47 = t33 * qJD(3);
t46 = 0.2e1 * t47;
t29 = pkin(3) * t48;
t28 = -cos(pkin(7)) * pkin(1) - pkin(2);
t43 = sin(pkin(7)) * pkin(1) + pkin(6);
t41 = qJ(4) + t43;
t20 = t41 * t33;
t37 = t41 * t32;
t10 = t30 * t20 + t50 * t37;
t11 = t50 * t20 - t30 * t37;
t36 = qJD(3) * t41;
t16 = t33 * qJD(4) - t32 * t36;
t34 = -t32 * qJD(4) - t33 * t36;
t7 = t30 * t16 - t50 * t34;
t8 = t50 * t16 + t30 * t34;
t45 = t10 * t7 + t11 * t8;
t44 = t50 * t33;
t23 = -t33 * pkin(3) + t28;
t18 = t22 * qJD(3);
t19 = qJD(3) * t44 - t30 * t48;
t21 = t30 * t32 - t44;
t42 = 0.2e1 * t21 * t18 + 0.2e1 * t22 * t19;
t40 = t10 * t18 + t11 * t19 + t7 * t21 + t8 * t22;
t39 = t43 * qJD(3);
t35 = 0.2e1 * t10 * t19 - 0.2e1 * t11 * t18 - 0.2e1 * t8 * t21 + 0.2e1 * t7 * t22;
t27 = -t50 * pkin(3) - pkin(4);
t25 = t30 * pkin(3) + qJ(5);
t9 = t21 * pkin(4) - t22 * qJ(5) + t23;
t4 = t18 * pkin(4) - t19 * qJ(5) + t29 - t49;
t1 = [0, 0, 0, 0, t32 * t46, 0.2e1 * (-t32 ^ 2 + t33 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * t28 * t48, t28 * t46, 0.2e1 * t23 * t18 + 0.2e1 * t21 * t29, 0.2e1 * t23 * t19 + 0.2e1 * t22 * t29, t35, 0.2e1 * t23 * t29 + 0.2e1 * t45, 0.2e1 * t9 * t18 + 0.2e1 * t4 * t21, t35, -0.2e1 * t9 * t19 - 0.2e1 * t4 * t22, 0.2e1 * t9 * t4 + 0.2e1 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, 0, 0, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, t42; 0, 0, 0, 0, 0, 0, t47, -t48, 0, -t33 * t39, t32 * t39, -t7, -t8, (-t18 * t30 - t50 * t19) * pkin(3), (t30 * t8 - t50 * t7) * pkin(3), -t7, -qJD(5) * t21 - t25 * t18 + t27 * t19, t8, t11 * qJD(5) + t8 * t25 + t7 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, -t47, -t18, -t19, 0, (-t50 * t18 + t19 * t30) * pkin(3), -t18, 0, t19, t18 * t27 + t19 * t25 + t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, t25 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t19, 0, t29, t18, 0, -t19, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
