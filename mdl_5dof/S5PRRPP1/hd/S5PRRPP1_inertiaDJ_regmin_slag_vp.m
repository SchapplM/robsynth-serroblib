% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x19]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:23
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRPP1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP1_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:22:40
% EndTime: 2021-01-15 15:22:41
% DurationCPUTime: 0.23s
% Computational Cost: add. (303->53), mult. (764->102), div. (0->0), fcn. (629->4), ass. (0->36)
t48 = 2 * qJD(5);
t29 = sin(pkin(8));
t30 = sin(qJ(3));
t47 = t29 * t30;
t46 = -qJ(4) - pkin(6);
t45 = cos(pkin(8));
t31 = cos(qJ(3));
t39 = t45 * t30;
t20 = t29 * t31 + t39;
t44 = t20 * qJD(5);
t43 = t30 * qJD(3);
t42 = t31 * qJD(3);
t41 = -0.2e1 * pkin(2) * qJD(3);
t28 = pkin(3) * t43;
t27 = -t31 * pkin(3) - pkin(2);
t22 = t46 * t31;
t13 = -t29 * t22 - t46 * t39;
t14 = -t45 * t22 + t46 * t47;
t37 = qJD(3) * t46;
t16 = t31 * qJD(4) + t30 * t37;
t33 = -t30 * qJD(4) + t31 * t37;
t7 = t29 * t16 - t45 * t33;
t8 = t45 * t16 + t29 * t33;
t40 = t13 * t7 + t14 * t8;
t38 = t45 * t31;
t17 = t20 * qJD(3);
t18 = qJD(3) * t38 - t29 * t43;
t19 = -t38 + t47;
t36 = 0.2e1 * t19 * t17 + 0.2e1 * t20 * t18;
t35 = t17 * t13 + t18 * t14 + t19 * t7 + t20 * t8;
t32 = 0.2e1 * t13 * t18 - 0.2e1 * t14 * t17 - 0.2e1 * t8 * t19 + 0.2e1 * t7 * t20;
t26 = -t45 * pkin(3) - pkin(4);
t24 = t29 * pkin(3) + qJ(5);
t9 = t19 * pkin(4) - t20 * qJ(5) + t27;
t2 = t17 * pkin(4) - t18 * qJ(5) + t28 - t44;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, 0, 0, t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, 0, 0, t35; 0, 0, 0, 0, 0.2e1 * t30 * t42, 0.2e1 * (-t30 ^ 2 + t31 ^ 2) * qJD(3), 0, 0, 0, t30 * t41, t31 * t41, 0.2e1 * t27 * t17 + 0.2e1 * t19 * t28, 0.2e1 * t27 * t18 + 0.2e1 * t20 * t28, t32, 0.2e1 * t27 * t28 + 0.2e1 * t40, 0.2e1 * t9 * t17 + 0.2e1 * t2 * t19, t32, -0.2e1 * t9 * t18 - 0.2e1 * t2 * t20, 0.2e1 * t9 * t2 + 0.2e1 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t42, -t17, -t18, 0, (-t45 * t17 + t18 * t29) * pkin(3), -t17, 0, t18, t17 * t26 + t18 * t24 + t44; 0, 0, 0, 0, 0, 0, t42, -t43, 0, -pkin(6) * t42, pkin(6) * t43, -t7, -t8, (-t17 * t29 - t45 * t18) * pkin(3), (t29 * t8 - t45 * t7) * pkin(3), -t7, -qJD(5) * t19 - t24 * t17 + t26 * t18, t8, t14 * qJD(5) + t8 * t24 + t7 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t24 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t18, 0, t28, t17, 0, -t18, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
