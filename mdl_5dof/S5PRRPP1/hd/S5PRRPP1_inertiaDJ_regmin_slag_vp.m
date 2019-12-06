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
% MMD_reg [((5+1)*5/2)x17]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 16:06:53
% EndTime: 2019-12-05 16:06:54
% DurationCPUTime: 0.22s
% Computational Cost: add. (281->47), mult. (698->94), div. (0->0), fcn. (577->4), ass. (0->36)
t48 = 2 * qJD(5);
t28 = sin(pkin(8));
t29 = sin(qJ(3));
t47 = t28 * t29;
t46 = -qJ(4) - pkin(6);
t45 = cos(pkin(8));
t30 = cos(qJ(3));
t38 = t45 * t29;
t20 = t28 * t30 + t38;
t44 = t20 * qJD(5);
t43 = t29 * qJD(3);
t42 = t30 * qJD(3);
t41 = -0.2e1 * pkin(2) * qJD(3);
t27 = pkin(3) * t43;
t40 = -t30 * pkin(3) - pkin(2);
t22 = t46 * t30;
t13 = -t28 * t22 - t46 * t38;
t14 = -t45 * t22 + t46 * t47;
t36 = qJD(3) * t46;
t16 = t30 * qJD(4) + t29 * t36;
t32 = -t29 * qJD(4) + t30 * t36;
t7 = t28 * t16 - t45 * t32;
t8 = t45 * t16 + t28 * t32;
t39 = t13 * t7 + t14 * t8;
t37 = t45 * t30;
t17 = t20 * qJD(3);
t18 = qJD(3) * t37 - t28 * t43;
t19 = -t37 + t47;
t35 = 0.2e1 * t19 * t17 + 0.2e1 * t20 * t18;
t34 = t17 * t13 + t18 * t14 + t19 * t7 + t20 * t8;
t31 = 0.2e1 * t13 * t18 - 0.2e1 * t14 * t17 - 0.2e1 * t8 * t19 + 0.2e1 * t7 * t20;
t26 = -t45 * pkin(3) - pkin(4);
t24 = t28 * pkin(3) + qJ(5);
t9 = t19 * pkin(4) - t20 * qJ(5) + t40;
t2 = t17 * pkin(4) - t18 * qJ(5) + t27 - t44;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, 0, 0, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, 0, 0, t34; 0, 0, 0, 0, 0.2e1 * t29 * t42, 0.2e1 * (-t29 ^ 2 + t30 ^ 2) * qJD(3), 0, 0, 0, t29 * t41, t30 * t41, t31, 0.2e1 * t27 * t40 + 0.2e1 * t39, 0.2e1 * t9 * t17 + 0.2e1 * t2 * t19, t31, -0.2e1 * t9 * t18 - 0.2e1 * t2 * t20, 0.2e1 * t9 * t2 + 0.2e1 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t42, 0, (-t45 * t17 + t18 * t28) * pkin(3), -t17, 0, t18, t17 * t26 + t18 * t24 + t44; 0, 0, 0, 0, 0, 0, t42, -t43, 0, -pkin(6) * t42, pkin(6) * t43, (-t17 * t28 - t45 * t18) * pkin(3), (t28 * t8 - t45 * t7) * pkin(3), -t7, -qJD(5) * t19 - t24 * t17 + t26 * t18, t8, t14 * qJD(5) + t8 * t24 + t7 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t24 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t17, 0, -t18, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
