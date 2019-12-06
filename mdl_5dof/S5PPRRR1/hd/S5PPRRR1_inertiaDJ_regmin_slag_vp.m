% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x15]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PPRRR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR1_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:12:50
% EndTime: 2019-12-05 15:12:51
% DurationCPUTime: 0.17s
% Computational Cost: add. (120->32), mult. (339->56), div. (0->0), fcn. (332->8), ass. (0->33)
t28 = cos(qJ(4));
t20 = -t28 * pkin(3) - pkin(4);
t27 = cos(qJ(5));
t21 = qJD(5) * t27;
t24 = sin(qJ(5));
t25 = sin(qJ(4));
t39 = qJD(4) * t25;
t34 = pkin(3) * t39;
t40 = t20 * t21 + t24 * t34;
t38 = qJD(4) * t28;
t37 = qJD(5) * t24;
t36 = pkin(4) * t37;
t35 = pkin(4) * t21;
t33 = pkin(3) * t38;
t22 = sin(pkin(9));
t23 = cos(pkin(9));
t26 = sin(qJ(3));
t29 = cos(qJ(3));
t13 = t29 * t22 + t26 * t23;
t32 = t26 * t22 - t29 * t23;
t6 = t28 * t13 - t25 * t32;
t31 = t20 * t37 - t27 * t34;
t30 = t13 * qJD(3);
t19 = t25 * pkin(3) + pkin(7);
t18 = 0.2e1 * t24 * t21;
t14 = 0.2e1 * (-t24 ^ 2 + t27 ^ 2) * qJD(5);
t11 = t32 * qJD(3);
t5 = t25 * t13 + t28 * t32;
t4 = t6 * qJD(4) - t25 * t11 + t28 * t30;
t3 = t28 * t11 + t13 * t39 + t25 * t30 + t32 * t38;
t2 = -t4 * t27 + t5 * t37;
t1 = t5 * t21 + t4 * t24;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t30, t11, 0, -t4, t3, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -0.2e1 * t34, -0.2e1 * t33, t18, t14, 0, 0, 0, 0.2e1 * t31, 0.2e1 * t40; 0, 0, 0, 0, 0, 0, -t4, t3, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -t34, -t33, t18, t14, 0, 0, 0, t31 - t36, -t35 + t40; 0, 0, 0, 0, 0, 0, 0, 0, t18, t14, 0, 0, 0, -0.2e1 * t36, -0.2e1 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6 * t21 + t24 * t3, t27 * t3 + t6 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, -t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t37, 0, -t19 * t21 - t24 * t33, t19 * t37 - t27 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t37, 0, -pkin(7) * t21, pkin(7) * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
