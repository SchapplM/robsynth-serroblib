% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x13]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PPRPR3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR3_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:05:22
% EndTime: 2019-12-05 15:05:23
% DurationCPUTime: 0.14s
% Computational Cost: add. (65->30), mult. (229->68), div. (0->0), fcn. (226->8), ass. (0->26)
t13 = sin(pkin(9));
t15 = cos(pkin(9));
t18 = sin(qJ(3));
t20 = cos(qJ(3));
t8 = t13 * t20 + t15 * t18;
t5 = t8 * qJD(3);
t17 = sin(qJ(5));
t27 = t17 * qJD(5);
t26 = t18 * qJD(3);
t19 = cos(qJ(5));
t25 = t19 * qJD(5);
t24 = t20 * qJD(3);
t23 = 0.2e1 * t25;
t14 = sin(pkin(8));
t22 = t14 * t26;
t21 = t14 * t24;
t7 = t13 * t18 - t15 * t20;
t16 = cos(pkin(8));
t12 = -t15 * pkin(3) - pkin(4);
t11 = t13 * pkin(3) + pkin(6);
t6 = t7 * qJD(3);
t4 = t7 * t14;
t3 = t8 * t14;
t2 = t14 * t5;
t1 = t13 * t22 - t15 * t21;
t9 = [0, 0, 0, 0, 0, -0.2e1 * t3 * t1 + 0.2e1 * t4 * t2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t1 * t7 - t2 * t8 + t3 * t5 + t4 * t6, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0.2e1 * t7 * t5 - 0.2e1 * t8 * t6, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t21, t22, (t1 * t15 - t13 * t2) * pkin(3), 0, 0, 0, 0, 0, t1 * t19 + t3 * t27, -t1 * t17 + t3 * t25; 0, 0, 0, -t26, -t24, (-t13 * t6 - t15 * t5) * pkin(3), 0, 0, 0, 0, 0, -t5 * t19 + t7 * t27, t5 * t17 + t7 * t25; 0, 0, 0, 0, 0, 0, t17 * t23, 0.2e1 * (-t17 ^ 2 + t19 ^ 2) * qJD(5), 0, 0, 0, 0.2e1 * t12 * t27, t12 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17 * t2 + (t16 * t17 + t19 * t4) * qJD(5), t19 * t2 + (t16 * t19 - t17 * t4) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17 * t6 - t8 * t25, t19 * t6 + t8 * t27; 0, 0, 0, 0, 0, 0, 0, 0, t25, -t27, 0, -t11 * t25, t11 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, -t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
