% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PPRPR1
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
% MMD_reg [((5+1)*5/2)x16]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PPRPR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR1_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:01:19
% EndTime: 2019-12-05 15:01:20
% DurationCPUTime: 0.15s
% Computational Cost: add. (79->26), mult. (263->55), div. (0->0), fcn. (244->8), ass. (0->26)
t16 = sin(pkin(9));
t18 = cos(pkin(9));
t19 = sin(qJ(5));
t21 = cos(qJ(5));
t5 = t19 * t16 - t21 * t18;
t1 = t5 * qJD(5);
t17 = sin(pkin(8));
t20 = sin(qJ(3));
t26 = cos(pkin(8));
t32 = cos(qJ(3));
t6 = t20 * t17 - t32 * t26;
t33 = -0.2e1 * t1;
t28 = pkin(6) + qJ(4);
t27 = t16 ^ 2 + t18 ^ 2;
t3 = t6 * qJD(3);
t25 = t27 * t3;
t24 = t27 * qJD(4);
t23 = 0.2e1 * t24;
t7 = t21 * t16 + t19 * t18;
t2 = t7 * qJD(5);
t8 = t32 * t17 + t20 * t26;
t13 = -t18 * pkin(4) - pkin(3);
t10 = t28 * t18;
t9 = t28 * t16;
t4 = t8 * qJD(3);
t11 = [0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t8 * t25 + 0.2e1 * t6 * t4, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t4, t3, -t4 * t18, t4 * t16, -t25, -t4 * pkin(3) - qJ(4) * t25 + t8 * t24, 0, 0, 0, 0, 0, t6 * t2 + t4 * t5, -t6 * t1 + t4 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t23, qJ(4) * t23, t7 * t33, 0.2e1 * t1 * t5 - 0.2e1 * t7 * t2, 0, 0, 0, 0.2e1 * t13 * t2, t13 * t33; 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8 * t1 + t7 * t3, t8 * t2 - t5 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, -t2, 0, (-t10 * t21 + t19 * t9) * qJD(5) - t7 * qJD(4), (t10 * t19 + t21 * t9) * qJD(5) + t5 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t11;
