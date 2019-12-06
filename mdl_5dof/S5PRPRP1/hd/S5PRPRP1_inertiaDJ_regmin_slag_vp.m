% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x19]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRPRP1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP1_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP1_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP1_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:28:39
% EndTime: 2019-12-05 15:28:41
% DurationCPUTime: 0.18s
% Computational Cost: add. (196->45), mult. (509->83), div. (0->0), fcn. (453->4), ass. (0->27)
t34 = 2 * qJD(5);
t21 = sin(pkin(8));
t22 = cos(pkin(8));
t23 = sin(qJ(4));
t31 = t23 * t22;
t32 = cos(qJ(4));
t13 = t32 * t21 + t31;
t27 = qJD(4) * t32;
t29 = qJD(4) * t23;
t9 = t21 * t29 - t22 * t27;
t33 = t13 * t9;
t30 = pkin(6) + qJ(3);
t18 = -t22 * pkin(3) - pkin(2);
t28 = t30 * t21;
t26 = t32 * qJD(3);
t25 = 0.2e1 * (t21 ^ 2 + t22 ^ 2) * qJD(3);
t24 = t32 * t28;
t10 = t13 * qJD(4);
t1 = t10 * pkin(4) + t9 * qJ(5) - t13 * qJD(5);
t14 = t30 * t22;
t12 = t23 * t21 - t32 * t22;
t8 = t32 * t14 - t23 * t28;
t7 = t23 * t14 + t24;
t4 = t12 * pkin(4) - t13 * qJ(5) + t18;
t3 = t14 * t27 + qJD(3) * t31 + (-t30 * t29 + t26) * t21;
t2 = qJD(4) * t24 - t22 * t26 + (qJD(3) * t21 + qJD(4) * t14) * t23;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t12 * t10 - 0.2e1 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 * t7 + t12 * t3 - t13 * t2 - t9 * t8; 0, 0, 0, 0, 0, 0, t25, qJ(3) * t25, -0.2e1 * t33, -0.2e1 * t10 * t13 + 0.2e1 * t12 * t9, 0, 0, 0, 0.2e1 * t18 * t10, -0.2e1 * t18 * t9, 0.2e1 * t1 * t12 + 0.2e1 * t4 * t10, -0.2e1 * t8 * t10 + 0.2e1 * t2 * t12 + 0.2e1 * t3 * t13 - 0.2e1 * t7 * t9, -0.2e1 * t1 * t13 + 0.2e1 * t4 * t9, 0.2e1 * t4 * t1 - 0.2e1 * t8 * t2 + 0.2e1 * t7 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, t10, 0, t9, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t9, -t10, 0, -t9, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, -t10, 0, -t3, t2, -t3, pkin(4) * t9 - t10 * qJ(5) - t12 * qJD(5), -t2, -t3 * pkin(4) - t2 * qJ(5) + t8 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, qJ(5) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;
