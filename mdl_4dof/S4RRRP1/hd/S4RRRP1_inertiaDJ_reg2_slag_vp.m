% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S4RRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRRP1_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP1_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP1_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP1_inertiaDJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:36:03
% EndTime: 2019-03-08 18:36:04
% DurationCPUTime: 0.15s
% Computational Cost: add. (124->29), mult. (380->48), div. (0->0), fcn. (229->4), ass. (0->31)
t23 = cos(qJ(2));
t19 = t23 * pkin(1) + pkin(2);
t20 = sin(qJ(3));
t21 = sin(qJ(2));
t22 = cos(qJ(3));
t33 = t21 * t22;
t12 = pkin(1) * t33 + t20 * t19;
t30 = qJD(3) * t22;
t25 = pkin(2) * t30;
t32 = pkin(1) * qJD(2);
t27 = t23 * t32;
t29 = t20 * t21 * pkin(1);
t7 = -t19 * t30 - t22 * t27 + (qJD(2) + qJD(3)) * t29;
t35 = -t7 * t20 * pkin(2) + t12 * t25;
t34 = t22 * pkin(2);
t31 = qJD(3) * t20;
t28 = t21 * t32;
t26 = pkin(2) * t31;
t11 = t22 * t19 - t29;
t24 = (-t21 * t30 + (-t20 * t23 - t33) * qJD(2)) * pkin(1);
t18 = pkin(3) + t34;
t17 = -0.2e1 * t25;
t16 = -0.2e1 * t26;
t10 = pkin(3) + t11;
t8 = -t19 * t31 + t24;
t6 = 0.2e1 * t8;
t5 = 0.2e1 * t7;
t3 = t7 - t25;
t2 = (-pkin(2) - t19) * t31 + t24;
t1 = t12 * t7;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t28, -0.2e1 * t27, 0, 0, 0, 0, 0, 0, 0, 0, t6, t5, 0, 0.2e1 * t11 * t8 - 0.2e1 * t1, 0, 0, 0, 0, 0, 0, t6, t5, 0, 0.2e1 * t10 * t8 - 0.2e1 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, -t27, 0, 0, 0, 0, 0, 0, 0, 0, t2, t3, 0 (-t11 * t31 + t22 * t8) * pkin(2) + t35, 0, 0, 0, 0, 0, 0, t2, t3, 0, -t10 * t26 + t8 * t18 + t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t17, 0, 0, 0, 0, 0, 0, 0, 0, t16, t17, 0, 0.2e1 * (-t18 + t34) * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t7, 0, 0, 0, 0, 0, 0, 0, 0, t8, t7, 0, t8 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, -t25, 0, 0, 0, 0, 0, 0, 0, 0, -t26, -t25, 0, -pkin(3) * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t4;
