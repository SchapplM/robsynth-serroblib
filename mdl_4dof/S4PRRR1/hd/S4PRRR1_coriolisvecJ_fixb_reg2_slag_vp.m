% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4PRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:45
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRRR1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:44:29
% EndTime: 2018-11-14 13:44:30
% DurationCPUTime: 0.16s
% Computational Cost: add. (151->38), mult. (388->64), div. (0->0), fcn. (198->4), ass. (0->34)
t17 = cos(qJ(4));
t18 = cos(qJ(3));
t34 = pkin(2) * qJD(2);
t30 = t18 * t34;
t28 = qJD(3) * t30;
t15 = sin(qJ(4));
t16 = sin(qJ(3));
t31 = t16 * t34;
t29 = t15 * t31;
t35 = (qJD(3) + qJD(4)) * t29;
t14 = qJD(2) + qJD(3);
t9 = t14 * pkin(3) + t30;
t1 = (qJD(4) * t9 + t28) * t17 - t35;
t37 = t15 * t16;
t36 = t16 * t17;
t33 = qJD(4) * t15;
t32 = qJD(4) * t17;
t27 = (-qJD(3) + t14) * t34;
t26 = pkin(2) * qJD(3) * (-qJD(2) - t14);
t13 = qJD(4) + t14;
t25 = (-pkin(3) * t13 - t9) * qJD(4);
t24 = -t15 * t18 - t36;
t23 = t17 * t18 - t37;
t20 = (t24 * qJD(3) - t16 * t32) * pkin(2);
t19 = qJD(2) * t20;
t2 = -t9 * t33 + t19;
t12 = t18 * pkin(2) + pkin(3);
t8 = t23 * t34;
t7 = t24 * t34;
t6 = t15 * t9 + t17 * t31;
t5 = t17 * t9 - t29;
t4 = -t12 * t33 + t20;
t3 = t12 * t32 + (t23 * qJD(3) - t16 * t33) * pkin(2);
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16 * t26, t18 * t26, 0, 0, 0, 0, 0, 0, 0, 0, t4 * t13 + t2, -t3 * t13 - t1, 0, t1 * (pkin(2) * t36 + t15 * t12) + t6 * t3 + t2 * (-pkin(2) * t37 + t17 * t12) + t5 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16 * t27, t18 * t27, 0, 0, 0, 0, 0, 0, 0, 0, -t7 * t13 + t15 * t25 + t19, t8 * t13 + (t25 - t28) * t17 + t35, 0, -t5 * t7 - t6 * t8 + (t1 * t15 + t2 * t17 + (-t15 * t5 + t17 * t6) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6 * t13 + t2, t5 * t13 - t1, 0, 0;];
tauc_reg  = t10;
