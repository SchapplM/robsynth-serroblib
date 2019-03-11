% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4PRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRRP2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRRP2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:24:10
% EndTime: 2019-03-08 18:24:11
% DurationCPUTime: 0.16s
% Computational Cost: add. (212->39), mult. (514->58), div. (0->0), fcn. (358->4), ass. (0->37)
t24 = qJD(2) + qJD(3);
t28 = cos(qJ(2));
t40 = t28 * qJD(1);
t23 = qJD(2) * pkin(2) + t40;
t27 = cos(qJ(3));
t37 = qJD(2) * t40;
t25 = sin(qJ(3));
t26 = sin(qJ(2));
t43 = qJD(1) * t26;
t38 = t25 * t43;
t44 = t24 * t38;
t9 = (qJD(3) * t23 + t37) * t27 - t44;
t34 = -t25 * t26 + t27 * t28;
t12 = t24 * t34;
t45 = t12 * t24;
t42 = qJD(3) * t25;
t41 = qJD(3) * t27;
t35 = t25 * t28 + t27 * t26;
t32 = t35 * qJD(2);
t30 = (-t26 * t41 - t32) * qJD(1);
t10 = -t23 * t42 + t30;
t17 = t25 * t23 + t27 * t43;
t39 = t10 * t34 + t17 * t12 + t9 * t35;
t19 = t34 * qJD(1);
t36 = -t17 * t19 + (t17 * t41 + t25 * t9) * pkin(2);
t16 = t27 * t23 - t38;
t33 = (-pkin(2) * t24 - t23) * qJD(3);
t29 = qJD(2) ^ 2;
t18 = t35 * qJD(1);
t15 = t24 * pkin(3) + t16;
t13 = -t35 * qJD(3) - t32;
t11 = t13 * t24;
t4 = t17 * t24 + t10;
t3 = t16 * t24 - t9;
t2 = t19 * t24 + (t33 - t37) * t27 + t44;
t1 = t18 * t24 + t25 * t33 + t30;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29 * t26, -t29 * t28, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t45, 0, t16 * t13 + t39, 0, 0, 0, 0, 0, 0, t11, -t45, 0, t15 * t13 + t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t16 * t18 + (t10 * t27 - t16 * t42) * pkin(2) + t36, 0, 0, 0, 0, 0, 0, t1, t2, 0, t10 * (t27 * pkin(2) + pkin(3)) + (-pkin(2) * t42 + t18) * t15 + t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, 0, t10 * pkin(3) + (t15 - t16) * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
tauc_reg  = t5;
