% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2,theta3]';
% 
% Output:
% tauc_reg [4x15]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPPR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPPR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:37:56
% EndTime: 2019-12-31 16:37:57
% DurationCPUTime: 0.20s
% Computational Cost: add. (132->33), mult. (363->61), div. (0->0), fcn. (248->6), ass. (0->31)
t44 = 2 * qJD(3);
t20 = sin(pkin(6)) * pkin(1) + qJ(3);
t25 = sin(pkin(7));
t27 = cos(pkin(7));
t39 = t25 ^ 2 + t27 ^ 2;
t36 = qJD(1) * t39;
t43 = t20 * t36;
t29 = sin(qJ(4));
t30 = cos(qJ(4));
t16 = t30 * t25 + t29 * t27;
t10 = t16 * qJD(1);
t42 = pkin(5) + t20;
t41 = t29 * t25;
t40 = t30 * t27;
t34 = -t40 + t41;
t11 = t34 * qJD(4);
t38 = t11 * qJD(4);
t37 = qJD(1) * t41;
t17 = -cos(pkin(6)) * pkin(1) - t27 * pkin(3) - pkin(2);
t33 = t16 * qJD(3);
t12 = t16 * qJD(4);
t32 = t34 * qJD(3);
t19 = qJD(4) * qJD(1) * t40;
t14 = t42 * t27;
t13 = t42 * t25;
t9 = t34 * qJD(1);
t8 = t17 * qJD(1) + qJD(3);
t7 = t12 * qJD(4);
t6 = qJD(1) * t12;
t5 = -qJD(4) * t37 + t19;
t1 = [0, 0, 0, 0, 0, 0, t36 * t44, t43 * t44, -t10 * t11 + t5 * t16, -t10 * t12 + t11 * t9 - t16 * t6 - t5 * t34, -t38, -t7, 0, t8 * t12 + t17 * t6 + ((t13 * t29 - t14 * t30) * qJD(4) - t33) * qJD(4), -t8 * t11 + t17 * t5 + ((t13 * t30 + t14 * t29) * qJD(4) + t32) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, t38; 0, 0, 0, 0, 0, 0, -t39 * qJD(1) ^ 2, -t43 * qJD(1), 0, 0, 0, 0, 0, 0.2e1 * t10 * qJD(4), t19 + (-t9 - t37) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, t10 * t9, t10 ^ 2 - t9 ^ 2, t19 + (t9 - t37) * qJD(4), 0, 0, -qJD(1) * t33 - t8 * t10, qJD(1) * t32 + t8 * t9;];
tauc_reg = t1;
