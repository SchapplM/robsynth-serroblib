% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4PRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% 
% Output:
% tauc_reg [4x15]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRPR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:20:54
% EndTime: 2019-12-31 16:20:55
% DurationCPUTime: 0.18s
% Computational Cost: add. (108->31), mult. (323->59), div. (0->0), fcn. (224->4), ass. (0->30)
t41 = 2 * qJD(3);
t23 = sin(pkin(7));
t24 = cos(pkin(7));
t36 = t23 ^ 2 + t24 ^ 2;
t32 = t36 * qJD(2);
t40 = qJ(3) * t32;
t25 = sin(qJ(4));
t26 = cos(qJ(4));
t10 = t26 * t23 + t25 * t24;
t7 = t10 * qJD(2);
t39 = t25 * t23;
t38 = t26 * t24;
t37 = pkin(5) + qJ(3);
t30 = -t38 + t39;
t8 = t30 * qJD(4);
t35 = t8 * qJD(4);
t33 = qJD(2) * t39;
t18 = -t24 * pkin(3) - pkin(2);
t29 = t10 * qJD(3);
t9 = t10 * qJD(4);
t28 = t30 * qJD(3);
t16 = qJD(4) * qJD(2) * t38;
t15 = t37 * t24;
t14 = t37 * t23;
t13 = t18 * qJD(2) + qJD(3);
t6 = t30 * qJD(2);
t3 = t9 * qJD(4);
t2 = qJD(2) * t9;
t1 = -qJD(4) * t33 + t16;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, t35; 0, 0, 0, 0, 0, 0, t32 * t41, t40 * t41, t1 * t10 - t7 * t8, -t1 * t30 - t10 * t2 + t8 * t6 - t7 * t9, -t35, -t3, 0, t13 * t9 + t18 * t2 + ((t14 * t25 - t15 * t26) * qJD(4) - t29) * qJD(4), t18 * t1 - t13 * t8 + ((t14 * t26 + t15 * t25) * qJD(4) + t28) * qJD(4); 0, 0, 0, 0, 0, 0, -t36 * qJD(2) ^ 2, -t40 * qJD(2), 0, 0, 0, 0, 0, 0.2e1 * t7 * qJD(4), t16 + (-t6 - t33) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, t7 * t6, -t6 ^ 2 + t7 ^ 2, t16 + (t6 - t33) * qJD(4), 0, 0, -qJD(2) * t29 - t13 * t7, qJD(2) * t28 + t13 * t6;];
tauc_reg = t4;
