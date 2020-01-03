% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPPR4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:38:54
% EndTime: 2019-12-31 16:38:55
% DurationCPUTime: 0.14s
% Computational Cost: add. (165->29), mult. (327->56), div. (0->0), fcn. (144->4), ass. (0->27)
t13 = sin(qJ(4));
t14 = cos(qJ(4));
t7 = -cos(pkin(6)) * pkin(1) - pkin(2) - pkin(5);
t5 = t7 * qJD(1) + qJD(3);
t3 = -t13 * qJD(2) + t14 * t5;
t1 = t3 * qJD(4);
t4 = t14 * qJD(2) + t13 * t5;
t2 = t4 * qJD(4);
t33 = (t13 * t3 - t14 * t4) * qJD(4) - t1 * t13 + t2 * t14;
t32 = 0.2e1 * qJD(3);
t31 = -t13 ^ 2 + t14 ^ 2;
t30 = t13 * t14;
t15 = qJD(4) ^ 2;
t29 = t15 * t13;
t28 = t15 * t14;
t16 = qJD(1) ^ 2;
t27 = -t15 - t16;
t8 = sin(pkin(6)) * pkin(1) + qJ(3);
t6 = qJD(1) * t8;
t26 = t6 * qJD(1);
t25 = qJD(4) * t13;
t24 = qJD(4) * t14;
t23 = qJD(1) * qJD(4);
t22 = t16 * t30;
t21 = t23 * t30;
t18 = t6 * t32;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(1) * t32, t18, -0.2e1 * t21, -0.2e1 * t31 * t23, -t29, 0.2e1 * t21, -t28, 0, t6 * t24 - t7 * t29 + (t13 * t32 + t8 * t24) * qJD(1), -t6 * t25 - t7 * t28 + (t14 * t32 - t8 * t25) * qJD(1), t33, -t33 * t7 + t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, t29, 0, t1 * t14 + t2 * t13 + (-t13 * t4 - t14 * t3) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t26, 0, 0, 0, 0, 0, 0, t27 * t13, t27 * t14, 0, -t33 - t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t31 * t16, 0, -t22, 0, 0, -t14 * t26, t13 * t26, 0, 0;];
tauc_reg = t9;
