% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4PRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
% 
% Output:
% tauc_reg [4x14]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRPR4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:22:00
% EndTime: 2019-12-31 16:22:00
% DurationCPUTime: 0.08s
% Computational Cost: add. (31->14), mult. (93->31), div. (0->0), fcn. (32->2), ass. (0->16)
t10 = 2 * qJD(2);
t4 = sin(qJ(4));
t5 = cos(qJ(4));
t18 = t4 * t5;
t7 = qJD(4) ^ 2;
t17 = t7 * t4;
t16 = t7 * t5;
t15 = t4 ^ 2 - t5 ^ 2;
t8 = qJD(2) ^ 2;
t14 = -t7 - t8;
t13 = t8 * qJ(3);
t12 = qJ(3) * qJD(4);
t11 = qJD(2) * qJD(4);
t9 = qJD(3) * t10;
t6 = -pkin(2) - pkin(5);
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, t17; 0, 0, 0, 0, 0, t9, qJ(3) * t9, -0.2e1 * t11 * t18, 0.2e1 * t15 * t11, -t17, -t16, 0, -t6 * t17 + (qJD(3) * t4 + t5 * t12) * t10, -t6 * t16 + (qJD(3) * t5 - t4 * t12) * t10; 0, 0, 0, 0, 0, -t8, -t13, 0, 0, 0, 0, 0, t14 * t4, t14 * t5; 0, 0, 0, 0, 0, 0, 0, t8 * t18, -t15 * t8, 0, 0, 0, -t5 * t13, t4 * t13;];
tauc_reg = t1;
