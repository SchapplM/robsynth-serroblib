% Calculate minimal parameter regressor of coriolis joint torque vector for
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
% tauc_reg [4x14]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPPR4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:38:54
% EndTime: 2019-12-31 16:38:54
% DurationCPUTime: 0.09s
% Computational Cost: add. (48->18), mult. (123->37), div. (0->0), fcn. (49->4), ass. (0->18)
t22 = 2 * qJD(3);
t10 = cos(qJ(4));
t9 = sin(qJ(4));
t21 = -t10 ^ 2 + t9 ^ 2;
t20 = t10 * t9;
t11 = qJD(4) ^ 2;
t19 = t11 * t9;
t18 = t11 * t10;
t12 = qJD(1) ^ 2;
t17 = -t11 - t12;
t4 = sin(pkin(6)) * pkin(1) + qJ(3);
t2 = qJD(1) * t4;
t16 = t2 * qJD(1);
t15 = t9 * qJD(4);
t14 = t10 * qJD(4);
t13 = qJD(1) * qJD(4);
t3 = -cos(pkin(6)) * pkin(1) - pkin(2) - pkin(5);
t1 = [0, 0, 0, 0, 0, qJD(1) * t22, t2 * t22, -0.2e1 * t13 * t20, 0.2e1 * t21 * t13, -t19, -t18, 0, t2 * t14 - t3 * t19 + (t4 * t14 + t9 * t22) * qJD(1), -t2 * t15 - t3 * t18 + (t10 * t22 - t4 * t15) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, t19; 0, 0, 0, 0, 0, -t12, -t16, 0, 0, 0, 0, 0, t17 * t9, t17 * t10; 0, 0, 0, 0, 0, 0, 0, t12 * t20, -t21 * t12, 0, 0, 0, -t10 * t16, t9 * t16;];
tauc_reg = t1;
