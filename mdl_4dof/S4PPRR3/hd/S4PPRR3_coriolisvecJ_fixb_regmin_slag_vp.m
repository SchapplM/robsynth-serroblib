% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4PPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% 
% Output:
% tauc_reg [4x12]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PPRR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:17:26
% EndTime: 2019-12-31 16:17:27
% DurationCPUTime: 0.09s
% Computational Cost: add. (34->11), mult. (115->30), div. (0->0), fcn. (58->4), ass. (0->18)
t10 = qJD(4) ^ 2;
t11 = qJD(3) ^ 2;
t7 = sin(qJ(3));
t24 = (t10 + t11) * t7;
t6 = sin(qJ(4));
t8 = cos(qJ(4));
t23 = t6 * t8;
t22 = t6 ^ 2 - t8 ^ 2;
t21 = t10 * t6;
t3 = t10 * t8;
t19 = (qJD(3) * pkin(3));
t17 = (qJD(3) * qJD(4));
t16 = 2 * t17;
t9 = cos(qJ(3));
t14 = -0.2e1 * t9 * t17;
t13 = t19 * qJD(3);
t12 = -2 * qJD(4) * t19;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t3; 0, 0, 0, -t11 * t7, -t11 * t9, 0, 0, 0, 0, 0, t6 * t14 - t8 * t24, t8 * t14 + t6 * t24; 0, 0, 0, 0, 0, t16 * t23, -t22 * t16, t3, -t21, 0, -pkin(5) * t3 + t6 * t12, pkin(5) * t21 + t8 * t12; 0, 0, 0, 0, 0, -t11 * t23, t22 * t11, 0, 0, 0, t6 * t13, t8 * t13;];
tauc_reg = t1;
