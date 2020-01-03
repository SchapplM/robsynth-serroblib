% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4PPRR5
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
% Datum: 2019-12-31 16:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PPRR5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:19:52
% EndTime: 2019-12-31 16:19:52
% DurationCPUTime: 0.11s
% Computational Cost: add. (49->14), mult. (155->38), div. (0->0), fcn. (86->4), ass. (0->21)
t10 = cos(qJ(4));
t8 = sin(qJ(4));
t28 = -t10 ^ 2 + t8 ^ 2;
t13 = qJD(3) ^ 2;
t9 = sin(qJ(3));
t27 = t9 * t13;
t11 = cos(qJ(3));
t26 = t11 * t13;
t12 = qJD(4) ^ 2;
t25 = t12 + t13;
t24 = (qJD(3) * pkin(3));
t23 = (qJD(3) * qJD(4));
t22 = 2 * t23;
t21 = t25 * t8;
t20 = t10 * t25;
t19 = t24 * qJD(3);
t18 = t10 * t22;
t17 = -0.2e1 * t11 * t23;
t15 = pkin(5) * t12;
t14 = -2 * qJD(4) * t24;
t1 = [0, 0, 0, -t26, t27, 0, 0, 0, 0, 0, t8 * t9 * t22 - t11 * t20, t11 * t21 + t9 * t18; 0, 0, 0, -t27, -t26, 0, 0, 0, 0, 0, t8 * t17 - t9 * t20, t10 * t17 + t9 * t21; 0, 0, 0, 0, 0, t8 * t18, -t28 * t22, t12 * t10, -t12 * t8, 0, -t15 * t10 + t8 * t14, t10 * t14 + t15 * t8; 0, 0, 0, 0, 0, -t8 * t13 * t10, t28 * t13, 0, 0, 0, t19 * t8, t19 * t10;];
tauc_reg = t1;
