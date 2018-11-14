% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4PPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta2]';
% 
% Output:
% taug_reg [4x8]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:00
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc_reg = S4PPRR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:59:33
% EndTime: 2018-11-14 13:59:34
% DurationCPUTime: 0.15s
% Computational Cost: add. (84->28), mult. (230->42), div. (0->0), fcn. (186->6), ass. (0->22)
t14 = sin(pkin(6));
t15 = cos(pkin(6));
t17 = sin(qJ(3));
t19 = cos(qJ(3));
t29 = -t17 * t14 + t19 * t15;
t6 = t29 * qJD(1);
t13 = qJD(3) + qJD(4);
t28 = qJD(4) - t13;
t11 = t19 * t14 + t17 * t15;
t7 = t11 * qJD(1);
t18 = cos(qJ(4));
t27 = t18 * t7;
t3 = qJD(3) * pkin(3) + t6;
t24 = -pkin(3) * t13 - t3;
t16 = sin(qJ(4));
t8 = t29 * qJD(3);
t4 = qJD(1) * t8;
t9 = t11 * qJD(3);
t5 = qJD(1) * t9;
t23 = -t16 * t4 - t18 * t5;
t20 = (t28 * t7 + t5) * t16;
t1 = [0, 0, 0, -t9 * qJD(3), -t8 * qJD(3), 0 (-t16 * t8 - t18 * t9 + (-t11 * t18 - t16 * t29) * qJD(4)) * t13 -(-t16 * t9 + t18 * t8 + (-t11 * t16 + t18 * t29) * qJD(4)) * t13; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0 -(-t16 * t6 - t27) * t13 + (t24 * t16 - t27) * qJD(4) + t23 (t24 * qJD(4) + t6 * t13 - t4) * t18 + t20; 0, 0, 0, 0, 0, 0, t23 + t28 * (-t16 * t3 - t27) (-t28 * t3 - t4) * t18 + t20;];
tauc_reg  = t1;
