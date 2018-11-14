% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% taug_reg [4x10]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:51
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPRR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:50:35
% EndTime: 2018-11-14 13:50:35
% DurationCPUTime: 0.16s
% Computational Cost: add. (174->41), mult. (430->53), div. (0->0), fcn. (234->6), ass. (0->28)
t15 = cos(pkin(7)) * pkin(1) + pkin(2);
t14 = t15 * qJD(1);
t21 = sin(qJ(3));
t23 = cos(qJ(3));
t32 = pkin(1) * sin(pkin(7));
t30 = qJD(1) * t32;
t6 = t23 * t14 - t21 * t30;
t17 = qJD(1) + qJD(3);
t16 = qJD(4) + t17;
t33 = qJD(4) - t16;
t22 = cos(qJ(4));
t7 = t21 * t14 + t23 * t30;
t31 = t22 * t7;
t2 = t17 * pkin(3) + t6;
t29 = -pkin(3) * t16 - t2;
t20 = sin(qJ(4));
t4 = t6 * qJD(3);
t5 = t7 * qJD(3);
t28 = -t20 * t4 - t22 * t5;
t26 = -t20 * t2 - t31;
t11 = t21 * t15 + t23 * t32;
t25 = t23 * t15 - t21 * t32;
t3 = qJD(4) * t20 * t7;
t24 = t3 + (-t7 * t16 + t5) * t20;
t10 = pkin(3) + t25;
t9 = t11 * qJD(3);
t8 = t25 * qJD(3);
t1 = [0, 0, 0, 0, 0, -t9 * t17 - t5, -t8 * t17 - t4, 0 (-t20 * t8 - t22 * t9) * t16 + ((-t10 * t20 - t11 * t22) * t16 + t26) * qJD(4) + t28, t3 + (-(-qJD(4) * t11 - t9) * t16 + t5) * t20 + (-(qJD(4) * t10 + t8) * t16 - t4 - qJD(4) * t2) * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t7 * t17 - t5, t6 * t17 - t4, 0 -(-t20 * t6 - t31) * t16 + (t29 * t20 - t31) * qJD(4) + t28 (t29 * qJD(4) + t6 * t16 - t4) * t22 + t24; 0, 0, 0, 0, 0, 0, 0, 0, t33 * t26 + t28 (-t33 * t2 - t4) * t22 + t24;];
tauc_reg  = t1;
