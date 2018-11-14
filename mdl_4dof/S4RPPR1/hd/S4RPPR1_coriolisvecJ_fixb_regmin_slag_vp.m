% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RPPR1
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
% taug_reg [4x10]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:47
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc_reg = S4RPPR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:46:46
% EndTime: 2018-11-14 13:46:46
% DurationCPUTime: 0.09s
% Computational Cost: add. (64->21), mult. (110->32), div. (0->0), fcn. (45->4), ass. (0->14)
t14 = qJD(1) - qJD(4);
t16 = qJD(4) + t14;
t4 = sin(pkin(6)) * pkin(1) + qJ(3);
t2 = qJD(1) * t4;
t15 = qJD(3) * qJD(1);
t3 = -cos(pkin(6)) * pkin(1) - pkin(2) - pkin(3);
t13 = qJD(3) * (qJD(1) + t14);
t12 = t14 ^ 2;
t1 = t3 * qJD(1) + qJD(3);
t8 = sin(qJ(4));
t9 = cos(qJ(4));
t11 = t9 * t1 - t8 * t2;
t10 = -t8 * t1 - t9 * t2;
t5 = [0, 0, 0, 0, 0, 0.2e1 * t15, 0.2e1 * qJD(3) * t2, 0, t8 * t13 + (-(-t3 * t8 - t4 * t9) * t14 - t10) * qJD(4), t9 * t13 + ((t3 * t9 - t4 * t8) * t14 + t11) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -qJD(1) ^ 2, -t2 * qJD(1), 0, -t8 * t12, -t9 * t12; 0, 0, 0, 0, 0, 0, 0, 0, t16 * t10 - t8 * t15, -t16 * t11 - t9 * t15;];
tauc_reg  = t5;
