% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4PRPR1
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
% taug_reg [4x10]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:43
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRPR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:42:23
% EndTime: 2018-11-14 13:42:24
% DurationCPUTime: 0.13s
% Computational Cost: add. (44->21), mult. (76->33), div. (0->0), fcn. (26->2), ass. (0->15)
t5 = -pkin(2) - pkin(3);
t1 = t5 * qJD(2) + qJD(3);
t3 = sin(qJ(4));
t15 = t3 * t1;
t14 = qJ(3) * t3;
t11 = qJD(2) - qJD(4);
t13 = qJD(4) + t11;
t12 = qJD(2) * qJ(3);
t10 = 0.2e1 * qJD(2) * qJD(3);
t4 = cos(qJ(4));
t9 = t13 * t4;
t8 = qJD(3) * (qJD(2) + t11);
t7 = t11 ^ 2;
t6 = qJD(2) ^ 2;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t10, qJ(3) * t10, 0, t3 * t8 + (-(-qJ(3) * t4 - t3 * t5) * t11 + t4 * t12 + t15) * qJD(4), t4 * t8 + ((t4 * t5 - t14) * t11 - t3 * t12 + t4 * t1) * qJD(4); 0, 0, 0, 0, 0, -t6, -t6 * qJ(3), 0, -t3 * t7, -t4 * t7; 0, 0, 0, 0, 0, 0, 0, 0, -t13 * t15 + (-qJ(3) * t9 - qJD(3) * t3) * qJD(2), -t1 * t9 + (-qJD(3) * t4 + t13 * t14) * qJD(2);];
tauc_reg  = t2;
