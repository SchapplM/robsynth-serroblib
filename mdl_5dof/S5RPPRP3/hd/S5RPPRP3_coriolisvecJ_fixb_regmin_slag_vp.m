% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% tauc_reg [5x16]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPRP3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:51:02
% EndTime: 2019-12-31 17:51:03
% DurationCPUTime: 0.24s
% Computational Cost: add. (270->53), mult. (516->90), div. (0->0), fcn. (231->4), ass. (0->41)
t14 = -cos(pkin(7)) * pkin(1) - pkin(2) - pkin(6);
t41 = qJ(5) - t14;
t23 = sin(qJ(4));
t24 = cos(qJ(4));
t52 = -t41 * qJD(1) + qJD(3);
t5 = t24 * qJD(2) + t52 * t23;
t4 = -t23 * qJD(2) + t52 * t24;
t33 = qJD(1) * qJD(5);
t1 = t4 * qJD(4) - t23 * t33;
t2 = -t5 * qJD(4) - t24 * t33;
t43 = qJD(4) * pkin(4);
t3 = t4 + t43;
t51 = -t1 * t23 - t2 * t24 + (t23 * t3 - t24 * t5) * qJD(4);
t50 = 0.2e1 * qJD(3);
t49 = t3 - t4;
t25 = qJD(4) ^ 2;
t48 = t25 * t23;
t47 = t25 * t24;
t18 = qJD(3) * qJD(1);
t34 = qJD(1) * qJD(4);
t32 = t24 * t34;
t46 = pkin(4) * t32 + t18;
t19 = t23 ^ 2;
t20 = t24 ^ 2;
t45 = t19 - t20;
t26 = qJD(1) ^ 2;
t44 = -t25 - t26;
t16 = sin(pkin(7)) * pkin(1) + qJ(3);
t27 = t23 * pkin(4) + t16;
t9 = t27 * qJD(1) + qJD(5);
t42 = t9 * qJD(1);
t13 = qJD(1) * t16;
t40 = t13 * qJD(1);
t38 = t23 * qJD(4);
t36 = t24 * qJD(4);
t11 = t41 * t24;
t29 = t23 * t5 + t24 * t3;
t10 = t41 * t23;
t7 = -qJD(4) * t11 - t23 * qJD(5);
t6 = -t24 * qJD(5) + t41 * t38;
t8 = [0, 0, 0, 0, 0, 0.2e1 * t18, t13 * t50, -0.2e1 * t23 * t32, 0.2e1 * t45 * t34, -t48, -t47, 0, t13 * t36 - t14 * t48 + (t16 * t36 + t23 * t50) * qJD(1), -t13 * t38 - t14 * t47 + (-t16 * t38 + t24 * t50) * qJD(1), (-t23 * t7 - t24 * t6 + (t10 * t24 - t11 * t23) * qJD(4)) * qJD(1) + t51, -t1 * t10 + t5 * t7 - t2 * t11 + t3 * t6 + t46 * t27 + t9 * (pkin(4) * t36 + qJD(3)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, t48, 0, -t29 * qJD(4) + t1 * t24 - t2 * t23; 0, 0, 0, 0, 0, -t26, -t40, 0, 0, 0, 0, 0, t44 * t23, t44 * t24, 0, -t42 - t51; 0, 0, 0, 0, 0, 0, 0, t24 * t26 * t23, -t45 * t26, 0, 0, 0, -t24 * t40, t23 * t40, (t43 - t49) * t23 * qJD(1), t49 * t5 + (-t24 * t42 + t2) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t19 - t20) * t26, t29 * qJD(1) + t46;];
tauc_reg = t8;
