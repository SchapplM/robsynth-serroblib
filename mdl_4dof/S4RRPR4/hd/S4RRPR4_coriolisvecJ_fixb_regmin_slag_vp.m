% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% tauc_reg [4x17]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRPR4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:33
% EndTime: 2019-12-31 17:02:34
% DurationCPUTime: 0.31s
% Computational Cost: add. (373->75), mult. (731->123), div. (0->0), fcn. (448->6), ass. (0->59)
t41 = qJD(1) + qJD(2);
t47 = cos(qJ(2));
t63 = pkin(1) * qJD(2);
t59 = qJD(1) * t63;
t23 = t41 * qJD(3) + t47 * t59;
t42 = sin(pkin(7));
t43 = cos(pkin(7));
t65 = t42 ^ 2 + t43 ^ 2;
t73 = t65 * t23;
t44 = sin(qJ(4));
t46 = cos(qJ(4));
t25 = t46 * t42 + t44 * t43;
t12 = t25 * t41;
t20 = t25 * qJD(4);
t45 = sin(qJ(2));
t72 = t45 * pkin(1);
t71 = t47 * pkin(1);
t67 = t46 * t43;
t68 = t44 * t42;
t24 = -t67 + t68;
t53 = t45 * t59;
t36 = -t43 * pkin(3) - pkin(2);
t64 = pkin(1) * qJD(1);
t50 = -t47 * t64 + qJD(3);
t9 = t36 * t41 + t50;
t70 = t9 * t20 + t24 * t53;
t19 = t24 * qJD(4);
t69 = -t9 * t19 + t25 * t53;
t62 = t41 * t68;
t61 = t41 * t67;
t60 = t45 * t63;
t58 = t65 * t47;
t34 = t47 * t63 + qJD(3);
t56 = t65 * t34;
t55 = t41 * t42 * t72;
t54 = t65 * qJD(3);
t52 = (-qJD(2) + t41) * t64;
t51 = (-qJD(1) - t41) * t63;
t49 = t45 * t51;
t48 = t45 * t52;
t37 = t43 * pkin(6);
t35 = qJ(3) + t72;
t32 = t42 * t53;
t31 = t43 * qJ(3) + t37;
t30 = (-pkin(6) - qJ(3)) * t42;
t29 = t36 - t71;
t28 = qJD(4) * t61;
t27 = t41 * qJ(3) + t45 * t64;
t26 = -t41 * pkin(2) + t50;
t22 = t43 * t35 + t37;
t21 = (-pkin(6) - t35) * t42;
t16 = t20 * qJD(4);
t15 = t19 * qJD(4);
t10 = -t61 + t62;
t8 = t41 * t20;
t7 = -qJD(4) * t62 + t28;
t2 = -t12 * t19 + t7 * t25;
t1 = t19 * t10 - t12 * t20 - t7 * t24 - t25 * t8;
t3 = [0, 0, 0, 0, t49, t47 * t51, t43 * t49, qJD(2) * t55 + t32, t41 * t56 + t73, t27 * t56 + t35 * t73 + (t26 + (-pkin(2) - t71) * qJD(1)) * t60, t2, t1, -t15, -t16, 0, t10 * t60 + t29 * t8 + ((-t21 * t44 - t22 * t46) * qJD(4) - t25 * t34) * qJD(4) + t70, t12 * t60 + t29 * t7 + ((-t21 * t46 + t22 * t44) * qJD(4) + t24 * t34) * qJD(4) + t69; 0, 0, 0, 0, t48, t47 * t52, t43 * t48, -qJD(1) * t55 + t32, (-t58 * t64 + t54) * t41 + t73, t27 * t54 + qJ(3) * t73 + ((-pkin(2) * qJD(2) - t26) * t45 - t27 * t58) * t64, t2, t1, -t15, -t16, 0, t36 * t8 + ((-t30 * t44 - t31 * t46) * qJD(4) - t25 * qJD(3)) * qJD(4) + (-t45 * t10 + t47 * t20) * t64 + t70, t36 * t7 + ((-t30 * t46 + t31 * t44) * qJD(4) + t24 * qJD(3)) * qJD(4) + (-t45 * t12 - t47 * t19) * t64 + t69; 0, 0, 0, 0, 0, 0, 0, 0, -t65 * t41 ^ 2, -t65 * t41 * t27 + t53, 0, 0, 0, 0, 0, 0.2e1 * t12 * qJD(4), t28 + (-t10 - t62) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12 * t10, -t10 ^ 2 + t12 ^ 2, t28 + (t10 - t62) * qJD(4), 0, 0, -t9 * t12 - t25 * t23, t9 * t10 + t24 * t23;];
tauc_reg = t3;
