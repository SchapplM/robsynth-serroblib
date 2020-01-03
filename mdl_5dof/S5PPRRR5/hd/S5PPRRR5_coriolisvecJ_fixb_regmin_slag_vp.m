% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
% 
% Output:
% tauc_reg [5x15]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PPRRR5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:35:45
% EndTime: 2019-12-31 17:35:46
% DurationCPUTime: 0.23s
% Computational Cost: add. (212->52), mult. (473->83), div. (0->0), fcn. (310->6), ass. (0->52)
t70 = 2 * qJD(5);
t33 = cos(qJ(3));
t52 = t33 * qJD(2);
t20 = qJD(3) * pkin(3) + t52;
t30 = sin(qJ(3));
t29 = sin(qJ(4));
t32 = cos(qJ(4));
t16 = t29 * t33 + t32 * t30;
t39 = t16 * qJD(3);
t53 = qJD(4) * t32;
t36 = (t30 * t53 + t39) * qJD(2);
t54 = qJD(4) * t29;
t3 = t20 * t54 + t36;
t25 = qJD(3) + qJD(4);
t15 = t29 * t30 - t32 * t33;
t69 = t15 * t25;
t13 = t16 * qJD(2);
t34 = qJD(5) ^ 2;
t67 = (pkin(3) * t54 - t13) * t25 + (t29 * pkin(3) + pkin(7)) * t34;
t46 = qJD(3) * t52;
t55 = qJD(2) * t30;
t48 = t29 * t55;
t58 = t25 * t48;
t66 = -(qJD(4) * t20 + t46) * t32 + t58;
t28 = sin(qJ(5));
t31 = cos(qJ(5));
t10 = t32 * t20 - t48;
t64 = t25 * pkin(4);
t8 = -t10 - t64;
t56 = qJD(5) * t8;
t65 = t3 * t28 + t31 * t56;
t63 = (t16 * qJD(4) + t39) * t25;
t62 = (t29 * t20 + t32 * t55) * t25;
t60 = t28 * t31;
t59 = t34 * t28;
t57 = t28 ^ 2 - t31 ^ 2;
t51 = t25 * t70;
t47 = -t8 * t25 + t66;
t44 = pkin(7) * t34 - t62;
t43 = t16 * t34 + t63;
t42 = qJD(5) * (t10 - t64);
t41 = t69 * t70;
t40 = (-pkin(3) * t25 - t20) * qJD(4);
t14 = t15 * qJD(2);
t37 = qJD(5) * (-pkin(3) * t53 + (-t32 * pkin(3) - pkin(4)) * t25 - t14);
t35 = qJD(3) ^ 2;
t24 = t25 ^ 2;
t23 = t34 * t31;
t17 = t51 * t60;
t12 = t57 * t51;
t6 = t28 * t56;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, t23; 0, 0, 0, -t35 * t30, -t35 * t33, 0, -t63, t69 * t25, 0, 0, 0, 0, 0, t28 * t41 - t43 * t31, t43 * t28 + t31 * t41; 0, 0, 0, 0, 0, 0, t13 * t25 + t29 * t40 - t36, -t14 * t25 + (t40 - t46) * t32 + t58, t17, -t12, t23, -t59, 0, t6 + t28 * t37 + (-t3 - t67) * t31, t67 * t28 + t31 * t37 + t65; 0, 0, 0, 0, 0, 0, t62 - t3, t10 * t25 + t66, t17, -t12, t23, -t59, 0, t6 + t28 * t42 + (-t3 - t44) * t31, t44 * t28 + t31 * t42 + t65; 0, 0, 0, 0, 0, 0, 0, 0, -t24 * t60, t57 * t24, 0, 0, 0, t47 * t28, t47 * t31;];
tauc_reg = t1;
