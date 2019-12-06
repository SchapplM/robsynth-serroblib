% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% tauc_reg [5x17]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:06:18
% EndTime: 2019-12-05 17:06:19
% DurationCPUTime: 0.25s
% Computational Cost: add. (348->71), mult. (715->98), div. (0->0), fcn. (366->6), ass. (0->60)
t31 = sin(qJ(4));
t35 = cos(qJ(3));
t32 = sin(qJ(3));
t34 = cos(qJ(4));
t68 = t32 * t34;
t47 = t31 * t35 + t68;
t64 = pkin(2) * qJD(2);
t15 = t47 * t64;
t27 = qJD(2) + qJD(3);
t25 = qJD(4) + t27;
t36 = qJD(5) ^ 2;
t62 = qJD(4) * t31;
t77 = (pkin(3) * t62 - t15) * t25 + (pkin(3) * t31 + pkin(8)) * t36;
t58 = t35 * t64;
t18 = pkin(3) * t27 + t58;
t53 = qJD(3) * t58;
t59 = t32 * t64;
t54 = t31 * t59;
t66 = (qJD(3) + qJD(4)) * t54;
t40 = -(qJD(4) * t18 + t53) * t34 + t66;
t61 = qJD(4) * t34;
t76 = t47 * qJD(3) + t32 * t61;
t39 = t76 * pkin(2);
t57 = t18 * t62;
t3 = qJD(2) * t39 + t57;
t30 = sin(qJ(5));
t33 = cos(qJ(5));
t10 = t18 * t34 - t54;
t74 = pkin(4) * t25;
t8 = -t10 - t74;
t63 = qJD(5) * t8;
t75 = t3 * t30 + t33 * t63;
t23 = pkin(2) * t35 + pkin(3);
t73 = t25 * (t23 * t62 + t39);
t72 = (t18 * t31 + t34 * t59) * t25;
t70 = t30 * t33;
t69 = t31 * t32;
t67 = t36 * t30;
t65 = t30 ^ 2 - t33 ^ 2;
t60 = 0.2e1 * qJD(5) * t25;
t55 = -t25 * t8 + t40;
t51 = (-qJD(3) + t27) * t64;
t50 = pkin(2) * qJD(3) * (-qJD(2) - t27);
t49 = pkin(8) * t36 - t72;
t48 = (pkin(2) * t68 + t23 * t31 + pkin(8)) * t36 + t73;
t46 = t34 * t35 - t69;
t45 = qJD(5) * (t10 - t74);
t4 = t23 * t61 + (t46 * qJD(3) - t32 * t62) * pkin(2);
t44 = qJD(5) * ((pkin(2) * t69 - t23 * t34 - pkin(4)) * t25 - t4);
t43 = (-pkin(3) * t25 - t18) * qJD(4);
t16 = t46 * t64;
t41 = qJD(5) * (-pkin(3) * t61 + (-pkin(3) * t34 - pkin(4)) * t25 + t16);
t38 = t76 * t64;
t37 = -t38 - t57;
t26 = t36 * t33;
t24 = t25 ^ 2;
t17 = t60 * t70;
t12 = t65 * t60;
t6 = t30 * t63;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, -t26; 0, 0, 0, 0, 0, t32 * t50, t35 * t50, 0, t37 - t73, -t25 * t4 + t40, t17, -t12, t26, -t67, 0, t6 + t30 * t44 + (-t3 - t48) * t33, t48 * t30 + t33 * t44 + t75; 0, 0, 0, 0, 0, t32 * t51, t35 * t51, 0, t15 * t25 + t31 * t43 - t38, t16 * t25 + (t43 - t53) * t34 + t66, t17, -t12, t26, -t67, 0, t6 + t30 * t41 + (-t3 - t77) * t33, t77 * t30 + t33 * t41 + t75; 0, 0, 0, 0, 0, 0, 0, 0, t37 + t72, t10 * t25 + t40, t17, -t12, t26, -t67, 0, t6 + t30 * t45 + (-t3 - t49) * t33, t49 * t30 + t33 * t45 + t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24 * t70, t65 * t24, 0, 0, 0, t55 * t30, t55 * t33;];
tauc_reg = t1;
