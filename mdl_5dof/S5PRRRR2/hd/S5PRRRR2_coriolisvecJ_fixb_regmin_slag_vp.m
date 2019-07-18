% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
% 
% Output:
% tauc_reg [5x17]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:30
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:30:18
% EndTime: 2019-07-18 13:30:19
% DurationCPUTime: 0.24s
% Computational Cost: add. (314->67), mult. (701->98), div. (0->0), fcn. (364->6), ass. (0->60)
t26 = qJD(2) + qJD(3);
t24 = qJD(4) + t26;
t34 = cos(qJ(3));
t62 = pkin(2) * qJD(2);
t55 = t34 * t62;
t18 = t26 * pkin(3) + t55;
t33 = cos(qJ(4));
t49 = qJD(3) * t55;
t30 = sin(qJ(4));
t31 = sin(qJ(3));
t56 = t31 * t62;
t50 = t30 * t56;
t64 = (qJD(3) + qJD(4)) * t50;
t39 = -(qJD(4) * t18 + t49) * t33 + t64;
t9 = -t33 * t18 + t50;
t51 = -t9 * t24 + t39;
t66 = t31 * t33;
t44 = t30 * t34 + t66;
t60 = qJD(4) * t33;
t75 = t44 * qJD(3) + t31 * t60;
t29 = sin(qJ(5));
t38 = t75 * pkin(2);
t61 = qJD(4) * t30;
t54 = t18 * t61;
t3 = qJD(2) * t38 + t54;
t1 = t3 * t29;
t32 = cos(qJ(5));
t58 = t32 * qJD(5);
t74 = t9 * t58 + t1;
t22 = t34 * pkin(2) + pkin(3);
t73 = (t22 * t61 + t38) * t24;
t71 = (t30 * t18 + t33 * t56) * t24;
t13 = t44 * t62;
t70 = t13 * t24;
t35 = qJD(5) ^ 2;
t69 = (t30 * pkin(3) + pkin(6)) * t35;
t68 = t29 * t32;
t67 = t30 * t31;
t65 = t35 * t29;
t63 = t29 ^ 2 - t32 ^ 2;
t59 = t29 * qJD(5);
t57 = 0.2e1 * qJD(5) * t24;
t7 = t9 * t59;
t52 = t33 * t58;
t48 = (-qJD(3) + t26) * t62;
t47 = pkin(2) * qJD(3) * (-qJD(2) - t26);
t46 = pkin(6) * t35 - t71;
t45 = (pkin(2) * t66 + t30 * t22 + pkin(6)) * t35 + t73;
t43 = t33 * t34 - t67;
t4 = t22 * t60 + (t43 * qJD(3) - t31 * t61) * pkin(2);
t42 = qJD(5) * ((pkin(2) * t67 - t33 * t22) * t24 - t4);
t41 = (-pkin(3) * t24 - t18) * qJD(4);
t37 = t75 * t62;
t36 = -t37 - t54;
t25 = t35 * t32;
t23 = t24 ^ 2;
t17 = t57 * t68;
t14 = t43 * t62;
t11 = t63 * t57;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, -t25; 0, 0, 0, 0, 0, t31 * t47, t34 * t47, 0, t36 - t73, -t4 * t24 + t39, t17, -t11, t25, -t65, 0, t7 + t29 * t42 + (-t3 - t45) * t32, t45 * t29 + t32 * t42 + t74; 0, 0, 0, 0, 0, t31 * t48, t34 * t48, 0, t30 * t41 - t37 + t70, t14 * t24 + (t41 - t49) * t33 + t64, t17, -t11, t25, -t65, 0, t7 + (t14 + (-qJD(4) - t24) * t33 * pkin(3)) * t59 + (-t69 - t3 + (-pkin(3) * t61 + t13) * t24) * t32, t14 * t58 + (t69 - t70) * t29 + (-t24 * t52 + (t24 * t29 * t30 - t52) * qJD(4)) * pkin(3) + t74; 0, 0, 0, 0, 0, 0, 0, 0, t36 + t71, t51, t17, -t11, t25, -t65, 0, (-t3 - t46) * t32, t46 * t29 + t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23 * t68, t63 * t23, 0, 0, 0, t51 * t29, t51 * t32;];
tauc_reg  = t2;
