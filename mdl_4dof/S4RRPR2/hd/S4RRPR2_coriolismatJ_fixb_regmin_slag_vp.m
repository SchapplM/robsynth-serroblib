% Calculate minimal parameter regressor of coriolis matrix for
% S4RRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x12]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRPR2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 18:16:35
% EndTime: 2019-07-18 18:16:36
% DurationCPUTime: 0.23s
% Computational Cost: add. (251->81), mult. (425->89), div. (0->0), fcn. (277->4), ass. (0->64)
t79 = -pkin(2) - pkin(3);
t46 = sin(qJ(2));
t78 = t46 * pkin(1);
t47 = cos(qJ(4));
t71 = t47 * t46;
t45 = sin(qJ(4));
t48 = cos(qJ(2));
t73 = t45 * t48;
t12 = (-t71 + t73) * pkin(1);
t39 = t45 * qJD(3);
t77 = t12 * qJD(2) + t39;
t70 = t47 * t48;
t74 = t45 * t46;
t13 = (t70 + t74) * pkin(1);
t40 = t47 * qJD(3);
t76 = t13 * qJD(2) + t40;
t32 = qJ(3) + t78;
t75 = t45 * t32;
t72 = t47 * t32;
t69 = pkin(1) * qJD(1);
t68 = pkin(1) * qJD(2);
t51 = -qJ(3) / 0.2e1 - t32 / 0.2e1;
t55 = t71 / 0.2e1;
t56 = t45 * t79;
t1 = pkin(1) * t55 + t51 * t47 - t56;
t67 = t1 * qJD(1);
t52 = -pkin(1) * t48 - pkin(2);
t50 = -pkin(3) + t52;
t23 = t47 * t50;
t35 = t47 * t79;
t53 = -t23 / 0.2e1 - t35 / 0.2e1;
t54 = -t70 / 0.2e1;
t2 = pkin(1) * t54 + (-t78 / 0.2e1 - t51) * t45 + t53;
t66 = t2 * qJD(1);
t65 = t45 * qJ(3);
t64 = t47 * qJ(3);
t5 = (t32 * t48 + t52 * t46) * pkin(1);
t63 = t5 * qJD(1);
t62 = qJD(1) * t12;
t61 = qJD(1) * t13;
t38 = t48 * t68;
t60 = t38 + qJD(3);
t43 = qJD(1) + qJD(2);
t59 = qJD(1) - qJD(4);
t58 = qJD(2) - qJD(4);
t57 = t46 * t68;
t21 = t43 * t45;
t22 = t43 * t47;
t49 = t45 * t50;
t44 = qJ(3) * qJD(3);
t37 = t48 * t69;
t36 = t46 * t69;
t29 = t43 * qJ(3);
t28 = t32 * qJD(3);
t18 = t56 + t64;
t17 = -t35 + t65;
t14 = -t36 - t57;
t11 = -qJD(4) * t47 + t22;
t10 = -qJD(4) * t45 + t21;
t7 = t49 + t72;
t6 = -t23 + t75;
t4 = t64 / 0.2e1 + t56 / 0.2e1 + t72 / 0.2e1 + t49 / 0.2e1 + (-t73 / 0.2e1 + t55) * pkin(1);
t3 = -t65 / 0.2e1 - t75 / 0.2e1 + (t54 - t74 / 0.2e1) * pkin(1) - t53;
t8 = [0, 0, 0, 0, -t57, -t38, -t57, t60, qJD(2) * t5 + t28, 0, qJD(4) * t7 + t77, -qJD(4) * t6 + t76; 0, 0, 0, 0, t14, -t37 - t38, t14, t37 + t60, t63 + t28 + (-pkin(2) * t46 + qJ(3) * t48) * t68, 0, qJD(4) * t4 + t62 + t77, qJD(4) * t3 + t61 + t76; 0, 0, 0, 0, 0, 0, 0, t43, t43 * t32, 0, t21, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4 * qJD(2) + t59 * t7, t3 * qJD(2) - t59 * t6; 0, 0, 0, 0, t36, t37, t36, -t37 + qJD(3), t44 - t63, 0, -qJD(4) * t1 + t39 - t62, -qJD(4) * t2 + t40 - t61; 0, 0, 0, 0, 0, 0, 0, qJD(3), t44, 0, qJD(4) * t18 + t39, -qJD(4) * t17 + t40; 0, 0, 0, 0, 0, 0, 0, t43, t29, 0, t21, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58 * t18 - t67, -t58 * t17 - t66; 0, 0, 0, 0, 0, 0, 0, -t43, -qJ(3) * qJD(2) - qJD(1) * t32, 0, -t10, -t11; 0, 0, 0, 0, 0, 0, 0, -t43, -t29, 0, -t10, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(1) * t7 + qJD(2) * t1 - t39, qJD(1) * t6 + qJD(2) * t2 - t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2) * t18 - t39 + t67, qJD(2) * t17 - t40 + t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t8;
