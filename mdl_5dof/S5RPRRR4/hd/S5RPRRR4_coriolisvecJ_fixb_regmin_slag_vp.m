% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% tauc_reg [5x17]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:35
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRR4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:34:46
% EndTime: 2022-01-23 09:34:47
% DurationCPUTime: 0.32s
% Computational Cost: add. (636->79), mult. (1413->101), div. (0->0), fcn. (798->8), ass. (0->62)
t32 = cos(pkin(9)) * pkin(1) + pkin(2);
t31 = t32 * qJD(1);
t45 = sin(qJ(3));
t48 = cos(qJ(3));
t81 = pkin(1) * sin(pkin(9));
t67 = qJD(1) * t81;
t22 = t48 * t31 - t45 * t67;
t44 = sin(qJ(4));
t23 = t45 * t31 + t48 * t67;
t47 = cos(qJ(4));
t73 = t47 * t23;
t12 = t44 * t22 + t73;
t38 = qJD(1) + qJD(3);
t36 = qJD(4) + t38;
t49 = qJD(5) ^ 2;
t69 = qJD(4) * t44;
t83 = (pkin(3) * t69 - t12) * t36 + (t44 * pkin(3) + pkin(8)) * t49;
t17 = t38 * pkin(3) + t22;
t19 = t22 * qJD(3);
t20 = t23 * qJD(3);
t63 = -t44 * t20 - t23 * t69;
t51 = -(qJD(4) * t17 + t19) * t47 - t63;
t11 = t44 * t17 + t73;
t64 = t44 * t19 + t47 * t20;
t3 = t11 * qJD(4) + t64;
t43 = sin(qJ(5));
t46 = cos(qJ(5));
t74 = t44 * t23;
t10 = t47 * t17 - t74;
t80 = t36 * pkin(4);
t8 = -t10 - t80;
t70 = qJD(5) * t8;
t82 = t3 * t43 + t46 * t70;
t79 = t47 * pkin(3);
t53 = t48 * t32 - t45 * t81;
t24 = t53 * qJD(3);
t27 = t45 * t32 + t48 * t81;
t25 = t27 * qJD(3);
t26 = pkin(3) + t53;
t56 = t44 * t26 + t47 * t27;
t78 = (t56 * qJD(4) + t44 * t24 + t47 * t25) * t36;
t77 = t11 * t36;
t75 = t43 * t46;
t72 = t49 * t43;
t71 = t43 ^ 2 - t46 ^ 2;
t68 = 0.2e1 * qJD(5) * t36;
t66 = -t8 * t36 + t51;
t65 = -pkin(3) * t36 - t17;
t59 = pkin(8) * t49 - t77;
t58 = (pkin(8) + t56) * t49 + t78;
t57 = t47 * t26 - t44 * t27;
t55 = qJD(5) * (t10 - t80);
t4 = t57 * qJD(4) + t47 * t24 - t44 * t25;
t54 = qJD(5) * ((-pkin(4) - t57) * t36 - t4);
t13 = t47 * t22 - t74;
t52 = qJD(5) * (-qJD(4) * t79 + (-pkin(4) - t79) * t36 + t13);
t37 = t49 * t46;
t35 = t36 ^ 2;
t28 = t68 * t75;
t21 = t71 * t68;
t6 = t43 * t70;
t1 = [0, 0, 0, 0, 0, -t25 * t38 - t20, -t24 * t38 - t19, 0, -t3 - t78, -t4 * t36 + t51, t28, -t21, t37, -t72, 0, t6 + t43 * t54 + (-t3 - t58) * t46, t58 * t43 + t46 * t54 + t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, -t37; 0, 0, 0, 0, 0, t23 * t38 - t20, t22 * t38 - t19, 0, t12 * t36 + (t65 * t44 - t73) * qJD(4) - t64, t13 * t36 + (t65 * qJD(4) - t19) * t47 - t63, t28, -t21, t37, -t72, 0, t6 + t43 * t52 + (-t3 - t83) * t46, t83 * t43 + t46 * t52 + t82; 0, 0, 0, 0, 0, 0, 0, 0, -t3 + t77, t10 * t36 + t51, t28, -t21, t37, -t72, 0, t6 + t43 * t55 + (-t3 - t59) * t46, t59 * t43 + t46 * t55 + t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35 * t75, t71 * t35, 0, 0, 0, t66 * t43, t66 * t46;];
tauc_reg = t1;
