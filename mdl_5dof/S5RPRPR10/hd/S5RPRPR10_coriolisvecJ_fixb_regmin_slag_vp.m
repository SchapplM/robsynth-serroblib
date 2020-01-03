% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% tauc_reg [5x17]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPR10_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR10_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:26:04
% EndTime: 2019-12-31 18:26:05
% DurationCPUTime: 0.32s
% Computational Cost: add. (564->71), mult. (965->115), div. (0->0), fcn. (494->6), ass. (0->61)
t67 = qJD(1) - qJD(3);
t85 = t67 ^ 2;
t43 = sin(pkin(8));
t44 = cos(pkin(8));
t46 = sin(qJ(3));
t48 = cos(qJ(3));
t29 = t43 * t46 - t44 * t48;
t75 = t29 * t67;
t30 = t43 * t48 + t44 * t46;
t50 = qJD(5) ^ 2;
t76 = t67 * t30;
t84 = t30 * t50 + t67 * t76;
t70 = t46 * qJD(2);
t72 = qJD(3) * t48;
t49 = -pkin(1) - pkin(2);
t35 = t49 * qJD(1) + qJD(2);
t79 = t46 * t35;
t52 = -(qJ(2) * t72 + t70) * qJD(1) - qJD(3) * t79;
t69 = qJ(2) * qJD(1);
t21 = t48 * t69 + t79;
t82 = t43 * t21;
t81 = t44 * t21;
t45 = sin(qJ(5));
t47 = cos(qJ(5));
t80 = t45 * t47;
t78 = t50 * t45;
t77 = t50 * t47;
t59 = -t46 * qJ(2) + t48 * t49;
t32 = -pkin(3) + t59;
t34 = t48 * qJ(2) + t46 * t49;
t74 = t43 * t32 + t44 * t34;
t73 = t45 ^ 2 - t47 ^ 2;
t71 = qJD(5) * t67;
t68 = qJD(1) * qJD(2);
t65 = 0.2e1 * t68;
t63 = t46 * t69;
t14 = -qJD(3) * t63 + t35 * t72 + t48 * t68;
t2 = t44 * t14 + t43 * t52;
t20 = t48 * t35 - t63;
t16 = -pkin(3) * t67 + t20;
t5 = t44 * t16 - t82;
t3 = pkin(4) * t67 - t5;
t64 = t3 * t67 - t2;
t61 = t71 * t80;
t60 = t44 * t32 - t43 * t34;
t1 = t43 * t14 - t44 * t52;
t18 = t48 * qJD(2) + t59 * qJD(3);
t19 = -t34 * qJD(3) - t70;
t7 = t43 * t18 - t44 * t19;
t58 = (-pkin(7) + t74) * t50 - t67 * t7 - t1;
t9 = t43 * t20 + t81;
t57 = (t43 * pkin(3) + pkin(7)) * t50 + t67 * t9 + t1;
t8 = t44 * t18 + t43 * t19;
t56 = qJD(5) * (-(pkin(4) - t60) * t67 - t3 - t8);
t10 = t44 * t20 - t82;
t55 = qJD(5) * (-(-t44 * pkin(3) - pkin(4)) * t67 + t10 + t3);
t54 = -0.2e1 * qJD(5) * t75;
t51 = qJD(1) ^ 2;
t22 = t73 * t71;
t6 = t43 * t16 + t81;
t4 = [0, 0, 0, 0, t65, qJ(2) * t65, 0, -t19 * t67 - t52, t18 * t67 + t14, -t1 * t60 + t2 * t74 - t5 * t7 + t6 * t8, 0.2e1 * t61, -0.2e1 * t22, -t77, t78, 0, t45 * t56 - t58 * t47, t58 * t45 + t47 * t56; 0, 0, 0, 0, -t51, -t51 * qJ(2), 0, -t46 * t85, -t48 * t85, t1 * t29 + t2 * t30 + t76 * t5 + t75 * t6, 0, 0, 0, 0, 0, t45 * t54 - t84 * t47, t84 * t45 + t47 * t54; 0, 0, 0, 0, 0, 0, 0, -t21 * t67 + t52, -t20 * t67 - t14, -t6 * t10 + t5 * t9 + (-t1 * t44 + t2 * t43) * pkin(3), -0.2e1 * t61, 0.2e1 * t22, t77, -t78, 0, t45 * t55 - t57 * t47, t57 * t45 + t47 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, -t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85 * t80, t73 * t85, 0, 0, 0, t64 * t45, t64 * t47;];
tauc_reg = t4;
