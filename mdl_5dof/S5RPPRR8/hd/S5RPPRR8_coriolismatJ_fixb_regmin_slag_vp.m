% Calculate minimal parameter regressor of coriolis matrix for
% S5RPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x19]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPPRR8_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR8_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:01:15
% EndTime: 2019-12-31 18:01:17
% DurationCPUTime: 0.43s
% Computational Cost: add. (545->77), mult. (860->99), div. (0->0), fcn. (820->6), ass. (0->60)
t46 = sin(pkin(8));
t47 = cos(pkin(8));
t51 = -pkin(1) - pkin(2);
t36 = t47 * qJ(2) + t46 * t51;
t49 = sin(qJ(4));
t61 = -t46 * qJ(2) + t47 * t51;
t59 = -pkin(3) + t61;
t81 = cos(qJ(4));
t17 = t49 * t36 - t81 * t59;
t15 = pkin(4) + t17;
t18 = t81 * t36 + t49 * t59;
t64 = qJD(1) - qJD(4);
t60 = t64 * t18;
t48 = sin(qJ(5));
t50 = cos(qJ(5));
t39 = -t48 ^ 2 + t50 ^ 2;
t87 = t64 * t39;
t35 = t81 * t46 + t49 * t47;
t69 = t35 * qJD(4);
t71 = t35 * qJD(1);
t86 = t71 - t69;
t34 = t49 * t46 - t81 * t47;
t74 = t34 * qJD(1);
t85 = -t34 * qJD(4) + t74;
t84 = t17 / 0.2e1;
t82 = pkin(4) * t48;
t80 = t15 * t48;
t79 = qJD(1) * t50;
t78 = qJD(4) * t50;
t77 = qJD(5) * t48;
t76 = qJD(5) * t50;
t19 = t36 * t47 - t61 * t46;
t75 = t19 * qJD(1);
t73 = t34 * qJD(2);
t70 = t35 * qJD(2);
t68 = t39 * qJD(5);
t67 = t46 * qJD(1);
t66 = t47 * qJD(1);
t65 = qJ(2) * qJD(1);
t63 = qJD(1) * t80;
t62 = t15 * t79;
t58 = t84 - pkin(4) / 0.2e1 - t15 / 0.2e1;
t57 = t18 * qJD(1) + t70;
t56 = t18 * qJD(4) + t70;
t1 = t58 * t48;
t55 = t1 * qJD(1) + qJD(4) * t82;
t2 = t58 * t50;
t54 = pkin(4) * t78 + t2 * qJD(1);
t12 = t34 * t48;
t27 = t50 * t71;
t53 = -t12 * qJD(5) + t50 * t69 - t27;
t14 = t34 * t50;
t26 = t48 * t71;
t52 = t14 * qJD(5) + t48 * t69 - t26;
t41 = t48 * t76;
t29 = (-t78 + t79) * t48;
t16 = -pkin(7) + t18;
t4 = t15 * t50;
t3 = t82 / 0.2e1 + t80 / 0.2e1 + t48 * t84;
t5 = [0, 0, 0, 0, qJD(2), qJ(2) * qJD(2), t46 * qJD(2), t47 * qJD(2), t19 * qJD(2), 0, t56, -t17 * qJD(4) - t73, t41, t68, 0, 0, 0, -t15 * t77 + t56 * t50, -t15 * t76 - t56 * t48; 0, 0, 0, 0, qJD(1), t65, t67, t66, t75, 0, t71, -t74, 0, 0, 0, 0, 0, t27, -t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, -t64 * t17, -t41, -t68, 0, 0, 0, t3 * qJD(5) + t50 * t60, t4 * qJD(5) - t48 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t87, -t76, t77, 0, t3 * qJD(4) - t16 * t76 - t63, t4 * qJD(4) + t16 * t77 - t62; 0, 0, 0, 0, -qJD(1), -t65, -t67, -t66, -t75, 0, -t86, t85, 0, 0, 0, 0, 0, t53, -t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, -t85, 0, 0, 0, 0, 0, -t53, t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64 * t12 - t35 * t76, -t64 * t14 + t35 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, -t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, t17 * qJD(1) + t73, -t41, -t68, 0, 0, 0, -t1 * qJD(5) - t57 * t50, -t2 * qJD(5) + t57 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, t74, 0, 0, 0, 0, 0, -t27, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t68, 0, 0, 0, -pkin(4) * t77, -pkin(4) * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t87, t76, -t77, 0, -pkin(7) * t76 - t55, pkin(7) * t77 - t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t87, 0, 0, 0, t12 * qJD(2) + t1 * qJD(4) + t63, t14 * qJD(2) + t2 * qJD(4) + t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12 * qJD(1), t14 * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t87, 0, 0, 0, t55, t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t5;
