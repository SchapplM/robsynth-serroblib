% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% 
% Output:
% tauc_reg [4x16]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRPP3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:57:51
% EndTime: 2019-12-31 16:57:52
% DurationCPUTime: 0.36s
% Computational Cost: add. (544->110), mult. (1482->160), div. (0->0), fcn. (915->4), ass. (0->73)
t76 = (qJD(1) * qJD(2));
t89 = -2 * t76;
t57 = sin(pkin(6));
t58 = cos(pkin(6));
t59 = sin(qJ(2));
t60 = cos(qJ(2));
t42 = t57 * t60 + t58 * t59;
t35 = t42 * qJD(1);
t31 = t35 ^ 2;
t85 = t58 * t60;
t74 = qJD(1) * t85;
t78 = t59 * qJD(1);
t32 = t57 * t78 - t74;
t88 = -t32 ^ 2 - t31;
t81 = -qJ(3) - pkin(5);
t47 = t81 * t60;
t72 = t81 * t59;
t19 = -t57 * t47 - t58 * t72;
t70 = qJD(2) * t81;
t29 = t60 * qJD(3) + t59 * t70;
t24 = t29 * qJD(1);
t65 = -t59 * qJD(3) + t60 * t70;
t63 = t65 * qJD(1);
t5 = t57 * t24 - t58 * t63;
t87 = t5 * t19;
t45 = qJD(1) * t47;
t86 = t57 * t45;
t38 = t58 * t45;
t62 = qJD(1) ^ 2;
t84 = t60 * t62;
t61 = qJD(2) ^ 2;
t83 = t61 * t59;
t82 = t61 * t60;
t6 = t58 * t24 + t57 * t63;
t44 = qJD(1) * t72;
t40 = qJD(2) * pkin(2) + t44;
t16 = t57 * t40 - t38;
t80 = t59 ^ 2 - t60 ^ 2;
t79 = qJD(2) * t59;
t18 = t58 * t44 + t86;
t77 = qJD(4) - t18;
t75 = pkin(2) * t79;
t73 = -t60 * pkin(2) - pkin(1);
t71 = t59 * t76;
t69 = pkin(1) * t89;
t68 = t73 * qJD(1);
t15 = t58 * t40 + t86;
t46 = qJD(3) + t68;
t8 = t32 * pkin(3) - t35 * qJ(4) + t46;
t67 = t8 * t35 + t5;
t34 = t42 * qJD(2);
t26 = qJD(1) * t34;
t48 = t57 * t71;
t27 = qJD(2) * t74 - t48;
t51 = pkin(2) * t71;
t66 = t26 * pkin(3) - t27 * qJ(4) + t51;
t10 = t57 * t29 - t58 * t65;
t11 = t58 * t29 + t57 * t65;
t20 = -t58 * t47 + t57 * t72;
t64 = t10 * t35 - t11 * t32 + t19 * t27 - t20 * t26 + t5 * t42;
t54 = -t58 * pkin(2) - pkin(3);
t52 = t57 * pkin(2) + qJ(4);
t41 = t57 * t59 - t85;
t37 = qJD(2) * t85 - t57 * t79;
t17 = t57 * t44 - t38;
t14 = t41 * pkin(3) - t42 * qJ(4) + t73;
t13 = qJD(2) * qJ(4) + t16;
t12 = -qJD(2) * pkin(3) + qJD(4) - t15;
t9 = pkin(2) * t78 + t35 * pkin(3) + t32 * qJ(4);
t4 = qJD(2) * qJD(4) + t6;
t3 = t34 * pkin(3) - t37 * qJ(4) - t42 * qJD(4) + t75;
t1 = -t35 * qJD(4) + t66;
t2 = [0, 0, 0, 0.2e1 * t60 * t71, t80 * t89, t82, -t83, 0, -pkin(5) * t82 + t59 * t69, pkin(5) * t83 + t60 * t69, -t15 * t37 - t16 * t34 - t6 * t41 + t64, -t15 * t10 + t16 * t11 + t87 + t6 * t20 + (t46 + t68) * t75, -t10 * qJD(2) + t1 * t41 + t14 * t26 + t3 * t32 + t8 * t34, t12 * t37 - t13 * t34 - t4 * t41 + t64, t11 * qJD(2) - t1 * t42 - t14 * t27 - t3 * t35 - t8 * t37, t1 * t14 + t12 * t10 + t13 * t11 + t4 * t20 + t8 * t3 + t87; 0, 0, 0, -t59 * t84, t80 * t62, 0, 0, 0, t62 * pkin(1) * t59, pkin(1) * t84, (t16 - t17) * t35 + (-t15 + t18) * t32 + (-t26 * t57 - t27 * t58) * pkin(2), t15 * t17 - t16 * t18 + (-t46 * t78 - t5 * t58 + t57 * t6) * pkin(2), t17 * qJD(2) - t9 * t32 - t67, -t52 * t26 + t54 * t27 + (t13 - t17) * t35 + (t12 - t77) * t32, -t8 * t32 + t9 * t35 + (0.2e1 * qJD(4) - t18) * qJD(2) + t6, -t12 * t17 + t77 * t13 + t4 * t52 + t5 * t54 - t8 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, t15 * t35 + t16 * t32 + t51, 0.2e1 * t35 * qJD(2), t88, t48 + (t32 - t74) * qJD(2), t13 * t32 + (-qJD(4) - t12) * t35 + t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35 * t32, -t48 + (t32 + t74) * qJD(2), -t31 - t61, -t13 * qJD(2) + t67;];
tauc_reg = t2;
