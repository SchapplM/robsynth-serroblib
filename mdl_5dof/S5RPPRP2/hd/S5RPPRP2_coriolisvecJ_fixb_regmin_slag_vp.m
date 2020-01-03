% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% tauc_reg [5x19]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPRP2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:49:24
% EndTime: 2019-12-31 17:49:25
% DurationCPUTime: 0.34s
% Computational Cost: add. (678->114), mult. (1648->148), div. (0->0), fcn. (1134->6), ass. (0->71)
t58 = sin(pkin(8));
t60 = cos(pkin(8));
t62 = sin(qJ(4));
t63 = cos(qJ(4));
t43 = t63 * t58 + t62 * t60;
t35 = t43 * qJD(1);
t92 = t35 ^ 2;
t52 = sin(pkin(7)) * pkin(1) + qJ(3);
t91 = pkin(6) + t52;
t44 = -cos(pkin(7)) * pkin(1) - pkin(2) - t60 * pkin(3);
t31 = t44 * qJD(1) + qJD(3);
t81 = qJD(1) * t60;
t76 = t63 * t81;
t87 = t62 * t58;
t77 = qJD(1) * t87;
t33 = -t76 + t77;
t10 = t33 * pkin(4) - t35 * qJ(5) + t31;
t90 = t10 * t35;
t89 = t35 * t33;
t45 = t52 * qJD(1);
t26 = t58 * qJD(2) + t60 * t45;
t24 = pkin(6) * t81 + t26;
t88 = t62 * t24;
t38 = t43 * qJD(4);
t28 = qJD(1) * t38;
t79 = qJD(4) * t63;
t80 = qJD(4) * t62;
t37 = t58 * t80 - t60 * t79;
t86 = -t43 * t28 + t37 * t33;
t85 = t58 ^ 2 + t60 ^ 2;
t42 = -t63 * t60 + t87;
t39 = t91 * t58;
t40 = t91 * t60;
t68 = -t63 * t39 - t62 * t40;
t8 = -t42 * qJD(3) + t68 * qJD(4);
t84 = t8 * qJD(4);
t14 = -t62 * t39 + t63 * t40;
t9 = t43 * qJD(3) + t14 * qJD(4);
t83 = t9 * qJD(4);
t55 = t60 * qJD(2);
t23 = t55 + (-pkin(6) * qJD(1) - t45) * t58;
t6 = t63 * t23 - t88;
t82 = qJD(5) - t6;
t30 = t38 * qJD(4);
t78 = qJD(1) * qJD(3);
t75 = t6 + t88;
t74 = t62 * t78;
t73 = t63 * t78;
t72 = qJD(1) * t85;
t2 = t23 * t80 + t24 * t79 + t58 * t73 + t60 * t74;
t49 = qJD(4) * t76;
t27 = qJD(4) * t77 - t49;
t71 = t28 * pkin(4) + t27 * qJ(5);
t7 = t62 * t23 + t63 * t24;
t70 = (-t58 * t45 + t55) * t58 - t26 * t60;
t69 = -t42 * t27 + t38 * t35;
t67 = t7 * qJD(4) - t2;
t66 = -t23 * t79 + t58 * t74 - t60 * t73;
t65 = 0.2e1 * t35 * qJD(4);
t32 = t33 ^ 2;
t29 = t37 * qJD(4);
t17 = t35 * pkin(4) + t33 * qJ(5);
t16 = t49 + (t33 - t77) * qJD(4);
t15 = -t49 + (t33 + t77) * qJD(4);
t12 = t42 * pkin(4) - t43 * qJ(5) + t44;
t11 = t38 * pkin(4) + t37 * qJ(5) - t43 * qJD(5);
t5 = -t35 * qJD(5) + t71;
t4 = qJD(4) * qJ(5) + t7;
t3 = -qJD(4) * pkin(4) + t82;
t1 = (qJD(5) - t88) * qJD(4) - t66;
t13 = [0, 0, 0, 0, 0, 0, 0.2e1 * qJD(3) * t72, (t52 * t72 - t70) * qJD(3), -t27 * t43 - t35 * t37, -t69 + t86, -t29, -t30, 0, t44 * t28 + t31 * t38 - t83, -t44 * t27 - t31 * t37 - t84, t10 * t38 + t11 * t33 + t12 * t28 + t5 * t42 - t83, -t1 * t42 - t14 * t28 + t2 * t43 + t27 * t68 - t3 * t37 - t8 * t33 + t9 * t35 - t4 * t38, t10 * t37 - t11 * t35 + t12 * t27 - t5 * t43 + t84, t1 * t14 + t10 * t11 + t5 * t12 - t2 * t68 + t3 * t9 + t4 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, t29, -t30, t69 + t86, -t29, t1 * t43 + t2 * t42 + t3 * t38 - t4 * t37; 0, 0, 0, 0, 0, 0, -t85 * qJD(1) ^ 2, t70 * qJD(1), 0, 0, 0, 0, 0, t65, -t15, t65, -t32 - t92, t15, t4 * t33 + (-qJD(5) - t3) * t35 + t71; 0, 0, 0, 0, 0, 0, 0, 0, t89, -t32 + t92, t16, 0, 0, -t31 * t35 + t67, t75 * qJD(4) + t31 * t33 + t66, -t17 * t33 + t67 - t90, pkin(4) * t27 - t28 * qJ(5) + (t4 - t7) * t35 + (t3 - t82) * t33, -t10 * t33 + t17 * t35 + (0.2e1 * qJD(5) - t75) * qJD(4) - t66, -t2 * pkin(4) + t1 * qJ(5) - t10 * t17 - t3 * t7 + t82 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, t16, -qJD(4) ^ 2 - t92, -t4 * qJD(4) + t2 + t90;];
tauc_reg = t13;
