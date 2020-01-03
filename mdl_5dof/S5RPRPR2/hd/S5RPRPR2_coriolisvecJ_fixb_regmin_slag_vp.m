% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tauc_reg [5x18]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:34:02
% EndTime: 2020-01-03 11:34:05
% DurationCPUTime: 0.44s
% Computational Cost: add. (572->82), mult. (1201->123), div. (0->0), fcn. (750->8), ass. (0->65)
t59 = sin(pkin(9));
t61 = cos(pkin(9));
t80 = t59 ^ 2 + t61 ^ 2;
t58 = qJD(1) + qJD(3);
t50 = cos(pkin(8)) * pkin(1) + pkin(2);
t47 = t50 * qJD(1);
t64 = sin(qJ(3));
t89 = pkin(1) * sin(pkin(8));
t75 = qJD(3) * t89;
t73 = qJD(1) * t75;
t66 = cos(qJ(3));
t79 = qJD(3) * t66;
t71 = t47 * t79 - t64 * t73;
t15 = t58 * qJD(4) + t71;
t96 = t80 * t15;
t63 = sin(qJ(5));
t65 = cos(qJ(5));
t39 = t65 * t59 + t63 * t61;
t30 = t39 * t58;
t95 = t58 * t80;
t83 = t64 * t47;
t21 = qJD(3) * t83 + t66 * t73;
t67 = t64 * t50 + t66 * t89;
t31 = t67 * qJD(3);
t94 = -t31 * t58 - t21;
t76 = qJD(1) * t89;
t26 = t66 * t76 + t83;
t93 = t26 * t58 - t21;
t25 = t66 * t47 - t64 * t76;
t92 = t25 - qJD(4);
t88 = t61 * pkin(4);
t51 = -pkin(3) - t88;
t11 = t51 * t58 - t92;
t82 = t65 * t61;
t84 = t63 * t59;
t38 = -t82 + t84;
t35 = t38 * qJD(5);
t91 = -t11 * t35 + t21 * t39;
t36 = t39 * qJD(5);
t90 = t11 * t36 + t21 * t38;
t85 = t58 * t59;
t78 = t58 * t84;
t77 = t58 * t82;
t72 = -t66 * t50 + t64 * t89 - pkin(3);
t69 = t80 * (t58 * qJ(4) + t26);
t68 = t50 * t79 - t64 * t75;
t54 = t61 * pkin(7);
t43 = t61 * qJ(4) + t54;
t42 = (-pkin(7) - qJ(4)) * t59;
t40 = qJD(5) * t77;
t34 = t36 * qJD(5);
t33 = t35 * qJD(5);
t32 = qJ(4) + t67;
t28 = -t77 + t78;
t27 = qJD(4) + t68;
t24 = t72 - t88;
t23 = t58 * t36;
t22 = -qJD(5) * t78 + t40;
t20 = t61 * t32 + t54;
t19 = (-pkin(7) - t32) * t59;
t17 = t21 * t59;
t16 = -t58 * pkin(3) - t92;
t2 = t22 * t39 - t30 * t35;
t1 = -t22 * t38 - t39 * t23 + t35 * t28 - t30 * t36;
t3 = [0, 0, 0, 0, 0, t94, -t68 * t58 - t71, t94 * t61, t31 * t85 + t17, t27 * t95 + t96, t16 * t31 + t21 * t72 + t69 * t27 + t32 * t96, t2, t1, -t33, -t34, 0, t24 * t23 + t31 * t28 + ((-t19 * t63 - t20 * t65) * qJD(5) - t39 * t27) * qJD(5) + t90, t24 * t22 + t31 * t30 + ((-t19 * t65 + t20 * t63) * qJD(5) + t38 * t27) * qJD(5) + t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, t33; 0, 0, 0, 0, 0, t93, t25 * t58 - t71, t93 * t61, -t26 * t85 + t17, -t92 * t95 + t96, -t21 * pkin(3) + qJ(4) * t96 - t16 * t26 - t92 * t69, t2, t1, -t33, -t34, 0, t51 * t23 - t26 * t28 + ((-t42 * t63 - t43 * t65) * qJD(5) + t92 * t39) * qJD(5) + t90, t51 * t22 - t26 * t30 + ((-t42 * t65 + t43 * t63) * qJD(5) - t92 * t38) * qJD(5) + t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80 * t58 ^ 2, -t69 * t58 + t21, 0, 0, 0, 0, 0, 0.2e1 * t30 * qJD(5), t40 + (-t28 - t78) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30 * t28, -t28 ^ 2 + t30 ^ 2, t40 + (t28 - t78) * qJD(5), 0, 0, -t11 * t30 - t39 * t15, t11 * t28 + t38 * t15;];
tauc_reg = t3;
