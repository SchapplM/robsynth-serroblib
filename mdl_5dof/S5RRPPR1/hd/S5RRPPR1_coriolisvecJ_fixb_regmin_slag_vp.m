% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% tauc_reg [5x17]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:52
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPPR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:51:42
% EndTime: 2022-01-20 09:51:44
% DurationCPUTime: 0.42s
% Computational Cost: add. (562->92), mult. (1200->136), div. (0->0), fcn. (773->8), ass. (0->74)
t64 = sin(pkin(9));
t66 = cos(pkin(9));
t88 = t64 ^ 2 + t66 ^ 2;
t65 = sin(pkin(8));
t69 = sin(qJ(2));
t87 = pkin(1) * qJD(1);
t83 = t69 * t87;
t50 = t65 * t83;
t67 = cos(pkin(8));
t71 = cos(qJ(2));
t82 = t71 * t87;
t36 = t67 * t82 - t50;
t26 = t36 * qJD(2);
t63 = qJD(1) + qJD(2);
t19 = t63 * qJD(4) + t26;
t100 = t19 * t88;
t68 = sin(qJ(5));
t70 = cos(qJ(5));
t43 = t70 * t64 + t68 * t66;
t29 = t43 * t63;
t99 = t63 * t88;
t98 = t36 - qJD(4);
t97 = t66 * pkin(4);
t96 = t67 * pkin(2);
t86 = pkin(1) * qJD(2);
t93 = t67 * t69;
t35 = (t65 * t71 + t93) * t86;
t25 = qJD(1) * t35;
t38 = t43 * qJD(5);
t91 = t70 * t66;
t92 = t68 * t64;
t42 = -t91 + t92;
t44 = t63 * pkin(2) + t82;
t20 = t67 * t44 - t50;
t75 = qJD(4) - t20;
t81 = -pkin(3) - t97;
t9 = t81 * t63 + t75;
t95 = t25 * t42 + t9 * t38;
t37 = t42 * qJD(5);
t94 = t25 * t43 - t9 * t37;
t51 = t67 * t83;
t21 = t65 * t44 + t51;
t56 = t71 * pkin(1) + pkin(2);
t89 = pkin(1) * t93 + t65 * t56;
t53 = t65 * t69 * pkin(1);
t85 = t63 * t92;
t84 = t63 * t91;
t79 = t67 * t56 - t53;
t77 = -pkin(3) - t79;
t76 = t88 * (t63 * qJ(4) + t21);
t74 = (-qJD(2) + t63) * t87;
t73 = (-qJD(1) - t63) * t86;
t72 = t67 * t71 * t86 - qJD(2) * t53;
t59 = t66 * pkin(7);
t55 = t65 * pkin(2) + qJ(4);
t46 = t81 - t96;
t45 = qJD(5) * t84;
t40 = t66 * t55 + t59;
t39 = (-pkin(7) - t55) * t64;
t34 = t65 * t82 + t51;
t33 = t38 * qJD(5);
t32 = t37 * qJD(5);
t31 = qJ(4) + t89;
t30 = qJD(4) + t72;
t27 = -t84 + t85;
t24 = t77 - t97;
t23 = t63 * t38;
t22 = -qJD(5) * t85 + t45;
t18 = t66 * t31 + t59;
t17 = (-pkin(7) - t31) * t64;
t13 = -t63 * pkin(3) + t75;
t2 = t22 * t43 - t29 * t37;
t1 = -t22 * t42 - t43 * t23 + t37 * t27 - t29 * t38;
t3 = [0, 0, 0, 0, t69 * t73, t71 * t73, -t20 * t35 + t21 * t72 - t25 * t79 + t26 * t89, (-t35 * t63 - t25) * t66, t30 * t99 + t100, t100 * t31 + t13 * t35 + t25 * t77 + t76 * t30, t2, t1, -t32, -t33, 0, t24 * t23 + t35 * t27 + ((-t17 * t68 - t18 * t70) * qJD(5) - t43 * t30) * qJD(5) + t95, t24 * t22 + t35 * t29 + ((-t17 * t70 + t18 * t68) * qJD(5) + t42 * t30) * qJD(5) + t94; 0, 0, 0, 0, t69 * t74, t71 * t74, t20 * t34 - t21 * t36 + (-t25 * t67 + t26 * t65) * pkin(2), (t34 * t63 - t25) * t66, -t98 * t99 + t100, t25 * (-pkin(3) - t96) - t13 * t34 + t55 * t100 - t98 * t76, t2, t1, -t32, -t33, 0, t46 * t23 - t34 * t27 + ((-t39 * t68 - t40 * t70) * qJD(5) + t98 * t43) * qJD(5) + t95, t46 * t22 - t34 * t29 + ((-t39 * t70 + t40 * t68) * qJD(5) - t98 * t42) * qJD(5) + t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t32; 0, 0, 0, 0, 0, 0, 0, 0, -t88 * t63 ^ 2, -t76 * t63 + t25, 0, 0, 0, 0, 0, 0.2e1 * t29 * qJD(5), t45 + (-t27 - t85) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29 * t27, -t27 ^ 2 + t29 ^ 2, t45 + (t27 - t85) * qJD(5), 0, 0, -t43 * t19 - t9 * t29, t42 * t19 + t9 * t27;];
tauc_reg = t3;
