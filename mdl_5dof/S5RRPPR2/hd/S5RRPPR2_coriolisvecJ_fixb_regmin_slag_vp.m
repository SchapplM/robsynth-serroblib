% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPPR2
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
% Datum: 2022-01-20 10:06
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPPR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:05:51
% EndTime: 2022-01-20 10:05:53
% DurationCPUTime: 0.45s
% Computational Cost: add. (590->99), mult. (1286->163), div. (0->0), fcn. (742->8), ass. (0->87)
t51 = sin(pkin(8));
t55 = sin(qJ(2));
t85 = pkin(1) * qJD(1);
t78 = t55 * t85;
t36 = t51 * t78;
t53 = cos(pkin(8));
t57 = cos(qJ(2));
t77 = t57 * t85;
t27 = t53 * t77 - t36;
t22 = t27 * qJD(2);
t47 = qJD(1) + qJD(2);
t16 = t47 * qJD(4) + t22;
t50 = sin(pkin(9));
t45 = t50 ^ 2;
t52 = cos(pkin(9));
t87 = t52 ^ 2 + t45;
t105 = t16 * t87;
t104 = qJD(4) - t27;
t84 = pkin(1) * qJD(2);
t91 = t53 * t55;
t26 = (t51 * t57 + t91) * t84;
t21 = qJD(1) * t26;
t60 = -t52 * pkin(4) - t50 * pkin(7) - pkin(3);
t32 = t47 * pkin(2) + t77;
t18 = t53 * t32 - t36;
t65 = qJD(4) - t18;
t3 = t60 * t47 + t65;
t54 = sin(qJ(5));
t56 = cos(qJ(5));
t37 = t53 * t78;
t19 = t51 * t32 + t37;
t14 = t47 * qJ(4) + t19;
t6 = t50 * qJD(3) + t52 * t14;
t92 = t52 * t56;
t96 = t45 * t56;
t103 = (t16 * t92 + t54 * t21 + (t56 * t3 - t54 * t6) * qJD(5)) * t52 + t16 * t96;
t5 = -t52 * qJD(3) + t50 * t14;
t102 = t5 * t50;
t101 = t53 * pkin(2);
t39 = t51 * t55 * pkin(1);
t61 = t53 * t57 * t84 - qJD(2) * t39;
t23 = qJD(4) + t61;
t100 = t23 * t47;
t44 = t47 ^ 2;
t99 = t45 * t44;
t98 = t45 * t47;
t97 = t45 * t54;
t95 = t47 * t50;
t94 = t52 * t47;
t93 = t52 * t54;
t35 = -qJD(5) + t94;
t90 = t56 * t35;
t42 = t57 * pkin(1) + pkin(2);
t88 = pkin(1) * t91 + t51 * t42;
t86 = t54 ^ 2 - t56 ^ 2;
t83 = qJD(5) * t54;
t82 = qJD(5) * t56;
t81 = qJD(5) + t35;
t80 = t5 * t95;
t79 = t47 * t96;
t76 = qJD(5) * t98;
t75 = t50 * t83;
t74 = t50 * t82;
t72 = t53 * t42 - t39;
t67 = -t54 * t3 - t56 * t6;
t2 = t67 * qJD(5) - t16 * t93 + t56 * t21;
t70 = t16 * t97 - t2 * t52 + t5 * t74;
t69 = t35 * t75;
t68 = t81 * t95;
t66 = t6 * t52 + t102;
t64 = (-qJD(2) + t47) * t85;
t63 = (-qJD(1) - t47) * t84;
t62 = t35 * t52 + t98;
t41 = t51 * pkin(2) + qJ(4);
t59 = t104 * t54 + t41 * t82;
t58 = -t35 ^ 2 - t99;
t31 = t75 * t94;
t30 = -0.2e1 * t56 * t54 * t76;
t28 = t60 - t101;
t25 = t51 * t77 + t37;
t24 = qJ(4) + t88;
t17 = 0.2e1 * t86 * t76;
t15 = t60 - t72;
t11 = -t47 * pkin(3) + t65;
t10 = (t35 + t94) * t74;
t9 = t31 + t69;
t1 = [0, 0, 0, 0, t55 * t63, t57 * t63, -t18 * t26 + t19 * t61 - t21 * t72 + t22 * t88, (-t26 * t47 - t21) * t52, t87 * t100 + t105, t21 * (-pkin(3) - t72) + t11 * t26 + t66 * t23 + t24 * t105, t30, t17, t9, t10, 0, -(-t23 * t93 + t56 * t26) * t35 + t97 * t100 + (-(-t15 * t54 - t24 * t92) * t35 + t24 * t79) * qJD(5) + t70, (t23 * t92 + t54 * t26) * t35 + t23 * t79 + (t15 * t90 + (-t62 * t24 - t102) * t54) * qJD(5) + t103; 0, 0, 0, 0, t55 * t64, t57 * t64, t18 * t25 - t19 * t27 + (-t21 * t53 + t22 * t51) * pkin(2), (t25 * t47 - t21) * t52, t104 * t47 * t87 + t105, t21 * (-pkin(3) - t101) - t11 * t25 + t41 * t105 + t104 * t66, t30, t17, t9, t10, 0, (t56 * t25 + t28 * t83 + t59 * t52) * t35 + t59 * t98 + t70, -t54 * t25 * t35 + t104 * t62 * t56 + (t28 * t90 + (-t62 * t41 - t102) * t54) * qJD(5) + t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t35 - t94) * t74, t31 - t69; 0, 0, 0, 0, 0, 0, 0, 0, -t87 * t44, -t66 * t47 + t21, 0, 0, 0, 0, 0, t58 * t54, t58 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44 * t54 * t96, -t86 * t99, -t54 * t68, -t56 * t68, 0, t67 * t35 - t56 * t80 + t2, (-t52 * t16 - t81 * t3) * t56 + (t81 * t6 - t21 + t80) * t54;];
tauc_reg = t1;
