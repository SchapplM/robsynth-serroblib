% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% tauc_reg [5x17]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:31
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:30:53
% EndTime: 2022-01-20 11:30:54
% DurationCPUTime: 0.40s
% Computational Cost: add. (769->100), mult. (1687->146), div. (0->0), fcn. (960->8), ass. (0->81)
t52 = qJD(1) + qJD(2);
t62 = cos(qJ(2));
t92 = pkin(1) * qJD(1);
t85 = t62 * t92;
t40 = t52 * pkin(2) + t85;
t59 = sin(qJ(2));
t61 = cos(qJ(3));
t100 = t59 * t61;
t58 = sin(qJ(3));
t75 = -t58 * t62 - t100;
t88 = qJD(3) * t61;
t66 = (t75 * qJD(2) - t59 * t88) * pkin(1);
t65 = qJD(1) * t66;
t89 = qJD(3) * t58;
t64 = -t40 * t89 + t65;
t81 = qJD(2) * t85;
t86 = t59 * t92;
t82 = t58 * t86;
t95 = (qJD(2) + qJD(3)) * t82;
t15 = (qJD(3) * t40 + t81) * t61 - t95;
t50 = qJD(3) + t52;
t63 = qJD(5) ^ 2;
t56 = cos(pkin(9));
t103 = t56 * t58;
t47 = t61 * pkin(2) + pkin(3);
t55 = sin(pkin(9));
t94 = pkin(2) * t103 + t55 * t47;
t35 = t75 * t92;
t101 = t58 * t59;
t74 = t61 * t62 - t101;
t36 = t74 * t92;
t91 = pkin(2) * qJD(3);
t97 = -t56 * t35 + t55 * t36 - (t55 * t61 + t103) * t91;
t109 = -t97 * t50 + (pkin(8) + t94) * t63;
t2 = t55 * t15 - t56 * t64;
t57 = sin(qJ(5));
t60 = cos(qJ(5));
t26 = t58 * t40 + t61 * t86;
t106 = t55 * t26;
t25 = t61 * t40 - t82;
t23 = t50 * pkin(3) + t25;
t10 = t56 * t23 - t106;
t8 = -t50 * pkin(4) - t10;
t90 = qJD(5) * t8;
t108 = t2 * t57 + t60 * t90;
t105 = t55 * t58;
t104 = t56 * t26;
t102 = t57 * t60;
t99 = t63 * t57;
t48 = t62 * pkin(1) + pkin(2);
t34 = -pkin(1) * t101 + t61 * t48 + pkin(3);
t37 = pkin(1) * t100 + t58 * t48;
t98 = t55 * t34 + t56 * t37;
t96 = -t55 * t35 - t56 * t36 + (t56 * t61 - t105) * t91;
t93 = t57 ^ 2 - t60 ^ 2;
t87 = 0.2e1 * qJD(5) * t50;
t3 = t56 * t15 + t64 * t55;
t83 = -t8 * t50 - t3;
t80 = (-qJD(2) + t52) * t92;
t79 = pkin(1) * qJD(2) * (-qJD(1) - t52);
t20 = t48 * t88 + (t74 * qJD(2) - t59 * t89) * pkin(1);
t21 = -t48 * t89 + t66;
t4 = t55 * t20 - t56 * t21;
t78 = (pkin(8) + t98) * t63 + t4 * t50;
t12 = t55 * t25 + t104;
t77 = -t12 * t50 + (t55 * pkin(3) + pkin(8)) * t63;
t76 = t56 * t34 - t55 * t37;
t5 = t56 * t20 + t55 * t21;
t73 = qJD(5) * ((-pkin(4) - t76) * t50 - t5);
t72 = (-pkin(2) * t50 - t40) * qJD(3);
t13 = t56 * t25 - t106;
t71 = qJD(5) * ((-t56 * pkin(3) - pkin(4)) * t50 + t13);
t70 = -pkin(2) * t105 + t56 * t47;
t69 = qJD(5) * ((-pkin(4) - t70) * t50 - t96);
t51 = t63 * t60;
t49 = t50 ^ 2;
t39 = t87 * t102;
t27 = t93 * t87;
t11 = t55 * t23 + t104;
t6 = t57 * t90;
t1 = [0, 0, 0, 0, t59 * t79, t62 * t79, 0, t21 * t50 + t64, -t20 * t50 - t15, -t10 * t4 + t11 * t5 - t2 * t76 + t3 * t98, t39, -t27, t51, -t99, 0, t6 + t57 * t73 + (-t2 - t78) * t60, t78 * t57 + t60 * t73 + t108; 0, 0, 0, 0, t59 * t80, t62 * t80, 0, -t35 * t50 + t58 * t72 + t65, t36 * t50 + (t72 - t81) * t61 + t95, t97 * t10 + t96 * t11 - t2 * t70 + t3 * t94, t39, -t27, t51, -t99, 0, t6 + t57 * t69 + (-t109 - t2) * t60, t109 * t57 + t60 * t69 + t108; 0, 0, 0, 0, 0, 0, 0, t26 * t50 + t64, t25 * t50 - t15, t10 * t12 - t11 * t13 + (-t2 * t56 + t3 * t55) * pkin(3), t39, -t27, t51, -t99, 0, t6 + t57 * t71 + (-t2 - t77) * t60, t77 * t57 + t60 * t71 + t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t99, -t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49 * t102, t93 * t49, 0, 0, 0, t83 * t57, t83 * t60;];
tauc_reg = t1;
