% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% tauc_reg [5x21]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:15
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPRR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:14:36
% EndTime: 2022-01-23 09:14:37
% DurationCPUTime: 0.52s
% Computational Cost: add. (791->116), mult. (2027->173), div. (0->0), fcn. (1564->8), ass. (0->81)
t64 = cos(pkin(9));
t69 = cos(qJ(4));
t94 = t69 * t64;
t62 = sin(pkin(9));
t67 = sin(qJ(4));
t96 = t67 * t62;
t77 = -t94 + t96;
t43 = t77 * qJD(1);
t68 = cos(qJ(5));
t36 = t68 * t43;
t89 = t69 * qJD(4);
t92 = qJD(1) * t64;
t54 = t89 * t92;
t88 = qJD(1) * t96;
t39 = -qJD(4) * t88 + t54;
t51 = t69 * t62 + t67 * t64;
t46 = t51 * qJD(4);
t40 = qJD(1) * t46;
t44 = t51 * qJD(1);
t66 = sin(qJ(5));
t91 = qJD(5) * t66;
t4 = -qJD(5) * t36 + t68 * t39 - t66 * t40 - t44 * t91;
t18 = t66 * t44 + t36;
t61 = qJD(4) + qJD(5);
t98 = t18 * t61;
t112 = t4 + t98;
t80 = -t66 * t43 + t68 * t44;
t111 = t80 * t18;
t5 = t80 * qJD(5) + t66 * t39 + t68 * t40;
t99 = t80 * t61;
t110 = -t5 + t99;
t109 = -t18 ^ 2 + t80 ^ 2;
t56 = sin(pkin(8)) * pkin(1) + qJ(3);
t53 = t56 * qJD(1);
t58 = t64 * qJD(2);
t27 = t58 + (-pkin(6) * qJD(1) - t53) * t62;
t35 = t62 * qJD(2) + t64 * t53;
t28 = pkin(6) * t92 + t35;
t82 = -t67 * t27 - t69 * t28;
t10 = -t43 * pkin(7) - t82;
t52 = -cos(pkin(8)) * pkin(1) - pkin(2) - t64 * pkin(3);
t42 = t52 * qJD(1) + qJD(3);
t22 = t43 * pkin(4) + t42;
t75 = t51 * qJD(3);
t74 = qJD(1) * t75;
t3 = -t39 * pkin(7) + t82 * qJD(4) - t74;
t108 = t22 * t18 + t10 * t91 + (-t10 * t61 - t3) * t66;
t106 = t69 * t27 - t67 * t28;
t105 = qJD(3) * t43;
t104 = qJD(5) - t61;
t2 = -t40 * pkin(7) + t106 * qJD(4) - t105;
t103 = -t66 * t2 - t22 * t80 + t68 * t3;
t102 = pkin(4) * t44;
t101 = pkin(6) + t56;
t45 = t77 * qJD(4);
t78 = -t66 * t51 - t68 * t77;
t11 = t78 * qJD(5) - t68 * t45 - t66 * t46;
t100 = t11 * t61;
t95 = t68 * t10;
t93 = t62 ^ 2 + t64 ^ 2;
t90 = t45 * qJD(4);
t9 = -t44 * pkin(7) + t106;
t8 = qJD(4) * pkin(4) + t9;
t87 = -pkin(4) * t61 - t8;
t84 = qJD(1) * t93;
t81 = (-t62 * t53 + t58) * t62 - t35 * t64;
t47 = t101 * t62;
t48 = t101 * t64;
t79 = t67 * t47 - t69 * t48;
t24 = t68 * t51 - t66 * t77;
t73 = -t47 * t89 + qJD(3) * t94 + (-qJD(3) * t62 - qJD(4) * t48) * t67;
t71 = t79 * qJD(4) - t75;
t41 = t46 * qJD(4);
t26 = pkin(4) * t77 + t52;
t16 = -pkin(7) * t77 - t79;
t15 = -t51 * pkin(7) - t69 * t47 - t67 * t48;
t14 = t45 * pkin(7) + t71;
t13 = -t46 * pkin(7) + t73;
t12 = t24 * qJD(5) - t66 * t45 + t68 * t46;
t7 = t12 * t61;
t1 = [0, 0, 0, 0, 0, 0.2e1 * qJD(3) * t84, (t56 * t84 - t81) * qJD(3), t39 * t51 - t44 * t45, -t39 * t77 - t51 * t40 + t45 * t43 - t44 * t46, -t90, -t41, 0, t71 * qJD(4) + t52 * t40 + t42 * t46, -t73 * qJD(4) + t52 * t39 - t42 * t45, t11 * t80 + t4 * t24, -t11 * t18 - t12 * t80 - t24 * t5 + t4 * t78, t100, -t7, 0, t26 * t5 + t22 * t12 + (-t66 * t13 + t68 * t14 + (-t15 * t66 - t16 * t68) * qJD(5)) * t61 + (t46 * t18 - t40 * t78) * pkin(4), t26 * t4 + t22 * t11 - (t68 * t13 + t66 * t14 + (t15 * t68 - t16 * t66) * qJD(5)) * t61 + (t40 * t24 + t46 * t80) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, t90, 0, 0, 0, 0, 0, -t7, -t100; 0, 0, 0, 0, 0, -t93 * qJD(1) ^ 2, t81 * qJD(1), 0, 0, 0, 0, 0, 0.2e1 * t44 * qJD(4), t54 + (-t43 - t88) * qJD(4), 0, 0, 0, 0, 0, t5 + t99, t4 - t98; 0, 0, 0, 0, 0, 0, 0, t44 * t43, -t43 ^ 2 + t44 ^ 2, t54 + (t43 - t88) * qJD(4), 0, 0, -t42 * t44 - t74, t42 * t43 + t105, t111, t109, t112, t110, 0, -t18 * t102 - (-t66 * t9 - t95) * t61 + (t87 * t66 - t95) * qJD(5) + t103, -t80 * t102 + (t87 * qJD(5) + t9 * t61 - t2) * t68 + t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t111, t109, t112, t110, 0, t104 * (-t66 * t8 - t95) + t103, (-t104 * t8 - t2) * t68 + t108;];
tauc_reg = t1;
