% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPPRR6
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
% tauc_reg [5x22]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPRR6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:58:07
% EndTime: 2019-12-31 17:58:09
% DurationCPUTime: 0.61s
% Computational Cost: add. (937->137), mult. (2371->200), div. (0->0), fcn. (1772->8), ass. (0->83)
t66 = cos(qJ(4));
t61 = cos(pkin(9));
t90 = qJD(1) * t61;
t50 = t66 * t90;
t59 = sin(pkin(9));
t64 = sin(qJ(4));
t92 = t64 * t59;
t84 = qJD(1) * t92;
t38 = t50 - t84;
t36 = qJD(5) - t38;
t45 = -t66 * t61 + t92;
t71 = t45 * qJD(3);
t108 = qJD(1) * t71;
t46 = t66 * t59 + t64 * t61;
t107 = t46 * qJD(1);
t106 = -qJD(5) + t36;
t72 = t46 * qJD(3);
t51 = sin(pkin(8)) * pkin(1) + qJ(3);
t48 = t51 * qJD(1);
t55 = t61 * qJD(2);
t24 = t55 + (-pkin(6) * qJD(1) - t48) * t59;
t32 = t59 * qJD(2) + t61 * t48;
t25 = pkin(6) * t90 + t32;
t9 = t64 * t24 + t66 * t25;
t4 = qJD(1) * t72 + t9 * qJD(4);
t105 = (t107 * pkin(4) + t36 * pkin(7)) * t36 + t4;
t41 = t46 * qJD(4);
t34 = qJD(1) * t41;
t65 = cos(qJ(5));
t30 = t65 * t34;
t40 = t45 * qJD(4);
t63 = sin(qJ(5));
t88 = qJD(5) * t63;
t74 = t65 * t40 + t46 * t88;
t104 = t46 * t30 - t74 * t36;
t28 = t63 * qJD(4) + t107 * t65;
t49 = qJD(4) * t50;
t33 = -qJD(4) * t84 + t49;
t14 = qJD(5) * t28 + t63 * t33;
t102 = pkin(6) + t51;
t42 = t102 * t59;
t43 = t102 * t61;
t76 = -t66 * t42 - t64 * t43;
t10 = t76 * qJD(4) - t71;
t47 = -cos(pkin(8)) * pkin(1) - t61 * pkin(3) - pkin(2);
t37 = t47 * qJD(1) + qJD(3);
t12 = -t38 * pkin(4) - pkin(7) * t107 + t37;
t18 = t45 * pkin(4) - t46 * pkin(7) + t47;
t20 = -t64 * t42 + t66 * t43;
t8 = t66 * t24 - t64 * t25;
t3 = t8 * qJD(4) - t108;
t6 = -qJD(4) * pkin(4) - t8;
t103 = -t20 * t34 - (qJD(5) * t18 + t10) * t36 + t4 * t46 - t6 * t40 - (qJD(5) * t12 + t3) * t45;
t86 = t65 * qJD(4);
t13 = qJD(5) * t86 - t107 * t88 + t65 * t33;
t101 = t13 * t45 + t28 * t41;
t100 = t13 * t63;
t99 = t18 * t34;
t26 = t107 * t63 - t86;
t98 = t26 * t36;
t97 = t28 * t36;
t96 = t28 * t107;
t95 = t107 * t26;
t93 = t63 * t34;
t91 = t59 ^ 2 + t61 ^ 2;
t89 = qJD(5) * t46;
t87 = t40 * qJD(4);
t81 = qJD(1) * t91;
t80 = t36 * t65;
t7 = qJD(4) * pkin(7) + t9;
t1 = t65 * t12 - t63 * t7;
t2 = t63 * t12 + t65 * t7;
t78 = -t45 * t14 - t41 * t26;
t77 = (-t59 * t48 + t55) * t59 - t32 * t61;
t75 = t30 + (t38 * t63 - t88) * t36;
t69 = -pkin(7) * t34 + (t6 + t8) * t36;
t68 = (t63 * t40 - t65 * t89) * t36 - t46 * t93;
t35 = t41 * qJD(4);
t23 = t41 * pkin(4) + t40 * pkin(7);
t16 = t34 * pkin(4) - t33 * pkin(7);
t15 = t65 * t16;
t11 = t20 * qJD(4) + t72;
t5 = [0, 0, 0, 0, 0, 0, 0.2e1 * qJD(3) * t81, (t51 * t81 - t77) * qJD(3), -t107 * t40 + t33 * t46, -t107 * t41 - t33 * t45 - t46 * t34 - t40 * t38, -t87, -t35, 0, -t11 * qJD(4) + t47 * t34 + t37 * t41, -t10 * qJD(4) + t47 * t33 - t37 * t40, t13 * t65 * t46 - t74 * t28, -(-t26 * t65 - t28 * t63) * t40 + (-t100 - t14 * t65 + (t26 * t63 - t28 * t65) * qJD(5)) * t46, t101 + t104, t68 + t78, t34 * t45 + t36 * t41, t1 * t41 + t11 * t26 - t76 * t14 + t15 * t45 + (t99 + t23 * t36 + (-t20 * t36 - t7 * t45 + t6 * t46) * qJD(5)) * t65 + t103 * t63, t11 * t28 - t76 * t13 - t2 * t41 + (-(-qJD(5) * t20 + t23) * t36 - t99 - (-qJD(5) * t7 + t16) * t45 - t6 * t89) * t63 + t103 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, t87, 0, 0, 0, 0, 0, t68 - t78, t101 - t104; 0, 0, 0, 0, 0, 0, -t91 * qJD(1) ^ 2, t77 * qJD(1), 0, 0, 0, 0, 0, 0.2e1 * t107 * qJD(4), t49 + (t38 - t84) * qJD(4), 0, 0, 0, 0, 0, t75 - t95, -t36 ^ 2 * t65 - t93 - t96; 0, 0, 0, 0, 0, 0, 0, 0, -t107 * t38, t107 ^ 2 - t38 ^ 2, t49 + (-t38 - t84) * qJD(4), 0, 0, (-qJD(3) - t37) * t107, -t37 * t38 + t108, t28 * t80 + t100, (t13 - t98) * t65 + (-t14 - t97) * t63, t36 * t80 + t93 - t96, t75 + t95, -t36 * t107, -pkin(4) * t14 - t1 * t107 - t105 * t65 - t9 * t26 + t69 * t63, -pkin(4) * t13 + t105 * t63 + t107 * t2 - t9 * t28 + t69 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28 * t26, -t26 ^ 2 + t28 ^ 2, t13 + t98, -t14 + t97, t34, t106 * t2 - t6 * t28 - t63 * t3 + t15, t106 * t1 - t63 * t16 + t6 * t26 - t65 * t3;];
tauc_reg = t5;
