% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRRP3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:14:15
% EndTime: 2019-12-31 17:14:16
% DurationCPUTime: 0.45s
% Computational Cost: add. (467->108), mult. (964->148), div. (0->0), fcn. (420->4), ass. (0->84)
t43 = qJD(1) + qJD(2);
t46 = sin(qJ(3));
t44 = t46 ^ 2;
t48 = cos(qJ(3));
t45 = t48 ^ 2;
t83 = t44 + t45;
t100 = t43 * t83;
t99 = 2 * qJD(3);
t49 = cos(qJ(2));
t81 = pkin(1) * qJD(2);
t66 = qJD(1) * t81;
t61 = t49 * t66;
t34 = t48 * t61;
t47 = sin(qJ(2));
t82 = pkin(1) * qJD(1);
t72 = t47 * t82;
t27 = t43 * pkin(6) + t72;
t90 = t46 * t27;
t4 = t34 + (qJD(4) - t90) * qJD(3);
t32 = t46 * t61;
t79 = qJD(3) * t48;
t7 = t27 * t79 + t32;
t98 = t4 * t48 + t7 * t46;
t97 = t49 * t83;
t84 = -t44 + t45;
t96 = t84 * t43 * t99;
t50 = qJD(3) ^ 2;
t95 = pkin(6) * t50;
t94 = t49 * pkin(1);
t29 = -t48 * pkin(3) - t46 * qJ(4) - pkin(2);
t93 = t29 * t43;
t92 = t43 * t46;
t91 = t43 * t48;
t89 = t48 * t27;
t40 = t50 * t46;
t41 = t50 * t48;
t71 = t49 * t82;
t28 = -t43 * pkin(2) - t71;
t62 = t47 * t66;
t88 = t28 * t79 + t46 * t62;
t73 = t49 * t81;
t87 = t73 * t100;
t78 = qJD(3) * t49;
t68 = t46 * t78;
t86 = t68 * t82 + t72 * t91;
t85 = t83 * t61;
t80 = qJD(3) * t46;
t74 = qJD(3) * qJ(4);
t13 = t74 + t89;
t77 = t13 * qJD(3);
t76 = t46 * qJD(4);
t75 = -qJD(1) - t43;
t70 = t43 * t80;
t69 = t43 * t79;
t58 = pkin(3) * t46 - qJ(4) * t48;
t1 = t62 + (t58 * qJD(3) - t76) * t43;
t67 = -t1 - t95;
t64 = -qJD(3) * pkin(3) + qJD(4);
t63 = t46 * t69;
t60 = (-qJD(2) + t43) * t82;
t59 = t75 * t81;
t10 = t64 + t90;
t57 = t10 * t46 + t13 * t48;
t37 = t47 * pkin(1) + pkin(6);
t15 = pkin(3) * t80 - t48 * t74 - t76;
t9 = t47 * t81 + t15;
t56 = -t37 * t50 - t43 * t9 - t1;
t20 = t29 - t94;
t55 = t20 * t43 - t73;
t54 = t10 * t79 - t46 * t77 + t98;
t53 = t71 * t100;
t52 = -t47 * t92 + t48 * t78;
t51 = (t10 * t48 - t13 * t46) * qJD(3) + t98;
t42 = t43 ^ 2;
t38 = -pkin(2) - t94;
t36 = t46 * t42 * t48;
t26 = -0.2e1 * t63;
t25 = 0.2e1 * t63;
t21 = t84 * t42;
t18 = t28 * t80;
t16 = t58 * t43;
t6 = -t71 + t93;
t3 = t6 * t80;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47 * t59, t49 * t59, 0, 0, t25, t96, t41, t26, -t40, 0, t38 * t70 - t37 * t41 + t18 + (t75 * t48 * t47 - t68) * t81, t37 * t40 + t38 * t69 - t52 * t81 + t88, t85 + t87, ((qJD(1) * t38 + t28) * t47 + (qJD(1) * t37 + t27) * t97) * t81, t25, t41, -t96, 0, t40, t26, t56 * t48 + t55 * t80 + t3, t54 + t87, t56 * t46 + (-t55 - t6) * t79, t1 * t20 + t51 * t37 + t57 * t73 + t6 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47 * t60, t49 * t60, 0, 0, t25, t96, t41, t26, -t40, 0, -pkin(2) * t70 + t18 + (-t62 - t95) * t48 + t86, -pkin(2) * t69 + pkin(6) * t40 + t52 * t82 + t88, -t53 + t85, ((-pkin(2) * qJD(2) - t28) * t47 + (pkin(6) * qJD(2) - t27) * t97) * t82, t25, t41, -t96, 0, t40, t26, t29 * t70 + t3 + (-t15 * t43 + t67) * t48 + t86, -t53 + t54, (-t6 - t71 - t93) * t79 + ((-t15 + t72) * t43 + t67) * t46, t1 * t29 + t6 * t15 + (-t47 * t6 - t57 * t49) * t82 + t51 * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t21, 0, t36, 0, 0, -t28 * t92 - t32, -t28 * t91 - t34, 0, 0, -t36, 0, t21, 0, 0, t36, -t32 + (t16 * t48 - t46 * t6) * t43, ((t13 - t74) * t46 + (-t10 + t64) * t48) * t43, qJD(4) * t99 + t34 + (t16 * t46 + t48 * t6) * t43, -t10 * t89 - t7 * pkin(3) + t4 * qJ(4) - t6 * t16 + (qJD(4) + t90) * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, 0, -t44 * t42 - t50, t6 * t92 + t7 - t77;];
tauc_reg = t2;
