% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRPR5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR5_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:03:35
% EndTime: 2019-12-31 17:03:36
% DurationCPUTime: 0.42s
% Computational Cost: add. (525->126), mult. (701->141), div. (0->0), fcn. (317->8), ass. (0->88)
t40 = sin(qJ(2));
t88 = pkin(1) * qJD(1);
t72 = t40 * t88;
t35 = qJD(1) + qJD(2);
t86 = t35 * qJ(3);
t12 = t72 + t86;
t38 = qJ(1) + qJ(2);
t30 = sin(t38);
t31 = cos(t38);
t93 = -g(1) * t30 + g(2) * t31;
t63 = -t12 * t35 + t93;
t79 = qJDD(1) * t40;
t43 = cos(qJ(2));
t85 = qJD(1) * t43;
t94 = g(1) * t31 + g(2) * t30;
t106 = -((qJD(2) - t35) * t85 + t79) * pkin(1) + t94;
t45 = -pkin(2) - pkin(6);
t60 = -pkin(1) * t85 + qJD(3);
t39 = sin(qJ(4));
t36 = t39 ^ 2;
t42 = cos(qJ(4));
t37 = t42 ^ 2;
t89 = t36 + t37;
t105 = t89 * (t45 * t35 + t60);
t34 = qJDD(1) + qJDD(2);
t104 = t34 * t45;
t33 = t35 ^ 2;
t41 = sin(qJ(1));
t102 = g(1) * t41;
t101 = t34 * pkin(2);
t100 = t43 * pkin(1);
t83 = qJD(2) * t43;
t50 = (qJD(1) * t83 + t79) * pkin(1);
t81 = t35 * qJD(3);
t87 = t34 * qJ(3);
t5 = t50 + t81 + t87;
t82 = qJD(4) * t42;
t99 = t12 * t82 + t5 * t39;
t97 = t12 * t43;
t96 = t42 * t34;
t95 = t31 * pkin(2) + t30 * qJ(3);
t92 = -qJD(2) * t72 + qJDD(1) * t100;
t46 = qJD(4) ^ 2;
t91 = -t33 - t46;
t90 = t36 - t37;
t84 = qJD(2) * t40;
t27 = -pkin(2) - t100;
t17 = -pkin(6) + t27;
t78 = qJDD(4) * t17;
t77 = qJDD(4) * t39;
t76 = qJDD(4) * t45;
t75 = t42 * t33 * t39;
t44 = cos(qJ(1));
t74 = t44 * pkin(1) + t95;
t73 = pkin(1) * t84;
t71 = t35 * t84;
t69 = qJDD(3) - t92;
t4 = t69 + t104;
t70 = t89 * t4;
t19 = t31 * qJ(3);
t68 = -t30 * pkin(2) + t19;
t67 = -t92 + t93;
t66 = pkin(1) * t71;
t65 = t35 * t72;
t64 = t35 * t39 * t82;
t61 = -g(2) * t44 + t102;
t15 = pkin(1) * t83 + qJD(3);
t20 = t40 * pkin(1) + qJ(3);
t59 = t12 * t15 + t5 * t20;
t58 = t5 * qJ(3) + t12 * qJD(3);
t6 = t69 - t101;
t57 = g(1) * (t45 * t30 + t19);
t56 = t20 * t35 + t73;
t55 = -t72 + t86;
t52 = -t4 - t63;
t51 = t65 - t67;
t49 = t15 * t35 - t17 * t46 + t20 * t34 - t94;
t48 = t60 * t35 - t45 * t46 + t87 - t94;
t29 = qJDD(4) * t42;
t21 = t31 * pkin(6);
t14 = -t46 * t39 + t29;
t13 = -t46 * t42 - t77;
t11 = -t35 * pkin(2) + t60;
t8 = t37 * t34 - 0.2e1 * t64;
t7 = t36 * t34 + 0.2e1 * t64;
t3 = t5 * t42;
t1 = 0.2e1 * t90 * t35 * qJD(4) - 0.2e1 * t39 * t96;
t2 = [0, 0, 0, 0, 0, qJDD(1), t61, g(1) * t44 + g(2) * t41, 0, 0, 0, 0, 0, 0, 0, t34, (t34 * t43 - t71) * pkin(1) - t67, ((-qJDD(1) - t34) * t40 + (-qJD(1) - t35) * t83) * pkin(1) + t94, 0, (t61 + (t40 ^ 2 + t43 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t34, 0, 0, 0, 0, 0, 0, t66 + qJDD(3) + (-pkin(2) + t27) * t34 + t67, (qJD(3) + t15) * t35 + (qJ(3) + t20) * t34 + t50 - t94, t6 * t27 + t11 * t73 - g(1) * (-t41 * pkin(1) + t68) - g(2) * t74 + t59, t8, t1, t14, t7, t13, 0, (t56 * qJD(4) + t78) * t42 + t49 * t39 + t99, t3 + (-t78 + (-t12 - t56) * qJD(4)) * t39 + t49 * t42, -t93 + t89 * (-t17 * t34 - t4 - t66), -t57 - g(2) * (t21 + t74) + t17 * t70 + (t84 * t105 + t102) * pkin(1) + t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t51, t106, 0, 0, t34, 0, 0, 0, 0, 0, 0, qJDD(3) - t51 - 0.2e1 * t101, -t106 + 0.2e1 * t81 + 0.2e1 * t87, -t6 * pkin(2) - g(1) * t68 - g(2) * t95 + (-t11 * t40 - t97) * t88 + t58, t8, t1, t14, t7, t13, 0, (t55 * qJD(4) + t76) * t42 + t48 * t39 + t99, t3 + (-t76 + (-t12 - t55) * qJD(4)) * t39 + t48 * t42, -t93 + t89 * (-t4 + t65 - t104), -t57 - g(2) * (t21 + t95) + t45 * t70 + (-t40 * t105 - t97) * t88 + t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t33, t6 + t63, 0, 0, 0, 0, 0, 0, t91 * t39 + t29, t91 * t42 - t77, -t89 * t34, t70 + t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, -t90 * t33, t96, -t75, -t39 * t34, qJDD(4), g(3) * t39 - t52 * t42, g(3) * t42 + t52 * t39, 0, 0;];
tau_reg = t2;
