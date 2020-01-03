% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4PRRR3
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PRRR3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR3_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:31:41
% EndTime: 2019-12-31 16:31:42
% DurationCPUTime: 0.38s
% Computational Cost: add. (471->100), mult. (737->136), div. (0->0), fcn. (381->8), ass. (0->75)
t45 = cos(qJ(3));
t82 = t45 * pkin(2);
t30 = qJDD(2) * t82;
t43 = sin(qJ(3));
t74 = pkin(2) * qJD(2);
t66 = t43 * t74;
t36 = qJDD(2) + qJDD(3);
t85 = t36 * pkin(3);
t8 = qJD(3) * t66 - t30 - t85;
t37 = pkin(7) + qJ(2);
t34 = qJ(3) + t37;
t27 = cos(t34);
t86 = g(2) * t27;
t89 = t8 + t86;
t42 = sin(qJ(4));
t44 = cos(qJ(4));
t38 = qJD(2) + qJD(3);
t14 = t38 * pkin(6) + t66;
t6 = t44 * qJD(1) - t42 * t14;
t73 = t6 * qJD(4);
t70 = qJDD(2) * t43;
t72 = qJD(3) * t45;
t9 = t36 * pkin(6) + (qJD(2) * t72 + t70) * pkin(2);
t2 = t42 * qJDD(1) + t44 * t9 + t73;
t33 = t44 * qJDD(1);
t80 = t44 * t14;
t7 = t42 * qJD(1) + t80;
t3 = -t7 * qJD(4) - t42 * t9 + t33;
t49 = t2 * t44 - t3 * t42 + (-t42 * t7 - t44 * t6) * qJD(4);
t39 = t42 ^ 2;
t40 = t44 ^ 2;
t75 = t39 + t40;
t88 = t38 * t75;
t26 = sin(t34);
t77 = g(1) * t27 + g(2) * t26;
t65 = t45 * t74;
t84 = t38 * pkin(3);
t15 = -t65 - t84;
t57 = t42 * t6 - t44 * t7;
t87 = -t15 * t43 + t57 * t45;
t24 = g(1) * t26;
t83 = t43 * pkin(2);
t79 = t44 * t36;
t78 = t15 * qJD(4) * t42 + t44 * t24;
t76 = t39 - t40;
t71 = qJD(4) * t44;
t41 = qJDD(1) - g(3);
t69 = t15 * t71 + t89 * t42;
t35 = t38 ^ 2;
t68 = t42 * t35 * t44;
t67 = pkin(2) * t72;
t63 = t75 * t36;
t62 = qJD(2) * (-qJD(3) + t38);
t61 = qJD(3) * (-qJD(2) - t38);
t60 = t42 * t38 * t71;
t59 = t24 + t30 - t86;
t31 = sin(t37);
t32 = cos(t37);
t58 = g(1) * t31 - g(2) * t32;
t56 = -g(1) * (-t26 * pkin(3) + t27 * pkin(6)) - g(2) * (t27 * pkin(3) + t26 * pkin(6));
t46 = qJD(4) ^ 2;
t55 = pkin(6) * t46 - t38 * t66 - t85;
t28 = pkin(6) + t83;
t29 = -pkin(3) - t82;
t53 = qJD(3) * t38 * t83 + t28 * t46 + t29 * t36;
t52 = -qJD(1) * qJD(4) - t15 * t38 + t77 - t9;
t51 = -pkin(6) * qJDD(4) + (t65 - t84) * qJD(4);
t50 = -qJDD(4) * t28 + (t29 * t38 - t67) * qJD(4);
t48 = -t77 + t49;
t17 = qJDD(4) * t44 - t46 * t42;
t16 = qJDD(4) * t42 + t46 * t44;
t11 = t40 * t36 - 0.2e1 * t60;
t10 = t39 * t36 + 0.2e1 * t60;
t4 = -0.2e1 * t76 * t38 * qJD(4) + 0.2e1 * t42 * t79;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, 0, 0, 0, 0, 0, t17, -t16, 0, -t57 * qJD(4) + t2 * t42 + t3 * t44 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t58, g(1) * t32 + g(2) * t31, 0, 0, 0, 0, 0, 0, 0, t36, (t36 * t45 + t43 * t61) * pkin(2) + t59, ((-qJDD(2) - t36) * t43 + t45 * t61) * pkin(2) + t77, 0, (t58 + (t43 ^ 2 + t45 ^ 2) * qJDD(2) * pkin(2)) * pkin(2), t10, t4, t16, t11, t17, 0, t50 * t42 + (-t53 - t89) * t44 + t78, t50 * t44 + (t53 - t24) * t42 + t69, t28 * t63 + t67 * t88 + t48, t8 * t29 + (-t87 * qJD(3) + t58) * pkin(2) + t49 * t28 + t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, t62 * t83 + t59, (t45 * t62 - t70) * pkin(2) + t77, 0, 0, t10, t4, t16, t11, t17, 0, t51 * t42 + (-t55 - t89) * t44 + t78, t51 * t44 + (t55 - t24) * t42 + t69, pkin(6) * t63 - t65 * t88 + t48, -t8 * pkin(3) + t49 * pkin(6) + t87 * t74 + t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, t76 * t35, t42 * t36, t68, t79, qJDD(4), -g(3) * t44 + t33 + (t7 - t80) * qJD(4) + t52 * t42, t73 + (qJD(4) * t14 - t41) * t42 + t52 * t44, 0, 0;];
tau_reg = t1;
