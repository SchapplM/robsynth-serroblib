% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4PRRR5
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
% tau_reg [4x14]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PRRR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:33:47
% EndTime: 2019-12-31 16:33:48
% DurationCPUTime: 0.36s
% Computational Cost: add. (350->85), mult. (625->126), div. (0->0), fcn. (458->10), ass. (0->69)
t41 = qJ(2) + qJ(3);
t34 = sin(t41);
t42 = sin(pkin(7));
t43 = cos(pkin(7));
t66 = g(1) * t43 + g(2) * t42;
t98 = t66 * t34;
t49 = cos(qJ(2));
t33 = t49 * qJDD(1);
t46 = sin(qJ(2));
t73 = qJD(1) * qJD(2);
t16 = qJDD(2) * pkin(2) - t46 * t73 + t33;
t45 = sin(qJ(3));
t48 = cos(qJ(3));
t29 = qJD(2) * pkin(2) + t49 * qJD(1);
t60 = -t46 * qJDD(1) - t49 * t73;
t55 = qJD(3) * t29 - t60;
t35 = cos(t41);
t76 = qJD(1) * t46;
t69 = qJD(3) * t76;
t68 = g(3) * t34 + t66 * t35 + t45 * t69;
t97 = -t45 * t16 - t55 * t48 + t68;
t37 = qJDD(2) + qJDD(3);
t83 = t37 * pkin(3);
t84 = g(3) * t35;
t13 = t48 * t16;
t89 = t55 * t45 + t48 * t69 - t13;
t96 = -t83 + t89 + t84;
t38 = qJD(2) + qJD(3);
t36 = t38 ^ 2;
t17 = t45 * t46 - t48 * t49;
t93 = t17 * t38;
t91 = -t84 + t98;
t18 = t45 * t49 + t48 * t46;
t14 = t18 * qJD(1);
t31 = t45 * pkin(2) + pkin(6);
t32 = -t48 * pkin(2) - pkin(3);
t50 = qJD(4) ^ 2;
t77 = pkin(2) * qJD(3);
t90 = t31 * t50 + t32 * t37 + (t45 * t77 - t14) * t38;
t87 = pkin(2) * t37;
t82 = t38 * pkin(3);
t81 = (t45 * t29 + t48 * t76) * t38;
t47 = cos(qJ(4));
t79 = t47 * t37;
t44 = sin(qJ(4));
t39 = t44 ^ 2;
t78 = -t47 ^ 2 + t39;
t75 = qJD(4) * t47;
t74 = qJDD(1) - g(3);
t10 = t48 * t29 - t45 * t76;
t9 = -t10 - t82;
t72 = t96 * t44 + t9 * t75;
t71 = t9 * qJD(4) * t44 + t47 * t98;
t65 = g(1) * t42 - g(2) * t43;
t64 = -t17 * t37 - t36 * t18;
t59 = pkin(6) * t50 - t81 - t83;
t58 = t18 * t50 - t64;
t57 = -pkin(6) * qJDD(4) + (t10 - t82) * qJD(4);
t56 = 0.2e1 * t93 * qJD(4) - qJDD(4) * t18;
t54 = -t37 * pkin(6) - t9 * t38 + t97;
t15 = t17 * qJD(1);
t53 = -qJDD(4) * t31 + (t32 * t38 - t48 * t77 - t15) * qJD(4);
t52 = (-pkin(2) * t38 - t29) * qJD(3) + t60;
t51 = qJD(2) ^ 2;
t22 = qJDD(4) * t47 - t50 * t44;
t21 = qJDD(4) * t44 + t50 * t47;
t12 = 0.2e1 * t44 * t38 * t75 + t39 * t37;
t8 = -0.2e1 * t78 * t38 * qJD(4) + 0.2e1 * t44 * t79;
t1 = [t74, 0, t49 * qJDD(2) - t51 * t46, -qJDD(2) * t46 - t51 * t49, 0, t64, -t18 * t37 + t38 * t93, 0, 0, 0, 0, 0, t56 * t44 - t58 * t47, t58 * t44 + t56 * t47; 0, qJDD(2), -g(3) * t49 + t66 * t46 + t33, -t74 * t46 + t66 * t49, t37, t14 * t38 + t13 + (-t69 + t87) * t48 + t52 * t45 + t91, -t15 * t38 + (-t16 - t87) * t45 + t52 * t48 + t68, t12, t8, t21, t22, 0, t53 * t44 + (-t96 - t90) * t47 + t71, t53 * t47 + (-t98 + t90) * t44 + t72; 0, 0, 0, 0, t37, t81 - t89 + t91, t10 * t38 + t97, t12, t8, t21, t22, 0, t57 * t44 + (-t59 - t96) * t47 + t71, t57 * t47 + (-t98 + t59) * t44 + t72; 0, 0, 0, 0, 0, 0, 0, -t44 * t36 * t47, t78 * t36, t44 * t37, t79, qJDD(4), t54 * t44 - t65 * t47, t65 * t44 + t54 * t47;];
tau_reg = t1;
