% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
% 
% Output:
% tau_reg [5x15]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PPRRR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:35:46
% EndTime: 2019-12-31 17:35:47
% DurationCPUTime: 0.37s
% Computational Cost: add. (401->90), mult. (677->130), div. (0->0), fcn. (514->10), ass. (0->72)
t37 = qJDD(3) + qJDD(4);
t83 = t37 * pkin(4);
t48 = cos(qJ(3));
t35 = t48 * qJDD(2);
t45 = sin(qJ(3));
t71 = qJD(2) * qJD(3);
t21 = qJDD(3) * pkin(3) - t45 * t71 + t35;
t47 = cos(qJ(4));
t16 = t47 * t21;
t44 = sin(qJ(4));
t29 = qJD(3) * pkin(3) + t48 * qJD(2);
t72 = t45 * qJDD(2);
t59 = -t48 * t71 - t72;
t53 = qJD(4) * t29 - t59;
t75 = qJD(2) * t45;
t70 = qJD(4) * t75;
t89 = -t53 * t44 - t47 * t70 + t16;
t3 = -t83 - t89;
t42 = cos(pkin(8));
t73 = qJ(3) + qJ(4);
t67 = sin(t73);
t68 = cos(t73);
t76 = sin(pkin(8));
t17 = -t42 * t68 - t76 * t67;
t18 = t42 * t67 - t76 * t68;
t65 = g(1) * t18 - g(2) * t17;
t94 = t3 - t65;
t61 = -g(1) * t17 - g(2) * t18 + t44 * t70;
t93 = -t44 * t21 - t47 * t53 + t61;
t38 = qJD(3) + qJD(4);
t36 = t38 ^ 2;
t24 = t44 * t45 - t47 * t48;
t92 = t24 * t38;
t49 = qJD(5) ^ 2;
t91 = pkin(7) * t49 - t83;
t25 = t44 * t48 + t47 * t45;
t19 = t25 * qJD(2);
t33 = t44 * pkin(3) + pkin(7);
t34 = -t47 * pkin(3) - pkin(4);
t77 = pkin(3) * qJD(4);
t90 = t33 * t49 + t34 * t37 + (t44 * t77 - t19) * t38;
t87 = pkin(3) * t37;
t82 = t38 * pkin(4);
t81 = (t44 * t29 + t47 * t75) * t38;
t46 = cos(qJ(5));
t79 = t46 * t37;
t43 = sin(qJ(5));
t39 = t43 ^ 2;
t78 = -t46 ^ 2 + t39;
t74 = t46 * qJD(5);
t41 = qJDD(1) - g(3);
t13 = t47 * t29 - t44 * t75;
t9 = -t13 - t82;
t69 = t94 * t43 + t9 * t74;
t63 = -t24 * t37 - t36 * t25;
t58 = t65 + t81;
t57 = t25 * t49 - t63;
t56 = -t37 * pkin(7) - t9 * t38 + t93;
t55 = -pkin(7) * qJDD(5) + (t13 - t82) * qJD(5);
t54 = 0.2e1 * t92 * qJD(5) - qJDD(5) * t25;
t20 = t24 * qJD(2);
t52 = -qJDD(5) * t33 + (t34 * t38 - t47 * t77 - t20) * qJD(5);
t51 = (-pkin(3) * t38 - t29) * qJD(4) + t59;
t50 = qJD(3) ^ 2;
t27 = qJDD(5) * t46 - t49 * t43;
t26 = qJDD(5) * t43 + t49 * t46;
t23 = t42 * t45 - t76 * t48;
t22 = -t42 * t48 - t76 * t45;
t15 = 0.2e1 * t43 * t38 * t74 + t39 * t37;
t8 = -0.2e1 * t78 * t38 * qJD(5) + 0.2e1 * t43 * t79;
t6 = t9 * qJD(5) * t43;
t1 = [t41, t41, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, t26; 0, -g(1) * t76 + g(2) * t42 + qJDD(2), 0, t48 * qJDD(3) - t50 * t45, -qJDD(3) * t45 - t50 * t48, 0, t63, -t25 * t37 + t38 * t92, 0, 0, 0, 0, 0, t43 * t54 - t46 * t57, t43 * t57 + t46 * t54; 0, 0, qJDD(3), g(1) * t23 - g(2) * t22 + t35, -g(1) * t22 - g(2) * t23 - t72, t37, t19 * t38 + t16 + (-t70 + t87) * t47 + t51 * t44 + t65, -t20 * t38 + (-t21 - t87) * t44 + t51 * t47 + t61, t15, t8, t26, t27, 0, t6 + t52 * t43 + (-t90 - t94) * t46, t90 * t43 + t52 * t46 + t69; 0, 0, 0, 0, 0, t37, t58 + t89, t13 * t38 + t93, t15, t8, t26, t27, 0, t6 + t55 * t43 + (-t3 + t58 - t91) * t46, t55 * t46 + (-t81 + t91) * t43 + t69; 0, 0, 0, 0, 0, 0, 0, 0, -t43 * t36 * t46, t78 * t36, t43 * t37, t79, qJDD(5), -t41 * t46 + t56 * t43, t41 * t43 + t56 * t46;];
tau_reg = t1;
