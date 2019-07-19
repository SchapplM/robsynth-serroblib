% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RRPR2
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% tau_reg [4x12]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRPR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 18:16:35
% EndTime: 2019-07-18 18:16:36
% DurationCPUTime: 0.30s
% Computational Cost: add. (442->108), mult. (484->134), div. (0->0), fcn. (246->8), ass. (0->61)
t41 = qJ(1) + qJ(2);
t37 = sin(t41);
t38 = cos(t41);
t80 = -g(1) * t38 - g(2) * t37;
t48 = -pkin(2) - pkin(3);
t39 = qJDD(1) + qJDD(2);
t79 = pkin(2) * t39;
t40 = qJD(1) + qJD(2);
t43 = sin(qJ(2));
t74 = pkin(1) * qJD(1);
t66 = t43 * t74;
t14 = qJ(3) * t40 + t66;
t42 = sin(qJ(4));
t78 = t42 * t14;
t46 = cos(qJ(2));
t62 = qJD(2) * t74;
t73 = pkin(1) * qJDD(1);
t77 = t43 * t73 + t46 * t62;
t76 = t38 * pkin(2) + t37 * qJ(3);
t75 = -t43 * t62 + t46 * t73;
t72 = qJD(2) * t43;
t71 = qJD(2) * t46;
t70 = qJD(4) * t42;
t45 = cos(qJ(4));
t69 = qJD(4) * t45;
t68 = -qJD(4) + t40;
t67 = t68 * t74;
t65 = t46 * t74;
t64 = pkin(1) * t72;
t63 = t40 * t72;
t30 = -pkin(1) * t46 - pkin(2);
t61 = -qJDD(3) + t75;
t60 = -pkin(2) * t37 + t38 * qJ(3);
t59 = t68 ^ 2;
t58 = t77 + t80;
t33 = t39 * qJ(3);
t34 = t40 * qJD(3);
t5 = t33 + t34 + t77;
t57 = qJD(3) - t65;
t56 = t45 * qJ(3) + t42 * t48;
t55 = t42 * qJ(3) - t45 * t48;
t54 = g(1) * t37 - g(2) * t38 + t75;
t53 = -qJDD(3) + t54;
t10 = -t37 * t42 - t38 * t45;
t11 = -t37 * t45 + t38 * t42;
t4 = t48 * t39 - t61;
t9 = t48 * t40 + t57;
t52 = g(1) * t10 + g(2) * t11 + t42 * t4 + t45 * t5 + t9 * t69;
t51 = g(1) * t11 - g(2) * t10 - t14 * t69 - t42 * t5 - t9 * t70;
t50 = t40 * t65 - t58;
t49 = t45 * t4 + t51;
t47 = cos(qJ(1));
t44 = sin(qJ(1));
t35 = -qJDD(4) + t39;
t25 = pkin(1) * t43 + qJ(3);
t22 = -pkin(3) + t30;
t19 = pkin(1) * t71 + qJD(3);
t15 = t40 * t66;
t13 = -pkin(2) * t40 + t57;
t8 = -t61 - t79;
t1 = [qJDD(1), g(1) * t44 - g(2) * t47, g(1) * t47 + g(2) * t44, t39, (t39 * t46 - t63) * pkin(1) + t54, (-t39 * t43 - t40 * t71) * pkin(1) - t58, -pkin(1) * t63 + (pkin(2) - t30) * t39 + t53, t19 * t40 + t25 * t39 + t5 + t80, t5 * t25 + t14 * t19 + t8 * t30 + t13 * t64 - g(1) * (-pkin(1) * t44 + t60) - g(2) * (pkin(1) * t47 + t76), t35, (-(-qJD(4) * t22 - t19) * t68 + t25 * t35) * t42 + (-(-qJD(4) * t25 + t64) * t68 - t22 * t35 - t4) * t45 - t51, (t19 * t45 + t42 * t64) * t68 + (t22 * t42 + t25 * t45) * t35 + ((t22 * t45 - t25 * t42) * t68 - t78) * qJD(4) + t52; 0, 0, 0, t39, t15 + t54, t50, t15 + t53 + 0.2e1 * t79, 0.2e1 * t33 + 0.2e1 * t34 - t50, t5 * qJ(3) + t14 * qJD(3) - t8 * pkin(2) - g(1) * t60 - g(2) * t76 + (-t13 * t43 - t14 * t46) * t74, t35, -(-t42 * qJD(3) - t56 * qJD(4)) * t68 + t55 * t35 + (-t42 * t46 + t43 * t45) * t67 - t49, t45 * qJD(3) * t68 + t56 * t35 + (-t55 * t68 - t78) * qJD(4) - (t42 * t43 + t45 * t46) * t67 + t52; 0, 0, 0, 0, 0, 0, -t39, -t40 ^ 2, -t14 * t40 - t53 - t79, 0, -t45 * t35 - t42 * t59, t42 * t35 - t45 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, (-t14 * t45 - t42 * t9) * t68 + t49, t14 * t70 - (t45 * t9 - t78) * t68 - t52;];
tau_reg  = t1;
