% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RRPP1
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
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% 
% Output:
% tau_reg [4x10]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRPP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:33:10
% EndTime: 2019-03-08 18:33:10
% DurationCPUTime: 0.23s
% Computational Cost: add. (320->92), mult. (487->107), div. (0->0), fcn. (265->10), ass. (0->67)
t57 = cos(qJ(2));
t76 = pkin(1) * qJD(2);
t69 = qJD(1) * t76;
t55 = sin(qJ(2));
t75 = qJDD(1) * t55;
t88 = pkin(1) * t75 + t57 * t69;
t52 = qJ(1) + qJ(2);
t46 = sin(t52);
t47 = cos(t52);
t87 = g(1) * t46 - g(2) * t47;
t86 = pkin(2) * t46;
t56 = sin(qJ(1));
t84 = t56 * pkin(1);
t83 = t57 * pkin(1);
t53 = sin(pkin(6));
t82 = t53 * t55;
t54 = cos(pkin(6));
t81 = t54 * t55;
t80 = t54 * t57;
t50 = qJD(1) + qJD(2);
t77 = pkin(1) * qJD(1);
t71 = t57 * t77;
t19 = t50 * pkin(2) + t71;
t72 = t55 * t77;
t25 = t54 * t72;
t8 = t53 * t19 + t25;
t40 = pkin(2) + t83;
t79 = pkin(1) * t81 + t53 * t40;
t78 = g(1) * t47 + g(2) * t46;
t74 = pkin(1) * t82;
t41 = qJDD(1) * t83;
t49 = qJDD(1) + qJDD(2);
t12 = t49 * pkin(2) - t55 * t69 + t41;
t4 = t53 * t12 + t88 * t54;
t3 = t54 * t12 - t88 * t53;
t45 = pkin(6) + t52;
t33 = sin(t45);
t34 = cos(t45);
t39 = pkin(2) * t47;
t73 = t34 * pkin(3) + t33 * qJ(4) + t39;
t68 = t49 * qJ(4) + t4;
t67 = qJD(1) * (-qJD(2) + t50);
t66 = qJD(2) * (-qJD(1) - t50);
t64 = t41 + t87;
t2 = -t49 * pkin(3) + qJDD(4) - t3;
t63 = -qJD(2) * t74 + t76 * t80;
t62 = -t33 * pkin(3) + t34 * qJ(4) - t86;
t61 = t54 * t40 - t74;
t7 = t54 * t19 - t53 * t72;
t60 = -g(1) * t34 - g(2) * t33 + t68;
t59 = -g(1) * t33 + g(2) * t34 + t2;
t58 = cos(qJ(1));
t51 = qJDD(3) - g(3);
t48 = t58 * pkin(1);
t43 = t50 * qJD(4);
t35 = -t54 * pkin(2) - pkin(3);
t32 = t53 * pkin(2) + qJ(4);
t17 = (t80 - t82) * t77;
t16 = (t53 * t57 + t81) * t76;
t15 = t53 * t71 + t25;
t14 = -pkin(3) - t61;
t13 = qJ(4) + t79;
t11 = qJD(4) + t63;
t6 = t50 * qJ(4) + t8;
t5 = -t50 * pkin(3) + qJD(4) - t7;
t1 = t43 + t68;
t9 = [qJDD(1), g(1) * t56 - g(2) * t58, g(1) * t58 + g(2) * t56, t49 (t49 * t57 + t55 * t66) * pkin(1) + t64 ((-qJDD(1) - t49) * t55 + t57 * t66) * pkin(1) + t78, t4 * t79 + t8 * t63 + t3 * t61 - t7 * t16 - g(1) * (-t84 - t86) - g(2) * (t39 + t48) -t14 * t49 - t16 * t50 - t59, t11 * t50 + t13 * t49 + t43 + t60, t1 * t13 + t6 * t11 + t2 * t14 + t5 * t16 - g(1) * (t62 - t84) - g(2) * (t48 + t73); 0, 0, 0, t49, t55 * pkin(1) * t67 + t64 (t57 * t67 - t75) * pkin(1) + t78, t7 * t15 - t8 * t17 + (t3 * t54 + t4 * t53 + t87) * pkin(2), t15 * t50 - t35 * t49 - t59, -t17 * t50 + t32 * t49 + 0.2e1 * t43 + t60, t1 * t32 + t2 * t35 - t5 * t15 - g(1) * t62 - g(2) * t73 + (qJD(4) - t17) * t6; 0, 0, 0, 0, 0, 0, t51, 0, 0, t51; 0, 0, 0, 0, 0, 0, 0, -t49, -t50 ^ 2, -t6 * t50 + t59;];
tau_reg  = t9;
