% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RPRP2
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
%   pkin=[a2,a3,a4,d1,d3]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPRP2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP2_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_invdynJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:30:55
% EndTime: 2019-03-08 18:30:55
% DurationCPUTime: 0.26s
% Computational Cost: add. (413->87), mult. (606->97), div. (0->0), fcn. (305->4), ass. (0->54)
t51 = -pkin(1) - pkin(2);
t67 = qJ(2) * qJD(1);
t83 = -qJD(3) * t67 + t51 * qJDD(1) + qJDD(2);
t65 = (qJD(1) * qJD(2));
t66 = (qJ(2) * qJDD(1));
t82 = t65 + t66;
t47 = sin(qJ(3));
t49 = cos(qJ(3));
t26 = -t47 * qJ(2) + t49 * t51;
t50 = cos(qJ(1));
t48 = sin(qJ(1));
t75 = t48 * t47;
t20 = -t50 * t49 - t75;
t74 = t50 * t47;
t21 = -t48 * t49 + t74;
t81 = g(1) * t21 - g(2) * t20;
t14 = t49 * qJD(2) + t26 * qJD(3);
t29 = t51 * qJD(1) + qJD(2);
t17 = t47 * t29 + t49 * t67;
t27 = t49 * qJ(2) + t47 * t51;
t68 = qJD(3) * t49;
t6 = t29 * t68 + t83 * t47 + t82 * t49;
t80 = t17 * t14 + t6 * t27;
t45 = qJDD(1) - qJDD(3);
t77 = t45 * pkin(3);
t76 = t17 * t49;
t73 = t50 * pkin(1) + t48 * qJ(2);
t72 = g(1) * t48 - g(2) * t50;
t70 = pkin(1) * qJDD(1);
t69 = qJD(3) * t47;
t64 = qJD(1) - qJD(3);
t63 = 2 * t65;
t61 = t17 * t68 + t6 * t47 - t72;
t60 = t64 ^ 2;
t59 = qJDD(2) - t70;
t58 = t29 * t69 + t82 * t47 - t83 * t49;
t57 = g(1) * t50 + g(2) * t48;
t16 = t49 * t29 - t47 * t67;
t56 = -t58 + t81;
t15 = -t47 * qJD(2) - t27 * qJD(3);
t55 = -t15 * t64 - t56;
t54 = -t17 * t64 + t56;
t53 = -g(1) * t20 - g(2) * t21 - t6;
t52 = qJD(1) ^ 2;
t41 = t50 * qJ(2);
t38 = t49 * pkin(3) + pkin(2);
t22 = -pkin(3) + t26;
t11 = -pkin(3) * t64 + t16;
t10 = t47 * t45 - t49 * t60;
t9 = -t49 * t45 - t47 * t60;
t4 = -t58 - t77;
t2 = -t16 * t64 + t53;
t1 = t14 * t64 + t27 * t45 - t53;
t3 = [0, 0, 0, 0, 0, qJDD(1), t72, t57, 0, 0, 0, 0, 0, qJDD(1), 0, 0, -qJDD(2) + 0.2e1 * t70 + t72, 0, -t57 + t63 + (2 * t66) -t59 * pkin(1) - g(1) * (-t48 * pkin(1) + t41) - g(2) * t73 + (t63 + t66) * qJ(2), 0, 0, 0, 0, 0, t45, -t26 * t45 + t55, t1, 0, -t58 * t26 + t16 * t15 - g(1) * (t51 * t48 + t41) - g(2) * (t50 * pkin(2) + t73) + t80, 0, 0, 0, 0, 0, t45 (pkin(3) - t22) * t45 + t55, t1, 0, t4 * t22 + t11 * t15 - g(1) * (pkin(3) * t74 + t41 + (-pkin(1) - t38) * t48) - g(2) * (pkin(3) * t75 + t50 * t38 + t73) + t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t52, -t52 * qJ(2) + t59 - t72, 0, 0, 0, 0, 0, 0, t9, t10, 0, -t16 * t69 - t58 * t49 + (t16 * t47 - t76) * qJD(1) + t61, 0, 0, 0, 0, 0, 0, t9, t10, 0, -t11 * t69 + t4 * t49 + (t11 * t47 - t76) * qJD(1) + t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, t54, t2, 0, 0, 0, 0, 0, 0, 0, -t45, t54 - 0.2e1 * t77, t2, 0 (-t16 + t11) * t17 + (t4 + t81) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(4) + g(3);];
tau_reg  = t3;
