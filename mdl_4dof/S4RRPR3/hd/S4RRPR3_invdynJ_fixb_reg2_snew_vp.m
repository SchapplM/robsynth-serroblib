% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4RRPR3
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RRPR3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR3_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:34
% EndTime: 2019-12-31 17:01:35
% DurationCPUTime: 0.35s
% Computational Cost: add. (1133->87), mult. (1634->129), div. (0->0), fcn. (944->8), ass. (0->68)
t76 = sin(qJ(1));
t79 = cos(qJ(1));
t88 = t76 * g(1) - t79 * g(2);
t52 = qJDD(1) * pkin(1) + t88;
t83 = t79 * g(1) + t76 * g(2);
t53 = -qJD(1) ^ 2 * pkin(1) - t83;
t75 = sin(qJ(2));
t78 = cos(qJ(2));
t30 = t75 * t52 + t78 * t53;
t69 = qJD(1) + qJD(2);
t67 = t69 ^ 2;
t28 = -t67 * pkin(2) + t30;
t72 = sin(pkin(7));
t73 = cos(pkin(7));
t29 = t78 * t52 - t75 * t53;
t68 = qJDD(1) + qJDD(2);
t81 = t68 * pkin(2) + t29;
t17 = t73 * t28 + t72 * t81;
t15 = -t67 * pkin(3) + t68 * pkin(6) + t17;
t74 = sin(qJ(4));
t77 = cos(qJ(4));
t91 = -g(3) + qJDD(3);
t10 = t77 * t15 + t74 * t91;
t9 = t74 * t15 - t77 * t91;
t4 = t77 * t10 + t74 * t9;
t58 = t77 * t67 * t74;
t54 = qJDD(4) + t58;
t95 = t74 * t54;
t94 = t74 * t68;
t55 = qJDD(4) - t58;
t93 = t77 * t55;
t92 = qJD(4) * t69;
t16 = -t72 * t28 + t73 * t81;
t14 = -t68 * pkin(3) - t67 * pkin(6) - t16;
t2 = -t73 * t14 + t72 * t4;
t90 = pkin(2) * t2 - pkin(3) * t14 + pkin(6) * t4;
t48 = -t73 * t67 - t72 * t68;
t89 = pkin(2) * t48 - t17;
t70 = t74 ^ 2;
t61 = t70 * t67;
t80 = qJD(4) ^ 2;
t56 = -t61 - t80;
t38 = -t74 * t56 - t93;
t44 = 0.2e1 * t77 * t92 + t94;
t21 = t72 * t38 - t73 * t44;
t87 = pkin(2) * t21 - pkin(3) * t44 + pkin(6) * t38 + t74 * t14;
t71 = t77 ^ 2;
t62 = t71 * t67;
t57 = -t62 - t80;
t37 = t77 * t57 - t95;
t60 = t77 * t68;
t45 = -0.2e1 * t74 * t92 + t60;
t20 = t72 * t37 + t73 * t45;
t86 = pkin(2) * t20 + pkin(3) * t45 + pkin(6) * t37 - t77 * t14;
t50 = (t70 + t71) * t68;
t51 = t61 + t62;
t27 = t72 * t50 + t73 * t51;
t85 = pkin(2) * t27 + pkin(3) * t51 + pkin(6) * t50 + t4;
t82 = t72 * t67 - t73 * t68;
t84 = -pkin(2) * t82 + t16;
t36 = t95 + t77 * (-t61 + t80);
t35 = t74 * (t62 - t80) + t93;
t32 = t44 * t74;
t31 = t45 * t77;
t25 = t77 * t44 + t74 * t45;
t6 = t73 * t16 + t72 * t17;
t5 = pkin(2) * t6;
t1 = [0, 0, 0, 0, 0, qJDD(1), t88, t83, 0, 0, 0, 0, 0, 0, 0, t68, pkin(1) * (-t75 * t67 + t78 * t68) + t29, pkin(1) * (-t78 * t67 - t75 * t68) - t30, 0, pkin(1) * (t78 * t29 + t75 * t30), 0, 0, 0, 0, 0, t68, pkin(1) * (t75 * t48 - t78 * t82) + t84, pkin(1) * (t78 * t48 + t75 * t82) + t89, 0, pkin(1) * (t75 * (-t72 * t16 + t73 * t17) + t78 * t6) + t5, t32, t25, t36, t31, t35, 0, pkin(1) * (t75 * (t73 * t37 - t72 * t45) + t78 * t20) + t86, pkin(1) * (t75 * (t73 * t38 + t72 * t44) + t78 * t21) + t87, pkin(1) * (t75 * (t73 * t50 - t72 * t51) + t78 * t27) + t85, pkin(1) * (t75 * (t72 * t14 + t73 * t4) + t78 * t2) + t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, t29, -t30, 0, 0, 0, 0, 0, 0, 0, t68, t84, t89, 0, t5, t32, t25, t36, t31, t35, 0, t86, t87, t85, t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, 0, 0, 0, 0, 0, 0, t77 * t54 + t74 * t57, -t74 * t55 + t77 * t56, 0, t74 * t10 - t77 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, t61 - t62, t94, t58, t60, qJDD(4), -t9, -t10, 0, 0;];
tauJ_reg = t1;
