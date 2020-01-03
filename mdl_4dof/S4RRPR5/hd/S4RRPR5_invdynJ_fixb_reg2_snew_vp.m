% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
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
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RRPR5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR5_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_invdynJ_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:03:33
% EndTime: 2019-12-31 17:03:34
% DurationCPUTime: 0.30s
% Computational Cost: add. (702->66), mult. (933->89), div. (0->0), fcn. (482->6), ass. (0->60)
t89 = pkin(2) + pkin(6);
t52 = (qJD(1) + qJD(2));
t50 = t52 ^ 2;
t51 = qJDD(1) + qJDD(2);
t56 = sin(qJ(2));
t59 = cos(qJ(2));
t88 = pkin(1) * (t56 * t50 - t59 * t51);
t87 = pkin(1) * (t59 * t50 + t56 * t51);
t55 = sin(qJ(4));
t58 = cos(qJ(4));
t57 = sin(qJ(1));
t60 = cos(qJ(1));
t73 = t57 * g(1) - t60 * g(2);
t36 = qJDD(1) * pkin(1) + t73;
t70 = t60 * g(1) + t57 * g(2);
t37 = -qJD(1) ^ 2 * pkin(1) - t70;
t15 = t59 * t36 - t56 * t37;
t86 = t51 * pkin(2);
t66 = qJDD(3) - t15 - t86;
t13 = -t50 * qJ(3) + t66;
t62 = -t51 * pkin(6) + t13;
t4 = t55 * g(3) + t58 * t62;
t5 = t58 * g(3) - t55 * t62;
t2 = t58 * t4 - t55 * t5;
t53 = t55 ^ 2;
t84 = t53 * t50;
t54 = t58 ^ 2;
t83 = t54 * t50;
t75 = t55 * t50 * t58;
t82 = t55 * (qJDD(4) + t75);
t81 = t55 * t51;
t80 = t58 * (qJDD(4) - t75);
t16 = t56 * t36 + t59 * t37;
t74 = (2 * qJD(3) * t52) + t16;
t77 = t51 * qJ(3);
t11 = -t50 * pkin(2) + t74 + t77;
t79 = -pkin(2) * t13 + qJ(3) * t11;
t78 = t53 + t54;
t76 = qJD(4) * t52;
t72 = 0.2e1 * t77 + t74;
t8 = -t50 * pkin(6) + t11;
t71 = qJ(3) * t8 - t89 * t2;
t40 = t58 * t51;
t30 = -0.2e1 * t55 * t76 + t40;
t61 = qJD(4) ^ 2;
t22 = -t82 + t58 * (-t61 - t83);
t67 = qJ(3) * t30 - t89 * t22 + t58 * t8;
t29 = 0.2e1 * t58 * t76 + t81;
t21 = t55 * (-t61 - t84) + t80;
t65 = qJ(3) * t29 - t89 * t21 + t55 * t8;
t64 = -t86 + t66;
t34 = t78 * t51;
t35 = t78 * t50;
t63 = -qJ(3) * t35 + t89 * t34 - t2;
t24 = t80 - t55 * (t61 - t83);
t23 = t58 * (-t61 + t84) - t82;
t18 = t30 * t58;
t17 = t29 * t55;
t14 = -t58 * t29 - t55 * t30;
t1 = [0, 0, 0, 0, 0, qJDD(1), t73, t70, 0, 0, 0, 0, 0, 0, 0, t51, t15 - t88, -t16 - t87, 0, pkin(1) * (t59 * t15 + t56 * t16), t51, 0, 0, 0, 0, 0, 0, t64 + t88, t72 + t87, pkin(1) * (t56 * t11 - t59 * t13) + t79, t18, t14, t24, t17, t23, 0, pkin(1) * (-t59 * t21 + t56 * t29) + t65, pkin(1) * (-t59 * t22 + t56 * t30) + t67, pkin(1) * (t59 * t34 - t56 * t35) + t63, pkin(1) * (-t59 * t2 + t56 * t8) + t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, t15, -t16, 0, 0, t51, 0, 0, 0, 0, 0, 0, t64, t72, t79, t18, t14, t24, t17, t23, 0, t65, t67, t63, t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, -t50, t13, 0, 0, 0, 0, 0, 0, t21, t22, -t34, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, (-t53 + t54) * t50, t40, -t75, -t81, qJDD(4), t4, t5, 0, 0;];
tauJ_reg = t1;
