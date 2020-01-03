% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPPRR5
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPPRR5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR5_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR5_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR5_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:56:38
% EndTime: 2019-12-31 17:56:40
% DurationCPUTime: 0.42s
% Computational Cost: add. (1648->115), mult. (2505->130), div. (0->0), fcn. (1196->8), ass. (0->81)
t104 = -pkin(2) - pkin(3);
t65 = -qJD(1) + qJD(4);
t63 = t65 ^ 2;
t74 = sin(qJ(5));
t77 = cos(qJ(5));
t51 = t77 * t63 * t74;
t45 = qJDD(5) + t51;
t103 = t74 * t45;
t64 = qJDD(1) - qJDD(4);
t102 = t74 * t64;
t46 = qJDD(5) - t51;
t101 = t77 * t46;
t81 = qJD(1) ^ 2;
t76 = sin(qJ(1));
t79 = cos(qJ(1));
t88 = g(1) * t79 + g(2) * t76;
t44 = -pkin(1) * t81 - t88;
t72 = sin(pkin(8));
t73 = cos(pkin(8));
t96 = t76 * g(1) - g(2) * t79;
t86 = qJDD(1) * pkin(1) + t96;
t100 = t73 * t44 + t72 * t86;
t61 = 2 * qJD(3) * qJD(1);
t66 = qJDD(1) * qJ(3);
t92 = t61 + t66 + t100;
t20 = t104 * t81 + t92;
t75 = sin(qJ(4));
t78 = cos(qJ(4));
t71 = qJDD(1) * pkin(2);
t93 = -t44 * t72 + t73 * t86;
t87 = -qJDD(3) + t93;
t22 = -t81 * qJ(3) - t71 - t87;
t82 = -qJDD(1) * pkin(3) + t22;
t12 = t78 * t20 + t75 * t82;
t99 = qJD(5) * t65;
t98 = t74 * t99;
t97 = t77 * t99;
t10 = -pkin(4) * t63 - pkin(7) * t64 + t12;
t70 = g(3) - qJDD(2);
t7 = t10 * t74 - t77 * t70;
t8 = t10 * t77 + t70 * t74;
t3 = t74 * t7 + t77 * t8;
t95 = pkin(1) * t72 + qJ(3);
t94 = t75 * t20 - t78 * t82;
t38 = -t78 * t63 + t64 * t75;
t91 = pkin(1) * t73 - t104;
t90 = -pkin(1) * (qJDD(1) * t72 + t73 * t81) - t100;
t9 = t64 * pkin(4) - t63 * pkin(7) + t94;
t89 = -pkin(4) * t9 + pkin(7) * t3;
t39 = t63 * t75 + t78 * t64;
t69 = t77 ^ 2;
t56 = t69 * t63;
t80 = qJD(5) ^ 2;
t50 = -t56 - t80;
t30 = t50 * t77 - t103;
t53 = t77 * t64;
t33 = -t53 - 0.2e1 * t98;
t85 = pkin(4) * t33 + pkin(7) * t30 - t77 * t9;
t68 = t74 ^ 2;
t55 = t68 * t63;
t49 = -t55 - t80;
t31 = -t49 * t74 - t101;
t32 = 0.2e1 * t97 - t102;
t84 = pkin(4) * t32 - pkin(7) * t31 - t74 * t9;
t37 = (-t68 - t69) * t64;
t42 = t55 + t56;
t83 = pkin(4) * t42 + pkin(7) * t37 + t3;
t48 = t74 * t97;
t43 = pkin(1) * (qJDD(1) * t73 - t72 * t81);
t29 = t103 + t77 * (-t55 + t80);
t28 = t74 * (t56 - t80) + t101;
t27 = -t48 + t77 * (-t53 - t98);
t26 = -t48 - t74 * (t97 - t102);
t24 = t37 * t75 + t42 * t78;
t23 = t32 * t77 + t33 * t74;
t21 = -pkin(2) * t81 + t92;
t14 = t31 * t75 - t32 * t78;
t13 = t30 * t75 + t33 * t78;
t4 = t12 * t75 - t78 * t94;
t1 = t3 * t75 - t78 * t9;
t2 = [0, 0, 0, 0, 0, qJDD(1), t96, t88, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t43 + t93, t90, 0, pkin(1) * (t100 * t72 + t73 * t93), 0, 0, 0, qJDD(1), 0, 0, t43 + 0.2e1 * t71 + t87, 0, t61 + 0.2e1 * t66 - t90, pkin(1) * (t21 * t72 - t22 * t73) + qJ(3) * t21 - pkin(2) * t22, 0, 0, 0, 0, 0, t64, t38 * t95 + t39 * t91 + t94, -t38 * t91 + t39 * t95 + t12, 0, t95 * (t12 * t78 + t75 * t94) - t91 * t4, t26, -t23, -t29, -t27, -t28, 0, t95 * (t30 * t78 - t33 * t75) - t91 * t13 - t85, t95 * (t31 * t78 + t32 * t75) - t91 * t14 + t84, t95 * (t37 * t78 - t42 * t75) - t91 * t24 - t83, t95 * (t3 * t78 + t75 * t9) - t91 * t1 - t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, 0, 0, 0, 0, 0, 0, -t45 * t77 - t50 * t74, t46 * t74 - t49 * t77, 0, t7 * t77 - t74 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t81, t22, 0, 0, 0, 0, 0, 0, -t39, t38, 0, t4, 0, 0, 0, 0, 0, 0, t13, t14, t24, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, -t94, -t12, 0, 0, -t26, t23, t29, t27, t28, 0, t85, -t84, t83, t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, t55 - t56, -t102, t51, -t53, qJDD(5), -t7, -t8, 0, 0;];
tauJ_reg = t2;
