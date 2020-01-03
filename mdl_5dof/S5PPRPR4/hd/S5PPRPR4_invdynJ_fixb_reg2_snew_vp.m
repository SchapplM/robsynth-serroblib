% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PPRPR4
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PPRPR4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR4_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:32:24
% EndTime: 2019-12-31 17:32:26
% DurationCPUTime: 0.66s
% Computational Cost: add. (1070->125), mult. (2296->189), div. (0->0), fcn. (1594->8), ass. (0->83)
t93 = sin(pkin(7));
t94 = cos(pkin(7));
t55 = -t94 * g(1) - t93 * g(2);
t74 = sin(qJ(3));
t76 = cos(qJ(3));
t82 = -t93 * g(1) + t94 * g(2) + qJDD(2);
t38 = t76 * t55 + t74 * t82;
t78 = qJD(3) ^ 2;
t108 = -t78 * pkin(3) + qJDD(3) * qJ(4) + (2 * qJD(3) * qJD(4)) + t38;
t71 = sin(pkin(8));
t72 = cos(pkin(8));
t73 = sin(qJ(5));
t75 = cos(qJ(5));
t47 = (t71 * t73 - t72 * t75) * qJD(3);
t84 = t71 * t75 + t72 * t73;
t49 = t84 * qJD(3);
t36 = t49 * t47;
t102 = qJDD(5) - t36;
t107 = t102 * t73;
t106 = t102 * t75;
t105 = t72 * t78;
t67 = t72 ^ 2;
t79 = t71 ^ 2;
t104 = t67 + t79;
t69 = g(3) - qJDD(1);
t61 = t72 * t69;
t103 = t61 + (pkin(4) * t105 - pkin(6) * qJDD(3) - t108) * t71;
t54 = t104 * t78;
t44 = t47 ^ 2;
t45 = t49 ^ 2;
t100 = t67 * t78;
t23 = t108 * t72 + t71 * t69;
t89 = t72 * qJDD(3);
t16 = -pkin(4) * t100 + pkin(6) * t89 + t23;
t7 = -t75 * t103 + t73 * t16;
t8 = t103 * t73 + t75 * t16;
t2 = -t75 * t7 + t73 * t8;
t101 = t71 * t2;
t37 = -t74 * t55 + t76 * t82;
t70 = qJDD(3) * pkin(3);
t26 = -t78 * qJ(4) + qJDD(4) - t37 - t70;
t20 = -pkin(4) * t89 + t26 + (-t79 * t78 - t100) * pkin(6);
t99 = t73 * t20;
t29 = qJDD(5) + t36;
t98 = t73 * t29;
t97 = t75 * t20;
t96 = t75 * t29;
t92 = t47 * qJD(5);
t91 = t49 * qJD(5);
t90 = t71 * qJDD(3);
t88 = t76 * qJDD(3);
t3 = t73 * t7 + t75 * t8;
t22 = t108 * t71 - t61;
t10 = t71 * t22 + t72 * t23;
t85 = -t26 + t70;
t21 = -t73 * t90 + t75 * t89;
t46 = t84 * qJDD(3);
t77 = qJD(5) ^ 2;
t66 = t67 * qJDD(3);
t65 = t79 * qJDD(3);
t53 = t66 + t65;
t52 = t104 * t105;
t51 = t71 * t54;
t41 = -t45 - t77;
t40 = -t45 + t77;
t39 = t44 - t77;
t35 = t46 - t92;
t34 = t46 - 0.2e1 * t92;
t33 = t21 - t91;
t32 = -t21 + 0.2e1 * t91;
t27 = -t77 - t44;
t24 = -t44 - t45;
t18 = -t73 * t41 - t96;
t17 = t75 * t41 - t98;
t14 = t75 * t21 + t73 * t46;
t13 = t73 * t21 - t75 * t46;
t12 = t75 * t27 - t107;
t11 = t73 * t27 + t106;
t9 = -t71 * t17 + t72 * t18;
t5 = -t71 * t13 + t72 * t14;
t4 = -t71 * t11 + t72 * t12;
t1 = t72 * t3 - t101;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72 * t22 - t71 * t23, 0, 0, 0, 0, 0, 0, -t72 * t11 - t71 * t12, -t72 * t17 - t71 * t18, -t72 * t13 - t71 * t14, -t72 * t2 - t71 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, 0, 0, 0, 0, 0, 0, -t74 * t78 + t88, -t74 * qJDD(3) - t76 * t78, 0, t76 * t37 + t74 * t38, 0, 0, 0, 0, 0, 0, -t74 * t52 + t72 * t88, t74 * t51 - t71 * t88, t74 * t53 + t76 * t54, t74 * t10 - t76 * t26, 0, 0, 0, 0, 0, 0, -t76 * t32 + t74 * t4, -t76 * t34 + t74 * t9, -t76 * t24 + t74 * t5, t74 * t1 - t76 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t37, -t38, 0, 0, t65, 0.2e1 * t71 * t89, 0, t66, 0, 0, -qJ(4) * t52 + t85 * t72, qJ(4) * t51 - t85 * t71, pkin(3) * t54 + qJ(4) * t53 + t10, -pkin(3) * t26 + qJ(4) * t10, t71 * (t75 * t35 - t73 * t91) + t72 * (t73 * t35 + t75 * t91), t71 * (-t75 * t32 - t73 * t34) + t72 * (-t73 * t32 + t75 * t34), t71 * (-t73 * t40 + t106) + t72 * (t75 * t40 + t107), t71 * (-t73 * t33 + t75 * t92) + t72 * (t75 * t33 + t73 * t92), t71 * (t75 * t39 - t98) + t72 * (t73 * t39 + t96), (t71 * (-t47 * t75 + t49 * t73) + t72 * (-t47 * t73 - t49 * t75)) * qJD(5), t71 * (-pkin(6) * t11 + t99) + t72 * (-pkin(4) * t32 + pkin(6) * t12 - t97) - pkin(3) * t32 + qJ(4) * t4, t71 * (-pkin(6) * t17 + t97) + t72 * (-pkin(4) * t34 + pkin(6) * t18 + t99) - pkin(3) * t34 + qJ(4) * t9, t71 * (-pkin(6) * t13 - t2) + t72 * (-pkin(4) * t24 + pkin(6) * t14 + t3) - pkin(3) * t24 + qJ(4) * t5, -pkin(6) * t101 + t72 * (-pkin(4) * t20 + pkin(6) * t3) - pkin(3) * t20 + qJ(4) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89, t90, -t54, t26, 0, 0, 0, 0, 0, 0, t32, t34, t24, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, t45 - t44, t46, -t36, t21, qJDD(5), -t7, -t8, 0, 0;];
tauJ_reg = t6;
