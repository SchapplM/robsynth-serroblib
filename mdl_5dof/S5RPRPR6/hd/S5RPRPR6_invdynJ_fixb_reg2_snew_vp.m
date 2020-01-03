% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRPR6
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRPR6_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR6_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR6_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR6_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:17:49
% EndTime: 2019-12-31 18:17:51
% DurationCPUTime: 0.51s
% Computational Cost: add. (1497->103), mult. (2176->146), div. (0->0), fcn. (1196->8), ass. (0->80)
t69 = (qJD(1) + qJD(3));
t67 = t69 ^ 2;
t68 = qJDD(1) + qJDD(3);
t76 = sin(qJ(3));
t79 = cos(qJ(3));
t46 = t79 * t67 + t76 * t68;
t49 = t76 * t67 - t79 * t68;
t73 = sin(pkin(8));
t74 = cos(pkin(8));
t114 = -pkin(1) * (t73 * t46 + t74 * t49) - pkin(2) * t49;
t113 = -pkin(1) * (t74 * t46 - t73 * t49) - pkin(2) * t46;
t112 = pkin(3) + pkin(7);
t72 = -g(3) + qJDD(2);
t75 = sin(qJ(5));
t78 = cos(qJ(5));
t107 = t68 * pkin(3);
t77 = sin(qJ(1));
t80 = cos(qJ(1));
t93 = t77 * g(1) - t80 * g(2);
t51 = qJDD(1) * pkin(1) + t93;
t82 = qJD(1) ^ 2;
t89 = t80 * g(1) + t77 * g(2);
t52 = -t82 * pkin(1) - t89;
t92 = t74 * t51 - t73 * t52;
t26 = qJDD(1) * pkin(2) + t92;
t99 = t73 * t51 + t74 * t52;
t27 = -t82 * pkin(2) + t99;
t17 = t79 * t26 - t76 * t27;
t87 = qJDD(4) - t17 - t107;
t16 = -t67 * qJ(4) + t87;
t83 = -t68 * pkin(7) + t16;
t7 = t75 * t72 - t78 * t83;
t8 = t78 * t72 + t75 * t83;
t3 = -t78 * t7 + t75 * t8;
t70 = t75 ^ 2;
t105 = t70 * t67;
t71 = t78 ^ 2;
t104 = t71 * t67;
t95 = t75 * t67 * t78;
t53 = qJDD(5) + t95;
t103 = t75 * t53;
t102 = t75 * t68;
t54 = qJDD(5) - t95;
t101 = t78 * t54;
t18 = t76 * t26 + t79 * t27;
t94 = (2 * qJD(4) * t69) + t18;
t97 = t68 * qJ(4);
t14 = -t67 * pkin(3) + t94 + t97;
t100 = -pkin(3) * t16 + qJ(4) * t14;
t98 = t70 + t71;
t96 = qJD(5) * t69;
t91 = 0.2e1 * t97 + t94;
t11 = -t67 * pkin(7) + t14;
t90 = qJ(4) * t11 - t112 * t3;
t57 = t78 * t68;
t41 = -0.2e1 * t75 * t96 + t57;
t81 = qJD(5) ^ 2;
t56 = -t81 - t104;
t33 = t78 * t56 - t103;
t88 = qJ(4) * t41 + t78 * t11 - t112 * t33;
t40 = 0.2e1 * t78 * t96 + t102;
t55 = -t81 - t105;
t32 = t75 * t55 + t101;
t86 = qJ(4) * t40 + t75 * t11 - t112 * t32;
t85 = -t107 + t87;
t45 = t98 * t68;
t50 = t98 * t67;
t84 = -qJ(4) * t50 + t112 * t45 - t3;
t35 = t101 - t75 * (t81 - t104);
t34 = t78 * (-t81 + t105) - t103;
t29 = t41 * t78;
t28 = t40 * t75;
t25 = t79 * t45 - t76 * t50;
t21 = -t78 * t40 - t75 * t41;
t20 = -t79 * t33 + t76 * t41;
t19 = -t79 * t32 + t76 * t40;
t5 = t79 * t17 + t76 * t18;
t4 = t76 * t14 - t79 * t16;
t1 = t76 * t11 - t79 * t3;
t2 = [0, 0, 0, 0, 0, qJDD(1), t93, t89, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (t74 * qJDD(1) - t73 * t82) + t92, pkin(1) * (-t73 * qJDD(1) - t74 * t82) - t99, 0, pkin(1) * (t73 * t99 + t74 * t92), 0, 0, 0, 0, 0, t68, t114 + t17, t113 - t18, 0, pkin(1) * (t73 * (-t76 * t17 + t79 * t18) + t74 * t5) + pkin(2) * t5, t68, 0, 0, 0, 0, 0, 0, -t114 + t85, -t113 + t91, pkin(1) * (t73 * (t79 * t14 + t76 * t16) + t74 * t4) + pkin(2) * t4 + t100, t29, t21, t35, t28, t34, 0, pkin(1) * (t73 * (t76 * t32 + t79 * t40) + t74 * t19) + pkin(2) * t19 + t86, pkin(1) * (t73 * (t76 * t33 + t79 * t41) + t74 * t20) + pkin(2) * t20 + t88, pkin(1) * (t73 * (-t76 * t45 - t79 * t50) + t74 * t25) + pkin(2) * t25 + t84, pkin(1) * (t73 * (t79 * t11 + t76 * t3) + t74 * t1) + pkin(2) * t1 + t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, 0, 0, 0, 0, 0, 0, -t75 * t54 + t78 * t55, -t78 * t53 - t75 * t56, 0, t75 * t7 + t78 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, t17, -t18, 0, 0, t68, 0, 0, 0, 0, 0, 0, t85, t91, t100, t29, t21, t35, t28, t34, 0, t86, t88, t84, t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, -t67, t16, 0, 0, 0, 0, 0, 0, t32, t33, -t45, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, (-t70 + t71) * t67, t57, -t95, -t102, qJDD(5), -t7, -t8, 0, 0;];
tauJ_reg = t2;
