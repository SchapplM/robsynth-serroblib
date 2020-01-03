% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PPRRP4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PPRRP4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP4_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:34:41
% EndTime: 2019-12-31 17:34:43
% DurationCPUTime: 0.41s
% Computational Cost: add. (521->100), mult. (1066->122), div. (0->0), fcn. (602->6), ass. (0->73)
t81 = (qJD(3) * qJD(5));
t96 = 2 * t81;
t69 = sin(qJ(4));
t71 = cos(qJ(4));
t74 = qJD(3) ^ 2;
t51 = t69 * t74 * t71;
t45 = qJDD(4) + t51;
t95 = pkin(4) * t45;
t60 = t71 * qJDD(3);
t82 = qJD(3) * qJD(4);
t79 = t69 * t82;
t36 = t60 - t79;
t84 = qJD(3) * t69;
t44 = qJD(4) * pkin(4) - qJ(5) * t84;
t66 = sin(pkin(7));
t67 = cos(pkin(7));
t39 = -t66 * g(1) + t67 * g(2) + qJDD(2);
t43 = -t67 * g(1) - t66 * g(2);
t70 = sin(qJ(3));
t72 = cos(qJ(3));
t16 = t70 * t39 + t72 * t43;
t13 = -t74 * pkin(3) + qJDD(3) * pkin(6) + t16;
t63 = g(3) - qJDD(1);
t8 = t71 * t13 + t69 * t63;
t75 = t36 * qJ(5) - qJD(4) * t44 + t71 * t96 + t8;
t61 = t69 ^ 2;
t94 = t61 * t74;
t62 = t71 ^ 2;
t93 = t62 * t74;
t92 = t69 * t13;
t91 = t69 * t45;
t46 = qJDD(4) - t51;
t90 = t71 * t46;
t73 = qJD(4) ^ 2;
t49 = -t73 - t93;
t26 = t71 * t49 - t91;
t37 = t60 - 0.2e1 * t79;
t89 = pkin(3) * t37 + pkin(6) * t26;
t48 = -t73 - t94;
t27 = -t69 * t48 - t90;
t59 = t69 * qJDD(3);
t78 = t71 * t82;
t34 = t59 + 0.2e1 * t78;
t88 = -pkin(3) * t34 + pkin(6) * t27;
t86 = t61 + t62;
t40 = t86 * qJDD(3);
t41 = t86 * t74;
t87 = pkin(3) * t41 + pkin(6) * t40;
t85 = qJ(5) * t69;
t83 = qJ(5) * qJDD(3);
t57 = t71 * t63;
t7 = -t57 + t92;
t2 = t69 * t7 + t71 * t8;
t15 = t72 * t39 - t70 * t43;
t12 = -qJDD(3) * pkin(3) - t74 * pkin(6) - t15;
t35 = t59 + t78;
t76 = -t57 + (t35 - t78) * qJ(5) - t95;
t3 = -0.2e1 * t69 * t81 - t76 - t92;
t5 = -t36 * pkin(4) - qJ(5) * t93 + t44 * t84 + qJDD(5) + t12;
t42 = (t61 - t62) * t74;
t25 = t69 * t46 - t71 * t48;
t24 = t91 + t71 * (t73 - t94);
t23 = -t71 * t45 - t69 * t49;
t22 = t69 * (-t73 + t93) + t90;
t21 = (t35 + t78) * t69;
t20 = (t36 - t79) * t71;
t17 = t70 * t40 + t72 * t41;
t14 = t71 * t34 + t69 * t37;
t10 = t70 * t27 - t72 * t34;
t9 = t70 * t26 + t72 * t37;
t4 = -pkin(4) * t93 + t75;
t1 = -t69 * t3 + t71 * t4;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, 0, 0, 0, 0, 0, 0, t23, t25, 0, -t69 * t8 + t71 * t7, 0, 0, 0, 0, 0, 0, t23, t25, 0, -t71 * t3 - t69 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, 0, 0, 0, 0, 0, t72 * qJDD(3) - t70 * t74, -t70 * qJDD(3) - t72 * t74, 0, t72 * t15 + t70 * t16, 0, 0, 0, 0, 0, 0, t9, t10, t17, -t72 * t12 + t70 * t2, 0, 0, 0, 0, 0, 0, t9, t10, t17, t70 * t1 - t72 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t15, -t16, 0, 0, t21, t14, t24, t20, t22, 0, -t71 * t12 + t89, t69 * t12 + t88, t2 + t87, -pkin(3) * t12 + pkin(6) * t2, t21, t14, t24, t20, t22, 0, -t45 * t85 + t71 * (pkin(4) * t37 + qJ(5) * t49 - t5) + t89, t69 * (-qJ(5) * t48 + t5) + t71 * (-pkin(4) * t34 - qJ(5) * t46) + t88, t71 * (t71 * t83 + (t41 - t93) * pkin(4) + t75) + ((t13 + t96 + t83) * t69 + t76) * t69 + t87, -t3 * t85 + t71 * (-pkin(4) * t5 + qJ(5) * t4) - pkin(3) * t5 + pkin(6) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, t42, t59, t51, t60, qJDD(4), -t7, -t8, 0, 0, -t51, t42, t59, t51, t60, qJDD(4), t3 + t95, (t48 + t93) * pkin(4) - t75, -pkin(4) * t59, pkin(4) * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, t34, -t41, t5;];
tauJ_reg = t6;
