% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PPRRP2
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PPRRP2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP2_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:09:10
% EndTime: 2019-12-05 15:09:13
% DurationCPUTime: 0.55s
% Computational Cost: add. (796->115), mult. (1503->154), div. (0->0), fcn. (980->8), ass. (0->75)
t65 = sin(qJ(4));
t58 = t65 ^ 2;
t70 = qJD(3) ^ 2;
t103 = t58 * t70;
t69 = qJD(4) ^ 2;
t47 = t69 + t103;
t67 = cos(qJ(4));
t51 = t65 * t70 * t67;
t46 = qJDD(4) - t51;
t97 = t67 * t46;
t22 = -t65 * t47 + t97;
t84 = qJD(3) * qJD(4);
t86 = t65 * qJDD(3);
t36 = 0.2e1 * t67 * t84 + t86;
t62 = sin(pkin(8));
t63 = cos(pkin(8));
t66 = sin(qJ(3));
t68 = cos(qJ(3));
t109 = t62 * (t68 * t22 - t66 * t36) + t63 * (t66 * t22 + t68 * t36);
t108 = pkin(6) * t22;
t59 = t67 ^ 2;
t102 = t59 * t70;
t105 = t97 + t65 * (-t69 + t102);
t90 = sin(pkin(7));
t91 = cos(pkin(7));
t74 = -t91 * g(1) - t90 * g(2);
t87 = -g(3) + qJDD(1);
t25 = t62 * t87 + t63 * t74;
t71 = -t62 * t74 + t63 * t87;
t15 = t68 * t25 + t66 * t71;
t13 = -t70 * pkin(3) + qJDD(3) * pkin(6) + t15;
t39 = -t90 * g(1) + t91 * g(2) + qJDD(2);
t28 = t67 * t39;
t92 = t65 * qJ(5);
t78 = -t67 * pkin(4) - t92;
t89 = t70 * t78;
t5 = -qJDD(4) * pkin(4) - t69 * qJ(5) + (t13 + t89) * t65 + qJDD(5) - t28;
t104 = 2 * qJD(5);
t45 = qJDD(4) + t51;
t101 = t65 * t45;
t9 = t67 * t13 + t65 * t39;
t49 = -t69 - t102;
t21 = t67 * t49 - t101;
t83 = t65 * t84;
t85 = t67 * qJDD(3);
t37 = -0.2e1 * t83 + t85;
t95 = pkin(3) * t37 + pkin(6) * t21;
t93 = t58 + t59;
t40 = t93 * qJDD(3);
t43 = t93 * t70;
t94 = pkin(3) * t43 + pkin(6) * t40;
t88 = qJD(3) * t65;
t8 = t65 * t13 - t28;
t2 = t65 * t8 + t67 * t9;
t14 = -t66 * t25 + t68 * t71;
t79 = qJDD(4) * qJ(5) + (qJD(4) * t104) + t67 * t89 + t9;
t77 = t67 * t36 + t65 * t37;
t76 = t65 * t46 + t67 * t47;
t12 = -qJDD(3) * pkin(3) - t70 * pkin(6) - t14;
t75 = pkin(3) - t78;
t73 = t12 - (-t83 + t85) * pkin(4) - qJ(5) * t36;
t72 = t88 * t104 - t73;
t44 = (t58 - t59) * t70;
t42 = -t66 * qJDD(3) - t68 * t70;
t41 = t68 * qJDD(3) - t66 * t70;
t20 = t101 + t67 * (t69 - t103);
t19 = t67 * t45 + t65 * t49;
t18 = t36 * t65;
t17 = t37 * t67;
t10 = t62 * (t68 * t40 - t66 * t43) + t63 * (t66 * t40 + t68 * t43);
t6 = t62 * (t68 * t21 - t66 * t37) + t63 * (t66 * t21 + t68 * t37);
t4 = -t69 * pkin(4) + t79;
t3 = (pkin(4) * qJD(4) - (2 * qJD(5))) * t88 + t73;
t1 = t67 * t4 + t65 * t5;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t87, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62 * t25 + t63 * t71, 0, 0, 0, 0, 0, 0, t63 * t41 + t62 * t42, -t62 * t41 + t63 * t42, 0, t62 * (-t66 * t14 + t68 * t15) + t63 * (t68 * t14 + t66 * t15), 0, 0, 0, 0, 0, 0, t6, -t109, t10, t62 * (t66 * t12 + t68 * t2) + t63 * (-t68 * t12 + t66 * t2), 0, 0, 0, 0, 0, 0, t6, t10, t109, t62 * (t68 * t1 + t66 * t3) + t63 * (t66 * t1 - t68 * t3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, 0, 0, 0, 0, 0, t19, -t76, 0, t65 * t9 - t67 * t8, 0, 0, 0, 0, 0, 0, t19, 0, t76, t65 * t4 - t67 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t14, -t15, 0, 0, t18, t77, t20, t17, t105, 0, -t67 * t12 + t95, -pkin(3) * t36 + t65 * t12 - t108, t2 + t94, -pkin(3) * t12 + pkin(6) * t2, t18, t20, -t77, 0, -t105, t17, t37 * t92 + t67 * ((t37 - t83) * pkin(4) + t72) + t95, t67 * ((t43 - t69) * pkin(4) + t79) + (qJ(5) * t43 + t5) * t65 + t94, t65 * (-pkin(4) * t83 + t72) + t108 + t75 * t36, pkin(6) * t1 - t75 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, t44, t86, t51, t85, qJDD(4), -t8, -t9, 0, 0, -t51, t86, -t44, qJDD(4), -t85, t51, pkin(4) * t45 + qJ(5) * t49 - t5, (-pkin(4) * t65 + qJ(5) * t67) * qJDD(3), qJ(5) * t46 + (t47 - t69) * pkin(4) + t79, -pkin(4) * t5 + qJ(5) * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, t86, -t47, t5;];
tauJ_reg = t7;
