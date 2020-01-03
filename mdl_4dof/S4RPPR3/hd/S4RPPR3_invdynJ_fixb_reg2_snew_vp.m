% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4RPPR3
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
%   pkin=[a2,a3,a4,d1,d4,theta2,theta3]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RPPR3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR3_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR3_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPPR3_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:37:56
% EndTime: 2019-12-31 16:37:58
% DurationCPUTime: 0.60s
% Computational Cost: add. (1147->130), mult. (2516->196), div. (0->0), fcn. (1596->8), ass. (0->79)
t78 = qJD(1) ^ 2;
t103 = sin(qJ(1));
t104 = cos(qJ(1));
t85 = t104 * g(1) + t103 * g(2);
t53 = -t78 * pkin(1) - t85;
t72 = sin(pkin(6));
t74 = cos(pkin(6));
t84 = t103 * g(1) - t104 * g(2);
t82 = qJDD(1) * pkin(1) + t84;
t97 = t74 * t53 + t72 * t82;
t113 = -t78 * pkin(2) + qJDD(1) * qJ(3) + (2 * qJD(1) * qJD(3)) + t97;
t73 = cos(pkin(7));
t67 = t73 ^ 2;
t71 = sin(pkin(7));
t79 = t71 ^ 2;
t55 = (t67 + t79) * t78;
t75 = sin(qJ(4));
t76 = cos(qJ(4));
t45 = (t71 * t75 - t73 * t76) * qJD(1);
t86 = t71 * t76 + t73 * t75;
t47 = t86 * qJD(1);
t36 = t47 * t45;
t107 = qJDD(4) - t36;
t112 = t107 * t75;
t111 = t107 * t76;
t88 = -t72 * t53 + t74 * t82;
t26 = -qJDD(1) * pkin(2) - t78 * qJ(3) + qJDD(3) - t88;
t110 = (pkin(1) * t72 + qJ(3)) * t55 + t26 - (pkin(1) * t74 + pkin(2)) * qJDD(1);
t69 = -g(3) + qJDD(2);
t61 = t73 * t69;
t108 = t61 + (pkin(3) * t73 * t78 - pkin(5) * qJDD(1) - t113) * t71;
t42 = t45 ^ 2;
t43 = t47 ^ 2;
t102 = t67 * t78;
t22 = t113 * t73 + t71 * t69;
t92 = t73 * qJDD(1);
t16 = -pkin(3) * t102 + pkin(5) * t92 + t22;
t7 = -t76 * t108 + t75 * t16;
t8 = t108 * t75 + t76 * t16;
t2 = -t76 * t7 + t75 * t8;
t105 = t71 * t2;
t20 = -pkin(3) * t92 + t26 + (-t78 * t79 - t102) * pkin(5);
t101 = t75 * t20;
t30 = qJDD(4) + t36;
t100 = t75 * t30;
t99 = t76 * t20;
t98 = t76 * t30;
t95 = t45 * qJD(4);
t94 = t47 * qJD(4);
t93 = t71 * qJDD(1);
t3 = t75 * t7 + t76 * t8;
t21 = t113 * t71 - t61;
t10 = t71 * t21 + t73 * t22;
t23 = -t75 * t93 + t76 * t92;
t44 = t86 * qJDD(1);
t77 = qJD(4) ^ 2;
t66 = t67 * qJDD(1);
t65 = t79 * qJDD(1);
t54 = t66 + t65;
t39 = -t43 - t77;
t38 = -t43 + t77;
t37 = t42 - t77;
t35 = t44 - t95;
t34 = t44 - 0.2e1 * t95;
t33 = t23 - t94;
t32 = -t23 + 0.2e1 * t94;
t28 = -t77 - t42;
t24 = -t42 - t43;
t18 = -t75 * t39 - t98;
t17 = t76 * t39 - t100;
t14 = t76 * t23 + t75 * t44;
t13 = t75 * t23 - t76 * t44;
t12 = t76 * t28 - t112;
t11 = t75 * t28 + t111;
t9 = -t71 * t17 + t73 * t18;
t5 = -t71 * t13 + t73 * t14;
t4 = -t71 * t11 + t73 * t12;
t1 = t73 * t3 - t105;
t6 = [0, 0, 0, 0, 0, qJDD(1), t84, t85, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (t74 * qJDD(1) - t72 * t78) + t88, pkin(1) * (-t72 * qJDD(1) - t74 * t78) - t97, 0, pkin(1) * (t72 * t97 + t74 * t88), t65, 0.2e1 * t71 * t92, 0, t66, 0, 0, -t110 * t73, t110 * t71, pkin(2) * t55 + qJ(3) * t54 + pkin(1) * (t72 * t54 + t74 * t55) + t10, -pkin(2) * t26 + qJ(3) * t10 + pkin(1) * (t72 * t10 - t74 * t26), t71 * (t76 * t35 - t75 * t94) + t73 * (t75 * t35 + t76 * t94), t71 * (-t76 * t32 - t75 * t34) + t73 * (-t75 * t32 + t76 * t34), t71 * (-t75 * t38 + t111) + t73 * (t76 * t38 + t112), t71 * (-t75 * t33 + t76 * t95) + t73 * (t76 * t33 + t75 * t95), t71 * (t76 * t37 - t100) + t73 * (t75 * t37 + t98), (t71 * (-t45 * t76 + t47 * t75) + t73 * (-t45 * t75 - t47 * t76)) * qJD(4), t71 * (-pkin(5) * t11 + t101) + t73 * (-pkin(3) * t32 + pkin(5) * t12 - t99) - pkin(2) * t32 + qJ(3) * t4 + pkin(1) * (-t74 * t32 + t72 * t4), t71 * (-pkin(5) * t17 + t99) + t73 * (-pkin(3) * t34 + pkin(5) * t18 + t101) - pkin(2) * t34 + qJ(3) * t9 + pkin(1) * (-t74 * t34 + t72 * t9), t71 * (-pkin(5) * t13 - t2) + t73 * (-pkin(3) * t24 + pkin(5) * t14 + t3) - pkin(2) * t24 + qJ(3) * t5 + pkin(1) * (-t74 * t24 + t72 * t5), -pkin(5) * t105 + t73 * (-pkin(3) * t20 + pkin(5) * t3) - pkin(2) * t20 + qJ(3) * t1 + pkin(1) * (t72 * t1 - t74 * t20); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73 * t21 + t71 * t22, 0, 0, 0, 0, 0, 0, t73 * t11 + t71 * t12, t73 * t17 + t71 * t18, t73 * t13 + t71 * t14, t73 * t2 + t71 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t92, t93, -t55, t26, 0, 0, 0, 0, 0, 0, t32, t34, t24, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, t43 - t42, t44, -t36, t23, qJDD(4), -t7, -t8, 0, 0;];
tauJ_reg = t6;
