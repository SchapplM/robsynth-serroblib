% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRRRR2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:30
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRRRR2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR2_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_invdynJ_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:30:18
% EndTime: 2019-07-18 13:30:20
% DurationCPUTime: 0.36s
% Computational Cost: add. (1373->84), mult. (1772->124), div. (0->0), fcn. (1068->8), ass. (0->75)
t73 = sin(qJ(2));
t77 = cos(qJ(2));
t84 = t73 * g(1) - t77 * g(2);
t49 = qJDD(2) * pkin(2) + t84;
t81 = t77 * g(1) + t73 * g(2);
t50 = -qJD(2) ^ 2 * pkin(2) - t81;
t72 = sin(qJ(3));
t76 = cos(qJ(3));
t31 = t76 * t49 - t72 * t50;
t65 = qJDD(2) + qJDD(3);
t25 = t65 * pkin(3) + t31;
t32 = t72 * t49 + t76 * t50;
t66 = qJD(2) + qJD(3);
t64 = t66 ^ 2;
t26 = -t64 * pkin(3) + t32;
t71 = sin(qJ(4));
t75 = cos(qJ(4));
t16 = t71 * t25 + t75 * t26;
t61 = qJDD(4) + t65;
t13 = t61 * pkin(6) + t16;
t69 = -g(3) + qJDD(1);
t70 = sin(qJ(5));
t74 = cos(qJ(5));
t10 = t74 * t13 + t70 * t69;
t9 = t70 * t13 - t74 * t69;
t4 = t74 * t10 + t70 * t9;
t15 = t75 * t25 - t71 * t26;
t62 = qJD(4) + t66;
t60 = t62 ^ 2;
t14 = -t60 * pkin(6) - t15;
t2 = -t75 * t14 + t71 * t4;
t3 = pkin(6) * t4;
t94 = pkin(3) * t2 + t3;
t53 = t70 * t60 * t74;
t47 = qJDD(5) + t53;
t93 = t70 * t47;
t92 = t70 * t61;
t48 = qJDD(5) - t53;
t91 = t74 * t48;
t67 = t70 ^ 2;
t56 = t67 * t60;
t78 = qJD(5) ^ 2;
t51 = -t56 - t78;
t36 = -t70 * t51 - t91;
t90 = pkin(6) * t36 + t70 * t14;
t68 = t74 ^ 2;
t57 = t68 * t60;
t52 = -t57 - t78;
t35 = t74 * t52 - t93;
t89 = pkin(6) * t35 - t74 * t14;
t88 = qJD(5) * t62;
t42 = (t67 + t68) * t61;
t87 = pkin(6) * t42 + t4;
t40 = 0.2e1 * t74 * t88 + t92;
t20 = t71 * t36 - t75 * t40;
t86 = pkin(3) * t20 + t90;
t55 = t74 * t61;
t41 = -0.2e1 * t70 * t88 + t55;
t19 = t71 * t35 + t75 * t41;
t85 = pkin(3) * t19 + t89;
t45 = t56 + t57;
t23 = t71 * t42 + t75 * t45;
t83 = pkin(3) * t23 + t87;
t80 = t71 * t60 - t75 * t61;
t82 = -pkin(3) * t80 + t15;
t43 = -t75 * t60 - t71 * t61;
t79 = pkin(3) * t43 - t16;
t34 = t93 + t74 * (-t56 + t78);
t33 = t70 * (t57 - t78) + t91;
t28 = t40 * t70;
t27 = t41 * t74;
t21 = t74 * t40 + t70 * t41;
t6 = t75 * t15 + t71 * t16;
t5 = pkin(3) * t6;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t69, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, 0, 0, 0, 0, 0, 0, t74 * t47 + t70 * t52, -t70 * t48 + t74 * t51, 0, t70 * t10 - t74 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t84, t81, 0, 0, 0, 0, 0, 0, 0, t65, pkin(2) * (-t72 * t64 + t76 * t65) + t31, pkin(2) * (-t76 * t64 - t72 * t65) - t32, 0, pkin(2) * (t76 * t31 + t72 * t32), 0, 0, 0, 0, 0, t61, pkin(2) * (t72 * t43 - t76 * t80) + t82, pkin(2) * (t76 * t43 + t72 * t80) + t79, 0, pkin(2) * (t72 * (-t71 * t15 + t75 * t16) + t76 * t6) + t5, t28, t21, t34, t27, t33, 0, pkin(2) * (t72 * (t75 * t35 - t71 * t41) + t76 * t19) + t85, pkin(2) * (t72 * (t75 * t36 + t71 * t40) + t76 * t20) + t86, pkin(2) * (t72 * (t75 * t42 - t71 * t45) + t76 * t23) + t83, pkin(2) * (t72 * (t71 * t14 + t75 * t4) + t76 * t2) + t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, t31, -t32, 0, 0, 0, 0, 0, 0, 0, t61, t82, t79, 0, t5, t28, t21, t34, t27, t33, 0, t85, t86, t83, t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t15, -t16, 0, 0, t28, t21, t34, t27, t33, 0, t89, t90, t87, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, t56 - t57, t92, t53, t55, qJDD(5), -t9, -t10, 0, 0;];
tauJ_reg  = t1;
