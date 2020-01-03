% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4RRRR1
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
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RRRR1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR1_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR1_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR1_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:22:14
% EndTime: 2019-12-31 17:22:15
% DurationCPUTime: 0.34s
% Computational Cost: add. (1508->86), mult. (1886->124), div. (0->0), fcn. (1084->8), ass. (0->70)
t74 = sin(qJ(1));
t78 = cos(qJ(1));
t88 = t74 * g(1) - t78 * g(2);
t53 = qJDD(1) * pkin(1) + t88;
t82 = t78 * g(1) + t74 * g(2);
t54 = -qJD(1) ^ 2 * pkin(1) - t82;
t73 = sin(qJ(2));
t77 = cos(qJ(2));
t32 = t77 * t53 - t73 * t54;
t67 = qJDD(1) + qJDD(2);
t26 = t67 * pkin(2) + t32;
t33 = t73 * t53 + t77 * t54;
t68 = qJD(1) + qJD(2);
t66 = t68 ^ 2;
t27 = -t66 * pkin(2) + t33;
t72 = sin(qJ(3));
t76 = cos(qJ(3));
t17 = t72 * t26 + t76 * t27;
t64 = qJD(3) + t68;
t62 = t64 ^ 2;
t63 = qJDD(3) + t67;
t15 = -t62 * pkin(3) + t63 * pkin(7) + t17;
t71 = sin(qJ(4));
t75 = cos(qJ(4));
t10 = -t71 * g(3) + t75 * t15;
t9 = t75 * g(3) + t71 * t15;
t4 = t75 * t10 + t71 * t9;
t16 = t76 * t26 - t72 * t27;
t14 = -t63 * pkin(3) - t62 * pkin(7) - t16;
t96 = -pkin(3) * t14 + pkin(7) * t4;
t55 = t75 * t62 * t71;
t95 = t71 * (qJDD(4) + t55);
t94 = t71 * t63;
t93 = t75 * (qJDD(4) - t55);
t92 = qJD(4) * t64;
t2 = -t76 * t14 + t72 * t4;
t91 = pkin(2) * t2 + t96;
t69 = t71 ^ 2;
t58 = t69 * t62;
t79 = qJD(4) ^ 2;
t37 = -t93 - t71 * (-t58 - t79);
t43 = 0.2e1 * t75 * t92 + t94;
t90 = -pkin(3) * t43 + pkin(7) * t37 + t71 * t14;
t70 = t75 ^ 2;
t59 = t70 * t62;
t36 = t75 * (-t59 - t79) - t95;
t57 = t75 * t63;
t44 = -0.2e1 * t71 * t92 + t57;
t89 = pkin(3) * t44 + pkin(7) * t36 - t75 * t14;
t46 = (t69 + t70) * t63;
t49 = t58 + t59;
t87 = pkin(3) * t49 + pkin(7) * t46 + t4;
t21 = t72 * t37 - t76 * t43;
t86 = pkin(2) * t21 + t90;
t20 = t72 * t36 + t76 * t44;
t85 = pkin(2) * t20 + t89;
t24 = t72 * t46 + t76 * t49;
t84 = pkin(2) * t24 + t87;
t81 = t72 * t62 - t76 * t63;
t83 = -pkin(2) * t81 + t16;
t47 = -t76 * t62 - t72 * t63;
t80 = pkin(2) * t47 - t17;
t35 = t95 + t75 * (-t58 + t79);
t34 = t71 * (t59 - t79) + t93;
t29 = t43 * t71;
t28 = t44 * t75;
t22 = t75 * t43 + t71 * t44;
t6 = t76 * t16 + t72 * t17;
t5 = pkin(2) * t6;
t1 = [0, 0, 0, 0, 0, qJDD(1), t88, t82, 0, 0, 0, 0, 0, 0, 0, t67, pkin(1) * (-t73 * t66 + t77 * t67) + t32, pkin(1) * (-t77 * t66 - t73 * t67) - t33, 0, pkin(1) * (t77 * t32 + t73 * t33), 0, 0, 0, 0, 0, t63, pkin(1) * (t73 * t47 - t77 * t81) + t83, pkin(1) * (t77 * t47 + t73 * t81) + t80, 0, pkin(1) * (t73 * (-t72 * t16 + t76 * t17) + t77 * t6) + t5, t29, t22, t35, t28, t34, 0, pkin(1) * (t73 * (t76 * t36 - t72 * t44) + t77 * t20) + t85, pkin(1) * (t73 * (t76 * t37 + t72 * t43) + t77 * t21) + t86, pkin(1) * (t73 * (t76 * t46 - t72 * t49) + t77 * t24) + t84, pkin(1) * (t73 * (t72 * t14 + t76 * t4) + t77 * t2) + t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, t32, -t33, 0, 0, 0, 0, 0, 0, 0, t63, t83, t80, 0, t5, t29, t22, t35, t28, t34, 0, t85, t86, t84, t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, t16, -t17, 0, 0, t29, t22, t35, t28, t34, 0, t89, t90, t87, t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, t58 - t59, t94, t55, t57, qJDD(4), -t9, -t10, 0, 0;];
tauJ_reg = t1;
