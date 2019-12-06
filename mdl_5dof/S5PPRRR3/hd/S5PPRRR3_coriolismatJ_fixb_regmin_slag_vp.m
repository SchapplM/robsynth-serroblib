% Calculate minimal parameter regressor of coriolis matrix for
% S5PPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x19]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PPRRR3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:17:02
% EndTime: 2019-12-05 15:17:04
% DurationCPUTime: 0.60s
% Computational Cost: add. (452->81), mult. (1285->137), div. (0->0), fcn. (1316->8), ass. (0->82)
t102 = qJD(4) + qJD(5);
t76 = sin(pkin(9));
t79 = sin(qJ(3));
t132 = t76 * t79;
t77 = sin(qJ(5));
t78 = sin(qJ(4));
t130 = t77 * t78;
t133 = cos(qJ(5));
t80 = cos(qJ(4));
t92 = t133 * t80;
t61 = -t92 + t130;
t84 = t92 / 0.2e1 - t130 / 0.2e1;
t82 = t61 / 0.2e1 + t84;
t14 = t82 * t132;
t81 = cos(qJ(3));
t33 = t82 * t81;
t128 = -t14 * qJD(1) + t33 * qJD(2);
t129 = t77 * t80;
t93 = t133 * t78;
t85 = t129 / 0.2e1 + t93 / 0.2e1;
t88 = t93 + t129;
t83 = -t88 / 0.2e1 + t85;
t13 = t83 * t132;
t32 = t83 * t81;
t127 = -t13 * qJD(1) + t32 * qJD(2);
t94 = -t132 / 0.2e1;
t15 = t84 * t132 + t61 * t94;
t122 = cos(pkin(9));
t131 = t76 * t81;
t51 = t122 * t80 + t78 * t131;
t52 = -t122 * t78 + t80 * t131;
t143 = t15 * qJD(3) + t102 * (t133 * t51 + t77 * t52);
t16 = t85 * t132 - t88 * t94;
t142 = t16 * qJD(3) + t102 * (-t133 * t52 + t77 * t51);
t138 = pkin(6) + pkin(7);
t68 = t138 * t78;
t69 = t138 * t80;
t141 = t127 + t102 * (-t133 * t69 + t77 * t68);
t140 = t128 + t102 * (t133 * t68 + t77 * t69);
t139 = t81 / 0.2e1;
t137 = pkin(4) * t77;
t136 = pkin(4) * t78;
t75 = -t80 * pkin(4) - pkin(3);
t121 = qJD(3) * t75;
t120 = qJD(3) * t80;
t119 = qJD(5) * t75;
t21 = t61 ^ 2 - t88 ^ 2;
t116 = t21 * qJD(3);
t30 = t61 * t136 + t75 * t88;
t111 = t30 * qJD(3);
t31 = t136 * t88 - t75 * t61;
t110 = t31 * qJD(3);
t70 = -t78 ^ 2 + t80 ^ 2;
t107 = t70 * qJD(3);
t106 = t78 * qJD(4);
t105 = t79 * qJD(3);
t104 = t80 * qJD(4);
t103 = t81 * qJD(3);
t101 = pkin(3) * t78 * qJD(3);
t100 = pkin(3) * t120;
t99 = t61 * t121;
t98 = t88 * t121;
t97 = t78 * t120;
t96 = t76 * t105;
t95 = t76 * t103;
t91 = t133 * qJD(4);
t90 = t133 * qJD(5);
t89 = t79 * t102;
t37 = t102 * t88;
t87 = -t80 * t103 + t79 * t106;
t86 = t78 * t103 + t79 * t104;
t38 = t88 * t61 * qJD(3);
t36 = t102 * t61;
t35 = -t139 * t88 - t85 * t81;
t34 = t61 * t139 - t84 * t81;
t26 = t33 * qJD(3);
t24 = t32 * qJD(3);
t7 = t14 * qJD(3);
t5 = t13 * qJD(3);
t2 = t35 * qJD(3) + t61 * t89;
t1 = t34 * qJD(3) + t88 * t89;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t95, t96, 0, 0, 0, 0, 0, t87 * t76, t86 * t76, 0, 0, 0, 0, 0, t102 * t16 + t61 * t95, t102 * t15 + t88 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52 * qJD(4) + t78 * t96, t51 * qJD(4) + t80 * t96, 0, 0, 0, 0, 0, t142, t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t142, t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t105, -t103, 0, 0, 0, 0, 0, -t80 * t105 - t81 * t106, -t81 * t104 + t78 * t105, 0, 0, 0, 0, 0, t102 * t35 + t61 * t105, t102 * t34 + t105 * t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, t87, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t102 * t13, -t102 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102 * t32, t102 * t33; 0, 0, 0, 0, 0, t78 * t104, t70 * qJD(4), 0, 0, 0, -pkin(3) * t106, -pkin(3) * t104, -t61 * t37, t102 * t21, 0, 0, 0, t30 * qJD(4) + t119 * t88, t31 * qJD(4) - t61 * t119; 0, 0, 0, 0, 0, t97, t107, t104, -t106, 0, -pkin(6) * t104 - t101, pkin(6) * t106 - t100, -t38, t116, -t36, -t37, 0, t111 + t141, t110 + t140; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, t116, -t36, -t37, 0, t98 + t141, -t99 + t140; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t26; 0, 0, 0, 0, 0, -t97, -t107, 0, 0, 0, t101, t100, t38, -t116, 0, 0, 0, -t111 - t127, -t110 - t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(5) * t137, -pkin(4) * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t102 * t137, (-t91 - t90) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t116, 0, 0, 0, -t98 - t127, t99 - t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4) * t137, pkin(4) * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t3;
