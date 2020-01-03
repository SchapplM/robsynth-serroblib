% Calculate minimal parameter regressor of coriolis matrix for
% S5RPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x19]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRPP4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:14:52
% EndTime: 2019-12-31 18:14:54
% DurationCPUTime: 0.71s
% Computational Cost: add. (1378->108), mult. (2297->141), div. (0->0), fcn. (2456->4), ass. (0->84)
t133 = sin(pkin(7));
t90 = -pkin(1) - pkin(6);
t153 = -qJ(4) + t90;
t89 = cos(qJ(3));
t78 = t153 * t89;
t104 = t133 * t78;
t134 = cos(pkin(7));
t88 = sin(qJ(3));
t77 = t153 * t88;
t70 = t134 * t77;
t152 = t104 + t70;
t43 = t133 * t77 - t134 * t78;
t74 = -t133 * t88 + t134 * t89;
t75 = -t133 * t89 - t134 * t88;
t99 = t152 * t75 + t43 * t74;
t145 = t99 * qJD(1);
t72 = t75 ^ 2;
t73 = t74 ^ 2;
t102 = -t73 / 0.2e1 - t72 / 0.2e1;
t28 = -0.1e1 / 0.2e1 + t102;
t159 = -t28 * qJD(2) - t145;
t27 = 0.1e1 / 0.2e1 + t102;
t158 = t27 * qJD(2) + t145;
t155 = t72 + t73;
t157 = qJD(4) * t155;
t156 = t155 * qJD(1);
t140 = qJD(3) * pkin(3);
t154 = t140 * (t133 * t74 + t134 * t75);
t149 = qJD(4) * t99;
t80 = pkin(3) * t133 + qJ(5);
t84 = -pkin(3) * t134 - pkin(4);
t144 = -(t74 * t80 - t84 * t75) * qJD(3) + qJD(5) * t75;
t101 = t70 / 0.2e1;
t143 = t75 / 0.2e1;
t142 = t74 * pkin(4);
t141 = t89 * pkin(3);
t135 = t75 * qJ(5);
t85 = t88 * pkin(3) + qJ(2);
t40 = -t75 * pkin(4) - t74 * qJ(5) + t85;
t41 = -t135 + t141 + t142;
t15 = t40 * t74 - t41 * t75;
t131 = t15 * qJD(1);
t16 = -t40 * t75 - t41 * t74;
t130 = t16 * qJD(1);
t87 = t141 / 0.2e1;
t19 = t87 + (-qJ(5) / 0.2e1 - t80 / 0.2e1) * t75 + (pkin(4) / 0.2e1 - t84 / 0.2e1) * t74;
t127 = t19 * qJD(1);
t126 = t28 * qJD(1);
t92 = t133 * t143 - t134 * t74 / 0.2e1;
t39 = (-t89 / 0.2e1 + t92) * pkin(3);
t125 = t39 * qJD(1);
t124 = t40 * qJD(1);
t121 = t73 * qJD(1);
t120 = t74 * qJD(1);
t119 = t74 * qJD(5);
t118 = t75 * qJD(1);
t117 = t75 * qJD(3);
t79 = t88 ^ 2 - t89 ^ 2;
t116 = t79 * qJD(1);
t115 = t85 * qJD(1);
t114 = t88 * qJD(1);
t113 = t88 * qJD(3);
t112 = t89 * qJD(1);
t111 = t89 * qJD(3);
t110 = qJ(2) * qJD(3);
t109 = qJD(1) * qJ(2);
t108 = t74 * t118;
t107 = t88 * t112;
t106 = t88 * t109;
t105 = t89 * t109;
t103 = -qJD(2) + t119;
t1 = t40 * t41;
t97 = t1 * qJD(1);
t8 = t85 * t141;
t96 = t8 * qJD(1);
t42 = t101 - t70 / 0.2e1;
t95 = t42 * qJD(1) + t80 * qJD(3);
t63 = t74 * qJD(3);
t38 = pkin(3) * t92 + t87;
t31 = 0.2e1 * t101 + t104;
t23 = t28 * qJD(4);
t21 = t27 * qJD(4);
t20 = t80 * t143 + t84 * t74 / 0.2e1 + t87 - t135 / 0.2e1 + t142 / 0.2e1;
t2 = [0, 0, 0, 0, qJD(2), qJ(2) * qJD(2), -t88 * t111, t79 * qJD(3), 0, 0, 0, qJD(2) * t88 + t110 * t89, qJD(2) * t89 - t110 * t88, t157, t85 * qJD(2) + t8 * qJD(3) + t149, t15 * qJD(3) + t103 * t75, t157, -t74 * qJD(2) + t16 * qJD(3) + t73 * qJD(5), t1 * qJD(3) - t103 * t40 + t149; 0, 0, 0, 0, qJD(1), t109, 0, 0, 0, 0, 0, t114, t112, 0, t21 + t115, -t118, 0, -t120, t21 + t124; 0, 0, 0, 0, 0, 0, -t107, t116, -t113, -t111, 0, -t113 * t90 + t105, -t111 * t90 - t106, -t154, (-t133 * t43 - t134 * t152) * t140 + t38 * qJD(4) + t96, -qJD(3) * t152 + t131, t144, -qJD(3) * t43 + t130, (t152 * t84 - t43 * t80) * qJD(3) + t20 * qJD(4) + t31 * qJD(5) + t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t156, t38 * qJD(3) + t158, 0, t156, 0, t20 * qJD(3) + t158; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t108, t117, t121, t31 * qJD(3) - t120 * t40; 0, 0, 0, 0, -qJD(1), -t109, 0, 0, 0, 0, 0, -t114, -t112, 0, t23 - t115, t118, 0, t120, t23 - t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t113, -t111, 0, t154, t117, 0, t63, -t144; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t126, 0, 0, 0, t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t117; 0, 0, 0, 0, 0, 0, t107, -t116, 0, 0, 0, -t105, t106, 0, t39 * qJD(4) - t96, -t74 * qJD(4) - t131, 0, t75 * qJD(4) - t130, -t19 * qJD(4) + t42 * qJD(5) - t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), t80 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t125, -t120, 0, t118, -t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t156, -t39 * qJD(3) + t159, t63, -t156, -t117, t19 * qJD(3) - t119 + t159; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t126, 0, 0, 0, -t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t125, t120, 0, -t118, t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108, 0, -t121, -t42 * qJD(3) + (qJD(4) + t124) * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), -t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t2;
