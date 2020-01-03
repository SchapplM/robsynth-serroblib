% Calculate minimal parameter regressor of coriolis matrix for
% S5RPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x22]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRPP3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:56
% EndTime: 2019-12-31 18:12:58
% DurationCPUTime: 0.66s
% Computational Cost: add. (969->116), mult. (1845->124), div. (0->0), fcn. (2041->4), ass. (0->83)
t116 = cos(qJ(3));
t114 = pkin(6) + qJ(2);
t74 = sin(pkin(7));
t65 = t114 * t74;
t75 = cos(pkin(7));
t66 = t114 * t75;
t77 = sin(qJ(3));
t29 = t116 * t65 + t77 * t66;
t62 = t116 * t74 + t77 * t75;
t18 = t62 * pkin(4) + t29;
t126 = qJD(3) * t18;
t119 = t62 ^ 2;
t60 = -t116 * t75 + t77 * t74;
t120 = t60 ^ 2;
t121 = t120 + t119;
t125 = t121 * qJD(1);
t124 = t121 * qJD(2);
t50 = t62 * qJD(4);
t123 = t60 * qJD(5) - t50;
t54 = qJ(4) * t60;
t122 = pkin(3) * t62 + t54;
t117 = t60 * pkin(3);
t76 = pkin(3) + qJ(5);
t115 = t76 * t60;
t67 = t74 ^ 2 + t75 ^ 2;
t107 = t62 * qJ(4);
t69 = -t75 * pkin(2) - pkin(1);
t78 = t69 - t107;
t15 = t78 + t115;
t17 = t76 * t62 + t54;
t1 = t15 * t17;
t112 = t1 * qJD(1);
t22 = t78 + t117;
t3 = t22 * t122;
t110 = t3 * qJD(1);
t5 = t15 * t62 + t17 * t60;
t109 = t5 * qJD(1);
t6 = t15 * t60 - t17 * t62;
t108 = t6 * qJD(1);
t30 = t116 * t66 - t77 * t65;
t20 = -t60 * pkin(4) + t30;
t7 = t18 * t62 - t20 * t60;
t106 = t7 * qJD(1);
t8 = -t122 * t60 - t22 * t62;
t105 = t8 * qJD(1);
t9 = -t122 * t62 + t22 * t60;
t104 = t9 * qJD(1);
t10 = t29 * t62 - t30 * t60;
t103 = t10 * qJD(1);
t79 = pkin(3) / 0.2e1 + qJ(5) / 0.2e1;
t11 = t54 + (t76 / 0.2e1 + t79) * t62;
t102 = t11 * qJD(1);
t101 = t122 * qJD(1);
t100 = t20 * qJD(3);
t21 = t120 - t119;
t99 = t21 * qJD(1);
t96 = t29 * qJD(3);
t25 = t30 * qJD(3);
t95 = t120 * qJD(1);
t94 = t60 * qJD(1);
t93 = t60 * qJD(3);
t91 = t62 * qJD(1);
t90 = t62 * qJD(2);
t89 = t62 * qJD(3);
t88 = t62 * qJD(5);
t64 = t67 * qJ(2);
t87 = t64 * qJD(1);
t86 = t67 * qJD(1);
t85 = t76 * qJD(3);
t84 = t15 * t94;
t83 = t15 * t91;
t82 = t22 * t91;
t81 = t60 * t50;
t28 = t60 * t91;
t80 = t69 * t91;
t73 = qJ(4) * qJD(4);
t72 = qJD(3) * qJ(4);
t49 = t119 * qJD(1);
t48 = t119 * qJD(4);
t43 = t60 * qJD(2);
t41 = t60 * qJD(4);
t16 = (-t76 / 0.2e1 + t79) * t62;
t2 = [0, 0, 0, 0, 0, t67 * qJD(2), t64 * qJD(2), -t60 * t89, t21 * qJD(3), 0, 0, 0, t69 * t89, -t69 * t93, t124, t8 * qJD(3) + t81, t9 * qJD(3) + t48, t10 * qJD(2) + t3 * qJD(3) - t22 * t50, t124, t6 * qJD(3) - t60 * t88 + t48, t5 * qJD(3) + qJD(5) * t120 - t81, t7 * qJD(2) + t1 * qJD(3) + t123 * t15; 0, 0, 0, 0, 0, t86, t87, 0, 0, 0, 0, 0, 0, 0, t125, 0, 0, t103, t125, 0, 0, t16 * qJD(3) + t106; 0, 0, 0, 0, 0, 0, 0, -t28, t99, -t93, -t89, 0, t80 - t25, -t69 * t94 + t96, (-t107 + t117) * qJD(3) - t41, t25 + t105, -t96 + t104, t110 + (-t30 * pkin(3) - t29 * qJ(4)) * qJD(3) + t30 * qJD(4), (-t107 + t115) * qJD(3) - t41 - t88, t108 - t126, -t100 + t109, t112 + t16 * qJD(2) + (-qJ(4) * t18 - t20 * t76) * qJD(3) + t20 * qJD(4) - t18 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93, t28, t49, t25 - t82, -t93, t49, -t28, -t83 + t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89, -t28, t95, t84 - t126; 0, 0, 0, 0, 0, -t86, -t87, 0, 0, 0, 0, 0, t89, -t93, -t125, -t89, t93, qJD(3) * t122 - t103 - t50, -t125, t93, t89, t11 * qJD(3) - t106 + t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, -t94, 0, -t91, t94, t101, 0, t94, t91, t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91, 0, 0, 0, -t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94; 0, 0, 0, 0, 0, 0, 0, t28, -t99, 0, 0, 0, -t90 - t80, (qJD(1) * t69 + qJD(2)) * t60, 0, t90 - t105, -t43 - t104, -qJD(2) * t122 - t110, 0, -t43 - t108, -t90 - t109, -t11 * qJD(2) - t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91, t94, 0, t91, -t94, -t101, 0, -t94, -t91, -t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t73, 0, qJD(4), qJD(5), t76 * qJD(5) + t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t72, 0, qJD(3), 0, t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, -t49, t90 + t82, 0, -t49, t28, t90 + t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, 0, 0, 0, t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), -t72, 0, -qJD(3), 0, -qJD(5) - t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t95, -t43 - t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), qJD(4) - t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t2;
