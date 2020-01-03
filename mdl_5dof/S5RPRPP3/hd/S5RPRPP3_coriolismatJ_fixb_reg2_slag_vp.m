% Calculate inertial parameters regressor of coriolis matrix for
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
% cmat_reg [(5*5)x(5*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRPP3_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP3_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:59
% EndTime: 2019-12-31 18:13:01
% DurationCPUTime: 1.28s
% Computational Cost: add. (1166->129), mult. (2271->124), div. (0->0), fcn. (2525->4), ass. (0->86)
t134 = cos(qJ(3));
t84 = sin(pkin(7));
t85 = cos(pkin(7));
t87 = sin(qJ(3));
t73 = t134 * t84 + t85 * t87;
t137 = t73 ^ 2;
t71 = -t134 * t85 + t84 * t87;
t138 = t71 ^ 2;
t140 = t138 - t137;
t150 = t140 * qJD(1);
t145 = t140 * qJD(3);
t132 = pkin(6) + qJ(2);
t76 = t132 * t85;
t90 = t132 * t84;
t40 = t134 * t90 + t76 * t87;
t21 = pkin(4) * t73 + t40;
t149 = qJD(3) * t21;
t139 = t138 + t137;
t148 = t139 * qJD(1);
t147 = t139 * qJD(2);
t61 = t73 * qJD(4);
t144 = qJD(5) * t71 - t61;
t65 = qJ(4) * t71;
t143 = pkin(3) * t73 + t65;
t41 = t134 * t76 - t87 * t90;
t89 = t40 * t73 - t41 * t71;
t142 = t89 * qJD(1);
t141 = t89 * qJD(2);
t135 = t71 * pkin(3);
t86 = pkin(3) + qJ(5);
t133 = t86 * t71;
t77 = t84 ^ 2 + t85 ^ 2;
t125 = t73 * qJ(4);
t79 = -t85 * pkin(2) - pkin(1);
t88 = t79 - t125;
t18 = t88 + t133;
t20 = t86 * t73 + t65;
t1 = t18 * t20;
t130 = t1 * qJD(1);
t27 = t88 + t135;
t3 = t27 * t143;
t128 = t3 * qJD(1);
t6 = t18 * t73 + t20 * t71;
t127 = t6 * qJD(1);
t7 = t18 * t71 - t20 * t73;
t126 = t7 * qJD(1);
t23 = -t71 * pkin(4) + t41;
t8 = t21 * t73 - t23 * t71;
t124 = t8 * qJD(1);
t10 = -t143 * t71 - t27 * t73;
t123 = t10 * qJD(1);
t11 = -t143 * t73 + t27 * t71;
t122 = t11 * qJD(1);
t91 = pkin(3) / 0.2e1 + qJ(5) / 0.2e1;
t14 = t65 + (t86 / 0.2e1 + t91) * t73;
t119 = t14 * qJD(1);
t118 = t143 * qJD(1);
t117 = t23 * qJD(3);
t110 = t40 * qJD(3);
t30 = t41 * qJD(3);
t109 = t138 * qJD(1);
t108 = t71 * qJD(1);
t107 = t71 * qJD(3);
t105 = t73 * qJD(1);
t104 = t73 * qJD(2);
t103 = t73 * qJD(3);
t102 = t73 * qJD(5);
t75 = t77 * qJ(2);
t101 = t75 * qJD(1);
t100 = t77 * qJD(1);
t99 = t86 * qJD(3);
t98 = t18 * t108;
t97 = t18 * t105;
t96 = t27 * t105;
t95 = t71 * t105;
t94 = t71 * t103;
t93 = t71 * t61;
t92 = t79 * t105;
t83 = qJ(4) * qJD(4);
t82 = qJD(3) * qJ(4);
t60 = t137 * qJD(1);
t59 = t137 * qJD(4);
t54 = t71 * qJD(2);
t52 = t71 * qJD(4);
t19 = (-t86 / 0.2e1 + t91) * t73;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77 * qJD(2), t75 * qJD(2), -t94, t145, 0, t94, 0, 0, t79 * t103, -t79 * t107, t147, t141, 0, 0, 0, -t94, t145, t94, t147, qJD(3) * t10 + t93, qJD(3) * t11 + t59, qJD(3) * t3 - t27 * t61 + t141, 0, 0, 0, t94, -t145, -t94, t147, qJD(3) * t7 - t102 * t71 + t59, qJD(3) * t6 + qJD(5) * t138 - t93, t8 * qJD(2) + t1 * qJD(3) + t144 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, t101, 0, 0, 0, 0, 0, 0, 0, 0, t148, t142, 0, 0, 0, 0, 0, 0, t148, 0, 0, t142, 0, 0, 0, 0, 0, 0, t148, 0, 0, qJD(3) * t19 + t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t95, t150, -t107, t95, -t103, 0, t92 - t30, -t108 * t79 + t110, 0, 0, 0, t107, t103, -t95, t150, t95, (-t125 + t135) * qJD(3) - t52, t30 + t123, -t110 + t122, t128 + (-pkin(3) * t41 - qJ(4) * t40) * qJD(3) + t41 * qJD(4), 0, t103, -t107, t95, -t150, -t95, (-t125 + t133) * qJD(3) - t52 - t102, t126 - t149, -t117 + t127, t130 + t19 * qJD(2) + (-qJ(4) * t21 - t23 * t86) * qJD(3) + t23 * qJD(4) - t21 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t107, t95, t60, t30 - t96, 0, 0, 0, 0, 0, 0, -t107, t60, -t95, -t97 + t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t103, -t95, t109, t98 - t149; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100, -t101, 0, 0, 0, 0, 0, 0, t103, -t107, -t148, -t142, 0, 0, 0, 0, 0, 0, -t148, -t103, t107, qJD(3) * t143 - t142 - t61, 0, 0, 0, 0, 0, 0, -t148, t107, t103, qJD(3) * t14 - t124 + t144; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, -t108, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t105, t108, t118, 0, 0, 0, 0, 0, 0, 0, t108, t105, t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t105, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, -t150, 0, -t95, 0, 0, -t104 - t92, (qJD(1) * t79 + qJD(2)) * t71, 0, 0, 0, 0, 0, t95, -t150, -t95, 0, t104 - t123, -t54 - t122, -qJD(2) * t143 - t128, 0, 0, 0, -t95, t150, t95, 0, -t54 - t126, -t104 - t127, -qJD(2) * t14 - t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t105, t108, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, -t108, -t118, 0, 0, 0, 0, 0, 0, 0, -t108, -t105, -t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t83, 0, 0, 0, 0, 0, 0, 0, qJD(4), qJD(5), qJD(5) * t86 + t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t82, 0, 0, 0, 0, 0, 0, 0, qJD(3), 0, t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t95, -t60, t104 + t96, 0, 0, 0, 0, 0, 0, 0, -t60, t95, t104 + t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), -t82, 0, 0, 0, 0, 0, 0, 0, -qJD(3), 0, -qJD(5) - t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, -t109, -t54 - t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), qJD(4) - t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t2;
