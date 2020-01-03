% Calculate minimal parameter regressor of coriolis matrix for
% S5RPPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x21]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPPRP6_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP6_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:55:18
% EndTime: 2019-12-31 17:55:20
% DurationCPUTime: 0.52s
% Computational Cost: add. (832->90), mult. (1433->99), div. (0->0), fcn. (1571->4), ass. (0->68)
t65 = sin(pkin(7));
t68 = sin(qJ(4));
t100 = cos(qJ(4));
t66 = cos(pkin(7));
t71 = t100 * t66;
t49 = t68 * t65 - t71;
t99 = t68 * t66;
t50 = t100 * t65 + t99;
t103 = t50 * pkin(4) + t49 * qJ(5);
t80 = t50 * qJD(5);
t104 = t103 * qJD(4) - t80;
t46 = t49 ^ 2;
t47 = t50 ^ 2;
t67 = -pkin(1) - qJ(3);
t101 = -pkin(6) + t67;
t61 = t65 ^ 2;
t62 = t66 ^ 2;
t56 = t61 + t62;
t13 = t49 * pkin(4) - t50 * qJ(5);
t57 = t65 * pkin(3) + qJ(2);
t19 = t57 + t103;
t1 = t19 * t13;
t98 = t1 * qJD(1);
t7 = -t13 * t49 + t19 * t50;
t97 = t7 * qJD(1);
t8 = -t13 * t50 - t19 * t49;
t96 = t8 * qJD(1);
t53 = t101 * t65;
t22 = -t101 * t71 + t68 * t53;
t23 = t100 * t53 + t101 * t99;
t9 = -t22 * t49 - t23 * t50;
t95 = t9 * qJD(1);
t69 = -t46 / 0.2e1 - t47 / 0.2e1;
t11 = -0.1e1 / 0.2e1 + t69;
t94 = t11 * qJD(1);
t93 = t13 * qJD(1);
t18 = t47 - t46;
t92 = t18 * qJD(1);
t91 = t19 * qJD(1);
t90 = t22 * qJD(4);
t20 = t23 * qJD(4);
t25 = t47 + t46;
t89 = t25 * qJD(1);
t88 = t46 * qJD(1);
t48 = t56 * t67;
t87 = t48 * qJD(1);
t86 = t49 * qJD(1);
t85 = t49 * qJD(2);
t84 = t49 * qJD(4);
t83 = t49 * qJD(5);
t82 = t50 * qJD(1);
t81 = t50 * qJD(4);
t72 = -t61 / 0.2e1 - t62 / 0.2e1;
t55 = -0.1e1 / 0.2e1 + t72;
t79 = t55 * qJD(1);
t78 = t56 * qJD(1);
t77 = t65 * qJD(1);
t76 = t66 * qJD(1);
t75 = qJD(4) * qJ(5);
t74 = t49 * t82;
t73 = t57 * t86;
t64 = qJ(2) * qJD(2);
t63 = qJD(1) * qJ(2);
t54 = 0.1e1 / 0.2e1 + t72;
t43 = qJD(2) * t50;
t35 = t49 * qJD(3);
t10 = 0.1e1 / 0.2e1 + t69;
t2 = [0, 0, 0, 0, qJD(2), t64, qJD(2) * t65, qJD(2) * t66, t56 * qJD(3), -t48 * qJD(3) + t64, t49 * t81, t18 * qJD(4), 0, 0, 0, -t57 * t84 + t43, -t57 * t81 - t85, t8 * qJD(4) + t49 * t80 + t43, t25 * qJD(3), t7 * qJD(4) + t46 * qJD(5) + t85, t9 * qJD(3) - t1 * qJD(4) + (qJD(2) + t83) * t19; 0, 0, 0, 0, qJD(1), t63, t77, t76, 0, t54 * qJD(3) + t63, 0, 0, 0, 0, 0, t82, -t86, t82, 0, t86, t10 * qJD(3) + t91; 0, 0, 0, 0, 0, 0, 0, 0, t78, t54 * qJD(2) - t87, 0, 0, 0, 0, 0, 0, 0, 0, t89, 0, t10 * qJD(2) + t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, t92, -t81, t84, 0, -t20 - t73, -t57 * t82 + t90, -t20 + t96, t104, -t90 + t97, -t98 + (-t23 * pkin(4) - t22 * qJ(5)) * qJD(4) + t23 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, -t81, t88, t19 * t86 + t20; 0, 0, 0, 0, -qJD(1), -t63, -t77, -t76, 0, t55 * qJD(3) - t63, 0, 0, 0, 0, 0, -t82, t86, -t82, 0, -t86, t11 * qJD(3) - t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, t84, -t81, 0, -t84, -t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81; 0, 0, 0, 0, 0, 0, 0, 0, -t78, -t55 * qJD(2) + t87, 0, 0, 0, 0, 0, -t84, -t81, -t84, -t89, t81, -t11 * qJD(2) - t13 * qJD(4) + t83 - t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, -t82, -t86, 0, t82, -t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, -t92, 0, 0, 0, t35 + t73, (qJD(1) * t57 + qJD(3)) * t50, t35 - t96, 0, -t50 * qJD(3) - t97, t13 * qJD(3) + t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, t82, t86, 0, -t82, t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), qJ(5) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, 0, -t88, (-qJD(3) - t91) * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4), -t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t2;
