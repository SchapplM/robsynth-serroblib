% Calculate inertial parameters regressor of coriolis matrix for
% S4PRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% 
% Output:
% cmat_reg [(4*4)x(4*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4PRPR6_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR6_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:44
% EndTime: 2019-12-31 16:24:45
% DurationCPUTime: 0.46s
% Computational Cost: add. (555->64), mult. (1429->114), div. (0->0), fcn. (1527->6), ass. (0->67)
t62 = sin(pkin(7));
t60 = t62 ^ 2;
t63 = cos(pkin(7));
t61 = t63 ^ 2;
t53 = t61 + t60;
t94 = cos(qJ(4));
t76 = t94 * t63;
t64 = sin(qJ(4));
t93 = t64 * t62;
t69 = t76 - t93;
t77 = t94 * t62;
t92 = t64 * t63;
t44 = t77 + t92;
t98 = t44 ^ 2;
t97 = -t69 / 0.2e1;
t96 = -t44 / 0.2e1;
t65 = sin(qJ(2));
t95 = t65 / 0.2e1;
t91 = pkin(5) + qJ(3);
t31 = t44 * t65;
t66 = cos(qJ(2));
t32 = t66 * t44;
t33 = t69 * t65;
t34 = t69 * t66;
t58 = t65 * t66;
t8 = t31 * t32 + t33 * t34 - t58;
t90 = t8 * qJD(1);
t40 = t69 ^ 2;
t16 = t40 - t98;
t89 = t16 * qJD(2);
t68 = -t92 / 0.2e1 - t77 / 0.2e1;
t17 = (t44 / 0.2e1 + t68) * t66;
t88 = t17 * qJD(1);
t67 = -t76 / 0.2e1 + t93 / 0.2e1;
t18 = (t69 / 0.2e1 + t67) * t66;
t87 = t18 * qJD(1);
t27 = t40 + t98;
t86 = t27 * qJD(2);
t74 = t53 * t66;
t30 = t65 * t74 - t58;
t85 = t30 * qJD(1);
t84 = t69 * qJD(2);
t39 = t69 * qJD(4);
t83 = t44 * qJD(2);
t82 = t44 * qJD(4);
t81 = t53 * qJD(2);
t80 = t65 * qJD(2);
t79 = t69 * t83;
t78 = t69 * t82;
t75 = t91 * t62;
t46 = t53 * qJ(3);
t57 = -pkin(3) * t63 - pkin(2);
t73 = qJD(2) * t57 + qJD(3);
t70 = t31 * t96 + t33 * t97;
t10 = t95 + t70;
t51 = t91 * t63;
t28 = t64 * t51 + t75 * t94;
t29 = t51 * t94 - t64 * t75;
t7 = t28 * t44 + t29 * t69;
t72 = qJD(1) * t10 - qJD(2) * t7;
t35 = (0.1e1 / 0.2e1 - t61 / 0.2e1 - t60 / 0.2e1) * t65;
t71 = qJD(1) * t35 - qJD(2) * t46;
t36 = (0.1e1 + t53) * t95;
t20 = (t68 + t96) * t66;
t19 = (t67 + t97) * t66;
t11 = t95 - t70;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t8 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80, -t66 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, -t63 * t80, t62 * t80, qJD(2) * t74, t85 + (-t65 * pkin(2) + t46 * t66) * qJD(2) + t36 * qJD(3), 0, 0, 0, 0, 0, 0, qJD(4) * t20 - t69 * t80, qJD(4) * t19 + t44 * t80, (t32 * t44 + t34 * t69) * qJD(2), t90 + (t28 * t32 + t29 * t34 + t57 * t65) * qJD(2) + t11 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t11 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t20 - qJD(4) * t33, qJD(2) * t19 + qJD(4) * t31, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t35 - t85, 0, 0, 0, 0, 0, 0, -t17 * qJD(4), -t18 * qJD(4), 0, -qJD(3) * t10 - t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53 * qJD(3), t46 * qJD(3), t78, t16 * qJD(4), 0, -t78, 0, 0, t57 * t82, t57 * t39, t27 * qJD(3), t7 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, -t71, 0, 0, 0, 0, 0, 0, 0, 0, t86, -t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, t89, t39, -t79, -t82, 0, -qJD(4) * t29 + t57 * t83 - t88, qJD(4) * t28 + t57 * t84 - t87, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, t71, 0, 0, 0, 0, 0, 0, t82, t39, -t86, t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, t84, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17 * qJD(2), t18 * qJD(2), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, -t89, 0, t79, 0, 0, -t44 * t73 + t88, -t69 * t73 + t87, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83, -t84, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
