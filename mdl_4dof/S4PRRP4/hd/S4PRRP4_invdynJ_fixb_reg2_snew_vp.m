% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4PRRP4
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4PRRP4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP4_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP4_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP4_invdynJ_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:27:58
% EndTime: 2019-12-31 16:27:59
% DurationCPUTime: 0.36s
% Computational Cost: add. (380->88), mult. (811->103), div. (0->0), fcn. (448->6), ass. (0->63)
t58 = cos(qJ(3));
t76 = qJD(2) * qJD(3);
t56 = sin(qJ(3));
t78 = t56 * qJDD(2);
t97 = 0.2e1 * t58 * t76 + t78;
t60 = qJD(3) ^ 2;
t48 = t56 ^ 2;
t61 = qJD(2) ^ 2;
t91 = t48 * t61;
t35 = t60 + t91;
t39 = t56 * t61 * t58;
t34 = qJDD(3) - t39;
t86 = t58 * t34;
t96 = pkin(5) * (-t56 * t35 + t86);
t49 = t58 ^ 2;
t90 = t49 * t61;
t94 = t86 + t56 * (-t60 + t90);
t50 = -g(3) + qJDD(1);
t44 = t58 * t50;
t82 = t56 * qJ(4);
t68 = -t58 * pkin(3) - t82;
t81 = t61 * t68;
t53 = sin(pkin(6));
t54 = cos(pkin(6));
t31 = t53 * g(1) - t54 * g(2);
t32 = -t54 * g(1) - t53 * g(2);
t57 = sin(qJ(2));
t59 = cos(qJ(2));
t66 = -t57 * t31 - t59 * t32;
t79 = qJDD(2) * pkin(5);
t9 = -t61 * pkin(2) - t66 + t79;
t3 = -qJDD(3) * pkin(3) - t60 * qJ(4) + (t9 + t81) * t56 + qJDD(4) - t44;
t92 = 2 * qJD(4);
t6 = t56 * t50 + t58 * t9;
t33 = qJDD(3) + t39;
t89 = t56 * t33;
t74 = t56 * t76;
t77 = t58 * qJDD(2);
t27 = -0.2e1 * t74 + t77;
t37 = -t60 - t90;
t85 = pkin(5) * (t58 * t37 - t89) + pkin(2) * t27;
t83 = t48 + t49;
t29 = t83 * t61;
t84 = pkin(2) * t29 + t83 * t79;
t80 = qJD(2) * t56;
t5 = t56 * t9 - t44;
t75 = t56 * t5 + t58 * t6;
t72 = t59 * t31 - t57 * t32;
t69 = qJDD(3) * qJ(4) + (qJD(3) * t92) + t58 * t81 + t6;
t67 = t56 * t27 + t58 * t97;
t65 = t56 * t34 + t58 * t35;
t8 = -qJDD(2) * pkin(2) - t61 * pkin(5) - t72;
t64 = pkin(2) - t68;
t26 = -t74 + t77;
t63 = -t26 * pkin(3) - qJ(4) * t97 + t8;
t62 = t80 * t92 - t63;
t30 = (t48 - t49) * t61;
t14 = t89 + t58 * (t60 - t91);
t13 = t58 * t33 + t56 * t37;
t12 = t97 * t56;
t11 = (t26 - t74) * t58;
t2 = -t60 * pkin(3) + t69;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, 0, 0, 0, 0, 0, t13, -t65, 0, -t58 * t5 + t56 * t6, 0, 0, 0, 0, 0, 0, t13, 0, t65, t56 * t2 - t58 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t72, t66, 0, 0, t12, t67, t14, t11, t94, 0, -t58 * t8 + t85, -pkin(2) * t97 + t56 * t8 - t96, t75 + t84, -pkin(2) * t8 + pkin(5) * t75, t12, t14, -t67, 0, -t94, t11, t27 * t82 + t58 * ((t27 - t74) * pkin(3) + t62) + t85, t58 * ((t29 - t60) * pkin(3) + t69) + (qJ(4) * t29 + t3) * t56 + t84, t56 * (-pkin(3) * t74 + t62) + t96 + t64 * t97, pkin(5) * (t58 * t2 + t56 * t3) - t64 * ((pkin(3) * qJD(3) - (2 * qJD(4))) * t80 + t63); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, t30, t78, t39, t77, qJDD(3), -t5, -t6, 0, 0, -t39, t78, -t30, qJDD(3), -t77, t39, pkin(3) * t33 + qJ(4) * t37 - t3, (-pkin(3) * t56 + qJ(4) * t58) * qJDD(2), qJ(4) * t34 + (t35 - t60) * pkin(3) + t69, -pkin(3) * t3 + qJ(4) * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t78, -t35, t3;];
tauJ_reg = t1;
