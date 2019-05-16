% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S3PRP1
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% qJDD [3x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2]';
%
% Output:
% m_new_reg [(3*4)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:24
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S3PRP1_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRP1_invdynm_fixb_reg2_snew_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRP1_invdynm_fixb_reg2_snew_vp: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3PRP1_invdynm_fixb_reg2_snew_vp: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3PRP1_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PRP1_invdynm_fixb_reg2_snew_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:23:50
% EndTime: 2019-05-04 18:23:50
% DurationCPUTime: 0.41s
% Computational Cost: add. (389->84), mult. (566->58), div. (0->0), fcn. (296->2), ass. (0->41)
t79 = (qJDD(2) * qJ(3));
t96 = (qJD(3) * qJD(2));
t105 = 2 * t79 + 2 * t96;
t80 = g(2) - qJDD(1);
t81 = sin(qJ(2));
t82 = cos(qJ(2));
t69 = g(1) * t81 - t80 * t82;
t70 = g(1) * t82 + t81 * t80;
t55 = t69 * t82 - t70 * t81;
t104 = pkin(1) * t55;
t84 = qJD(2) ^ 2;
t71 = qJDD(2) * t81 + t82 * t84;
t103 = pkin(1) * t71;
t102 = pkin(3) * t55;
t101 = g(3) * t81;
t100 = qJ(1) * g(3);
t99 = qJ(3) * g(3);
t87 = t70 - (2 * t96);
t60 = pkin(2) * t84 - t79 + t87;
t88 = -qJDD(3) + t69;
t97 = qJDD(2) * pkin(2);
t62 = qJ(3) * t84 + t88 + t97;
t98 = pkin(2) * t62 - qJ(3) * t60;
t51 = t60 * t81 - t62 * t82;
t95 = pkin(3) * t51 + t82 * t99;
t94 = -t82 * t60 - t62 * t81;
t93 = -t69 * t81 - t82 * t70;
t92 = pkin(1) * t51 - t98;
t91 = -pkin(3) * t71 + g(3) * t82;
t72 = t82 * qJDD(2) - t81 * t84;
t63 = pkin(3) * t72 + t101;
t90 = -t70 + t103;
t68 = pkin(1) * t72;
t89 = t68 + t69;
t86 = qJ(1) * t72 - t90;
t85 = qJ(1) * t71 + t89;
t83 = pkin(1) * g(3);
t78 = 0.2e1 * t97;
t53 = pkin(3) * t93 + t83;
t48 = pkin(3) * t94 + t83 + (pkin(2) * t82 + qJ(3) * t81) * g(3);
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, 0, -t80, -g(3), -t100, 0, 0, t72, 0, -t71, 0, -t63, -t91, -t55, -t100 - t102, 0, t72, 0, 0, t71, 0, -t63, t51, t91, (-pkin(2) * t81 - qJ(1)) * g(3) + t95; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, g(3), g(1), 0, 0, 0, 0, t71, 0, t72, 0, t91, -t63, t93, t53, 0, t71, 0, 0, -t72, 0, t91, t94, t63, t48; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, -t80, 0, g(1), qJ(1) * g(1), 0, 0, 0, 0, 0, qJDD(2), t85, t86, 0, -qJ(1) * t93 + t104, 0, 0, 0, qJDD(2), 0, 0, -qJDD(3) + t78 + t85, 0, -t86 + t105, -qJ(1) * t94 - t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80, -g(3), 0, 0, 0, t72, 0, -t71, 0, -t63, -t91, -t55, -t102, 0, t72, 0, 0, t71, 0, -t63, t51, t91, -pkin(2) * t101 + t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, 0, -g(1), 0, 0, 0, 0, 0, 0, -qJDD(2), -t89, t90, 0, -t104, 0, 0, 0, -qJDD(2), 0, 0, -t68 - t88 - 0.2e1 * t97, 0, -(2 * t79) + t87 - t103, t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), g(1), 0, 0, 0, 0, t71, 0, t72, 0, t91, -t63, t93, t53, 0, t71, 0, 0, -t72, 0, t91, t94, t63, t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), 0, -t84, 0, 0, -g(3), -t69, 0, 0, qJDD(2), 0, 0, t84, 0, 0, -t62, g(3), t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, 0, qJDD(2), 0, g(3), 0, -t70, 0, 0, t84, 0, 0, -qJDD(2), 0, g(3), -t60, 0, pkin(2) * g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t69, t70, 0, 0, 0, 0, 0, qJDD(2), 0, 0, t78 + t88, 0, -t70 + t105, t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), 0, 0, t84, 0, 0, -t62, g(3), 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), 0, 0, t62, 0, -t60, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84, 0, 0, qJDD(2), 0, -g(3), t60, 0, 0;];
m_new_reg  = t1;
