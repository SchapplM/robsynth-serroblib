% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S3PPR1
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
%   pkin=[a2,a3,d3]';
%
% Output:
% m_new_reg [(3*4)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:21
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S3PPR1_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PPR1_invdynm_fixb_reg2_snew_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PPR1_invdynm_fixb_reg2_snew_vp: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3PPR1_invdynm_fixb_reg2_snew_vp: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3PPR1_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PPR1_invdynm_fixb_reg2_snew_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:21:23
% EndTime: 2019-05-04 18:21:23
% DurationCPUTime: 0.31s
% Computational Cost: add. (318->68), mult. (364->41), div. (0->0), fcn. (252->2), ass. (0->32)
t63 = sin(qJ(3));
t64 = cos(qJ(3));
t81 = qJD(3) ^ 2;
t53 = qJDD(3) * t64 - t63 * t81;
t82 = -pkin(3) * t53 + t63 * g(3);
t43 = pkin(1) * t53 - t82;
t80 = pkin(2) * g(3);
t60 = g(2) - qJDD(1);
t79 = pkin(1) * t60;
t61 = g(1) - qJDD(2);
t47 = -t60 * t63 + t64 * t61;
t48 = -t64 * t60 - t63 * t61;
t40 = t47 * t64 - t63 * t48;
t78 = pkin(2) * t40;
t52 = t63 * qJDD(3) + t81 * t64;
t77 = pkin(2) * t52;
t76 = pkin(2) * t53;
t75 = qJ(1) * g(3);
t74 = pkin(2) + qJ(1);
t67 = t63 * t47 + t48 * t64;
t72 = qJ(2) * t67;
t71 = qJ(2) * t60;
t37 = pkin(3) * t67;
t70 = pkin(1) * t67 + t37;
t46 = pkin(3) * t52 + g(3) * t64;
t68 = -qJ(2) * t53 + t48;
t66 = pkin(1) * t52 + t46;
t65 = -qJ(2) * t52 + t47;
t62 = qJ(2) * g(3);
t54 = pkin(1) * t61 + t62;
t36 = t62 - (-pkin(1) - pkin(3)) * t40;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t60, 0, -g(3), -t75, 0, 0, 0, 0, 0, 0, -g(3), t60, 0, -t75 - t79, 0, 0, t52, 0, t53, 0, -t66, -t43, t67, -t74 * g(3) + t70; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, g(1), -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, g(3), t54, 0, 0, t53, 0, -t52, 0, -t43, t66, t40, t36; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, g(1), qJ(1) * g(1), 0, 0, 0, 0, 0, 0, t61, 0, -t60, qJ(1) * t61 - t71, 0, 0, 0, 0, 0, -qJDD(3), -t74 * t53 + t65, t74 * t52 + t68, 0, t40 * t74 + t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, -g(1), 0, 0, 0, 0, 0, 0, 0, -t61, 0, t60, t71, 0, 0, 0, 0, 0, qJDD(3), -t65 + t76, -t68 - t77, 0, -t72 - t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, g(3), 0, 0, 0, 0, 0, 0, 0, g(3), -t60, 0, t79, 0, 0, -t52, 0, -t53, 0, t66, t43, -t67, -t70 + t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1), -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, g(3), t54, 0, 0, t53, 0, -t52, 0, -t43, t66, t40, t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, g(3), 0, 0, 0, t53, 0, -t52, 0, t82, t46, t40, pkin(3) * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, 0, -t60, 0, 0, 0, 0, 0, 0, -qJDD(3), t47 - t76, t48 + t77, 0, t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), t60, 0, 0, 0, 0, t52, 0, t53, 0, -t46, t82, t67, t37 - t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), 0, -t81, 0, 0, g(3), t47, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, 0, qJDD(3), 0, -g(3), 0, t48, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), -t47, -t48, 0, 0;];
m_new_reg  = t1;
