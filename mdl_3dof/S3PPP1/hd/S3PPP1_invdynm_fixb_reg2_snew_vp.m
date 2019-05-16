% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S3PPP1
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
%   pkin=[a2,a3,theta1]';
%
% Output:
% m_new_reg [(3*4)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:20
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S3PPP1_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PPP1_invdynm_fixb_reg2_snew_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PPP1_invdynm_fixb_reg2_snew_vp: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3PPP1_invdynm_fixb_reg2_snew_vp: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3PPP1_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PPP1_invdynm_fixb_reg2_snew_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:20:09
% EndTime: 2019-05-04 18:20:10
% DurationCPUTime: 0.25s
% Computational Cost: add. (208->60), mult. (231->40), div. (0->0), fcn. (190->2), ass. (0->29)
t55 = sin(pkin(3));
t56 = cos(pkin(3));
t48 = t55 * g(1) - t56 * g(2);
t46 = -qJDD(2) + t48;
t39 = t55 * t46;
t49 = t56 * g(1) + t55 * g(2);
t47 = -qJDD(3) + t49;
t69 = -t56 * t47 - t39;
t44 = t56 * t49;
t68 = -t39 - t44;
t43 = t55 * t49;
t62 = t56 * t46;
t67 = -t62 + t43;
t66 = -t55 * t47 + t62;
t65 = pkin(2) * t46;
t64 = pkin(2) * t47;
t54 = g(3) - qJDD(1);
t50 = t55 * t54;
t61 = t56 * t54;
t60 = qJ(2) * t54;
t59 = qJ(3) * t46;
t58 = t56 * t48 - t43;
t57 = -t55 * t48 - t44;
t45 = pkin(1) * t46;
t38 = t60 - t65;
t37 = -t64 + (pkin(1) + qJ(3)) * t54;
t36 = -qJ(2) * t49 + t45;
t35 = -qJ(2) * t47 + t45 + t59;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t50, -t61, -t58, -qJ(1) * t58, 0, 0, 0, 0, 0, 0, t67, t50, t61, qJ(1) * t67 + (-t55 * pkin(1) + t56 * qJ(2)) * t54, 0, 0, 0, 0, 0, 0, -t66, t61, -t50, -qJ(1) * t66 - t55 * t37 + t56 * t38; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t61, -t50, t57, qJ(1) * t57, 0, 0, 0, 0, 0, 0, t68, -t61, t50, qJ(1) * t68 + (t56 * pkin(1) + t55 * qJ(2)) * t54, 0, 0, 0, 0, 0, 0, t69, t50, t61, qJ(1) * t69 + t56 * t37 + t55 * t38; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t48, t49, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, -t49, t36, 0, 0, 0, 0, 0, 0, 0, -t47, t46, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, -t48, 0, 0, 0, 0, 0, 0, 0, -t46, 0, t54, t60, 0, 0, 0, 0, 0, 0, -t46, t54, 0, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, -t49, 0, 0, 0, 0, 0, 0, 0, -t49, -t54, 0, pkin(1) * t54, 0, 0, 0, 0, 0, 0, -t47, 0, t54, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t49, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, -t49, t36, 0, 0, 0, 0, 0, 0, 0, -t47, t46, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, -t49, 0, 0, 0, 0, 0, 0, 0, 0, -t47, t46, t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, -t54, 0, 0, 0, 0, 0, 0, 0, t46, -t54, 0, t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t54, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, -t54, -qJ(3) * t54 + t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, t46, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, -t54, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, t54, 0, 0;];
m_new_reg  = t1;
