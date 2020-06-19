% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S2RR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% qJDD [2x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,d1,d2]';
% 
% Output:
% tauB_reg [6x(3*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:14
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S2RR3_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(2,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR3_invdynB_fixb_reg2_snew_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR3_invdynB_fixb_reg2_snew_vp: qJD has to be [2x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [2 1]), ...
  'S2RR3_invdynB_fixb_reg2_snew_vp: qJDD has to be [2x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S2RR3_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S2RR3_invdynB_fixb_reg2_snew_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 09:14:28
% EndTime: 2020-06-19 09:14:28
% DurationCPUTime: 0.47s
% Computational Cost: add. (387->59), mult. (625->63), div. (0->0), fcn. (442->4), ass. (0->36)
t54 = qJD(1) + qJD(2);
t52 = t54 ^ 2;
t53 = qJDD(1) + qJDD(2);
t55 = sin(qJ(2));
t57 = cos(qJ(2));
t39 = t57 * t52 + t55 * t53;
t42 = t55 * t52 - t57 * t53;
t56 = sin(qJ(1));
t58 = cos(qJ(1));
t25 = t56 * t39 + t58 * t42;
t35 = pkin(3) * t39 - t57 * g(3);
t69 = pkin(3) * t42 - t55 * g(3);
t74 = pkin(2) * t25 + t56 * t35 + t58 * t69;
t60 = t58 * t39 - t56 * t42;
t73 = pkin(2) * t60 + t58 * t35 - t56 * t69;
t48 = t56 * g(1) - t58 * g(2);
t44 = qJDD(1) * pkin(1) + t48;
t49 = t58 * g(1) + t56 * g(2);
t59 = qJD(1) ^ 2;
t45 = -t59 * pkin(1) - t49;
t28 = -t57 * t44 + t55 * t45;
t29 = t55 * t44 + t57 * t45;
t62 = t55 * t28 + t57 * t29;
t21 = t57 * t28 - t55 * t29;
t63 = t58 * t21;
t70 = -t56 * t62 + t63;
t64 = t56 * t21;
t17 = t58 * t62 + t64;
t31 = -t56 * t48 - t58 * t49;
t47 = t58 * qJDD(1) - t56 * t59;
t61 = -pkin(2) * t47 - t56 * g(3);
t30 = t58 * t48 - t56 * t49;
t46 = t56 * qJDD(1) + t58 * t59;
t37 = -pkin(2) * t46 + t58 * g(3);
t18 = pkin(1) * g(3) + pkin(3) * t62;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t46, -t47, 0, t31, 0, 0, 0, 0, 0, 0, -t60, t25, 0, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t47, -t46, 0, t30, 0, 0, 0, 0, 0, 0, -t25, -t60, 0, -t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t47, 0, -t46, 0, t61, -t37, -t30, -pkin(2) * t30, 0, 0, -t25, 0, -t60, 0, t74, t73, t70, pkin(2) * t70 + pkin(3) * t63 - t56 * t18; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t46, 0, t47, 0, t37, t61, t31, pkin(2) * t31, 0, 0, t60, 0, -t25, 0, -t73, t74, t17, pkin(2) * t17 + pkin(3) * t64 + t58 * t18; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t48, t49, 0, 0, 0, 0, 0, 0, 0, t53, -pkin(1) * t42 - t28, -pkin(1) * t39 - t29, 0, -pkin(1) * t21;];
tauB_reg = t1;
