% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S1R1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% qJDD [1x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1]';
% 
% Output:
% tauB_reg [6x(2*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:13
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S1R1_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(1,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'S1R1_invdynB_fixb_reg2_snew_vp: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'S1R1_invdynB_fixb_reg2_snew_vp: qJD has to be [1x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [1 1]), ...
  'S1R1_invdynB_fixb_reg2_snew_vp: qJDD has to be [1x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S1R1_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S1R1_invdynB_fixb_reg2_snew_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 09:12:57
% EndTime: 2020-06-19 09:12:57
% DurationCPUTime: 0.11s
% Computational Cost: add. (44->21), mult. (86->19), div. (0->0), fcn. (68->2), ass. (0->12)
t17 = qJD(1) ^ 2;
t13 = sin(qJ(1));
t14 = cos(qJ(1));
t10 = t14 * g(1) + t13 * g(2);
t9 = t13 * g(1) - t14 * g(2);
t4 = -t14 * t10 - t13 * t9;
t8 = t14 * qJDD(1) - t13 * t17;
t16 = -pkin(1) * t8 - t13 * g(3);
t15 = t13 * t10 - t14 * t9;
t7 = t13 * qJDD(1) + t14 * t17;
t5 = -pkin(1) * t7 + t14 * g(3);
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t7, -t8, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t8, -t7, 0, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t8, 0, -t7, 0, t16, -t5, t15, pkin(1) * t15; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t7, 0, t8, 0, t5, t16, t4, pkin(1) * t4; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t9, t10, 0, 0;];
tauB_reg = t1;
