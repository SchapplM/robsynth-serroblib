% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S1P1
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
%   pkin=[dummy]';
% 
% Output:
% tauJ_reg [1x(1*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-14 12:22
% Revision: 96facaeb42edba38506bd76ea342a8981e82f256 (2020-11-21)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S1P1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(1,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'S1P1_invdynJ_fixb_reg2_snew_vp: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'S1P1_invdynJ_fixb_reg2_snew_vp: qJD has to be [1x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [1 1]), ...
  'S1P1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [1x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S1P1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S1P1_invdynJ_fixb_reg2_snew_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-14 12:21:56
% EndTime: 2021-01-14 12:21:56
% DurationCPUTime: 0.10s
% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) + qJDD(1);];
tauJ_reg = t1;
