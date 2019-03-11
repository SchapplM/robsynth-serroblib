% Calculate minimal parameter regressor of coriolis matrix for
% S2RR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
% 
% Output:
% cmat_reg [(2*%NQJ)%x10]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S2RR2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:01:01
% EndTime: 2019-03-08 18:01:01
% DurationCPUTime: 0.05s
% Computational Cost: add. (6->5), mult. (21->10), div. (0->0), fcn. (16->2), ass. (0->8)
t2 = sin(qJ(2));
t3 = cos(qJ(2));
t1 = -t2 ^ 2 + t3 ^ 2;
t7 = t1 * qJD(1);
t6 = t2 * qJD(2);
t5 = t3 * qJD(2);
t4 = t2 * t3 * qJD(1);
t8 = [0, 0, 0, t2 * t5, t1 * qJD(2), 0, 0, 0, 0, 0; 0, 0, 0, t4, t7, t5, -t6, 0, -pkin(1) * t5, pkin(1) * t6; 0, 0, 0, -t4, -t7, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t8;
