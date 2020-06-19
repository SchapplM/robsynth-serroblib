% Calculate minimal parameter regressor of coriolis matrix for
% S2RR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,d1,d2]';
% 
% Output:
% cmat_reg [(2*%NQJ)%x6]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:14
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S2RR3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S2RR3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 09:14:27
% EndTime: 2020-06-19 09:14:28
% DurationCPUTime: 0.09s
% Computational Cost: add. (4->3), mult. (16->9), div. (0->0), fcn. (8->2), ass. (0->6)
t5 = pkin(1) * qJD(1);
t4 = pkin(1) * qJD(2);
t3 = pkin(1) * (-qJD(1) - qJD(2));
t2 = cos(qJ(2));
t1 = sin(qJ(2));
t6 = [0, 0, 0, 0, -t1 * t4, -t2 * t4; 0, 0, 0, 0, t1 * t3, t2 * t3; 0, 0, 0, 0, t1 * t5, t2 * t5; 0, 0, 0, 0, 0, 0;];
cmat_reg = t6;
