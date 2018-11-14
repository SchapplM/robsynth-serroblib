% Calculate minimal parameter regressor of coriolis matrix for
% S3PPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d3]';
% 
% Output:
% cmat_reg [(3*%NQJ)%x5]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:11
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function cmat_reg = S3PPR2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PPR2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PPR2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PPR2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:11:11
% EndTime: 2018-11-14 10:11:11
% DurationCPUTime: 0.04s
% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
t4 = sin(qJ(3)) * qJD(3);
t3 = cos(qJ(3)) * qJD(3);
t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, -t3, t4; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, -t4, -t3; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
cmat_reg  = t1;
