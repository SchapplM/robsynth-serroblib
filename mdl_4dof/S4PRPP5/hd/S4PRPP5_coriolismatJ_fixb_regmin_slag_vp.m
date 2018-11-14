% Calculate minimal parameter regressor of coriolis matrix for
% S4PRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x10]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:11
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function cmat_reg = S4PRPP5_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP5_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP5_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PRPP5_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:10:23
% EndTime: 2018-11-14 14:10:23
% DurationCPUTime: 0.09s
% Computational Cost: add. (13->13), mult. (22->10), div. (0->0), fcn. (14->2), ass. (0->9)
t8 = cos(qJ(2));
t7 = sin(qJ(2));
t6 = qJ(3) * qJD(3);
t5 = qJD(2) * qJ(3);
t4 = t8 * qJ(3);
t3 = t8 * qJD(2);
t2 = t7 * qJD(2);
t1 = t7 * qJD(3);
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t2, -t3, -t2, t3 (-t7 * pkin(2) + t4) * qJD(2) + t1, -t2, t3 (t4 + t7 * (-pkin(2) - pkin(3))) * qJD(2) + t1; 0, 0, 0, 0, 0, 0, t2, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, qJD(3), t6, 0, qJD(3), t6; 0, 0, 0, 0, 0, qJD(2), t5, 0, qJD(2), t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -qJD(2), -t5, 0, -qJD(2), -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t9;
