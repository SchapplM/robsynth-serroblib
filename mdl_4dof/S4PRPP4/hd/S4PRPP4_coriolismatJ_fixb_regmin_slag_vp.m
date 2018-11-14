% Calculate minimal parameter regressor of coriolis matrix for
% S4PRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta3]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x8]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:10
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function cmat_reg = S4PRPP4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:09:30
% EndTime: 2018-11-14 14:09:30
% DurationCPUTime: 0.07s
% Computational Cost: add. (51->13), mult. (104->20), div. (0->0), fcn. (137->4), ass. (0->10)
t14 = cos(qJ(2));
t13 = sin(qJ(2));
t10 = cos(pkin(5));
t9 = sin(pkin(5));
t4 = -t10 * t13 - t9 * t14;
t12 = t4 * qJD(2);
t5 = t9 * pkin(2) + qJ(4);
t11 = t5 * qJD(2);
t3 = -t10 * t14 + t9 * t13;
t1 = [0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t13 * qJD(2), -t14 * qJD(2) (t10 * t4 - t9 * t3) * qJD(2) * pkin(2), t12, -t3 * qJD(2) (-t3 * t5 - t4 * (-t10 * pkin(2) - pkin(3))) * qJD(2) - t4 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, -t12; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, qJD(4), t5 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, qJD(2), t11; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -qJD(2), -t11; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t1;
