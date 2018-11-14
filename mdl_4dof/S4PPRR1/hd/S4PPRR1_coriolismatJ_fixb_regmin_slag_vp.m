% Calculate minimal parameter regressor of coriolis matrix for
% S4PPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x8]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:40
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function cmat_reg = S4PPRR1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:40:19
% EndTime: 2018-11-14 13:40:19
% DurationCPUTime: 0.07s
% Computational Cost: add. (18->9), mult. (42->17), div. (0->0), fcn. (42->4), ass. (0->11)
t12 = pkin(3) * qJD(4);
t11 = qJD(3) * pkin(3);
t10 = qJD(3) + qJD(4);
t9 = pkin(3) * t10;
t8 = cos(qJ(3));
t7 = cos(qJ(4));
t6 = sin(qJ(3));
t5 = sin(qJ(4));
t2 = t10 * (t5 * t6 - t7 * t8);
t1 = t10 * (-t5 * t8 - t7 * t6);
t3 = [0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t6 * qJD(3), -t8 * qJD(3), 0, t1, t2; 0, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -t5 * t12, -t7 * t12; 0, 0, 0, 0, 0, 0, -t5 * t9, -t7 * t9; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t5 * t11, t7 * t11; 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t3;
