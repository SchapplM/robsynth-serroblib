% Calculate inertial parameters regressor of fixed base kinetic energy for
% S3RRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% 
% Output:
% T_reg [1x(3*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:16
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_reg = S3RRR1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_energykin_fixb_reg2_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRR1_energykin_fixb_reg2_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_energykin_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:15:59
% EndTime: 2018-11-14 10:15:59
% DurationCPUTime: 0.07s
% Computational Cost: add. (25->9), mult. (58->29), div. (0->0), fcn. (20->4), ass. (0->14)
t14 = pkin(1) * qJD(1);
t5 = qJD(1) + qJD(2);
t7 = sin(qJ(2));
t13 = t7 * t14;
t9 = cos(qJ(2));
t12 = t9 * t14;
t10 = qJD(1) ^ 2;
t8 = cos(qJ(3));
t6 = sin(qJ(3));
t4 = qJD(3) + t5;
t3 = t5 * pkin(2) + t12;
t2 = t8 * t13 + t6 * t3;
t1 = -t6 * t13 + t8 * t3;
t11 = [0, 0, 0, 0, 0, t10 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 ^ 2 / 0.2e1, t5 * t12, -t5 * t13, 0 (t7 ^ 2 / 0.2e1 + t9 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t10, 0, 0, 0, 0, 0, t4 ^ 2 / 0.2e1, t1 * t4, -t2 * t4, 0, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t11;
