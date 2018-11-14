% Calculate inertial parameters regressor of fixed base kinetic energy for
% S3PRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';
% 
% Output:
% T_reg [1x(3*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:13
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_reg = S3PRR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRR2_energykin_fixb_reg2_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRR2_energykin_fixb_reg2_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3PRR2_energykin_fixb_reg2_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:12:50
% EndTime: 2018-11-14 10:12:50
% DurationCPUTime: 0.06s
% Computational Cost: add. (15->8), mult. (44->27), div. (0->0), fcn. (20->4), ass. (0->12)
t6 = sin(qJ(2));
t11 = qJD(1) * t6;
t10 = qJD(1) * qJD(2);
t9 = qJD(1) ^ 2;
t8 = cos(qJ(2));
t7 = cos(qJ(3));
t5 = sin(qJ(3));
t4 = qJD(2) + qJD(3);
t3 = qJD(2) * pkin(2) + t8 * qJD(1);
t2 = t7 * t11 + t5 * t3;
t1 = -t5 * t11 + t7 * t3;
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t9 / 0.2e1, 0, 0, 0, 0, 0, qJD(2) ^ 2 / 0.2e1, t8 * t10, -t6 * t10, 0 (t6 ^ 2 / 0.2e1 + t8 ^ 2 / 0.2e1) * t9, 0, 0, 0, 0, 0, t4 ^ 2 / 0.2e1, t1 * t4, -t2 * t4, 0, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t12;
