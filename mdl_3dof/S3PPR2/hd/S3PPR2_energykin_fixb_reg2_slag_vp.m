% Calculate inertial parameters regressor of fixed base kinetic energy for
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
% T_reg [1x(3*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:11
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_reg = S3PPR2_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PPR2_energykin_fixb_reg2_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PPR2_energykin_fixb_reg2_slag_vp: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PPR2_energykin_fixb_reg2_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:11:10
% EndTime: 2018-11-14 10:11:10
% DurationCPUTime: 0.05s
% Computational Cost: add. (7->5), mult. (22->15), div. (0->0), fcn. (8->2), ass. (0->6)
t5 = cos(qJ(3));
t4 = sin(qJ(3));
t3 = qJD(1) ^ 2 / 0.2e1;
t2 = qJD(1) * t5 + qJD(2) * t4;
t1 = -qJD(1) * t4 + qJD(2) * t5;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3 + qJD(2) ^ 2 / 0.2e1, 0, 0, 0, 0, 0, qJD(3) ^ 2 / 0.2e1, t1 * qJD(3), -t2 * qJD(3), 0, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1;];
T_reg  = t6;
