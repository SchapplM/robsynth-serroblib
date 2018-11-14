% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4PPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3]';
% 
% Output:
% T_reg [1x(4*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:59
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_reg = S4PPRP3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP3_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP3_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PPRP3_energykin_fixb_reg2_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:58:25
% EndTime: 2018-11-14 13:58:25
% DurationCPUTime: 0.06s
% Computational Cost: add. (16->9), mult. (42->19), div. (0->0), fcn. (16->2), ass. (0->10)
t8 = sin(qJ(3));
t9 = cos(qJ(3));
t4 = t9 * qJD(1) + t8 * qJD(2);
t10 = t4 * qJD(3);
t3 = -t8 * qJD(1) + t9 * qJD(2);
t7 = qJD(1) ^ 2 / 0.2e1;
t6 = qJD(3) ^ 2 / 0.2e1;
t2 = t4 ^ 2 / 0.2e1;
t1 = qJD(3) * pkin(3) + t3;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 + qJD(2) ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t6, t3 * qJD(3), -t10, 0, t2 + t3 ^ 2 / 0.2e1, 0, 0, 0, 0, 0, t6, t1 * qJD(3), -t10, 0, t2 + t1 ^ 2 / 0.2e1 + qJD(4) ^ 2 / 0.2e1;];
T_reg  = t5;
