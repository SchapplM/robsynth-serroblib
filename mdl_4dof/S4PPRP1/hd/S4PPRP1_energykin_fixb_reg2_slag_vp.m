% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4PPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,theta1]';
% 
% Output:
% T_reg [1x(4*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:40
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_reg = S4PPRP1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP1_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP1_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP1_energykin_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:39:14
% EndTime: 2018-11-14 13:39:14
% DurationCPUTime: 0.06s
% Computational Cost: add. (13->10), mult. (40->22), div. (0->0), fcn. (8->2), ass. (0->9)
t8 = qJD(2) * qJD(3);
t7 = qJD(2) ^ 2;
t6 = cos(qJ(3));
t5 = sin(qJ(3));
t4 = qJD(1) ^ 2 / 0.2e1;
t3 = qJD(3) ^ 2 / 0.2e1;
t2 = qJD(3) * qJ(4) + qJD(2) * t5;
t1 = -qJD(3) * pkin(3) - qJD(2) * t6 + qJD(4);
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4 + t7 / 0.2e1, 0, 0, 0, 0, 0, t3, t6 * t8, -t5 * t8, 0, t4 + (t5 ^ 2 / 0.2e1 + t6 ^ 2 / 0.2e1) * t7, 0, 0, 0, t3, 0, 0, -t1 * qJD(3), 0, t2 * qJD(3), t2 ^ 2 / 0.2e1 + t4 + t1 ^ 2 / 0.2e1;];
T_reg  = t9;
