% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4PRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% T_reg [1x10]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:45
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_reg = S4PRRR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR1_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR1_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR1_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:44:29
% EndTime: 2018-11-14 13:44:29
% DurationCPUTime: 0.02s
% Computational Cost: add. (17->7), mult. (30->20), div. (0->0), fcn. (10->4), ass. (0->9)
t46 = pkin(2) * qJD(2);
t39 = qJD(2) + qJD(3);
t45 = sin(qJ(3)) * t46;
t44 = cos(qJ(3)) * t46;
t42 = cos(qJ(4));
t40 = sin(qJ(4));
t38 = qJD(4) + t39;
t37 = t39 * pkin(3) + t44;
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, 0, 0, t39 ^ 2 / 0.2e1, t39 * t44, -t39 * t45, t38 ^ 2 / 0.2e1 (t42 * t37 - t40 * t45) * t38 -(t40 * t37 + t42 * t45) * t38;];
T_reg  = t1;
