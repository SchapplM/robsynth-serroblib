% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% T_reg [1x12]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:53
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_reg = S4RRPP2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP2_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP2_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP2_energykin_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:52:30
% EndTime: 2018-11-14 13:52:30
% DurationCPUTime: 0.02s
% Computational Cost: add. (37->12), mult. (50->22), div. (0->0), fcn. (10->2), ass. (0->11)
t39 = qJD(1) + qJD(2);
t45 = pkin(1) * qJD(1);
t44 = sin(qJ(2)) * t45;
t38 = t39 * qJ(3) + t44;
t46 = t38 * t39;
t43 = cos(qJ(2)) * t45;
t42 = qJD(3) - t43;
t37 = -t39 * pkin(2) + t42;
t36 = t38 ^ 2 / 0.2e1;
t35 = (-pkin(2) - pkin(3)) * t39 + t42;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t39 ^ 2 / 0.2e1, t39 * t43, -t39 * t44, -t37 * t39, t46, t36 + t37 ^ 2 / 0.2e1, -t35 * t39, t46, t36 + t35 ^ 2 / 0.2e1 + qJD(4) ^ 2 / 0.2e1;];
T_reg  = t1;
