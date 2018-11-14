% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4PPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d4,theta2]';
% 
% Output:
% T_reg [1x6]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:05
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_reg = S4PPPR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR5_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR5_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR5_energykin_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:05:11
% EndTime: 2018-11-14 14:05:11
% DurationCPUTime: 0.02s
% Computational Cost: add. (10->7), mult. (33->19), div. (0->0), fcn. (12->4), ass. (0->9)
t39 = qJD(1) ^ 2 / 0.2e1;
t32 = sin(pkin(5));
t38 = t32 ^ 2 * t39 + qJD(2) ^ 2 / 0.2e1;
t37 = qJD(1) * t32;
t35 = cos(qJ(4));
t34 = sin(qJ(4));
t33 = cos(pkin(5));
t29 = -t33 * qJD(1) + qJD(3);
t1 = [t39, t33 ^ 2 * t39 + t38, t29 ^ 2 / 0.2e1 + t38, qJD(4) ^ 2 / 0.2e1 (t35 * t29 - t34 * t37) * qJD(4) -(t34 * t29 + t35 * t37) * qJD(4);];
T_reg  = t1;
