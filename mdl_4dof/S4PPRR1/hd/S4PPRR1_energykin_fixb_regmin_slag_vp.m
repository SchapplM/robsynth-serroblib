% Calculate minimal parameter regressor of fixed base kinetic energy for
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
% T_reg [1x8]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4PPRR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR1_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR1_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR1_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:15:48
% EndTime: 2019-03-08 18:15:48
% DurationCPUTime: 0.02s
% Computational Cost: add. (10->7), mult. (26->19), div. (0->0), fcn. (10->4), ass. (0->10)
t39 = sin(qJ(3));
t43 = qJD(2) * t39;
t42 = qJD(2) * qJD(3);
t41 = cos(qJ(3));
t40 = cos(qJ(4));
t38 = sin(qJ(4));
t37 = qJD(1) ^ 2 / 0.2e1;
t36 = qJD(3) + qJD(4);
t35 = qJD(3) * pkin(3) + t41 * qJD(2);
t1 = [t37, t37 + qJD(2) ^ 2 / 0.2e1, qJD(3) ^ 2 / 0.2e1, t41 * t42, -t39 * t42, t36 ^ 2 / 0.2e1 (t40 * t35 - t38 * t43) * t36 -(t38 * t35 + t40 * t43) * t36;];
T_reg  = t1;
