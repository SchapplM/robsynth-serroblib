% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% T_reg [1x10]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRRP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP1_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP1_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP1_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:36:03
% EndTime: 2019-03-08 18:36:03
% DurationCPUTime: 0.02s
% Computational Cost: add. (28->10), mult. (49->25), div. (0->0), fcn. (18->4), ass. (0->12)
t53 = pkin(1) * qJD(1);
t45 = qJD(1) + qJD(2);
t52 = sin(qJ(2)) * t53;
t51 = cos(qJ(2)) * t53;
t43 = t45 * pkin(2) + t51;
t46 = sin(qJ(3));
t48 = cos(qJ(3));
t50 = t48 * t43 - t46 * t52;
t44 = qJD(3) + t45;
t41 = t46 * t43 + t48 * t52;
t40 = t44 * pkin(3) + t50;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t45 ^ 2 / 0.2e1, t45 * t51, -t45 * t52, t44 ^ 2 / 0.2e1, t50 * t44, -t41 * t44, t41 ^ 2 / 0.2e1 + t40 ^ 2 / 0.2e1 + qJD(4) ^ 2 / 0.2e1;];
T_reg  = t1;
