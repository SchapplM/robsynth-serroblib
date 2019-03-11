% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4PRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta3]';
% 
% Output:
% T_reg [1x8]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:22
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4PRPR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR2_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR2_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR2_energykin_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:21:49
% EndTime: 2019-03-08 18:21:49
% DurationCPUTime: 0.02s
% Computational Cost: add. (23->11), mult. (54->29), div. (0->0), fcn. (30->6), ass. (0->14)
t47 = sin(qJ(2));
t51 = qJD(1) * t47;
t50 = qJD(1) * qJD(2);
t49 = cos(qJ(2));
t42 = qJD(2) * pkin(2) + t49 * qJD(1);
t44 = sin(pkin(6));
t45 = cos(pkin(6));
t39 = t45 * t42 - t44 * t51;
t48 = cos(qJ(4));
t46 = sin(qJ(4));
t43 = qJD(2) + qJD(4);
t40 = t44 * t42 + t45 * t51;
t38 = qJD(2) * pkin(3) + t39;
t1 = [qJD(1) ^ 2 / 0.2e1, qJD(2) ^ 2 / 0.2e1, t49 * t50, -t47 * t50, t40 ^ 2 / 0.2e1 + t39 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t43 ^ 2 / 0.2e1 (t48 * t38 - t46 * t40) * t43 -(t46 * t38 + t48 * t40) * t43;];
T_reg  = t1;
