% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RRPR2
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRPR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR2_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_energykin_fixb_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 18:16:34
% EndTime: 2019-07-18 18:16:34
% DurationCPUTime: 0.03s
% Computational Cost: add. (41->13), mult. (50->26), div. (0->0), fcn. (14->4), ass. (0->12)
t52 = pkin(1) * qJD(1);
t44 = qJD(1) + qJD(2);
t51 = sin(qJ(2)) * t52;
t50 = cos(qJ(2)) * t52;
t49 = qJD(3) - t50;
t47 = cos(qJ(4));
t45 = sin(qJ(4));
t43 = qJD(4) - t44;
t42 = t44 * qJ(3) + t51;
t41 = -t44 * pkin(2) + t49;
t40 = (-pkin(2) - pkin(3)) * t44 + t49;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t44 ^ 2 / 0.2e1, t44 * t50, -t44 * t51, -t41 * t44, t42 * t44, t42 ^ 2 / 0.2e1 + t41 ^ 2 / 0.2e1, t43 ^ 2 / 0.2e1, (t47 * t40 - t45 * t42) * t43, -(t45 * t40 + t47 * t42) * t43;];
T_reg  = t1;
