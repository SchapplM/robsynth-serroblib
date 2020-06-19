% Calculate inertial parameters regressor of fixed base kinetic energy for
% S2RR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,d1,d2]';
% 
% Output:
% T_reg [1x(2*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:14
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S2RR3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR3_energykin_fixb_reg2_slag_vp: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR3_energykin_fixb_reg2_slag_vp: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S2RR3_energykin_fixb_reg2_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 09:14:26
% EndTime: 2020-06-19 09:14:26
% DurationCPUTime: 0.12s
% Computational Cost: add. (5->3), mult. (22->14), div. (0->0), fcn. (4->2), ass. (0->6)
t1 = qJD(1) + qJD(2);
t6 = pkin(1) * qJD(1) * t1;
t4 = qJD(1) ^ 2;
t3 = cos(qJ(2));
t2 = sin(qJ(2));
t5 = [0, 0, 0, 0, 0, t4 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 ^ 2 / 0.2e1, t3 * t6, -t2 * t6, 0, (t2 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t4;];
T_reg = t5;
