% Calculate inertial parameters regressor of fixed base kinetic energy for
% S1P1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% 
% Output:
% T_reg [1x(1*10)]
%   inertial parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-14 12:22
% Revision: 96facaeb42edba38506bd76ea342a8981e82f256 (2020-11-21)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S1P1_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'S1P1_energykin_fixb_reg2_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'S1P1_energykin_fixb_reg2_slag_vp: qJD has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S1P1_energykin_fixb_reg2_slag_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-14 12:21:54
% EndTime: 2021-01-14 12:21:54
% DurationCPUTime: 0.10s
% Computational Cost: add. (0->0), mult. (2->2), div. (0->0), fcn. (0->0), ass. (0->1)
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(1) ^ 2 / 0.2e1;];
T_reg = t1;
