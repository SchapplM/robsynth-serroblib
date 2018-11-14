% Calculate inertial parameters regressor of fixed base kinetic energy for
% S4PPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d4,theta3]';
% 
% Output:
% T_reg [1x(4*10)]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:57
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_reg = S4PPPR3_energykin_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR3_energykin_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR3_energykin_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR3_energykin_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:56:32
% EndTime: 2018-11-14 13:56:33
% DurationCPUTime: 0.06s
% Computational Cost: add. (20->10), mult. (50->24), div. (0->0), fcn. (28->4), ass. (0->11)
t10 = cos(qJ(4));
t9 = sin(qJ(4));
t8 = cos(pkin(5));
t7 = sin(pkin(5));
t6 = qJD(1) ^ 2 / 0.2e1;
t5 = qJD(3) ^ 2 / 0.2e1;
t4 = qJD(1) * t8 + qJD(2) * t7;
t3 = -qJD(1) * t7 + qJD(2) * t8;
t2 = t10 * t4 + t3 * t9;
t1 = t10 * t3 - t4 * t9;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6 + qJD(2) ^ 2 / 0.2e1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4 ^ 2 / 0.2e1 + t3 ^ 2 / 0.2e1 + t5, 0, 0, 0, 0, 0, qJD(4) ^ 2 / 0.2e1, t1 * qJD(4), -t2 * qJD(4), 0, t2 ^ 2 / 0.2e1 + t1 ^ 2 / 0.2e1 + t5;];
T_reg  = t11;
