% Calculate minimal parameter regressor of coriolis matrix for
% S5PPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x14]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PPRRP4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:34:43
% EndTime: 2019-12-31 17:34:43
% DurationCPUTime: 0.20s
% Computational Cost: add. (140->36), mult. (377->66), div. (0->0), fcn. (280->4), ass. (0->36)
t27 = sin(qJ(4));
t25 = t27 ^ 2;
t29 = cos(qJ(4));
t26 = t29 ^ 2;
t17 = t25 + t26;
t47 = pkin(4) * t27;
t46 = -qJ(5) - pkin(6);
t21 = -t29 * pkin(4) - pkin(3);
t1 = t21 * t47;
t45 = t1 * qJD(3);
t28 = sin(qJ(3));
t30 = cos(qJ(3));
t8 = (-0.1e1 + t17) * t30 * t28;
t44 = t8 * qJD(2);
t43 = qJD(3) * t27;
t42 = qJD(3) * t29;
t41 = t17 * qJD(3);
t18 = t26 - t25;
t40 = t18 * qJD(3);
t39 = t27 * qJD(4);
t38 = t28 * qJD(3);
t23 = t29 * qJD(4);
t37 = t30 * qJD(3);
t36 = pkin(3) * t43;
t35 = pkin(3) * t42;
t34 = pkin(4) * t43;
t33 = t27 * t42;
t32 = t28 * t23;
t16 = t46 * t29;
t7 = -t16 * t29 - t25 * t46;
t9 = (0.1e1 / 0.2e1 - t26 / 0.2e1 - t25 / 0.2e1) * t28;
t31 = -t9 * qJD(2) + t7 * qJD(3);
t22 = pkin(4) * t39;
t10 = (0.1e1 + t17) * t28 / 0.2e1;
t2 = t30 * t47;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t23, 0, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8 * qJD(3); 0, 0, 0, -t38, -t37, 0, 0, 0, 0, 0, -t29 * t38 - t30 * t39, -t30 * t23 + t27 * t38, t17 * t37, t44 + (t28 * t21 + t7 * t30) * qJD(3) - t2 * qJD(4) + t10 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27 * t37 - t32, t28 * t39 - t29 * t37, 0, -pkin(4) * t32 - t2 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9 * qJD(5) - t44; 0, 0, 0, 0, 0, t27 * t23, t18 * qJD(4), 0, 0, 0, -pkin(3) * t39, -pkin(3) * t23, t17 * qJD(5), t1 * qJD(4) + t7 * qJD(5); 0, 0, 0, 0, 0, t33, t40, t23, -t39, 0, -pkin(6) * t23 - t36, pkin(6) * t39 - t35, -pkin(4) * t23, t16 * pkin(4) * qJD(4) + t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t33, -t40, 0, 0, 0, t36, t35, 0, -qJD(5) * t47 - t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, t22 - t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t3;
