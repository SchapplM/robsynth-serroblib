% Calculate time derivative of joint inertia matrix for
% S4RPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR2_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR2_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR2_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR2_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR2_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR2_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:48:07
% EndTime: 2019-12-31 16:48:08
% DurationCPUTime: 0.13s
% Computational Cost: add. (123->36), mult. (317->61), div. (0->0), fcn. (184->6), ass. (0->25)
t37 = Ifges(5,1) - Ifges(5,2);
t20 = sin(qJ(4));
t22 = cos(qJ(4));
t32 = t22 * mrSges(5,1);
t9 = t20 * mrSges(5,2) - t32;
t36 = -mrSges(4,1) + t9;
t13 = cos(pkin(7)) * pkin(1) + pkin(2);
t21 = sin(qJ(3));
t23 = cos(qJ(3));
t35 = pkin(1) * sin(pkin(7));
t31 = t21 * t13 + t23 * t35;
t30 = t20 ^ 2 + t22 ^ 2;
t29 = qJD(4) * t20;
t28 = qJD(4) * t22;
t25 = t23 * t13 - t21 * t35;
t1 = t25 * qJD(3);
t27 = t30 * t1;
t26 = mrSges(5,3) * t30;
t24 = (-0.2e1 * Ifges(5,4) * t20 + t37 * t22) * t29 + (0.2e1 * Ifges(5,4) * t22 + t37 * t20) * t28;
t16 = Ifges(5,5) * t28;
t8 = -mrSges(5,1) * t29 - mrSges(5,2) * t28;
t4 = pkin(6) + t31;
t3 = -pkin(3) - t25;
t2 = t31 * qJD(3);
t5 = [-0.2e1 * t3 * t8 - 0.2e1 * t1 * mrSges(4,2) + 0.2e1 * m(5) * (t3 * t2 + t4 * t27) + 0.2e1 * m(4) * (t31 * t1 - t25 * t2) + t24 + 0.2e1 * t36 * t2 + 0.2e1 * t1 * t26; 0; 0; m(5) * pkin(6) * t27 + (pkin(3) - t3) * t8 + t24 + (-m(5) * pkin(3) + t36) * t2 + (-mrSges(4,2) + t26) * t1; 0; 0.2e1 * pkin(3) * t8 + t24; t16 + (-mrSges(5,1) * t20 - mrSges(5,2) * t22) * t1 + (-t4 * t32 + (mrSges(5,2) * t4 - Ifges(5,6)) * t20) * qJD(4); t8; t16 + (-Ifges(5,6) * t20 + t9 * pkin(6)) * qJD(4); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t5(1), t5(2), t5(4), t5(7); t5(2), t5(3), t5(5), t5(8); t5(4), t5(5), t5(6), t5(9); t5(7), t5(8), t5(9), t5(10);];
Mq = res;
