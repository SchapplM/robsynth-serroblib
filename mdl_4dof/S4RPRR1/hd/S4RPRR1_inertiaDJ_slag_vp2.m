% Calculate time derivative of joint inertia matrix for
% S4RPRR1
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
% Datum: 2019-03-08 18:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR1_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR1_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR1_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR1_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR1_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR1_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:31:39
% EndTime: 2019-03-08 18:31:39
% DurationCPUTime: 0.09s
% Computational Cost: add. (179->24), mult. (406->46), div. (0->0), fcn. (282->6), ass. (0->21)
t12 = cos(pkin(7)) * pkin(1) + pkin(2);
t15 = sin(qJ(3));
t17 = cos(qJ(3));
t24 = pkin(1) * sin(pkin(7));
t19 = t17 * t12 - t15 * t24;
t14 = sin(qJ(4));
t23 = qJD(4) * t14;
t16 = cos(qJ(4));
t22 = qJD(4) * t16;
t8 = pkin(3) + t19;
t9 = t15 * t12 + t17 * t24;
t4 = -t14 * t9 + t16 * t8;
t6 = t19 * qJD(3);
t7 = t9 * qJD(3);
t2 = qJD(4) * t4 - t14 * t7 + t16 * t6;
t5 = t14 * t8 + t16 * t9;
t3 = -qJD(4) * t5 - t14 * t6 - t16 * t7;
t20 = t3 * mrSges(5,1) - t2 * mrSges(5,2);
t18 = -t7 * mrSges(4,1) - t6 * mrSges(4,2) + t20;
t10 = (-t14 * mrSges(5,1) - t16 * mrSges(5,2)) * qJD(4) * pkin(3);
t1 = [0.2e1 * m(5) * (t5 * t2 + t4 * t3) + 0.2e1 * m(4) * (-t19 * t7 + t9 * t6) + 0.2e1 * t18; 0; 0; (m(5) * (t14 * t2 + t16 * t3 + t5 * t22 - t4 * t23) - mrSges(5,2) * t22 - mrSges(5,1) * t23) * pkin(3) + t18; 0; 0.2e1 * t10; t20; 0; t10; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1) t1(2) t1(4) t1(7); t1(2) t1(3) t1(5) t1(8); t1(4) t1(5) t1(6) t1(9); t1(7) t1(8) t1(9) t1(10);];
Mq  = res;
