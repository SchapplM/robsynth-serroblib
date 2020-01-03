% Calculate time derivative of joint inertia matrix for
% S5RPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR5_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR5_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR5_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR5_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR5_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR5_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:56:29
% EndTime: 2019-12-31 17:56:30
% DurationCPUTime: 0.28s
% Computational Cost: add. (292->61), mult. (545->99), div. (0->0), fcn. (321->6), ass. (0->34)
t25 = cos(qJ(4));
t22 = sin(qJ(5));
t24 = cos(qJ(5));
t36 = t22 ^ 2 + t24 ^ 2;
t48 = mrSges(6,3) * t36;
t50 = mrSges(5,2) - t48;
t51 = t50 * t25;
t49 = m(6) * pkin(7);
t12 = -t24 * mrSges(6,1) + t22 * mrSges(6,2);
t37 = mrSges(5,1) - t12;
t15 = -cos(pkin(8)) * pkin(1) - pkin(2) - pkin(3);
t16 = sin(pkin(8)) * pkin(1) + qJ(3);
t23 = sin(qJ(4));
t6 = t23 * t15 + t25 * t16;
t46 = t36 * t25;
t33 = qJD(5) * t24;
t34 = qJD(5) * t22;
t26 = -(t22 * (Ifges(6,4) * t22 + Ifges(6,2) * t24) - t24 * (Ifges(6,1) * t22 + Ifges(6,4) * t24)) * qJD(5) + t22 * (Ifges(6,1) * t33 - Ifges(6,4) * t34) + t24 * (Ifges(6,4) * t33 - Ifges(6,2) * t34);
t45 = -m(6) * pkin(4) - t37;
t44 = 0.2e1 * m(6);
t28 = mrSges(6,1) * t22 + mrSges(6,2) * t24;
t7 = t28 * qJD(5);
t43 = -0.2e1 * t7;
t5 = t25 * t15 - t23 * t16;
t1 = t25 * qJD(3) + qJD(4) * t5;
t42 = t23 * t1;
t2 = t23 * qJD(3) + qJD(4) * t6;
t41 = t25 * t2;
t40 = t25 * t7;
t35 = qJD(4) * t25;
t32 = t36 * t1;
t4 = -pkin(7) + t6;
t3 = pkin(4) - t5;
t8 = [t3 * t43 + 0.2e1 * t1 * mrSges(5,2) + (t3 * t2 + t4 * t32) * t44 + 0.2e1 * m(5) * (t6 * t1 - t5 * t2) + t26 + 0.2e1 * t37 * t2 + 0.2e1 * (m(4) * t16 + mrSges(4,3)) * qJD(3) - 0.2e1 * t48 * t1; 0; 0; t40 + m(6) * (t36 * t42 - t41) + m(5) * (-t41 + t42) + (t37 * t23 + t51 + m(6) * (t23 * t3 + t4 * t46) + m(5) * (-t23 * t5 + t25 * t6)) * qJD(4); 0; (-0.1e1 + t36) * t23 * t35 * t44; t32 * t49 + (t3 + pkin(4)) * t7 + t45 * t2 - t50 * t1 - t26; 0; -t40 + (t23 * t45 + t46 * t49 - t51) * qJD(4); pkin(4) * t43 + t26; -t28 * t1 + ((-mrSges(6,1) * t4 - Ifges(6,5)) * t24 + (mrSges(6,2) * t4 + Ifges(6,6)) * t22) * qJD(5); t7; (t23 * t34 - t24 * t35) * mrSges(6,2) + (-t22 * t35 - t23 * t33) * mrSges(6,1); (Ifges(6,5) * t24 - Ifges(6,6) * t22 + t12 * pkin(7)) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t8(1), t8(2), t8(4), t8(7), t8(11); t8(2), t8(3), t8(5), t8(8), t8(12); t8(4), t8(5), t8(6), t8(9), t8(13); t8(7), t8(8), t8(9), t8(10), t8(14); t8(11), t8(12), t8(13), t8(14), t8(15);];
Mq = res;
