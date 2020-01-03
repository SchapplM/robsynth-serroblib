% Calculate time derivative of joint inertia matrix for
% S4PRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-31 16:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRPR5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR5_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR5_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR5_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR5_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR5_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR5_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:22:58
% EndTime: 2019-12-31 16:22:59
% DurationCPUTime: 0.10s
% Computational Cost: add. (65->31), mult. (196->59), div. (0->0), fcn. (144->6), ass. (0->19)
t14 = cos(qJ(4));
t9 = t14 ^ 2;
t10 = sin(pkin(7));
t11 = cos(pkin(7));
t13 = sin(qJ(2));
t15 = cos(qJ(2));
t4 = t10 * t15 + t11 * t13;
t1 = t4 * qJD(2);
t3 = t10 * t13 - t11 * t15;
t19 = t3 * t1;
t12 = sin(qJ(4));
t18 = t12 ^ 2 + t9;
t17 = qJD(4) * t12;
t2 = t3 * qJD(2);
t16 = t2 * t18;
t7 = -t11 * pkin(2) - pkin(3);
t6 = t10 * pkin(2) + pkin(5);
t5 = (mrSges(5,1) * t12 + mrSges(5,2) * t14) * qJD(4);
t8 = [0.2e1 * m(4) * (-t4 * t2 + t19) + 0.2e1 * m(5) * (-t4 * t16 + t19); t3 * t5 + (-t14 * mrSges(5,1) + t12 * mrSges(5,2) - mrSges(4,1)) * t1 + (-t13 * mrSges(3,1) - t15 * mrSges(3,2)) * qJD(2) - (t18 * mrSges(5,3) - mrSges(4,2)) * t2 + m(5) * (t7 * t1 - t6 * t16) + m(4) * (-t1 * t11 - t10 * t2) * pkin(2); 0.2e1 * qJD(4) * t9 * Ifges(5,4) + 0.2e1 * t7 * t5 + 0.2e1 * (-Ifges(5,4) * t12 + (Ifges(5,1) - Ifges(5,2)) * t14) * t17; 0; 0; 0; (t14 * t2 + t4 * t17) * mrSges(5,2) + (-qJD(4) * t14 * t4 + t12 * t2) * mrSges(5,1); ((-mrSges(5,1) * t6 + Ifges(5,5)) * t14 + (mrSges(5,2) * t6 - Ifges(5,6)) * t12) * qJD(4); -t5; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t8(1), t8(2), t8(4), t8(7); t8(2), t8(3), t8(5), t8(8); t8(4), t8(5), t8(6), t8(9); t8(7), t8(8), t8(9), t8(10);];
Mq = res;
