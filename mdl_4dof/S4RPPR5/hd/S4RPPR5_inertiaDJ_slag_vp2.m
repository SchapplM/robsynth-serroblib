% Calculate time derivative of joint inertia matrix for
% S4RPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPPR5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR5_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR5_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR5_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR5_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:39:43
% EndTime: 2019-12-31 16:39:43
% DurationCPUTime: 0.14s
% Computational Cost: add. (55->29), mult. (135->48), div. (0->0), fcn. (71->4), ass. (0->18)
t10 = sin(qJ(4));
t29 = t10 ^ 2;
t11 = cos(qJ(4));
t28 = t11 ^ 2;
t24 = t28 + t29;
t12 = -pkin(1) - pkin(2);
t8 = sin(pkin(6));
t9 = cos(pkin(6));
t25 = t9 * qJ(2) + t8 * t12;
t21 = qJD(2) * t8;
t20 = qJD(2) * t9;
t15 = -t8 * qJ(2) + t9 * t12;
t14 = t11 * mrSges(5,1) - t10 * mrSges(5,2);
t13 = -mrSges(5,1) * t10 - mrSges(5,2) * t11;
t3 = t13 * qJD(4);
t2 = -pkin(5) + t25;
t1 = pkin(3) - t15;
t4 = [0.2e1 * t1 * t3 + 0.2e1 * (t14 + mrSges(4,1)) * t21 + 0.2e1 * (m(5) * (t24 * t9 * t2 + t1 * t8) + m(4) * (-t15 * t8 + t25 * t9) + m(3) * qJ(2) + mrSges(3,3)) * qJD(2) + 0.2e1 * (Ifges(5,1) - Ifges(5,2)) * t11 * t10 * qJD(4) + 0.2e1 * (-t24 * mrSges(5,3) + mrSges(4,2)) * t20 + 0.2e1 * (t28 - t29) * Ifges(5,4) * qJD(4); (-t3 + m(5) * (-0.1e1 + t24) * t21) * t9; 0; 0; 0; 0; t13 * t20 + ((-mrSges(5,1) * t2 - Ifges(5,5)) * t11 + (mrSges(5,2) * t2 + Ifges(5,6)) * t10) * qJD(4); -t14 * t8 * qJD(4); t3; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t4(1), t4(2), t4(4), t4(7); t4(2), t4(3), t4(5), t4(8); t4(4), t4(5), t4(6), t4(9); t4(7), t4(8), t4(9), t4(10);];
Mq = res;
