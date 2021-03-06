% Calculate time derivative of joint inertia matrix for
% S4PPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-31 16:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PPRR4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR4_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR4_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR4_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR4_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:18:27
% EndTime: 2019-12-31 16:18:28
% DurationCPUTime: 0.09s
% Computational Cost: add. (53->24), mult. (170->48), div. (0->0), fcn. (126->6), ass. (0->17)
t12 = cos(qJ(4));
t7 = t12 ^ 2;
t11 = sin(qJ(3));
t13 = cos(qJ(3));
t8 = sin(pkin(7));
t9 = cos(pkin(7));
t4 = t11 * t9 + t13 * t8;
t2 = t4 * qJD(3);
t3 = t11 * t8 - t13 * t9;
t17 = t3 * t2;
t10 = sin(qJ(4));
t16 = qJD(4) * t10;
t1 = t3 * qJD(3);
t15 = (t10 ^ 2 + t7) * t1;
t14 = -t12 * mrSges(5,1) + t10 * mrSges(5,2);
t5 = (mrSges(5,1) * t10 + mrSges(5,2) * t12) * qJD(4);
t6 = [0.2e1 * m(4) * (-t4 * t1 + t17) + 0.2e1 * m(5) * (-t4 * t15 + t17); 0; 0; t1 * mrSges(4,2) + t3 * t5 - (m(5) * pkin(5) + mrSges(5,3)) * t15 + (-m(5) * pkin(3) - mrSges(4,1) + t14) * t2; 0; 0.2e1 * qJD(4) * t7 * Ifges(5,4) - 0.2e1 * pkin(3) * t5 + 0.2e1 * (-Ifges(5,4) * t10 + (Ifges(5,1) - Ifges(5,2)) * t12) * t16; (t1 * t12 + t4 * t16) * mrSges(5,2) + (-qJD(4) * t12 * t4 + t1 * t10) * mrSges(5,1); -t5; (Ifges(5,5) * t12 - Ifges(5,6) * t10 + t14 * pkin(5)) * qJD(4); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t6(1), t6(2), t6(4), t6(7); t6(2), t6(3), t6(5), t6(8); t6(4), t6(5), t6(6), t6(9); t6(7), t6(8), t6(9), t6(10);];
Mq = res;
