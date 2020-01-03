% Calculate time derivative of joint inertia matrix for
% S4PPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PPRR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR3_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR3_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR3_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR3_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:17:22
% EndTime: 2019-12-31 16:17:22
% DurationCPUTime: 0.09s
% Computational Cost: add. (28->19), mult. (101->40), div. (0->0), fcn. (53->4), ass. (0->11)
t6 = cos(qJ(4));
t3 = t6 ^ 2;
t4 = sin(qJ(4));
t11 = t4 ^ 2 + t3;
t7 = cos(qJ(3));
t10 = qJD(3) * t7;
t5 = sin(qJ(3));
t9 = qJD(4) * t5;
t8 = -t6 * mrSges(5,1) + t4 * mrSges(5,2);
t1 = (mrSges(5,1) * t4 + mrSges(5,2) * t6) * qJD(4);
t2 = [0; 0; 0.2e1 * m(5) * (-0.1e1 + t11) * t5 * t10; 0; (m(5) * t11 * pkin(5) * t7 + (-m(5) * pkin(3) - mrSges(4,1) + t8) * t5) * qJD(3) + (-t1 + (t11 * mrSges(5,3) - mrSges(4,2)) * qJD(3)) * t7; -0.2e1 * pkin(3) * t1 + 0.2e1 * (t3 * Ifges(5,4) + (-Ifges(5,4) * t4 + (Ifges(5,1) - Ifges(5,2)) * t6) * t4) * qJD(4); t1; (-t6 * t10 + t4 * t9) * mrSges(5,2) + (-t4 * t10 - t6 * t9) * mrSges(5,1); (Ifges(5,5) * t6 - Ifges(5,6) * t4 + t8 * pkin(5)) * qJD(4); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t2(1), t2(2), t2(4), t2(7); t2(2), t2(3), t2(5), t2(8); t2(4), t2(5), t2(6), t2(9); t2(7), t2(8), t2(9), t2(10);];
Mq = res;
