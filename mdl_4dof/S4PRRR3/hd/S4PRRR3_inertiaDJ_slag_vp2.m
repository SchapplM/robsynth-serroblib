% Calculate time derivative of joint inertia matrix for
% S4PRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR3_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR3_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR3_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR3_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:31:34
% EndTime: 2019-12-31 16:31:35
% DurationCPUTime: 0.13s
% Computational Cost: add. (63->29), mult. (198->49), div. (0->0), fcn. (93->4), ass. (0->19)
t14 = sin(qJ(4));
t16 = cos(qJ(4));
t34 = t14 ^ 2 + t16 ^ 2;
t33 = Ifges(5,1) - Ifges(5,2);
t17 = cos(qJ(3));
t32 = t34 * t17;
t15 = sin(qJ(3));
t24 = t16 * mrSges(5,1);
t5 = t14 * mrSges(5,2) - t24;
t31 = (t34 * mrSges(5,3) - mrSges(4,2)) * t17 + (t5 - mrSges(4,1)) * t15;
t22 = pkin(2) * qJD(3);
t21 = qJD(4) * t14;
t20 = qJD(4) * t16;
t18 = (-0.2e1 * Ifges(5,4) * t14 + t33 * t16) * t21 + (0.2e1 * Ifges(5,4) * t16 + t33 * t14) * t20;
t11 = Ifges(5,5) * t20;
t8 = -t17 * pkin(2) - pkin(3);
t7 = t15 * pkin(2) + pkin(6);
t4 = -mrSges(5,1) * t21 - mrSges(5,2) * t20;
t1 = [0; 0; -0.2e1 * t8 * t4 + 0.2e1 * (m(5) * (t15 * t8 + t32 * t7) + t31) * t22 + t18; 0; (pkin(3) - t8) * t4 + (m(5) * (-pkin(3) * t15 + t32 * pkin(6)) + t31) * t22 + t18; 0.2e1 * pkin(3) * t4 + t18; t4; t11 + (-mrSges(5,1) * t14 - mrSges(5,2) * t16) * t17 * t22 + (-t7 * t24 + (mrSges(5,2) * t7 - Ifges(5,6)) * t14) * qJD(4); t11 + (-Ifges(5,6) * t14 + t5 * pkin(6)) * qJD(4); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
