% Calculate time derivative of joint inertia matrix for
% S4RRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
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
% Datum: 2019-03-08 18:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPP1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP1_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP1_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP1_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP1_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP1_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP1_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:32:59
% EndTime: 2019-03-08 18:32:59
% DurationCPUTime: 0.09s
% Computational Cost: add. (58->26), mult. (189->43), div. (0->0), fcn. (97->4), ass. (0->18)
t12 = cos(pkin(6));
t14 = cos(qJ(2));
t11 = sin(pkin(6));
t13 = sin(qJ(2));
t18 = t11 * t13 * pkin(1);
t19 = pkin(1) * qJD(2);
t4 = t12 * t14 * t19 - qJD(2) * t18;
t1 = qJD(4) + t4;
t20 = t12 * t13;
t3 = (t11 * t14 + t20) * t19;
t24 = (-mrSges(3,1) * t13 - mrSges(3,2) * t14) * t19 + (-mrSges(4,1) - mrSges(5,1)) * t3 + t1 * mrSges(5,3) - t4 * mrSges(4,2);
t9 = t14 * pkin(1) + pkin(2);
t23 = pkin(1) * t20 + t11 * t9;
t16 = t12 * t9 - t18;
t10 = qJD(4) * mrSges(5,3);
t8 = t11 * pkin(2) + qJ(4);
t2 = qJ(4) + t23;
t5 = [0.2e1 * m(4) * (-t16 * t3 + t23 * t4) + 0.2e1 * m(5) * (t2 * t1 + (-pkin(3) - t16) * t3) + 0.2e1 * t24; t10 + m(4) * (t11 * t4 - t12 * t3) * pkin(2) + m(5) * (qJD(4) * t2 + t8 * t1 + (-t12 * pkin(2) - pkin(3)) * t3) + t24; 0.2e1 * m(5) * t8 * qJD(4) + 0.2e1 * t10; 0; 0; 0; m(5) * t3; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t5(1) t5(2) t5(4) t5(7); t5(2) t5(3) t5(5) t5(8); t5(4) t5(5) t5(6) t5(9); t5(7) t5(8) t5(9) t5(10);];
Mq  = res;
