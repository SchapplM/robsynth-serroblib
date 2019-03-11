% Calculate time derivative of joint inertia matrix for
% S4PRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% Datum: 2019-03-08 18:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRPR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR1_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR1_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR1_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR1_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR1_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR1_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:20:52
% EndTime: 2019-03-08 18:20:52
% DurationCPUTime: 0.06s
% Computational Cost: add. (44->15), mult. (85->28), div. (0->0), fcn. (40->2), ass. (0->10)
t5 = sin(qJ(4));
t6 = cos(qJ(4));
t13 = -t5 * mrSges(5,1) - t6 * mrSges(5,2);
t7 = -pkin(2) - pkin(3);
t3 = -t5 * qJ(3) + t6 * t7;
t1 = t6 * qJD(3) + qJD(4) * t3;
t4 = t6 * qJ(3) + t5 * t7;
t2 = -t5 * qJD(3) - qJD(4) * t4;
t12 = t2 * mrSges(5,1) - t1 * mrSges(5,2);
t8 = [0; 0; 0.2e1 * m(5) * (t4 * t1 + t3 * t2) + 0.2e1 * (m(4) * qJ(3) + mrSges(4,3)) * qJD(3) - 0.2e1 * t12; 0; m(5) * (t5 * t1 + t6 * t2) + (m(5) * (-t3 * t5 + t4 * t6) - t13) * qJD(4); 0; 0; t12; t13 * qJD(4); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t8(1) t8(2) t8(4) t8(7); t8(2) t8(3) t8(5) t8(8); t8(4) t8(5) t8(6) t8(9); t8(7) t8(8) t8(9) t8(10);];
Mq  = res;
