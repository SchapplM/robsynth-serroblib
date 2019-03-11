% Calculate time derivative of joint inertia matrix for
% S4PRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-03-08 18:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRP1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP1_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP1_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP1_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP1_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP1_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP1_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:22:49
% EndTime: 2019-03-08 18:22:49
% DurationCPUTime: 0.05s
% Computational Cost: add. (20->13), mult. (62->23), div. (0->0), fcn. (16->2), ass. (0->10)
t9 = pkin(2) * qJD(3);
t4 = sin(qJ(3));
t8 = t4 * t9;
t5 = cos(qJ(3));
t7 = t5 * t9;
t1 = qJD(4) + t7;
t6 = -mrSges(4,2) * t7 + t1 * mrSges(5,3) + (-mrSges(4,1) - mrSges(5,1)) * t8;
t3 = qJD(4) * mrSges(5,3);
t2 = t4 * pkin(2) + qJ(4);
t10 = [0; 0; 0.2e1 * m(5) * (t2 * t1 + (-t5 * pkin(2) - pkin(3)) * t8) + 0.2e1 * t6; 0; t3 + m(5) * (-pkin(3) * t8 + qJ(4) * t1 + qJD(4) * t2) + t6; 0.2e1 * m(5) * qJ(4) * qJD(4) + 0.2e1 * t3; 0; m(5) * t8; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t10(1) t10(2) t10(4) t10(7); t10(2) t10(3) t10(5) t10(8); t10(4) t10(5) t10(6) t10(9); t10(7) t10(8) t10(9) t10(10);];
Mq  = res;
