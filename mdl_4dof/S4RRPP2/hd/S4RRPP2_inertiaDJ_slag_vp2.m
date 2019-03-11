% Calculate time derivative of joint inertia matrix for
% S4RRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2019-03-08 18:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPP2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP2_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP2_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP2_inertiaDJ_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP2_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP2_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP2_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:34:04
% EndTime: 2019-03-08 18:34:04
% DurationCPUTime: 0.07s
% Computational Cost: add. (42->22), mult. (110->29), div. (0->0), fcn. (28->2), ass. (0->14)
t20 = mrSges(4,3) + mrSges(5,2);
t14 = pkin(1) * qJD(2);
t9 = cos(qJ(2));
t4 = t9 * t14 + qJD(3);
t8 = sin(qJ(2));
t19 = t20 * t4 + (-mrSges(3,2) * t9 + (-mrSges(3,1) - mrSges(4,1) - mrSges(5,1)) * t8) * t14;
t18 = m(4) + m(5);
t5 = t8 * pkin(1) + qJ(3);
t17 = qJ(3) * t4 + qJD(3) * t5;
t16 = t20 * qJD(3);
t13 = t8 * t14;
t12 = -t9 * pkin(1) - pkin(2);
t1 = t5 * t4;
t2 = [0.2e1 * m(4) * (t12 * t13 + t1) + 0.2e1 * m(5) * (t1 + (-pkin(3) + t12) * t13) + 0.2e1 * t19; m(4) * (-pkin(2) * t13 + t17) + m(5) * ((-pkin(2) - pkin(3)) * t13 + t17) + t16 + t19; 0.2e1 * t18 * qJD(3) * qJ(3) + 0.2e1 * t16; t18 * t13; 0; 0; 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t2(1) t2(2) t2(4) t2(7); t2(2) t2(3) t2(5) t2(8); t2(4) t2(5) t2(6) t2(9); t2(7) t2(8) t2(9) t2(10);];
Mq  = res;
