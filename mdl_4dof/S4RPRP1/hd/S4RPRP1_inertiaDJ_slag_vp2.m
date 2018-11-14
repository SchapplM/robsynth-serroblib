% Calculate time derivative of joint inertia matrix for
% S4RPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:49
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S4RPRP1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP1_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP1_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP1_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP1_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP1_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP1_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:48:28
% EndTime: 2018-11-14 13:48:28
% DurationCPUTime: 0.08s
% Computational Cost: add. (62->19), mult. (152->32), div. (0->0), fcn. (80->4), ass. (0->14)
t11 = sin(qJ(3));
t12 = cos(qJ(3));
t18 = sin(pkin(6)) * pkin(1);
t8 = cos(pkin(6)) * pkin(1) + pkin(2);
t13 = -t11 * t18 + t12 * t8;
t21 = -mrSges(4,1) - mrSges(5,1);
t2 = t13 * qJD(3);
t1 = qJD(4) + t2;
t20 = -t2 * mrSges(4,2) + t1 * mrSges(5,3);
t19 = t11 * t8 + t12 * t18;
t9 = qJD(4) * mrSges(5,3);
t4 = qJ(4) + t19;
t3 = t19 * qJD(3);
t5 = [0.2e1 * t21 * t3 + 0.2e1 * m(4) * (-t13 * t3 + t19 * t2) + 0.2e1 * m(5) * (t4 * t1 + (-pkin(3) - t13) * t3) + 0.2e1 * t20; 0; 0; t9 + m(5) * (qJ(4) * t1 + qJD(4) * t4) + (-m(5) * pkin(3) + t21) * t3 + t20; 0; 0.2e1 * m(5) * qJ(4) * qJD(4) + 0.2e1 * t9; m(5) * t3; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t5(1) t5(2) t5(4) t5(7); t5(2) t5(3) t5(5) t5(8); t5(4) t5(5) t5(6) t5(9); t5(7) t5(8) t5(9) t5(10);];
Mq  = res;
