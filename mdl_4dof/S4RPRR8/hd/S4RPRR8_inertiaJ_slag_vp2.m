% Calculate joint inertia matrix for
% S4RPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR8_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR8_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR8_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR8_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:55:05
% EndTime: 2019-12-31 16:55:06
% DurationCPUTime: 0.18s
% Computational Cost: add. (166->61), mult. (293->89), div. (0->0), fcn. (226->4), ass. (0->25)
t22 = cos(qJ(3));
t37 = t22 ^ 2;
t23 = -pkin(1) - pkin(5);
t35 = -pkin(6) + t23;
t20 = sin(qJ(3));
t34 = t20 ^ 2 + t37;
t19 = sin(qJ(4));
t21 = cos(qJ(4));
t10 = -t19 * t20 + t21 * t22;
t8 = -t19 * t22 - t21 * t20;
t33 = t10 ^ 2 + t8 ^ 2;
t32 = t10 * mrSges(5,1) + t8 * mrSges(5,2);
t31 = m(4) * t34;
t30 = t34 * mrSges(4,3);
t11 = t35 * t20;
t12 = t35 * t22;
t2 = -t19 * t11 + t21 * t12;
t3 = t21 * t11 + t19 * t12;
t29 = t10 * t2 - t8 * t3;
t28 = t2 * mrSges(5,1) - t3 * mrSges(5,2) + Ifges(5,5) * t10 + Ifges(5,6) * t8;
t27 = t10 * t21 - t19 * t8;
t26 = (t21 * mrSges(5,1) - t19 * mrSges(5,2)) * pkin(3);
t24 = qJ(2) ^ 2;
t13 = t20 * pkin(3) + qJ(2);
t1 = [Ifges(3,1) + Ifges(2,3) + 0.2e1 * t13 * (-t8 * mrSges(5,1) + t10 * mrSges(5,2)) + t10 * (Ifges(5,1) * t10 + Ifges(5,4) * t8) + t8 * (Ifges(5,4) * t10 + Ifges(5,2) * t8) + Ifges(4,1) * t37 - (2 * pkin(1) * mrSges(3,2)) + m(5) * (t13 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(4) * (t34 * t23 ^ 2 + t24) + m(3) * ((pkin(1) ^ 2) + t24) + (-0.2e1 * Ifges(4,4) * t22 + Ifges(4,2) * t20) * t20 + 0.2e1 * (t20 * mrSges(4,1) + t22 * mrSges(4,2) + mrSges(3,3)) * qJ(2) - 0.2e1 * t29 * mrSges(5,3) - 0.2e1 * t23 * t30; -m(3) * pkin(1) + m(5) * t29 - t33 * mrSges(5,3) + t23 * t31 + mrSges(3,2) - t30; m(5) * t33 + m(3) + t31; (t23 * mrSges(4,1) + Ifges(4,5)) * t22 + (-t23 * mrSges(4,2) - Ifges(4,6)) * t20 + (m(5) * (t19 * t3 + t2 * t21) - t27 * mrSges(5,3)) * pkin(3) + t28; m(5) * t27 * pkin(3) + t22 * mrSges(4,1) - t20 * mrSges(4,2) + t32; Ifges(4,3) + Ifges(5,3) + m(5) * (t19 ^ 2 + t21 ^ 2) * pkin(3) ^ 2 + 0.2e1 * t26; t28; t32; Ifges(5,3) + t26; Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
