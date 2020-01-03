% Calculate joint inertia matrix for
% S4RPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
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
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRP3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP3_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP3_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP3_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:42:32
% EndTime: 2019-12-31 16:42:32
% DurationCPUTime: 0.14s
% Computational Cost: add. (82->48), mult. (160->65), div. (0->0), fcn. (87->4), ass. (0->17)
t21 = 2 * mrSges(5,3);
t13 = sin(qJ(3));
t14 = cos(qJ(3));
t18 = t13 ^ 2 + t14 ^ 2;
t20 = 0.2e1 * t18;
t19 = m(5) * pkin(3);
t11 = sin(pkin(6));
t6 = t11 * pkin(1) + pkin(5);
t17 = qJ(4) + t6;
t16 = mrSges(5,1) + t19;
t12 = cos(pkin(6));
t7 = -t12 * pkin(1) - pkin(2);
t8 = t13 * mrSges(5,2);
t3 = -t14 * pkin(3) + t7;
t2 = t17 * t14;
t1 = t17 * t13;
t4 = [0.2e1 * t3 * t8 + Ifges(2,3) + Ifges(3,3) + m(5) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(4) * (t18 * t6 ^ 2 + t7 ^ 2) + m(3) * (t11 ^ 2 + t12 ^ 2) * pkin(1) ^ 2 + (-0.2e1 * t7 * mrSges(4,1) - 0.2e1 * t3 * mrSges(5,1) + t2 * t21 + (Ifges(5,2) + Ifges(4,2)) * t14) * t14 + 0.2e1 * (t12 * mrSges(3,1) - t11 * mrSges(3,2)) * pkin(1) + t6 * mrSges(4,3) * t20 + (0.2e1 * t7 * mrSges(4,2) + t1 * t21 + 0.2e1 * (Ifges(4,4) + Ifges(5,4)) * t14 + (Ifges(5,1) + Ifges(4,1)) * t13) * t13; m(5) * (-t1 * t14 + t2 * t13); m(3) + (m(4) / 0.2e1 + m(5) / 0.2e1) * t20; -t2 * mrSges(5,2) - t16 * t1 + (-mrSges(4,2) * t6 + Ifges(4,6) + Ifges(5,6)) * t14 + (-mrSges(4,1) * t6 - mrSges(5,3) * pkin(3) + Ifges(4,5) + Ifges(5,5)) * t13; -t13 * mrSges(4,2) - t8 + (mrSges(4,1) + t16) * t14; Ifges(4,3) + Ifges(5,3) + (0.2e1 * mrSges(5,1) + t19) * pkin(3); m(5) * t3 - t14 * mrSges(5,1) + t8; 0; 0; m(5);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t4(1), t4(2), t4(4), t4(7); t4(2), t4(3), t4(5), t4(8); t4(4), t4(5), t4(6), t4(9); t4(7), t4(8), t4(9), t4(10);];
Mq = res;
