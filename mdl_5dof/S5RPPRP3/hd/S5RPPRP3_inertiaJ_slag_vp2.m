% Calculate joint inertia matrix for
% S5RPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRP3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP3_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP3_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP3_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:50:44
% EndTime: 2019-12-31 17:50:45
% DurationCPUTime: 0.18s
% Computational Cost: add. (135->63), mult. (227->74), div. (0->0), fcn. (117->4), ass. (0->22)
t16 = sin(pkin(7));
t9 = t16 * pkin(1) + qJ(3);
t29 = t9 ^ 2;
t28 = 0.2e1 * t9;
t27 = -2 * mrSges(6,3);
t18 = sin(qJ(4));
t19 = cos(qJ(4));
t7 = t18 ^ 2 + t19 ^ 2;
t4 = m(5) * t7;
t26 = m(6) * pkin(4);
t17 = cos(pkin(7));
t11 = -t17 * pkin(1) - pkin(2);
t8 = -pkin(6) + t11;
t25 = -qJ(5) + t8;
t24 = t18 * mrSges(6,1) + t19 * mrSges(6,2);
t23 = m(6) * t7 + m(4) + t4;
t22 = mrSges(6,1) + t26;
t21 = t7 * mrSges(5,3);
t5 = t18 * pkin(4) + t9;
t2 = t25 * t19;
t1 = t25 * t18;
t3 = [Ifges(4,1) + Ifges(2,3) + Ifges(3,3) + 0.2e1 * t5 * t24 + 0.2e1 * t11 * mrSges(4,2) + mrSges(4,3) * t28 + (mrSges(5,2) * t28 + t2 * t27 + (Ifges(6,1) + Ifges(5,1)) * t19) * t19 + m(6) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + m(5) * (t7 * t8 ^ 2 + t29) + m(4) * (t11 ^ 2 + t29) + m(3) * (t16 ^ 2 + t17 ^ 2) * pkin(1) ^ 2 + 0.2e1 * (t17 * mrSges(3,1) - t16 * mrSges(3,2)) * pkin(1) - 0.2e1 * t8 * t21 + (mrSges(5,1) * t28 + t1 * t27 + 0.2e1 * (-Ifges(5,4) - Ifges(6,4)) * t19 + (Ifges(6,2) + Ifges(5,2)) * t18) * t18; m(6) * (t1 * t19 - t2 * t18); m(3) + t23; mrSges(4,2) - t7 * mrSges(6,3) - t21 + m(6) * (t18 * t1 + t19 * t2) + m(4) * t11 + t8 * t4; 0; t23; -t1 * mrSges(6,2) + t22 * t2 + (-mrSges(5,2) * t8 - Ifges(5,6) - Ifges(6,6)) * t18 + (mrSges(5,1) * t8 - mrSges(6,3) * pkin(4) + Ifges(5,5) + Ifges(6,5)) * t19; -t19 * mrSges(5,2) + (-mrSges(5,1) - t26) * t18 - t24; (-mrSges(5,2) - mrSges(6,2)) * t18 + (mrSges(5,1) + t22) * t19; Ifges(5,3) + Ifges(6,3) + (0.2e1 * mrSges(6,1) + t26) * pkin(4); m(6) * t5 + t24; 0; 0; 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t3(1), t3(2), t3(4), t3(7), t3(11); t3(2), t3(3), t3(5), t3(8), t3(12); t3(4), t3(5), t3(6), t3(9), t3(13); t3(7), t3(8), t3(9), t3(10), t3(14); t3(11), t3(12), t3(13), t3(14), t3(15);];
Mq = res;
