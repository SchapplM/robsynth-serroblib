% Calculate joint inertia matrix for
% S5RPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
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
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRP4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP4_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP4_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP4_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:52:04
% EndTime: 2019-12-31 17:52:04
% DurationCPUTime: 0.21s
% Computational Cost: add. (192->77), mult. (305->93), div. (0->0), fcn. (183->4), ass. (0->28)
t36 = 2 * mrSges(6,3);
t35 = m(5) + m(6);
t34 = m(6) * pkin(4);
t25 = cos(qJ(4));
t33 = t25 * pkin(4);
t32 = -mrSges(5,2) - mrSges(6,2);
t22 = sin(pkin(7));
t23 = cos(pkin(7));
t26 = -pkin(1) - pkin(2);
t9 = t23 * qJ(2) + t22 * t26;
t7 = -pkin(6) + t9;
t31 = qJ(5) - t7;
t24 = sin(qJ(4));
t30 = t24 ^ 2 + t25 ^ 2;
t28 = mrSges(6,1) + t34;
t27 = t30 * mrSges(5,3);
t16 = t25 * mrSges(6,1);
t10 = -t24 * mrSges(6,2) + t16;
t8 = -t22 * qJ(2) + t23 * t26;
t6 = pkin(3) - t8;
t19 = t23 ^ 2;
t18 = t22 ^ 2;
t17 = t25 * mrSges(5,1);
t11 = -t24 * mrSges(5,2) + t17;
t4 = t6 + t33;
t2 = t31 * t25;
t1 = t31 * t24;
t3 = [(2 * pkin(1) * mrSges(3,1)) - 0.2e1 * t8 * mrSges(4,1) + 0.2e1 * t9 * mrSges(4,2) + 0.2e1 * qJ(2) * mrSges(3,3) + 0.2e1 * t4 * t10 + 0.2e1 * t6 * t11 + Ifges(3,2) + Ifges(2,3) + Ifges(4,3) + (t2 * t36 + (Ifges(6,2) + Ifges(5,2)) * t25) * t25 + m(6) * (t1 ^ 2 + t2 ^ 2 + t4 ^ 2) + m(5) * (t30 * t7 ^ 2 + t6 ^ 2) + m(4) * (t8 ^ 2 + t9 ^ 2) + m(3) * ((pkin(1) ^ 2) + qJ(2) ^ 2) - 0.2e1 * t7 * t27 + (t1 * t36 + 0.2e1 * (Ifges(5,4) + Ifges(6,4)) * t25 + (Ifges(6,1) + Ifges(5,1)) * t24) * t24; -m(3) * pkin(1) - mrSges(3,1) + (-mrSges(4,1) - t10 - t11) * t23 + (-t30 * mrSges(6,3) + mrSges(4,2) - t27) * t22 + m(6) * (-t23 * t4 + (-t1 * t24 - t2 * t25) * t22) + m(5) * (t30 * t7 * t22 - t23 * t6) + m(4) * (t22 * t9 + t23 * t8); m(3) + m(4) * (t18 + t19) + (t30 * t18 + t19) * t35; m(6) * (t1 * t25 - t2 * t24); 0; t30 * t35 + m(4); t2 * mrSges(6,2) + t28 * t1 + (-mrSges(5,2) * t7 - Ifges(5,6) - Ifges(6,6)) * t25 + (-mrSges(5,1) * t7 + mrSges(6,3) * pkin(4) - Ifges(5,5) - Ifges(6,5)) * t24; (t32 * t25 + (-mrSges(5,1) - t28) * t24) * t22; m(6) * t33 + t32 * t24 + t16 + t17; Ifges(5,3) + Ifges(6,3) + (0.2e1 * mrSges(6,1) + t34) * pkin(4); m(6) * t4 + t10; -m(6) * t23; 0; 0; m(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t3(1), t3(2), t3(4), t3(7), t3(11); t3(2), t3(3), t3(5), t3(8), t3(12); t3(4), t3(5), t3(6), t3(9), t3(13); t3(7), t3(8), t3(9), t3(10), t3(14); t3(11), t3(12), t3(13), t3(14), t3(15);];
Mq = res;
