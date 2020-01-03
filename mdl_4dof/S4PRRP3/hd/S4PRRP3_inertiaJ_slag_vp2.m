% Calculate joint inertia matrix for
% S4PRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
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
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:27
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRP3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP3_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP3_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP3_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:26:44
% EndTime: 2019-12-31 16:26:44
% DurationCPUTime: 0.13s
% Computational Cost: add. (58->42), mult. (123->54), div. (0->0), fcn. (63->2), ass. (0->13)
t16 = -2 * mrSges(5,3);
t8 = sin(qJ(3));
t9 = cos(qJ(3));
t13 = t8 ^ 2 + t9 ^ 2;
t15 = 0.2e1 * t13;
t14 = m(5) * pkin(3);
t12 = -qJ(4) - pkin(5);
t11 = mrSges(5,1) + t14;
t5 = t8 * mrSges(5,2);
t4 = -t9 * pkin(3) - pkin(2);
t2 = t12 * t9;
t1 = t12 * t8;
t3 = [m(2) + m(3) + (m(4) / 0.2e1 + m(5) / 0.2e1) * t15; m(5) * (t1 * t9 - t2 * t8); 0.2e1 * t4 * t5 + Ifges(3,3) + m(5) * (t1 ^ 2 + t2 ^ 2 + t4 ^ 2) + m(4) * (t13 * pkin(5) ^ 2 + pkin(2) ^ 2) + (0.2e1 * pkin(2) * mrSges(4,1) - 0.2e1 * t4 * mrSges(5,1) + t2 * t16 + (Ifges(5,2) + Ifges(4,2)) * t9) * t9 + pkin(5) * mrSges(4,3) * t15 + (-0.2e1 * pkin(2) * mrSges(4,2) + t1 * t16 + 0.2e1 * (Ifges(4,4) + Ifges(5,4)) * t9 + (Ifges(5,1) + Ifges(4,1)) * t8) * t8; -t8 * mrSges(4,2) - t5 + (mrSges(4,1) + t11) * t9; t2 * mrSges(5,2) + t11 * t1 + (-mrSges(4,2) * pkin(5) + Ifges(4,6) + Ifges(5,6)) * t9 + (-mrSges(4,1) * pkin(5) - mrSges(5,3) * pkin(3) + Ifges(4,5) + Ifges(5,5)) * t8; Ifges(4,3) + Ifges(5,3) + (0.2e1 * mrSges(5,1) + t14) * pkin(3); 0; m(5) * t4 - t9 * mrSges(5,1) + t5; 0; m(5);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t3(1), t3(2), t3(4), t3(7); t3(2), t3(3), t3(5), t3(8); t3(4), t3(5), t3(6), t3(9); t3(7), t3(8), t3(9), t3(10);];
Mq = res;
