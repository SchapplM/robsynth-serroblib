% Calculate joint inertia matrix for
% S5RPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3]';
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
% Datum: 2019-12-31 17:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPPR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_inertiaJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR3_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR3_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR3_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:43:51
% EndTime: 2019-12-31 17:43:52
% DurationCPUTime: 0.22s
% Computational Cost: add. (176->71), mult. (342->101), div. (0->0), fcn. (258->6), ass. (0->25)
t20 = sin(pkin(8));
t22 = cos(pkin(8));
t34 = t20 ^ 2 + t22 ^ 2;
t33 = 0.2e1 * t34;
t23 = cos(pkin(7));
t16 = -t23 * pkin(1) - pkin(2);
t27 = t20 * qJ(4) - t16;
t5 = -t22 * pkin(3) - t27;
t32 = -0.2e1 * t5;
t21 = sin(pkin(7));
t15 = t21 * pkin(1) + qJ(3);
t31 = -pkin(6) + t15;
t30 = t34 * t15 ^ 2;
t24 = sin(qJ(5));
t25 = cos(qJ(5));
t10 = t20 * t25 - t22 * t24;
t9 = -t20 * t24 - t22 * t25;
t3 = -t9 * mrSges(6,1) + t10 * mrSges(6,2);
t17 = t20 * mrSges(4,2);
t7 = t31 * t22;
t6 = t31 * t20;
t4 = (pkin(3) + pkin(4)) * t22 + t27;
t2 = t24 * t6 + t25 * t7;
t1 = -t24 * t7 + t25 * t6;
t8 = [0.2e1 * t16 * t17 + 0.2e1 * t4 * t3 + Ifges(2,3) + Ifges(3,3) + (0.2e1 * t2 * mrSges(6,3) + Ifges(6,2) * t9) * t9 + (-0.2e1 * t1 * mrSges(6,3) + Ifges(6,1) * t10 + 0.2e1 * Ifges(6,4) * t9) * t10 + (-0.2e1 * t16 * mrSges(4,1) + mrSges(5,1) * t32 + (Ifges(5,3) + Ifges(4,2)) * t22) * t22 + (mrSges(5,3) * t32 + (Ifges(4,1) + Ifges(5,1)) * t20 + 0.2e1 * (Ifges(4,4) - Ifges(5,5)) * t22) * t20 + m(6) * (t1 ^ 2 + t2 ^ 2 + t4 ^ 2) + m(5) * (t5 ^ 2 + t30) + m(4) * (t16 ^ 2 + t30) + (mrSges(5,2) + mrSges(4,3)) * t15 * t33 + (0.2e1 * t23 * mrSges(3,1) - 0.2e1 * t21 * mrSges(3,2) + m(3) * (t21 ^ 2 + t23 ^ 2) * pkin(1)) * pkin(1); m(6) * (t1 * t9 + t2 * t10); m(3) + m(6) * (t10 ^ 2 + t9 ^ 2) + (m(4) / 0.2e1 + m(5) / 0.2e1) * t33; -t20 * mrSges(5,3) + t17 + (-mrSges(5,1) - mrSges(4,1)) * t22 + m(4) * t16 + m(5) * t5 - m(6) * t4 - t3; 0; m(4) + m(5) + m(6); m(6) * (t25 * t1 + t24 * t2) + (m(5) * t15 + mrSges(5,2)) * t20 + (-t25 * t10 + t24 * t9) * mrSges(6,3); -m(5) * t22 + m(6) * (t24 * t10 + t25 * t9); 0; m(5) + m(6) * (t24 ^ 2 + t25 ^ 2); t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t10 + Ifges(6,6) * t9; -t3; 0; t25 * mrSges(6,1) - t24 * mrSges(6,2); Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t8(1), t8(2), t8(4), t8(7), t8(11); t8(2), t8(3), t8(5), t8(8), t8(12); t8(4), t8(5), t8(6), t8(9), t8(13); t8(7), t8(8), t8(9), t8(10), t8(14); t8(11), t8(12), t8(13), t8(14), t8(15);];
Mq = res;
