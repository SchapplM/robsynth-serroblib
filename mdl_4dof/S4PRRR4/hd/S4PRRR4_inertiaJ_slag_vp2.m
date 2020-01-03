% Calculate joint inertia matrix for
% S4PRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRR4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR4_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR4_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR4_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR4_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR4_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:32:30
% EndTime: 2019-12-31 16:32:31
% DurationCPUTime: 0.16s
% Computational Cost: add. (119->51), mult. (260->83), div. (0->0), fcn. (210->4), ass. (0->19)
t18 = cos(qJ(3));
t26 = t18 ^ 2;
t25 = -pkin(6) - pkin(5);
t16 = sin(qJ(3));
t24 = t16 ^ 2 + t26;
t15 = sin(qJ(4));
t17 = cos(qJ(4));
t8 = -t15 * t16 + t17 * t18;
t9 = t15 * t18 + t17 * t16;
t1 = -t8 * mrSges(5,1) + t9 * mrSges(5,2);
t10 = t25 * t16;
t11 = t25 * t18;
t3 = t17 * t10 + t15 * t11;
t4 = t15 * t10 - t17 * t11;
t23 = t3 * mrSges(5,1) - t4 * mrSges(5,2) + Ifges(5,5) * t9 + Ifges(5,6) * t8;
t22 = -t18 * mrSges(4,1) + t16 * mrSges(4,2);
t21 = (t17 * mrSges(5,1) - t15 * mrSges(5,2)) * pkin(3);
t12 = -t18 * pkin(3) - pkin(2);
t2 = [m(2) + m(3) + m(4) * t24 + m(5) * (t8 ^ 2 + t9 ^ 2); m(5) * (t3 * t8 + t4 * t9); 0.2e1 * t12 * t1 - 0.2e1 * pkin(2) * t22 + Ifges(4,2) * t26 + Ifges(3,3) + (-0.2e1 * t3 * mrSges(5,3) + Ifges(5,1) * t9) * t9 + 0.2e1 * t24 * pkin(5) * mrSges(4,3) + m(5) * (t12 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(4) * (t24 * pkin(5) ^ 2 + pkin(2) ^ 2) + (0.2e1 * t4 * mrSges(5,3) + 0.2e1 * Ifges(5,4) * t9 + Ifges(5,2) * t8) * t8 + (Ifges(4,1) * t16 + 0.2e1 * Ifges(4,4) * t18) * t16; m(5) * (t15 * t9 + t17 * t8) * pkin(3) - t22 - t1; Ifges(4,5) * t16 + Ifges(4,6) * t18 + (-t16 * mrSges(4,1) - t18 * mrSges(4,2)) * pkin(5) + (m(5) * (t15 * t4 + t17 * t3) + (t15 * t8 - t17 * t9) * mrSges(5,3)) * pkin(3) + t23; Ifges(4,3) + Ifges(5,3) + m(5) * (t15 ^ 2 + t17 ^ 2) * pkin(3) ^ 2 + 0.2e1 * t21; -t1; t23; Ifges(5,3) + t21; Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t2(1), t2(2), t2(4), t2(7); t2(2), t2(3), t2(5), t2(8); t2(4), t2(5), t2(6), t2(9); t2(7), t2(8), t2(9), t2(10);];
Mq = res;
