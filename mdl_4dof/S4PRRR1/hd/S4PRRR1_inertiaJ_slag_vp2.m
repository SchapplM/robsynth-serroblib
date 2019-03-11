% Calculate joint inertia matrix for
% S4PRRR1
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
% Datum: 2019-03-08 18:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR1_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR1_inertiaJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR1_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR1_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR1_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:25:19
% EndTime: 2019-03-08 18:25:19
% DurationCPUTime: 0.07s
% Computational Cost: add. (51->29), mult. (106->36), div. (0->0), fcn. (56->4), ass. (0->17)
t7 = sin(qJ(3));
t18 = pkin(2) * t7;
t9 = cos(qJ(3));
t4 = pkin(2) * t9 + pkin(3);
t6 = sin(qJ(4));
t8 = cos(qJ(4));
t3 = t8 * t18 + t4 * t6;
t17 = t3 * mrSges(5,2);
t16 = t6 * mrSges(5,2);
t15 = Ifges(4,3) + Ifges(5,3);
t14 = pkin(3) * t16;
t2 = -t6 * t18 + t4 * t8;
t1 = t2 * mrSges(5,1);
t13 = Ifges(5,3) + t1 - t17;
t12 = (t9 * mrSges(4,1) - t7 * mrSges(4,2)) * pkin(2);
t5 = t8 * pkin(3) * mrSges(5,1);
t10 = [m(2) + m(3) + m(4) + m(5); 0; -0.2e1 * t17 + Ifges(3,3) + 0.2e1 * t1 + 0.2e1 * t12 + m(5) * (t2 ^ 2 + t3 ^ 2) + m(4) * (t7 ^ 2 + t9 ^ 2) * pkin(2) ^ 2 + t15; 0; Ifges(4,3) + t5 + t12 + (-t16 + m(5) * (t2 * t8 + t3 * t6)) * pkin(3) + t13; -0.2e1 * t14 + 0.2e1 * t5 + m(5) * (t6 ^ 2 + t8 ^ 2) * pkin(3) ^ 2 + t15; 0; t13; Ifges(5,3) + t5 - t14; Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t10(1) t10(2) t10(4) t10(7); t10(2) t10(3) t10(5) t10(8); t10(4) t10(5) t10(6) t10(9); t10(7) t10(8) t10(9) t10(10);];
Mq  = res;
