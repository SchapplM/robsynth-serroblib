% Calculate joint inertia matrix for
% S4RPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
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
% Datum: 2019-03-08 18:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPPR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR1_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR1_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR1_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR1_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR1_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:27:25
% EndTime: 2019-03-08 18:27:26
% DurationCPUTime: 0.08s
% Computational Cost: add. (64->31), mult. (95->38), div. (0->0), fcn. (49->4), ass. (0->13)
t7 = cos(pkin(6));
t5 = -pkin(1) * t7 - pkin(2);
t3 = -pkin(3) + t5;
t6 = sin(pkin(6));
t4 = pkin(1) * t6 + qJ(3);
t8 = sin(qJ(4));
t9 = cos(qJ(4));
t1 = t3 * t9 - t4 * t8;
t13 = t1 * mrSges(5,1);
t2 = t3 * t8 + t4 * t9;
t12 = t2 * mrSges(5,2);
t11 = t9 * mrSges(5,1) - t8 * mrSges(5,2);
t10 = [-0.2e1 * t5 * mrSges(4,1) - 0.2e1 * t13 + 0.2e1 * t12 + 0.2e1 * t4 * mrSges(4,3) + Ifges(4,2) + Ifges(2,3) + Ifges(3,3) + Ifges(5,3) + m(4) * (t4 ^ 2 + t5 ^ 2) + m(5) * (t1 ^ 2 + t2 ^ 2) + (0.2e1 * t7 * mrSges(3,1) - 0.2e1 * t6 * mrSges(3,2) + m(3) * (t6 ^ 2 + t7 ^ 2) * pkin(1)) * pkin(1); 0; m(3) + m(4) + m(5); -mrSges(4,1) + m(4) * t5 + m(5) * (t1 * t9 + t2 * t8) - t11; 0; m(4) + m(5) * (t8 ^ 2 + t9 ^ 2); -Ifges(5,3) - t12 + t13; 0; t11; Ifges(5,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t10(1) t10(2) t10(4) t10(7); t10(2) t10(3) t10(5) t10(8); t10(4) t10(5) t10(6) t10(9); t10(7) t10(8) t10(9) t10(10);];
Mq  = res;
