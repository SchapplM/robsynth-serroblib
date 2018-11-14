% Calculate joint inertia matrix for
% S4PPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,theta2]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:08
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S4PPRP5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP5_inertiaJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP5_inertiaJ_slag_vp2: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP5_inertiaJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRP5_inertiaJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRP5_inertiaJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:07:13
% EndTime: 2018-11-14 14:07:13
% DurationCPUTime: 0.07s
% Computational Cost: add. (30->18), mult. (65->22), div. (0->0), fcn. (46->4), ass. (0->9)
t12 = m(4) + m(5);
t11 = -m(5) * pkin(3) - mrSges(5,1);
t10 = cos(qJ(3));
t8 = sin(qJ(3));
t7 = cos(pkin(5));
t6 = sin(pkin(5));
t4 = t10 * t6 + t8 * t7;
t2 = -t10 * t7 + t8 * t6;
t1 = [m(2) + m(3) * (t6 ^ 2 + t7 ^ 2) + t12 * (t2 ^ 2 + t4 ^ 2); 0; m(3) + t12; (m(5) * qJ(4) - mrSges(4,2) + mrSges(5,3)) * t4 + (-mrSges(4,1) + t11) * t2; 0; Ifges(5,2) + Ifges(4,3) + (2 * pkin(3) * mrSges(5,1)) + 0.2e1 * qJ(4) * mrSges(5,3) + m(5) * ((pkin(3) ^ 2) + qJ(4) ^ 2); m(5) * t2; 0; t11; m(5);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1) t1(2) t1(4) t1(7); t1(2) t1(3) t1(5) t1(8); t1(4) t1(5) t1(6) t1(9); t1(7) t1(8) t1(9) t1(10);];
Mq  = res;
