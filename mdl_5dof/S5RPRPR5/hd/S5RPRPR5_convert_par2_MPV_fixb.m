% Return the minimum parameter vector for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% m [6x1]
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
% MPV [24x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:56
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RPRPR5_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_convert_par2_MPV_fixb: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR5_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR5_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t18 = (pkin(7) * m(6));
t15 = pkin(4) ^ 2 * m(6);
t2 = -Ifges(5,1) + Ifges(5,2) + t15;
t6 = cos(pkin(9));
t3 = t6 ^ 2;
t17 = t2 * t3;
t4 = sin(pkin(9));
t16 = t4 * t6;
t14 = Ifges(5,4) * t16;
t13 = m(4) * pkin(6) + mrSges(4,3);
t12 = -mrSges(6,3) - t18;
t11 = -(2 * mrSges(4,3) * pkin(6)) - Ifges(3,1) - Ifges(5,1) - Ifges(4,2) - Ifges(6,2) - 0.2e1 * t14 - t17 + ((-2 * mrSges(6,3) - t18) * pkin(7));
t9 = pkin(6) ^ 2;
t7 = cos(pkin(8));
t5 = sin(pkin(8));
t1 = t12 * pkin(4) + Ifges(5,5);
t8 = [(m(4) * t9) + Ifges(2,3) + (0.2e1 * t5 * (t13 * pkin(2) + Ifges(3,4)) + (Ifges(3,2) + (pkin(2) ^ 2 - t9) * m(4) + t11) * t7) * t7 - t11; mrSges(2,1); mrSges(2,2); ((m(4) * pkin(2) + mrSges(3,1)) * t7 + t5 * (-mrSges(3,2) + t13)) / t7; mrSges(3,3); m(3) + m(4); Ifges(4,1) - Ifges(4,2) - 0.4e1 * t14 + t2 - 0.2e1 * t17; 0.2e1 * t3 * Ifges(5,4) - t2 * t16 + Ifges(4,4) - Ifges(5,4); -t4 * Ifges(5,6) + t1 * t6 + Ifges(4,5); t6 * Ifges(5,6) + t1 * t4 + Ifges(4,6); Ifges(4,3) + Ifges(5,3) + t15; mrSges(4,1); mrSges(4,2); pkin(4) * m(6) + mrSges(5,1); mrSges(5,2); mrSges(5,3) - t12; m(5) + m(6); Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t8;
