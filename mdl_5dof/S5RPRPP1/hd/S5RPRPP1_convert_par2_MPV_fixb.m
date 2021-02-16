% Return the minimum parameter vector for
% S5RPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
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
% MPV [19x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:14
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RPRPP1_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_convert_par2_MPV_fixb: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP1_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP1_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP1_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t14 = -Ifges(5,1) - Ifges(6,1);
t2 = Ifges(5,2) + Ifges(6,3) + t14;
t7 = cos(pkin(8));
t4 = t7 ^ 2;
t16 = t2 * t4;
t5 = sin(pkin(8));
t15 = t5 * t7;
t11 = Ifges(5,4) - Ifges(6,5);
t13 = t11 * t15;
t1 = m(4) * pkin(6) - mrSges(3,2) + mrSges(4,3);
t3 = m(4) * pkin(2) + mrSges(3,1);
t6 = sin(pkin(7));
t8 = cos(pkin(7));
t12 = t1 * t6 + t3 * t8;
t10 = Ifges(6,4) + Ifges(5,5);
t9 = Ifges(5,6) - Ifges(6,6);
t17 = [t16 + 0.2e1 * t13 + ((pkin(2) ^ 2 + pkin(6) ^ 2) * m(4)) + (2 * mrSges(4,3) * pkin(6)) + Ifges(4,2) + Ifges(2,3) + Ifges(3,3) + 0.2e1 * t12 * pkin(1) - t14; mrSges(2,1) + t12; -t1 * t8 + t3 * t6 + mrSges(2,2); m(3) + m(4); Ifges(4,1) - Ifges(4,2) - 0.4e1 * t13 + t2 - 0.2e1 * t16; 0.2e1 * t11 * t4 - t2 * t15 + Ifges(4,4) - t11; t10 * t7 - t9 * t5 + Ifges(4,5); t10 * t5 + t9 * t7 + Ifges(4,6); Ifges(4,3) + Ifges(5,3) + Ifges(6,2); mrSges(4,1); mrSges(4,2); mrSges(5,1); mrSges(5,2); mrSges(5,3); m(5); mrSges(6,1); mrSges(6,2); mrSges(6,3); m(6);];
MPV = t17;
