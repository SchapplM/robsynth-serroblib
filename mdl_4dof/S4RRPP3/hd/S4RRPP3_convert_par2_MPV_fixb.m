% Return the minimum parameter vector for
% S4RRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% m [5x1]
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
% MPV [18x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:36
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S4RRPP3_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_convert_par2_MPV_fixb: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP3_convert_par2_MPV_fixb: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP3_convert_par2_MPV_fixb: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP3_convert_par2_MPV_fixb: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t9 = -Ifges(4,1) - Ifges(5,1);
t1 = Ifges(4,2) + Ifges(5,3) + t9;
t4 = cos(pkin(6));
t2 = t4 ^ 2;
t11 = t1 * t2;
t3 = sin(pkin(6));
t10 = t3 * t4;
t7 = Ifges(4,4) - Ifges(5,5);
t8 = t7 * t10;
t6 = Ifges(5,4) + Ifges(4,5);
t5 = Ifges(4,6) - Ifges(5,6);
t12 = [t11 + 0.2e1 * t8 + ((pkin(1) ^ 2 + pkin(5) ^ 2) * m(3)) + (2 * pkin(5) * mrSges(3,3)) + Ifges(2,3) + Ifges(3,2) - t9; pkin(1) * m(3) + mrSges(2,1); -pkin(5) * m(3) + mrSges(2,2) - mrSges(3,3); Ifges(3,1) - Ifges(3,2) + t1 - 0.4e1 * t8 - 0.2e1 * t11; -t1 * t10 + 0.2e1 * t2 * t7 + Ifges(3,4) - t7; -t3 * t5 + t4 * t6 + Ifges(3,5); t3 * t6 + t4 * t5 + Ifges(3,6); Ifges(3,3) + Ifges(4,3) + Ifges(5,2); mrSges(3,1); mrSges(3,2); mrSges(4,1); mrSges(4,2); mrSges(4,3); m(4); mrSges(5,1); mrSges(5,2); mrSges(5,3); m(5);];
MPV = t12;
