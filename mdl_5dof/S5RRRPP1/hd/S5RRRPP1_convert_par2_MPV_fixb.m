% Return the minimum parameter vector for
% S5RRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% MPV [21x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:15
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RRRPP1_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_convert_par2_MPV_fixb: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP1_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP1_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP1_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t10 = -Ifges(5,1) - Ifges(6,1);
t1 = Ifges(5,2) + Ifges(6,3) + t10;
t4 = cos(pkin(8));
t2 = t4 ^ 2;
t12 = t1 * t2;
t3 = sin(pkin(8));
t11 = t3 * t4;
t7 = Ifges(5,4) - Ifges(6,5);
t9 = t7 * t11;
t8 = m(3) + m(4);
t6 = Ifges(6,4) + Ifges(5,5);
t5 = Ifges(5,6) - Ifges(6,6);
t13 = [pkin(1) ^ 2 * t8 + Ifges(2,3); pkin(1) * t8 + mrSges(2,1); mrSges(2,2); t12 + 0.2e1 * t9 + ((pkin(2) ^ 2 + pkin(7) ^ 2) * m(4)) + (2 * pkin(7) * mrSges(4,3)) + Ifges(3,3) + Ifges(4,2) - t10; pkin(2) * m(4) + mrSges(3,1); -pkin(7) * m(4) + mrSges(3,2) - mrSges(4,3); Ifges(4,1) - Ifges(4,2) + t1 - 0.4e1 * t9 - 0.2e1 * t12; -t1 * t11 + 0.2e1 * t2 * t7 + Ifges(4,4) - t7; -t3 * t5 + t4 * t6 + Ifges(4,5); t3 * t6 + t4 * t5 + Ifges(4,6); Ifges(4,3) + Ifges(5,3) + Ifges(6,2); mrSges(4,1); mrSges(4,2); mrSges(5,1); mrSges(5,2); mrSges(5,3); m(5); mrSges(6,1); mrSges(6,2); mrSges(6,3); m(6);];
MPV = t13;
