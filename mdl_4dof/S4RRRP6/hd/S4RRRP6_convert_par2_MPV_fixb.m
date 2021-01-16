% Return the minimum parameter vector for
% S4RRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% MPV [21x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:39
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S4RRRP6_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_convert_par2_MPV_fixb: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP6_convert_par2_MPV_fixb: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP6_convert_par2_MPV_fixb: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP6_convert_par2_MPV_fixb: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t10 = (pkin(1) ^ 2 + pkin(5) ^ 2);
t9 = (-Ifges(4,2) - Ifges(5,2));
t8 = 2 * pkin(6) * mrSges(4,3) - t9;
t7 = pkin(6) * m(4) + mrSges(4,3);
t5 = (pkin(2) ^ 2);
t3 = pkin(6) ^ 2;
t2 = m(3) + m(4);
t1 = [(t5 + t10) * m(4) + 2 * pkin(5) * mrSges(3,3) + Ifges(3,2) + Ifges(2,3) + t10 * m(3); t2 * pkin(1) + mrSges(2,1); -t2 * pkin(5) + mrSges(2,2) - mrSges(3,3); Ifges(3,1) - Ifges(3,2) + (t3 - t5) * m(4) + t8; t7 * pkin(2) + Ifges(3,4); Ifges(3,5); Ifges(3,6); Ifges(3,3) + (t3 + t5) * m(4) + t8; pkin(2) * m(4) + mrSges(3,1); mrSges(3,2) - t7; Ifges(4,1) + Ifges(5,1) + t9; Ifges(4,4) + Ifges(5,4); Ifges(4,5) + Ifges(5,5); Ifges(4,6) + Ifges(5,6); Ifges(4,3) + Ifges(5,3); mrSges(4,1); mrSges(4,2); mrSges(5,1); mrSges(5,2); mrSges(5,3); m(5);];
MPV = t1;
