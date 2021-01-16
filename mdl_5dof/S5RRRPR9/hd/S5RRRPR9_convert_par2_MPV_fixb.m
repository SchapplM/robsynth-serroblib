% Return the minimum parameter vector for
% S5RRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% MPV [28x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:25
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RRRPR9_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_convert_par2_MPV_fixb: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR9_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR9_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR9_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t24 = (pkin(8) * m(6));
t11 = pkin(4) ^ 2 * m(6);
t4 = -Ifges(5,1) + Ifges(5,2) + t11;
t9 = cos(pkin(9));
t5 = t9 ^ 2;
t1 = t4 * t5;
t8 = sin(pkin(9));
t23 = t8 * t9;
t22 = (pkin(1) ^ 2 + pkin(6) ^ 2);
t21 = Ifges(5,4) * t23;
t20 = pkin(7) * m(4) + mrSges(4,3);
t19 = -mrSges(6,3) - t24;
t18 = (2 * pkin(7) * mrSges(4,3)) + Ifges(5,1) + Ifges(4,2) + Ifges(6,2) + t1 + 0.2e1 * t21 + ((2 * mrSges(6,3) + t24) * pkin(8));
t16 = (pkin(2) ^ 2);
t13 = pkin(7) ^ 2;
t10 = (m(3) + m(4));
t3 = t19 * pkin(4) + Ifges(5,5);
t2 = [(t16 + t22) * m(4) + 2 * pkin(6) * mrSges(3,3) + Ifges(3,2) + Ifges(2,3) + t22 * m(3); t10 * pkin(1) + mrSges(2,1); -t10 * pkin(6) + mrSges(2,2) - mrSges(3,3); ((t13 - t16) * m(4)) - Ifges(3,2) + Ifges(3,1) + t18; t20 * pkin(2) + Ifges(3,4); Ifges(3,5); Ifges(3,6); ((t13 + t16) * m(4)) + Ifges(3,3) + t18; pkin(2) * m(4) + mrSges(3,1); mrSges(3,2) - t20; Ifges(4,1) - Ifges(4,2) - 0.4e1 * t21 + t4 - 0.2e1 * t1; 0.2e1 * t5 * Ifges(5,4) - t4 * t23 + Ifges(4,4) - Ifges(5,4); -t8 * Ifges(5,6) + t3 * t9 + Ifges(4,5); t9 * Ifges(5,6) + t3 * t8 + Ifges(4,6); Ifges(4,3) + Ifges(5,3) + t11; mrSges(4,1); mrSges(4,2); pkin(4) * m(6) + mrSges(5,1); mrSges(5,2); mrSges(5,3) - t19; m(5) + m(6); Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t2;
