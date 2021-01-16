% Return the minimum parameter vector for
% S5RRRPP4
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
% MPV [25x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:26
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RRRPP4_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_convert_par2_MPV_fixb: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP4_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP4_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP4_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t18 = 2 * pkin(6);
t12 = Ifges(5,1) + Ifges(6,1);
t1 = -Ifges(5,2) - Ifges(6,3) + t12;
t4 = cos(pkin(8));
t2 = t4 ^ 2;
t17 = t1 * t2;
t3 = sin(pkin(8));
t16 = t3 * t4;
t9 = (pkin(2) ^ 2);
t15 = t9 * m(4);
t14 = pkin(1) ^ 2 + pkin(6) ^ 2;
t13 = (mrSges(3,3) + mrSges(4,3));
t7 = Ifges(5,4) - Ifges(6,5);
t11 = t7 * t16;
t6 = Ifges(6,4) + Ifges(5,5);
t5 = Ifges(5,6) - Ifges(6,6);
t8 = [-t17 + 0.2e1 * t11 + ((t9 + (t18 + pkin(7)) * pkin(7) + t14) * m(4)) + (t13 * t18) + (2 * mrSges(4,3) * pkin(7)) + Ifges(3,2) + Ifges(4,2) + Ifges(2,3) + (t14 * m(3)) + t12; (m(3) + m(4)) * pkin(1) + mrSges(2,1); mrSges(2,2) - pkin(6) * m(3) + (-pkin(7) - pkin(6)) * m(4) - t13; Ifges(3,1) - Ifges(3,2) - t15; Ifges(3,4); (-m(4) * pkin(7) - mrSges(4,3)) * pkin(2) + Ifges(3,5); Ifges(3,6); Ifges(3,3) + t15; pkin(2) * m(4) + mrSges(3,1); mrSges(3,2); Ifges(4,1) - Ifges(4,2) - t1 - 0.4e1 * t11 + 0.2e1 * t17; t1 * t16 + 0.2e1 * t2 * t7 + Ifges(4,4) - t7; -t5 * t3 + t6 * t4 + Ifges(4,5); t6 * t3 + t5 * t4 + Ifges(4,6); Ifges(4,3) + Ifges(5,3) + Ifges(6,2); mrSges(4,1); mrSges(4,2); mrSges(5,1); mrSges(5,2); mrSges(5,3); m(5); mrSges(6,1); mrSges(6,2); mrSges(6,3); m(6);];
MPV = t8;
