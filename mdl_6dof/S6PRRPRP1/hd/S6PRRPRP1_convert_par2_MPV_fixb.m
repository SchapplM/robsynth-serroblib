% Return the minimum parameter vector for
% S6PRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% m [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MPV [26x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:40
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6PRRPRP1_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP1_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP1_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP1_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t10 = pkin(9) ^ 2;
t11 = pkin(4) ^ 2;
t18 = Ifges(6,2) + Ifges(7,2);
t15 = Ifges(5,1) + t18;
t12 = (Ifges(5,2) - t15 + (-t10 + t11) * m(6));
t19 = (pkin(9) * mrSges(6,3));
t1 = -2 * t19 + t12;
t3 = pkin(9) * m(6) + mrSges(6,3);
t14 = t3 * pkin(4) + Ifges(5,4);
t8 = sin(pkin(11));
t9 = cos(pkin(11));
t20 = t8 * t9;
t17 = t14 * t20;
t6 = 2 * t19;
t5 = t9 ^ 2;
t16 = 0.2e1 * t5;
t2 = [m(2) + m(3) + m(4); t1 * t5 + 0.2e1 * t17 + m(6) * t10 + t6 + ((pkin(2) ^ 2 + pkin(8) ^ 2) * m(4)) + (2 * pkin(8) * mrSges(4,3)) + Ifges(4,2) + Ifges(3,3) + t15; pkin(2) * m(4) + mrSges(3,1); -pkin(8) * m(4) + mrSges(3,2) - mrSges(4,3); (-t12 + t6) * t16 - 0.4e1 * t17 + Ifges(4,1) - Ifges(4,2) + t1; -t1 * t20 + t14 * t16 + Ifges(4,4) - t14; t9 * Ifges(5,5) - t8 * Ifges(5,6) + Ifges(4,5); t8 * Ifges(5,5) + t9 * Ifges(5,6) + Ifges(4,6); Ifges(4,3) + Ifges(5,3) + t6 + (t10 + t11) * m(6) + t18; mrSges(4,1); mrSges(4,2); pkin(4) * m(6) + mrSges(5,1); mrSges(5,2) - t3; mrSges(5,3); m(5) + m(6); Ifges(6,1) + Ifges(7,1) - t18; Ifges(6,4) + Ifges(7,4); Ifges(6,5) + Ifges(7,5); Ifges(6,6) + Ifges(7,6); Ifges(6,3) + Ifges(7,3); mrSges(6,1); mrSges(6,2); mrSges(7,1); mrSges(7,2); mrSges(7,3); m(7);];
MPV = t2;
