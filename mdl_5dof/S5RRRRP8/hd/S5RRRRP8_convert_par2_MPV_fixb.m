% Return the minimum parameter vector for
% S5RRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2021-01-16 00:22
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RRRRP8_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_convert_par2_MPV_fixb: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP8_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP8_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP8_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t22 = 2 * pkin(8);
t20 = (m(4) + m(5));
t3 = m(3) + t20;
t21 = (t3 * pkin(6));
t11 = (pkin(2) ^ 2);
t19 = (m(4) * t11);
t10 = (pkin(3) ^ 2);
t18 = (t10 * m(5));
t17 = (mrSges(4,3) + mrSges(5,3));
t16 = (-Ifges(5,2) - Ifges(6,2));
t9 = (pkin(7) ^ 2);
t15 = pkin(7) * t22 + pkin(8) ^ 2 + t10 + t9;
t14 = m(4) * t9 + mrSges(5,3) * t22 + 2 * t17 * pkin(7) + Ifges(4,2) - t16;
t13 = pkin(7) * m(4) + (pkin(8) + pkin(7)) * m(5) + t17;
t12 = (pkin(1) ^ 2);
t1 = [m(3) * t12 + Ifges(3,2) + Ifges(2,3) + t20 * (t11 + t12) + (2 * mrSges(3,3) + t21) * pkin(6); t3 * pkin(1) + mrSges(2,1); mrSges(2,2) - mrSges(3,3) - t21; (-t11 + t15) * m(5) - t19 - Ifges(3,2) + Ifges(3,1) + t14; t13 * pkin(2) + Ifges(3,4); Ifges(3,5); Ifges(3,6); (t11 + t15) * m(5) + t19 + Ifges(3,3) + t14; t20 * pkin(2) + mrSges(3,1); mrSges(3,2) - t13; Ifges(4,1) - Ifges(4,2) - t18; Ifges(4,4); (-m(5) * pkin(8) - mrSges(5,3)) * pkin(3) + Ifges(4,5); Ifges(4,6); Ifges(4,3) + t18; pkin(3) * m(5) + mrSges(4,1); mrSges(4,2); Ifges(5,1) + Ifges(6,1) + t16; Ifges(5,4) + Ifges(6,4); Ifges(5,5) + Ifges(6,5); Ifges(5,6) + Ifges(6,6); Ifges(5,3) + Ifges(6,3); mrSges(5,1); mrSges(5,2); mrSges(6,1); mrSges(6,2); mrSges(6,3); m(6);];
MPV = t1;
