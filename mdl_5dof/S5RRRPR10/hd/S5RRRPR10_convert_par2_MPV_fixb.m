% Return the minimum parameter vector for
% S5RRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2021-01-15 23:43
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RRRPR10_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR10_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR10_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR10_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t14 = m(3) + m(4);
t34 = (t14 * pkin(7));
t5 = m(6) * pkin(9) + mrSges(6,3);
t25 = t5 * pkin(4) + Ifges(5,4);
t32 = (pkin(9) * mrSges(6,3));
t11 = sin(pkin(10));
t12 = cos(pkin(10));
t31 = t11 * t12;
t30 = Ifges(5,1) + Ifges(6,2);
t8 = 2 * t32;
t7 = t12 ^ 2;
t29 = 0.2e1 * t7;
t28 = t25 * t31;
t16 = (pkin(9) ^ 2);
t27 = (m(6) * t16 + t30);
t26 = m(4) * pkin(8) + mrSges(4,3);
t19 = (pkin(4) ^ 2);
t24 = ((-t16 + t19) * m(6) + Ifges(5,2) - t30);
t23 = Ifges(3,2) + (2 * mrSges(3,3) + t34) * pkin(7);
t9 = -2 * t32;
t3 = m(6) * t19 + Ifges(5,2) - t27 + t9;
t22 = (2 * pkin(8) * mrSges(4,3)) + t3 * t7 + Ifges(4,2) + t27 + 0.2e1 * t28 + t8;
t21 = (pkin(1) ^ 2);
t20 = (pkin(2) ^ 2);
t17 = pkin(8) ^ 2;
t13 = cos(pkin(5));
t1 = [(-m(4) * t20 - t23) * t13 ^ 2 + ((t20 + t21) * m(4)) + (m(3) * t21) + Ifges(2,3) + t23; t14 * pkin(1) + mrSges(2,1); (-mrSges(3,3) - t34) * sin(pkin(5)) + mrSges(2,2); ((t17 - t20) * m(4)) - Ifges(3,2) + Ifges(3,1) + t22; t26 * pkin(2) + Ifges(3,4); Ifges(3,5); Ifges(3,6); ((t17 + t20) * m(4)) + Ifges(3,3) + t22; m(4) * pkin(2) + mrSges(3,1); mrSges(3,2) - t26; (-t24 + t8) * t29 - 0.4e1 * t28 + t9 + Ifges(4,1) - Ifges(4,2) + t24; t25 * t29 - t3 * t31 + Ifges(4,4) - t25; t12 * Ifges(5,5) - t11 * Ifges(5,6) + Ifges(4,5); t11 * Ifges(5,5) + t12 * Ifges(5,6) + Ifges(4,6); Ifges(4,3) + Ifges(5,3) + Ifges(6,2) + t8 + (t16 + t19) * m(6); mrSges(4,1); mrSges(4,2); m(6) * pkin(4) + mrSges(5,1); mrSges(5,2) - t5; mrSges(5,3); m(5) + m(6); Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t1;
