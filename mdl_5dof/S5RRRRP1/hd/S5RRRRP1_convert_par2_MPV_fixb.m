% Return the minimum parameter vector for
% S5RRRRP1
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
% Datum: 2021-01-15 23:53
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RRRRP1_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_convert_par2_MPV_fixb: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP1_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP1_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP1_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t16 = 2 * pkin(8);
t2 = (m(4) + m(5));
t7 = (pkin(2) ^ 2);
t15 = (t2 * t7);
t6 = (pkin(3) ^ 2);
t14 = (t6 * m(5));
t1 = pkin(7) + pkin(8);
t13 = (mrSges(4,3) + mrSges(5,3));
t12 = (-Ifges(5,2) - Ifges(6,2));
t8 = (pkin(1) ^ 2);
t11 = pkin(7) ^ 2 + t7 + t8;
t10 = -m(4) * pkin(7) - t13;
t9 = -t1 * m(5) + t10;
t5 = (pkin(6) ^ 2);
t3 = m(3) + m(4);
t4 = [(pkin(7) * t16 + pkin(8) ^ 2 + t11 + t5 + t6) * m(5) + t3 * t5 + m(3) * t8 + mrSges(5,3) * t16 + Ifges(3,2) + Ifges(4,2) + Ifges(2,3) + t11 * m(4) + 2 * t13 * pkin(7) + 2 * (mrSges(3,3) - t9) * pkin(6) - t12; (m(3) + t2) * pkin(1) + mrSges(2,1); (-pkin(6) - t1) * m(5) - t3 * pkin(6) + mrSges(2,2) - mrSges(3,3) + t10; Ifges(3,1) - Ifges(3,2) - t15; Ifges(3,4); t9 * pkin(2) + Ifges(3,5); Ifges(3,6); Ifges(3,3) + t15; t2 * pkin(2) + mrSges(3,1); mrSges(3,2); Ifges(4,1) - Ifges(4,2) - t14; Ifges(4,4); (-m(5) * pkin(8) - mrSges(5,3)) * pkin(3) + Ifges(4,5); Ifges(4,6); Ifges(4,3) + t14; pkin(3) * m(5) + mrSges(4,1); mrSges(4,2); Ifges(5,1) + Ifges(6,1) + t12; Ifges(5,4) + Ifges(6,4); Ifges(5,5) + Ifges(6,5); Ifges(5,6) + Ifges(6,6); Ifges(5,3) + Ifges(6,3); mrSges(5,1); mrSges(5,2); mrSges(6,1); mrSges(6,2); mrSges(6,3); m(6);];
MPV = t4;
