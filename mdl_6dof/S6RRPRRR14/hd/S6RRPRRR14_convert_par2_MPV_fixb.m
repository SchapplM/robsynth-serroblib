% Return the minimum parameter vector for
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
% m_mdh [7x1]
%   mass of all robot links (leg links until cut joint, platform)
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
% MPV [35x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2019-01-03 10:26
% Revision: 5fdbc45bcf2cc60deefd7ac2d71d743ed41bf7e4 (2018-12-21)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function MPV = S6RRPRRR14_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_convert_par2_MPV_fixb: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR14_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR14_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t241 = sin(pkin(8));
t245 = cos(pkin(8));
t247 = (m(6) + m(7));
t253 = (pkin(4) ^ 2);
t228 = (t253 * t247 + Ifges(5,2));
t237 = m(5) + t247;
t280 = (pkin(11) * t237);
t256 = t228 + (2 * mrSges(5,3) + t280) * pkin(11);
t219 = t256 * t245 * t241 + Ifges(4,6);
t262 = mrSges(5,3) + t280;
t278 = t262 * t245;
t223 = -pkin(3) * t278 + Ifges(4,5);
t240 = sin(pkin(14));
t244 = cos(pkin(14));
t257 = t219 * t244 + t223 * t240;
t242 = sin(pkin(7));
t246 = cos(pkin(7));
t269 = t242 * t246;
t281 = t257 * t269;
t279 = t262 * t241;
t268 = pkin(3) ^ 2 * t237;
t220 = t256 * t241 ^ 2 + Ifges(4,3) + t268;
t232 = t242 ^ 2;
t236 = t246 ^ 2;
t221 = t256 * t245 ^ 2 + Ifges(4,2) + t268;
t225 = Ifges(4,1) + t256;
t230 = t240 ^ 2;
t234 = t244 ^ 2;
t224 = pkin(3) * t279 + Ifges(4,4);
t263 = t244 * t240 * t224;
t255 = t221 * t234 + t225 * t230 + 0.2e1 * t263;
t277 = t232 * t220 + t255 * t236 + Ifges(3,2);
t276 = (m(3) * pkin(10));
t252 = (pkin(5) ^ 2);
t275 = (t252 * m(7) + Ifges(6,2));
t274 = 2 * pkin(13) * mrSges(7,3) + Ifges(7,2);
t271 = t223 * t244;
t270 = t224 * (t234 - t230);
t267 = 0.2e1 * t281;
t265 = 2 * pkin(12) * mrSges(6,3) + t275;
t264 = pkin(13) * m(7) + mrSges(7,3);
t261 = pkin(12) * t247 + mrSges(6,3);
t259 = (-t221 + t225) * t244;
t250 = pkin(12) ^ 2;
t249 = pkin(13) ^ 2;
t243 = sin(pkin(6));
t1 = [pkin(1) ^ 2 * m(3) + Ifges(2,3) + (-0.2e1 * t281 + ((2 * mrSges(3,3) + t276) * pkin(10)) + t277) * t243 ^ 2; m(3) * pkin(1) + mrSges(2,1); mrSges(2,2) + (-mrSges(3,3) - t276) * t243; t230 * t221 + t234 * t225 + Ifges(3,1) - 0.2e1 * t263 + t267 - t277; -t242 * t271 + t246 * t270 + Ifges(3,4) + (t219 * t242 + t246 * t259) * t240; t246 * t271 + t242 * t270 + Ifges(3,5) + (-t219 * t246 + t242 * t259) * t240; Ifges(3,6) + t257 * (t236 - t232) + (-t220 + t255) * t269; t236 * t220 + t255 * t232 + Ifges(3,3) + t267; mrSges(3,1); mrSges(3,2); pkin(3) * t237 + mrSges(4,1); mrSges(4,2) - t279; mrSges(4,3) + t278; m(4) + t237; t250 * t247 + Ifges(5,1) - t228 + t265; t261 * pkin(4) + Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3) + (t250 + t253) * t247 + t265; pkin(4) * t247 + mrSges(5,1); mrSges(5,2) - t261; m(7) * t249 + Ifges(6,1) + t274 - t275; t264 * pkin(5) + Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3) + (t249 + t252) * m(7) + t274; m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) - t264; Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
