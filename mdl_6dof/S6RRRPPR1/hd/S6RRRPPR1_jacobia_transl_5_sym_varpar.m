% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:03
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPPR1_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR1_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:03:28
% EndTime: 2019-02-26 22:03:28
% DurationCPUTime: 0.12s
% Computational Cost: add. (153->24), mult. (113->25), div. (0->0), fcn. (122->10), ass. (0->21)
t45 = r_i_i_C(3) + qJ(5);
t20 = qJ(2) + qJ(3);
t16 = pkin(10) + t20;
t13 = sin(t16);
t14 = cos(t16);
t22 = cos(pkin(11));
t31 = -r_i_i_C(1) * t22 - pkin(4);
t21 = sin(pkin(11));
t38 = r_i_i_C(2) * t21;
t25 = t14 * (-t31 - t38) + pkin(3) * cos(t20) + t45 * t13;
t44 = cos(qJ(2)) * pkin(2) + t25;
t43 = t13 * t38 + t45 * t14;
t26 = t31 * t13 - pkin(3) * sin(t20);
t41 = pkin(1) + t44;
t23 = sin(qJ(1));
t34 = t43 * t23;
t24 = cos(qJ(1));
t33 = t43 * t24;
t28 = t21 * r_i_i_C(1) + t22 * r_i_i_C(2) + pkin(7) + pkin(8) + qJ(4);
t27 = -sin(qJ(2)) * pkin(2) + t26;
t1 = [-t41 * t23 + t28 * t24, t24 * t27 + t33, t24 * t26 + t33, t23, t24 * t13, 0; t28 * t23 + t41 * t24, t23 * t27 + t34, t23 * t26 + t34, -t24, t23 * t13, 0; 0, t44, t25, 0, -t14, 0;];
Ja_transl  = t1;
