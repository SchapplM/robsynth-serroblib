% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR8
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPPR8_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR8_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:07:39
% EndTime: 2019-02-26 22:07:39
% DurationCPUTime: 0.18s
% Computational Cost: add. (223->48), mult. (556->78), div. (0->0), fcn. (707->10), ass. (0->34)
t21 = sin(qJ(3));
t25 = cos(qJ(3));
t20 = sin(qJ(6));
t24 = cos(qJ(6));
t37 = pkin(5) + qJ(4);
t30 = t24 * r_i_i_C(1) - t20 * r_i_i_C(2) + t37;
t34 = pkin(3) + pkin(4) + pkin(10) + r_i_i_C(3);
t42 = t30 * t21 + t34 * t25 + pkin(2);
t41 = cos(qJ(1));
t19 = sin(pkin(6));
t23 = sin(qJ(1));
t40 = t19 * t23;
t39 = t19 * t25;
t26 = cos(qJ(2));
t38 = t19 * t26;
t36 = pkin(9) - qJ(5);
t35 = cos(pkin(6));
t33 = t19 * t41;
t32 = t23 * t35;
t31 = t35 * t41;
t29 = t20 * r_i_i_C(1) + t24 * r_i_i_C(2) - t36;
t22 = sin(qJ(2));
t12 = t22 * t31 + t23 * t26;
t3 = t12 * t21 + t25 * t33;
t28 = t12 * t25 - t21 * t33;
t14 = -t22 * t32 + t41 * t26;
t13 = t41 * t22 + t26 * t32;
t11 = t23 * t22 - t26 * t31;
t9 = t19 * t22 * t21 - t35 * t25;
t8 = t14 * t25 + t21 * t40;
t7 = t14 * t21 - t23 * t39;
t2 = -t13 * t20 + t7 * t24;
t1 = -t13 * t24 - t7 * t20;
t4 = [-t23 * pkin(1) - t12 * pkin(2) + pkin(8) * t33 + t29 * t11 - t34 * t28 - t30 * t3, -t13 * t42 - t29 * t14, t30 * t8 - t34 * t7, t7, -t13, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t41 * pkin(1) + t14 * pkin(2) + pkin(8) * t40 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t36 * t13 + t34 * t8 + t37 * t7, -t11 * t42 - t29 * t12, t30 * t28 - t34 * t3, t3, -t11 (-t11 * t24 - t3 * t20) * r_i_i_C(1) + (t11 * t20 - t3 * t24) * r_i_i_C(2); 0 (-t29 * t22 + t42 * t26) * t19, -t34 * t9 + t30 * (t35 * t21 + t22 * t39) t9, t38 (-t9 * t20 + t24 * t38) * r_i_i_C(1) + (-t20 * t38 - t9 * t24) * r_i_i_C(2);];
Ja_transl  = t4;
