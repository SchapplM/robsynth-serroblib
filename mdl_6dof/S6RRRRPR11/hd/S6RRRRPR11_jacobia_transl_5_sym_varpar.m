% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR11
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:36
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPR11_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR11_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR11_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR11_jacobia_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:36:21
% EndTime: 2019-02-26 22:36:21
% DurationCPUTime: 0.20s
% Computational Cost: add. (231->53), mult. (461->89), div. (0->0), fcn. (581->12), ass. (0->38)
t26 = sin(qJ(3));
t30 = cos(qJ(3));
t29 = cos(qJ(4));
t19 = t29 * pkin(4) + pkin(3);
t22 = qJ(4) + pkin(12);
t20 = sin(t22);
t21 = cos(t22);
t34 = t21 * r_i_i_C(1) - t20 * r_i_i_C(2) + t19;
t44 = r_i_i_C(3) + qJ(5) + pkin(10);
t45 = t44 * t26 + t34 * t30 + pkin(2);
t43 = cos(qJ(1));
t23 = sin(pkin(6));
t28 = sin(qJ(1));
t42 = t23 * t28;
t41 = t23 * t30;
t31 = cos(qJ(2));
t40 = t23 * t31;
t39 = cos(pkin(6));
t25 = sin(qJ(4));
t38 = t25 * pkin(4) + pkin(9);
t37 = t23 * t43;
t27 = sin(qJ(2));
t35 = t39 * t43;
t12 = t27 * t35 + t28 * t31;
t4 = t12 * t30 - t26 * t37;
t36 = t28 * t39;
t3 = t12 * t26 + t30 * t37;
t33 = t20 * r_i_i_C(1) + t21 * r_i_i_C(2) + t38;
t14 = -t27 * t36 + t43 * t31;
t13 = t43 * t27 + t31 * t36;
t11 = t28 * t27 - t31 * t35;
t10 = t39 * t26 + t27 * t41;
t9 = t23 * t27 * t26 - t39 * t30;
t8 = t14 * t30 + t26 * t42;
t7 = t14 * t26 - t28 * t41;
t2 = t13 * t20 + t8 * t21;
t1 = t13 * t21 - t8 * t20;
t5 = [-t28 * pkin(1) - t12 * pkin(2) + pkin(8) * t37 - t33 * t11 - t44 * t3 - t34 * t4, -t13 * t45 + t33 * t14, -t34 * t7 + t44 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2) + (t13 * t29 - t25 * t8) * pkin(4), t7, 0; t43 * pkin(1) + t14 * pkin(2) + pkin(8) * t42 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t38 * t13 + t8 * t19 + t44 * t7, -t11 * t45 + t33 * t12, -t34 * t3 + t44 * t4 (t11 * t21 - t4 * t20) * r_i_i_C(1) + (-t11 * t20 - t4 * t21) * r_i_i_C(2) + (t11 * t29 - t4 * t25) * pkin(4), t3, 0; 0 (t33 * t27 + t31 * t45) * t23, t44 * t10 - t34 * t9 (-t10 * t20 - t21 * t40) * r_i_i_C(1) + (-t10 * t21 + t20 * t40) * r_i_i_C(2) + (-t10 * t25 - t29 * t40) * pkin(4), t9, 0;];
Ja_transl  = t5;
