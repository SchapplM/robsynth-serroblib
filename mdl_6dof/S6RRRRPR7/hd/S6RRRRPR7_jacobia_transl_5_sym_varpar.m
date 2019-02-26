% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR7
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
% Datum: 2019-02-26 22:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPR7_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR7_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_jacobia_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:34:06
% EndTime: 2019-02-26 22:34:06
% DurationCPUTime: 0.15s
% Computational Cost: add. (218->48), mult. (259->70), div. (0->0), fcn. (312->12), ass. (0->36)
t28 = cos(pkin(6));
t29 = sin(qJ(2));
t32 = cos(qJ(1));
t35 = t32 * t29;
t30 = sin(qJ(1));
t31 = cos(qJ(2));
t36 = t30 * t31;
t10 = t28 * t35 + t36;
t26 = qJ(3) + qJ(4);
t21 = pkin(12) + t26;
t18 = sin(t21);
t27 = sin(pkin(6));
t38 = t27 * t32;
t13 = t18 * t38;
t19 = cos(t21);
t45 = (-t10 * t18 - t19 * t38) * r_i_i_C(1) + (-t10 * t19 + t13) * r_i_i_C(2);
t34 = t32 * t31;
t37 = t30 * t29;
t12 = -t28 * t37 + t34;
t39 = t27 * t30;
t5 = -t12 * t18 + t19 * t39;
t6 = t12 * t19 + t18 * t39;
t44 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
t40 = t27 * t29;
t43 = (-t18 * t40 + t28 * t19) * r_i_i_C(1) + (-t28 * t18 - t19 * t40) * r_i_i_C(2);
t22 = sin(t26);
t15 = pkin(4) * t22 + sin(qJ(3)) * pkin(3);
t42 = pkin(8) + t15;
t41 = r_i_i_C(3) + qJ(5) + pkin(10) + pkin(9);
t23 = cos(t26);
t16 = pkin(4) * t23 + cos(qJ(3)) * pkin(3);
t14 = pkin(2) + t16;
t33 = t19 * r_i_i_C(1) - t18 * r_i_i_C(2) + t14;
t11 = t28 * t36 + t35;
t9 = -t28 * t34 + t37;
t1 = [-t30 * pkin(1) + t13 * r_i_i_C(1) - t41 * t9 - t33 * t10 + (r_i_i_C(2) * t19 + t42) * t38, -t33 * t11 + t41 * t12, -t12 * t15 + t16 * t39 + t44 (-t12 * t22 + t23 * t39) * pkin(4) + t44, t11, 0; t32 * pkin(1) + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t41 * t11 + t12 * t14 + t42 * t39, t41 * t10 - t33 * t9, -t10 * t15 - t16 * t38 + t45 (-t10 * t22 - t23 * t38) * pkin(4) + t45, t9, 0; 0 (t41 * t29 + t33 * t31) * t27, -t15 * t40 + t28 * t16 + t43 (-t22 * t40 + t23 * t28) * pkin(4) + t43, -t27 * t31, 0;];
Ja_transl  = t1;
