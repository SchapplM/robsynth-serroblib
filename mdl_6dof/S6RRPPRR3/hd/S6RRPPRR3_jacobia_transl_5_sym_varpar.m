% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPRR3_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR3_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_jacobia_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:29:33
% EndTime: 2019-02-26 21:29:33
% DurationCPUTime: 0.14s
% Computational Cost: add. (173->46), mult. (344->72), div. (0->0), fcn. (446->12), ass. (0->35)
t25 = sin(pkin(11));
t27 = cos(pkin(11));
t30 = sin(qJ(2));
t32 = cos(qJ(2));
t14 = t30 * t25 - t32 * t27;
t45 = pkin(4) * sin(pkin(12));
t44 = t32 * pkin(2);
t43 = r_i_i_C(3) + pkin(9) + qJ(4);
t28 = cos(pkin(6));
t42 = t28 * t32;
t26 = sin(pkin(6));
t31 = sin(qJ(1));
t40 = t31 * t26;
t33 = cos(qJ(1));
t38 = t33 * t26;
t36 = t32 * t25 + t30 * t27;
t12 = t36 * t28;
t5 = t33 * t12 - t31 * t14;
t37 = t31 * t12 + t33 * t14;
t19 = cos(pkin(12)) * pkin(4) + pkin(3);
t23 = pkin(12) + qJ(5);
t21 = sin(t23);
t22 = cos(t23);
t35 = t22 * r_i_i_C(1) - t21 * r_i_i_C(2) + t19;
t34 = t14 * t28;
t20 = pkin(1) + t44;
t16 = t21 * t38;
t13 = t28 * t30 * pkin(2) + (-pkin(8) - qJ(3)) * t26;
t11 = t36 * t26;
t10 = t14 * t26;
t7 = t31 * t34 - t33 * t36;
t4 = -t31 * t36 - t33 * t34;
t2 = t21 * t40 - t22 * t37;
t1 = t21 * t37 + t22 * t40;
t3 = [t16 * r_i_i_C(1) - t31 * t20 - t35 * t5 + t43 * t4 + (-t13 + (t22 * r_i_i_C(2) + t45) * t26) * t33, -t43 * t37 + (-t30 * t33 - t31 * t42) * pkin(2) + t35 * t7, t40, -t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t37 * t19 + t33 * t20 - t43 * t7 + (t26 * t45 - t13) * t31, t43 * t5 + (-t30 * t31 + t33 * t42) * pkin(2) + t35 * t4, -t38, -t4 (-t5 * t21 - t22 * t38) * r_i_i_C(1) + (-t5 * t22 + t16) * r_i_i_C(2), 0; 0, -t35 * t10 + t43 * t11 + t26 * t44, t28, t10 (-t11 * t21 + t28 * t22) * r_i_i_C(1) + (-t11 * t22 - t28 * t21) * r_i_i_C(2), 0;];
Ja_transl  = t3;
