% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRR9_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR9_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_jacobia_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:58:50
% EndTime: 2019-02-26 21:58:50
% DurationCPUTime: 0.14s
% Computational Cost: add. (198->42), mult. (232->64), div. (0->0), fcn. (283->12), ass. (0->34)
t25 = cos(pkin(6));
t26 = sin(qJ(2));
t29 = cos(qJ(1));
t32 = t29 * t26;
t27 = sin(qJ(1));
t28 = cos(qJ(2));
t33 = t27 * t28;
t10 = t25 * t32 + t33;
t23 = pkin(12) + qJ(4);
t21 = qJ(5) + t23;
t17 = sin(t21);
t24 = sin(pkin(6));
t35 = t24 * t29;
t14 = t17 * t35;
t18 = cos(t21);
t42 = (-t10 * t17 - t18 * t35) * r_i_i_C(1) + (-t10 * t18 + t14) * r_i_i_C(2);
t31 = t29 * t28;
t34 = t27 * t26;
t12 = -t25 * t34 + t31;
t36 = t24 * t27;
t5 = -t12 * t17 + t18 * t36;
t6 = t12 * t18 + t17 * t36;
t41 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
t37 = t24 * t26;
t40 = (-t17 * t37 + t25 * t18) * r_i_i_C(1) + (-t25 * t17 - t18 * t37) * r_i_i_C(2);
t19 = sin(t23);
t39 = pkin(8) + pkin(4) * t19 + sin(pkin(12)) * pkin(3);
t38 = r_i_i_C(3) + pkin(10) + pkin(9) + qJ(3);
t20 = cos(t23);
t13 = pkin(4) * t20 + cos(pkin(12)) * pkin(3) + pkin(2);
t30 = t18 * r_i_i_C(1) - t17 * r_i_i_C(2) + t13;
t11 = t25 * t33 + t32;
t9 = -t25 * t31 + t34;
t1 = [-t27 * pkin(1) + t14 * r_i_i_C(1) - t38 * t9 - t30 * t10 + (r_i_i_C(2) * t18 + t39) * t35, -t30 * t11 + t12 * t38, t11 (-t12 * t19 + t20 * t36) * pkin(4) + t41, t41, 0; t29 * pkin(1) + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t38 * t11 + t12 * t13 + t39 * t36, t10 * t38 - t30 * t9, t9 (-t10 * t19 - t20 * t35) * pkin(4) + t42, t42, 0; 0 (t26 * t38 + t30 * t28) * t24, -t24 * t28 (-t19 * t37 + t20 * t25) * pkin(4) + t40, t40, 0;];
Ja_transl  = t1;
